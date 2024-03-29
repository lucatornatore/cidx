
#include "cidx.h"



int main( int argc, char **argv)
{
  
  double tstart, telapsed;
  double tbegin;
  FILE *timings = fopen("timings", "w");
  details = fopen("details", "w");
   
 #if defined(USE_MPI)
  {
   #error "MPI support not yet implemented"
    int mpi_provided_thread_level;
    MPI_Init_threads( &argc, &argv, MPI_THREAD_FUNNELED, &mpi_provided_thread_level);
    if ( mpi_provided_thread_level < MPI_THREAD_FUNNELED )
      {
 	printf("a problem arise when asking for MPI_THREAD_FUNNELED level\n");
 	MPI_Finalize();
 	exit( 1 );
      }
    MPI_Comm_size(MPI_COMM_WORLD, &Ntasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &MyPid);
  }
 #endif

 #pragma omp parallel
  {
    me = omp_get_thread_num();
   #pragma omp master
    {
      Nthreads = omp_get_num_threads();
      printf("using %d threads\n", Nthreads);
    }    
  }

 #if defined(DEBUG)
 #pragma omp single
  printf(">> you're running a version with DEBUG code\n");
 #endif


  // ------------------------------------
  // get the command-line arguments
  //

  // prepare some defulat values of variables
  //
  verbose = 0;
  int ret;
  action = UNDEFINED;
  
  sprintf( snapnum, "---");
  sprintf( working_dir, "%s", working_dir_def );
  sprintf( working_dir_snap_base, "%s", working_dir_snap_base_def );
  sprintf( working_dir_subf_base, "%s", working_dir_subf_base_def );
  sprintf( snap_base, "%s", snap_base_def );
  sprintf( subf_base, "%s", subf_base_def );
  n_stars_generations = n_stars_generations_def;

  // process the command line
  //
  
  {
    int point = 1;
    while ( point < argc )
      {
	int a = 0;
	while( (a < n_args) && (strcasecmp(*(argv+point), args[a].argument) != 0) )
	  a++;
	if( a < n_args )
	  point = args[a].action( argv, point, argc );
	else
	  printf("unknown argument %d: %s\n", point, *(argv+point));
	point++;
      }
  }

  // final processing of variables
  //
  
  if ( ((action == UNDEFINED) || (strcmp(snapnum, "---")==0)) ) {
    if( !help_given )
      printf("Either you did not specify an action or "
	     "you did not specify a snapshot number.\n"
	     "I would not know what else to do more than suggesting to try -h.\n\n" );
    exit(1); }

  if( (action & SEARCH_PARTICLES) && (Nlists == 0) ) {
    printf("you did not specify any list of particles to look for in catalogs\n");
    exit(2); }

 #if defined(DEBUG)
  action |= BUILD_FOF_TABLE;
 #endif

  
  if( (strcmp( working_dir_snap_base, ".") == 0) ||
      (strcmp( working_dir_snap_base, "./") == 0) )
    snprintf( working_dir_snap, DIR_SIZE, "%s",
	      working_dir_snap_base);
  else
    snprintf( working_dir_snap, DIR_SIZE+NUM_SIZE, "%s%s",
	      working_dir_snap_base, snapnum);

  if( (strcmp( working_dir_subf_base, ".") == 0) ||
      (strcmp( working_dir_subf_base, "./") == 0) )
    snprintf( working_dir_subf, DIR_SIZE, "%s",
	      working_dir_subf_base);
  else
    snprintf( working_dir_subf, DIR_SIZE+NUM_SIZE, "%s%s",
	      working_dir_subf_base, snapnum);
  
  snprintf( snap_name, NAME_SIZE+NUM_SIZE, "%s%s", snap_base, snapnum );
  snprintf( subf_name, NAME_SIZE+NUM_SIZE, "%s%s", subf_base, snapnum );
  snprintf( catalog_name, CATALOG_NAME_SIZE, "%s/%s_%s_type_",
	    working_dir_subf, snapnum, catalog_base );
  snprintf( catalog_table_name, CATALOG_NAME_SIZE, "%s/%s_%s_table_fof",
	    working_dir_subf, snapnum, catalog_base );

  DPRINT(3,0, "working dir is %s\n"
	 "snap dir is    %s\n"
	 "subf dir is    %s\n"
	 "snap name is   %s\n"
	 "subf name is   %s\n",
	 working_dir, working_dir_snap, working_dir_subf, snap_name, subf_name );

  id_mask = get_stellargenerations_mask( n_stars_generations, &id_bitshift );      
  for( int t = 0; t < NTYPES; t++ )
    P_all[t] = (particle_t**)calloc( Nthreads, sizeof(particle_t*));

  tbegin = CPU_RTIME;
  

  
  // if at command line it was specified to create catalogs,
  // then let's do it
  //
  
  if ( (action & CREATE_CATALOGS ) ||
       ( action & WRITE_MASK_FILES ) ||
       ( action & BUILD_FOF_TABLE ) )
                                                              /* ------------------------ *
                                                               *                          *
                                                               *    create catalogs       *
                                                               *                          *
							       * ------------------------ */
    {

     #define Fof   0
     #define SubH  1
     #define FofP  2
     #define SubHP 3
      num_t HowMany[4];
      
      // ------------------------------------
      // get the subfind data
      //

      printf("loading subfind data..\n");
      tstart = CPU_TIME;
      ret = get_subfind_data( working_dir_subf, subf_name, HowMany );
      telapsed = CPU_TIME - tstart;
      PRINT_TIMINGS( "getting subfind data", "s", telapsed );      
      
      if( ret!= 0 )
	{
	  printf("[err] some error occurred while reading sub-find data, "
		 "let's clean up and terminate\n");
	 #pragma omp parallel
	  {
	    free( P[0] );
	  }
	  exit(1);
	}

      memset(diagnostic, 0, sizeof(num_t)*2);
     #pragma omp parallel
      get_howmany_in_subhaloes ( &diagnostic[0] );

      DPRINT(0, 0, "***** [D] %llu in shaloes, %llu in fof %d\n", diagnostic[0], diagnostic[1], FOFDBG );


      
      if ( !( ( action & CREATE_CATALOGS ) ||
	      ( action & WRITE_MASK_FILES ) ||
	      ( action & BUILD_FOF_TABLE) ) )
	return 0;

      
      /* if( action & WRITE_MASK_FILES ) */
      /* 	{ */
      /* 	  double tmask = CPU_TIME; */

      /* 	  printf("Writing ids for fof 0\n" ); */

      /* 	  write_fofgal_ids( working_dir_snap, snap_name, 0 ); */

      /* 	  tmask = CPU_TIME - tmask; */
      /* 	  PRINT_TIMINGS( "masking", "s", tmask ); */
      /* 	} */

      
      // note: at this point the particles are evenly distributed among the
      //       threads in the same order in which they appear in the subfind
      //       files, i.e. clustered by fof and sub-haloes

      
      // ------------------------------------  
      // re-distribute the data among threads
      // (particles are re-distributed by masked ids)

      if( Nthreads > 1 )
	printf("\tre-distributing subfind data among threads..\n");
      tstart = CPU_TIME;
      ret = distribute_particles( );
      telapsed = CPU_TIME - tstart;
      PRINT_TIMINGS( "distributing subfind data", "s", telapsed );
      if( ret != 0 )
	{
	  printf("[err] some error occurred while redistributing subfind particles, "
		 "let's clean up and terminate\n");
	 #pragma omp parallel
	  {
	    free( P[0] );
	  }
	  exit(2);
	}
      
      
      // ------------------------------------  
      // load the IDs from the snapshot so to reconstruct the types
      //

      printf("loading IDs and type data..\n");
      tstart = CPU_TIME;

      ret = get_id_data( working_dir_snap, snap_name );
      telapsed = CPU_TIME - tstart;
      PRINT_TIMINGS( "loading ids and type data", "s", telapsed );
  
      if ( ret != 0 )
	{
	  printf("[err] some error occurred while getting ids, "
		 "let's clean up and terminate\n");
	 #pragma omp parallel
	  {
	    free( P[0] );
	    free( IDs );
	  }
	  exit(3);
	}


      // ------------------------------------  
      // re-distribute the ID block in the snapshot
      // among the threads
      //
      if ( Nthreads > 1 )
	printf("\tre-distributing IDs and type data among threads..\n");
      tstart = CPU_TIME;
      ret = distribute_ids( );
      telapsed = CPU_TIME - tstart;
      PRINT_TIMINGS( "distributing ids and type data", "s", telapsed );
      
      if ( ret != 0 )
	{
	  printf("[err] some error occurred while re-distributing ids, "
		 "let's clean up and terminate\n");
	 #pragma omp parallel
	  {
	    free( P[0] );
	    free( IDs );
	  }
	  exit(4);
	}
      
      tstart = CPU_TIME;
     #pragma omp parallel
      {
	// ------------------------------------  
	//  assign to each subfind particle its type,
	//  the number of snapshot file to which it
	//  belongs and its position in the ID array
	//  of that file
	//
    
       #pragma omp single
	printf("\tsorting ids..\n");
    
	// sort ids-type structure for each thread
       #pragma omp master
	tstart = CPU_TIME;
	
	sort_thread_idtype();
       #pragma omp master
	{
	  telapsed = CPU_TIME - tstart;       
	  PRINT_TIMINGS( "sorting thread's ids", "s", telapsed );
	}
       #if defined(DEBUG)
	dprint(1, -1, "\tIDs limit for thread %d is: %llu\n",
	       me, IDranges[me] );
       #pragma omp barrier
       #pragma omp single
	dprint(0, -1, "\tchecking snap particles ids order..\n");
	
	{
	  num_t check = check_sorting( 0, myNID, 2);
	  if( check )
	    printf("[err] thread %d :: sorting of snap particles "
		   "by IDs is broken at %llu\n",
		   me, check);
	}
       #endif

       #pragma omp single
	printf("\tassigning type to subfind particles..\n");

	num_t OoR, rng_Fls, gen_Flrs;
	num_t fails = assign_type_to_subfind_particles( &OoR, &rng_Fls, &gen_Flrs );

	if( fails )
	  printf("\t[err] thread %d has got %llu out-of-range "
		 "%llu range failures and %llu generation failures over %llu particles\n",
		 me, OoR, rng_Fls, gen_Flrs, myNp);

       #pragma omp single
	memset(diagnostic, 0, sizeof(num_t)*2);
       #pragma omp barrier
	get_howmany_in_subhaloes ( &diagnostic[0] );
       #pragma omp barrier
	//#pragma omp single
	//printf("***** [D] %llu in shaloes, %llu in fof %d\n", diagnostic[0], diagnostic[1], FOFDBG );

	free(IDs);
       #pragma omp single
	{
	  free(IDranges);
	  free(all_NID);
	  free(all_IDs);
	}

      }

      
      // ------------------------------------  
      //  create a catalog table if required;
      //  that is a table which summarize how
      //  many particles of each type belong
      //  to every fof and galaxy
      //
      
      if( action & BUILD_FOF_TABLE )
	{	    
	  printf("building fof table..\n");
	  fof_table = (fof_table_t*)calloc(HowMany[Fof] , sizeof(fof_table_t));

	 #pragma omp parallel
	  make_fof_table( fof_table, HowMany[Fof], 1 );	  	 

	  FILE *file = fopen( catalog_table_name, "w" );
	  
	  for( num_t i = 0; i < HowMany[Fof]; i++ )
	    {
	      fprintf(file, "%10llu\t%lld\t%llu\t", i, (long long int)fof_table[i].nsubhaloes, fof_table[i].TotN );
	      
	      num_t subh_occupancy = 0;
	      for( int s = 0; s < fof_table[i].nsubhaloes; s++ ) {
		subh_occupancy += fof_table[i].subh_occupancy[s];
		if (fof_table[i].subh_occupancy[s] == 0)
		  printf(">>>> fof %llu has got sub-halo %u with 0 occupancy\n", i,s ); }
	      free( fof_table[i].subh_occupancy );
	      
	      fprintf(file, "%llu\t*\t", subh_occupancy);
	      
	      for( int j = 0; j < NTYPES; j++ )
		fprintf(file, "\t%10llu", fof_table[i].Nparts[j] );
	      fprintf(file, "\n");
	    }	    
	  
	  fclose(file);
	  free( fof_table );
	  
	}

      
       
      // ------------------------------------  
      //  write files for masked particles
      //
            
      if( action & WRITE_MASK_FILES )
	{
	  double tmask = CPU_TIME;

	  printf("Writing files for %d masked fofs\n", mask_fofgal_n );

	  write_fofgal_snapshots( working_dir_snap, snap_name, mask_fofgal, mask_fofgal_n );

	  tmask = CPU_TIME - tmask;
	  PRINT_TIMINGS( "masking", "s", tmask );
	}
      

      if( action & CREATE_CATALOGS )
	{
	  
	 #pragma omp parallel
	  {
    
	    /* if ( action & CREATE_CATALOGS ) */
	    /*   { */
	    
	    // ------------------------------------  
	    //  partition the data in each thread by particle type
	    //
	
	   #pragma omp single
	    printf("partitioning particles by type in each thread..\n");
        
	    {
	      memset( type_positions, 0, sizeof(num_t)*NTYPES );
	      num_t start = 0;
	      num_t stop = myNp;

	      telapsed = 0;
	      for( int t = 0; t < NTYPES; t++ ) {
		type_positions[t] = partition_P_by_type( start, stop, t );
	       #if defined(DEBUG)
		{
		  num_t check = check_partition(start, stop, type_positions[t], (PID_t)t, 1);
		  if ( check )
		    printf("[err] %d partitioning by type is broken at %llu\n",
			   me, check);

		 #if defined(MASKED_ID_DBG)
		  for( num_t i = 0; i < myNp; i++ )
		    if( P[0][i].pid == MASKED_ID_DBG )
		      dprint(0, me, "[ID DBG][P][th %d] found particle with masked id %llu, "
			     "type %d, gen %d at pos %llu\n", me,
			     (ull_t)P[0][i].pid, P[0][i].type, P[0][i].gen, i);	      
		 #endif
	      
		}
	       #endif	
		start = type_positions[t]; }
	  
	     #if defined(DEBUG)
	     #pragma omp single
	      dprint(0, -1, "\tcheck of type partitioning done..\n");
	     #endif

	    }
    
    
	    // ------------------------------------  
	    //  sort the data to create a catolog sorted by particle id
	    //

	   #pragma omp single
	    printf("sorting particles by id for each type in each thread..\n");

	    sort_thread_particles( type_positions );

	   #if defined(DEBUG)
	   #pragma omp single
	    dprint(0, -1, "\tcheck of subfind particles sorting by id..\n");

	    num_t ret = 0;
	    num_t start = 0;    
	    for( int t = 0; t < NTYPES; t++ ) {
	      num_t stop = type_positions[t];
	      ret += check_sorting( start, stop, 0);
	      start = stop; }	
	    if( ret )
	      printf("[err] thread %d :: sorting of particles by IDs is broken (%llu)\n",
		     me, ret);

	   #if defined(MASKED_ID_DBG)
	    for( num_t i = 0; i < myNp; i++ )
	      if( P[0][i].pid == MASKED_ID_DBG )
		dprint(0, me, "[ID DBG][S][th %d] found particle with masked id %llu, "
		       "type %d, gen %d at pos %llu\n", me,
		       (ull_t)P[0][i].pid, P[0][i].type, P[0][i].gen, i);	
	   #endif
	
	   #endif

	   #pragma omp master
	    {
	      Nparts_all[0] = (num_t*)malloc( (NTYPES+1)*Nthreads*sizeof(num_t));
	      for ( int i = 1; i <= NTYPES; i++ )
		Nparts_all[i] = Nparts_all[i-1] + Nthreads;
	    }
	   #pragma omp barrier
	
	    P_all[0][me] = P[0];
	    for( int t = 0; t < NTYPES; t++ )
	      {
		Nparts[t] = ( t>0 ? type_positions[t] - type_positions[t-1] : type_positions[t]);
		if( t > 0 ) {
		  P[t] = P[0] + type_positions[t-1];
		  P_all[t][me] = P[t];
		}
		Nparts_all[t][me] = Nparts[t];
	       #pragma omp atomic update
		Np_all[t] += Nparts[t];
	      }
	  }
	  telapsed = CPU_TIME - tstart;
	  PRINT_TIMINGS( "sorting", "s", telapsed );

      
	  // ------------------------------------  
	  // write the catalog files
	  //
	  
	  printf("creating catalogs..\n");
	  
	  FILE **catalog_files = malloc( NTYPES * sizeof(FILE*) );
	  
	  tstart = CPU_TIME;
	  
	  for( int t = 0; t < NTYPES; t++ )
	    {
	      char name[CATALOG_NAME_SIZE+10];
	      sprintf( name, "%s%d", catalog_name, t);
	      catalog_files[t] = fopen( name, "w" );
      
	      if( catalog_files[t] == NULL )
		{
		  for( int j = t-1; j>=0; j-- )
		    fclose(catalog_files[j]);
	  
		  printf("got a problem while creating catalog file %s, aborting\n",
			 name );
		  exit(11);
		}

	      printf("\tcatalog %s will have %llu particles\n",
		     name, Np_all[t]);

	      catalog_header_t catalog_header = {0};
	      catalog_header.id_size = (int)sizeof(PID_t);
	      catalog_header.nfiles = 1;
	      catalog_header.particle_t_size = (int)sizeof(particle_t);
	      catalog_header.Nparts = Np_all[t];                               // particles in this file
	      // writing of multiple files has yet to be
	      // implemented
	      catalog_header.Nparts_total = Np_all[t];                         // total number of particles

	      fwrite( &catalog_header, sizeof(catalog_header_t), 1, catalog_files[t] );        

	  
	     #pragma omp parallel for ordered
	      for( int th = 0; th < Nthreads; th++ )
	       #pragma omp ordered                  // force a pid-ordered writing
		write_particles( catalog_files[t], t );

	      fclose( catalog_files[t] );
	    }

	  free(catalog_files);

	  telapsed = CPU_TIME - tstart;
	  PRINT_TIMINGS( "writing catalogs", "s", telapsed );	  
	}

      
      
  
  
      if( !( action & SEARCH_PARTICLES) ){
	for( int t = 0; t < NTYPES; t++ )
	  free( P_all[t] );
       #pragma omp parallel
	{
	  free(P[0]);
	}}
    }

  

  if ( action & SEARCH_PARTICLES )                            /* ------------------------ *
                                                               *                          *
                                                               *    search in catalogs    *
                                                               *                          *
							       * ------------------------ */

    {

      int input_types[NTYPES] = {0};
      
      // explore which catalog we should load in
      //
      for( int ll = 0; ll < Nlists; ll++ )
	{
	  if( list_types[ll] >= 0 )
	    input_types[list_types[ll]] = 1;
	  else {
	    // one of the lists is either of type -1,
	    // i.e. "unknown", or of type -2, i.e. the
	    // file contains pairs (ID, type) and
	    // potentially all the types are present.
	    // then, let's set all the type to be present
	    // and to end the loop
	    for( int tt = 0; tt < NTYPES; tt++ )
	      input_types[tt] = 1;
	    break; }
	}
      
      if( !(action & CREATE_CATALOGS ) ) {
	// In case we just created catalogs, we assume to reuse
	//   P, Nparts and myNp
	// Otherwise, we read in the required catalog
	//
	dprint( 0, 0, "getting catalogs\n");
	tstart = CPU_TIME;
	int ret = get_catalog_data( catalog_name, &input_types[0] );
	telapsed = CPU_TIME - tstart;
	PRINT_TIMINGS( "getting catalogs data", "s", telapsed );
	if( ret ) {
	  printf( "a problem arose when reading catalog\n"
		 " I'm stopping here\n" );
	  exit(1); } }

      IDdomains[0] = (num_t*)malloc( NTYPES * Nthreads * sizeof(num_t) );
      for ( int i = 1; i < NTYPES; i++ )
	IDdomains[i] = IDdomains[i-1] + Nthreads;
      
     #pragma omp parallel
      {
	for( int t = 0; t < NTYPES; t++ )
	  if( Nparts[t] > 0 ) {
	    IDdomains[t][me] = P[t][Nparts[t]-1].pid;}
	  else
	    IDdomains[t][me] = 0;
      }
      
      // get the ids from the list
      //
      _Alignas(32) double search_timings[Nthreads];
      memset( search_timings, 0, Nthreads*sizeof(double));

      for ( int ll = 0; ll < Nlists; ll++ )
	{
	  dprint( 0, 0, "getting ids from list file %s, list type %d..\n", list_names[ll], list_types[ll]);
	  tstart = CPU_TIME;
	  get_list_ids( list_names[ll], list_types[ll] );
	  telapsed = CPU_TIME - tstart;
	  PRINT_TIMINGS( "getting ids from list", "s", telapsed );
	  dprint( 0, 0, "\t%llu ids found in the list file\n", Nl);
	  
	  int type_start = 0;
	  int type_end = NTYPES;
	  if( list_types[ll] >= 0 ) {
	    type_start = list_types[ll]; type_end = list_types[ll]+1;}

	  num_t  out_of_range = 0;
	  num_t  id_fails     = 0;
	  num_t  g_fails      = 0;
	  num_t  found        = 0;
	  double search_avg   = 0;

	  char  *name_out = (char*)malloc( strlen(list_names[ll])+5 );
	  snprintf( name_out, strlen(catalog_name)+5 , "%s.fh", list_names[ll] );
	  //FILE *list_out = fopen( name_out, "w" );
	  FILE *list_out = fopen( "list.out", "w" );
	  if( list_out == NULL )
	    {
	      printf("[search] it was impossible to create output file %s for list %s, I skip it\n",
		     name_out, list_names[ll]);
	      #pragma omp parallel
	      free(List);
	      continue;
	    }
	  free(name_out);
	  
	  dprint(0, 0, "searching for ids..\n" ); fflush(stdout);
	  
	  tstart = CPU_TIME;
	 #pragma omp parallel reduction(+:out_of_range, id_fails, g_fails, found, search_avg)
	  {

	    num_t my_out_of_range = 0;
	    num_t my_id_fails     = 0;
	    num_t my_g_fails      = 0;
	    num_t my_found        = 0;
	    int   bs              = sizeof(PID_t)*8 - id_bitshift;
	    
	    double mytstart = CPU_TIME;
	    for( num_t j = 0; j < myNl; j++ )
	      {

		// find masked ids and generation
		PID_t masked_id;
		int   generation;

		masked_id  = List[j].pid & id_mask;
		generation = List[j].pid >> bs;
		// in case the particles have different types and the
		// type is known, let's set-up the type boundaries
		// individually for each particle
		//
		type_start = ( list_types[ll] ==-2 ? Types[j] : type_start );
		type_end = ( list_types[ll] ==-2 ? Types[j]+1 : type_end );

		// search among all the catalogs
		// note: if list is typed, only the interested type is
		//       active in the search
		//
		for( int t = type_start; t < type_end; t++ )
		  {
		    // find where the particle is located among the threads
		    //
		    int target_thread = 0;
		    while( (target_thread < Nthreads) && (masked_id > IDdomains[t][target_thread]) )
		      target_thread++;
		    
		    if( target_thread < Nthreads )
		      // we found it, let's search
		      {
			particle_t *res;
		       #if !defined(USE_LIBC_BSEARCH)
			res = (particle_t*)mybsearch_in_P(P_all[t][target_thread], Nparts_all[t][target_thread], masked_id );
			
		       #else
			res = bsearch( (void*)&masked_id, P_all[t][target_thread], Nparts_all[t][target_thread],
				       sizeof(particle_t), compare_pid_with_particle_t );
		       #endif

			
			if( res != NULL )
			  {
			    // particle's been found

			    particle_t *stop = P_all[t][target_thread];
			    while( (res>stop) && ((res-1)-> pid == masked_id) && (res-> gen > generation) ) --res;
			    stop += Nparts_all[t][target_thread];
			    while( res < stop && ((res+1)-> pid == masked_id) && (res-> gen < generation) ) ++res;

			    if(res-> gen != generation) {
			     #if defined(DEBUG) && defined(MASKED_ID_DBG)
			      if( masked_id == MASKED_ID_DBG)
				printf("### %d %d\n", res->type, res->gen);
			     #endif

			      my_g_fails++; List[j].fofid = -1; List[j].gid = -1; }
			    else {
			     #if defined(DEBUG) && defined(MASKED_ID_DBG)
			      if( masked_id == MASKED_ID_DBG)
				printf("### %d %d\n", res->type, res->gen);
			     #endif

			      my_found++;
			      List[j].fofid = res->fofid;
			      List[j].gid   = res->gid; }
			  }
			else {
			  		       #if defined(DEBUG) && defined(MASKED_ID_DBG)
			if( masked_id == MASKED_ID_DBG)
			  printf("### %d %d\n", res->type, res->gen);
		       #endif

			  my_id_fails++; List[j].fofid = -1; List[j].gid = -1; }
		      }
		    else {
		      // oops, it seems that an invalid id
		      // has been provided
		      my_out_of_range++;
		      List[j].fofid = -1; List[j].gid = -1;
		    }
		  }
	      }
	    double mytelapsed = CPU_TIME - mytstart;
	   #pragma omp atomic update
	    search_timings[me] += mytelapsed;   // a little bit of false-sharing,
						// but that has a negligible impact on performance.
	    
	   #pragma omp single
	    dprint(1, me, "waiting for all the threads..\n");
	    
	   #pragma omp single
	    {
	      dprint(0, me, "writing output file..\n");
	      dprint(0, me, "writing %llu particles, with PID size %lu, item size %lu\n",
		     (unsigned long long)Nl, sizeof(PID_t), sizeof(list_t));
	      long long int N = Nl;
	      fwrite( &N, sizeof(long long), 1, list_out );
	      int PID_size = sizeof( PID_t );
	      fwrite( &PID_size, sizeof(int), 1, list_out );
	    }
	   #pragma omp for ordered
	    for( int th = 0; th < Nthreads; th++ )
	     #pragma omp ordered
		fwrite( List, sizeof(list_t), myNl, list_out);

	   #pragma omp barrier
	   #pragma omp master
	    fclose(list_out);
	    
	    if ( List != NULL ) { free( List ); List = NULL; }
	    if ( Types != NULL ) { free( Types ); Types = NULL; }

	   #pragma omp atomic update
	    search_avg += mytelapsed;
	   #pragma omp atomic update
	    found += my_found;
	   #pragma omp atomic update
	    out_of_range += my_out_of_range;
	   #pragma omp atomic update
	    g_fails += my_g_fails;
	   #pragma omp atomic update
	    id_fails += my_id_fails;


	   #pragma omp barrier                  // wait for everybody to finish
						// its own search
	    if( ll == Nlists - 1 )
	      free( Pbase );                    // no more lists to process
						// let's free the particles array

	  } // close parallel region

	  telapsed = CPU_TIME - tstart;	  
	  PRINT_TIMINGS( "searching for lists & writing files", "s", telapsed );
	  search_avg /= (Nthreads * Nlists);
	  PRINT_TIMINGS( "searching for lists (avg)", "s", search_avg );
	  {
	    double mint = DBL_MAX;
	    double maxt = 0;
	    for ( int t = 0; t < Nthreads; t++ ) {
	      mint = (mint > search_timings[t] ? search_timings[t] : mint );
	      maxt = (maxt < search_timings[t] ? search_timings[t] : maxt ); }
	    PRINT_TIMINGS( "searching for lists (imb)", "%", (maxt-mint)/search_avg*100 );
	    PRINT_TIMINGS( "searching for lists (min)", "s", mint );
	    PRINT_TIMINGS( "searching for lists (max)", "s", maxt );
	  }

	  printf("%llu / %llu particles found\n",
		 found, Nl);
	  
	  if ( out_of_range || g_fails || id_fails )
	    dprint(1, 0, "\t%llu ids were out-of-range\n"
		   "\t%llu ids have not been found in catalogs\n"
		   "\t%llu ids have not been found with exact generation\n",
		   out_of_range, id_fails, g_fails );
  
	}

      free( IDdomains[0] );
      free( Nparts_all[0] );
      for( int t = 0; t < NTYPES; t++ )
	free(P_all[t]);
    }



  telapsed = CPU_RTIME - tbegin;
  PRINT_TIMINGS( "total time", "s", telapsed );

  fclose( details );
  fclose( timings );
      
  // un-process arguments 
  //
  {
    int point = 1;
    while ( point < argc )
      {
	int a = 0;
	while( (a < n_args) && (strcasecmp(*(argv+point), args[a].argument) != 0) )
	  a++;
	if( a < n_args )
	  point = args[a].action( argv, point, -1 );
	point++;
      }
  }
  
 #if defined(USE_MPI)
  MPI_Finalize();
 #endif
  return 0;
}



