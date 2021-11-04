
#include "cidx.h"


#define CPU_TIME (clock_gettime( CLOCK_THREAD_CPUTIME_ID, &ts ), \
		  (double)ts.tv_sec +				  \
		  (double)ts.tv_nsec * 1e-9)


#define UNDEFINED        0
#define CREATE_CATALOGS  1
#define SEARCH_PARTICLES 2

int main( int argc, char **argv)
{
  struct timespec ts;
  double tstart, telapsed;
  double tbegin;
  FILE *timings = fopen("timings", "w");
  details = fopen("details", "w");
  int input_types[NTYPES];
  
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


  // ------------------------------------
  // get the command-line arguments
  //

  // prepare some defulat values of variables
  //
  verbose = 0;
  int ret;
  int action = UNDEFINED;
  
  sprintf( snapnum, "---");
  sprintf( working_dir, "%s", working_dir_def );
  sprintf( working_dir_snap_base, "%s", working_dir_snap_base_def );
  sprintf( working_dir_subf_base, "%s", working_dir_subf_base_def );
  sprintf( snap_base, "%s", snap_base_def );
  sprintf( subf_base, "%s", subf_base_def );


  // process the command line
  //
  n_stars_generations = n_stars_generations_def;
  {
    int point = 1;
    while ( point < argc )
      {
	if ( (strcasecmp(*(argv+point), "-c") == 0 ) ||
	     (strcasecmp(*(argv+point), "-cc") == 0 ) ||
	     (strcasecmp(*(argv+point), "-create") == 0 ) ||
	     (strcasecmp(*(argv+point), "-catalogs") == 0 ) )
	  action = CREATE_CATALOGS;
	else if ( (strcasecmp(*(argv+point), "-s") == 0 ) ||
		  (strcasecmp(*(argv+point), "-sp") == 0 ) ||
		  (strcasecmp(*(argv+point), "-search") == 0 ) )
	  action = SEARCH_PARTICLES;
	else if( strcasecmp(*(argv+point), "-snapf") == 0 )
	  snprintf( snap_base, NAME_SIZE, "%s", *(argv + ++point) );
	else if ( strcasecmp(*(argv+point), "-g") == 0 )
	  n_stars_generations = atoi( *(argv + ++point) );
	else if ( strcasecmp(*(argv+point), "-subf") == 0 )
	  snprintf( subf_base, NAME_SIZE, "%s", *(argv + ++point) );
	else if ( strcasecmp(*(argv+point), "-wdir") == 0 )
	  snprintf( working_dir, DIR_SIZE, "%s", *(argv + ++point) );
	else if ( strcasecmp(*(argv+point), "-snapdir") == 0 )
	  snprintf( working_dir_snap_base, DIR_SIZE, "%s", *(argv + ++point) );	    
	else if ( strcasecmp(*(argv+point), "-subfdir") == 0 )
	  snprintf( working_dir_subf_base, DIR_SIZE, "%s", *(argv + ++point) );
	else if ( strcasecmp(*(argv+point), "-num") == 0 )
	  snprintf( snapnum, NUM_SIZE, "%s", *(argv + ++point) );
	else if ( strcasecmp(*(argv+point), "-list") == 0 )
	  snprintf( list_name, CATALOG_NAME_SIZE, "%s", *(argv + ++point) );	
	else if ( strcasecmp(*(argv+point), "-v") == 0 )
	  verbose = atoi( *(argv + ++point) );
	else if ( strcasecmp(*(argv+point), "-ask") == 0 )
	  {
	    char qsort_v[12];
	    char bsearch_v[12];
	    char debug_checks[12];
	   #if !defined(USE_LIBC_QSORT)
	    sprintf(qsort_v, "custom");
	   #else
	    sprintf(qsort_v, "libc");
	   #endif
	   #if !defined(USE_LIBC_BSEARCH)
	    sprintf(bsearch_v, "custom");
	   #else
	    sprintf(bsearch_v, "libc");
	   #endif
	   #if defined(DEBUG)
	    sprintf(debug_checks, "active");
	   #else
	    sprintf(debug_checks, "inactive");
	   #endif

	    printf("here they go some infos on me:\n"
		   "  ID size is %lu\n"
		   "  I'm using %s qsort\n"
		   "  I'm using %s bsearch\n"
		   "  DEBUG checks are %s\n"
		   "\n",
		   sizeof(PID_t), qsort_v, bsearch_v, debug_checks );
	  }
	else if ( ( strcasecmp(*(argv+point), "-h") == 0 ) ||
		  ( strcasecmp(*(argv+point), "-help") == 0 ) )
	  {
	    printf( "Here it goes a little help from my side:\n\n"
		    " * how to specify what you want to do\n"
		    "   - read subfind files and create the inverse catalogs\n"
		    "     < -c | -cc | -create | -catalogs >\n"
		    "   - search a list of IDs from snapshot in inverse catalogs\n"
		    "     <  >\n"
		    "\n"
		    " * specify the snap/subfind number\n"
		    "   -num [num]  // insert also the trailing zeros\n"
		    "\n"
		    " * specify the basename of snapshot files\n"
		    "   -snapf [basename]  --> default is %s+snapnum\n"
		    "\n"
		    " * specify the basename of subfind files\n"
		    "   -subf [basename]  --> default is %s+snapnum\n"
		    "\n"
		    " * specify the subfind folder\n"
		    "   -subfdir [dirname]  --> default is %s%s+snapnum\n"
		    "\n"
		    " * specify the snapshot subdir\n"
		    "   -snapdir [dirname]  --> default is %s%s+snapnum\n"
		    "\n"
		    " * specify the working dir\n"
		    "   -wdir [dirname]  --> default is %s\n"
		    "\n"
		    " * specify the number of stellar generations\n"
		    "   -g [num]  --> default is %d\n"
		    "\n"
		    " * specify the verbose level\n"
		    "   -v [num]  --> default is %d\n"
		    "\n\n",
		    snap_base_def,
		    subf_base_def,
		    working_dir_def, working_dir_subf_base_def,
		    working_dir_def, working_dir_snap_base_def,
		    working_dir,
		    n_stars_generations_def,
		    verbose);
	  }

	else
	  printf("unknown argument %d: %s\n", point, *(argv+point));
	point++;
      }
  }

  // final processing of variables
  //
  
  if ( (action == UNDEFINED) || (strcmp(snapnum, "---")==0) ) {
    printf("Either you did not specify an action or you did not"
	   "specify a snapshot number.\n"
	   "I would not know what else to do more than saying this.\n\n" );
    exit(1); }
  
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

  DPRINT(3,0, "working dir is %s\n"
	 "snap dir is    %s\n"
	 "subf dir is    %s\n"
	 "snap name is   %s\n"
	 "subf name is   %s\n",
	 working_dir, working_dir_snap, working_dir_subf, snap_name, subf_name );

  id_mask = get_stellargenerations_mask( n_stars_generations, &id_bitshift );
  for( int t = 0; t < NTYPES; t++ )
    P_all[t] = (particle_t**)calloc( Nthreads, sizeof(particle_t*));
      

  tbegin = CPU_TIME;
  

  
  // if at command line it was specified to create catalogs,
  // then let's do it
  //
  
  if ( action | CREATE_CATALOGS )                             /* ------------------------ *
                                                               *                          *
                                                               *    create catalogs       *
                                                               *                          *
							       * ------------------------ */
    {

      // ------------------------------------
      // get the subfind data
      //

      printf("loading subfind data..\n");
      tstart = CPU_TIME;
      ret = get_subfind_data( working_dir_subf, subf_name );
      telapsed = CPU_TIME - tstart;
      fprintf( timings, "%35s%6.2g s\n", "getting subfind data", telapsed );      
      
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

      // ------------------------------------  
      // re-distribute the data among threads
      //

      if( Nthreads > 1 )
	printf("\tre-distributing subfind data among threads..\n");
      for( int t = 0; t < NTYPES; t++ )
	P_all[t] = (particle_t**)calloc( Nthreads, sizeof(particle_t*));
      tstart = CPU_TIME;
      ret = distribute_particles( );
      telapsed = CPU_TIME - tstart;
      fprintf( timings, "%35s%6.2g s\n", "distributing subfind data", telapsed );
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
      // load the IDs from the snapshot so that to reconstruct the types
      //

      printf("loading IDs and type data..\n");
      tstart = CPU_TIME;
      ret = get_id_data( working_dir_snap, snap_name );
      telapsed = CPU_TIME - tstart;
      fprintf( timings, "%35s%6.2g s\n", "loading ids and type data",
	       telapsed );
  
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
      // re-distribute the data among threads
      //
      if ( Nthreads > 1 )
	printf("\tre-distributing IDs and type data among threads..\n");
      tstart = CPU_TIME;
      ret = distribute_ids( );
      telapsed = CPU_TIME - tstart;
      fprintf( timings, "%35s%6.2g s\n", "distributing ids and type data",
	       telapsed );
      
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
	//  assigne a type to each subfind particle
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
	  fprintf( timings, "%35s%6.2g s\n", "sorting thread's ids",
		   telapsed );
	}
       #if defined(DEBUG)
       #pragma omp single
	dprint(0, -1, "\tchecking snap particles ids order..\n");
	
	{
	  ull_t check = check_sorting( 0, myNID, 2);
	  if( check )
	    printf("[err] thread %d :: sorting of snap particles "
		   "by IDs is broken at %llu\n",
		   me, check);
	}
       #endif

       #pragma omp single
	printf("\tassigning type to subfind particles..\n");

	ull_t OoR, Fls;
	ull_t fails = assign_type_to_subfind_particles( &OoR, &Fls );

	if( fails )
	  printf("[err] thread %d has got %llu out-of-range "
		 "and %llu failures over %llu particles\n",
		 me, OoR, Fls, myNp);

       #pragma omp barrier

	free(IDs);
	free(IDranges);
	free(all_NID);
	free(all_IDs);

	all_IDs[me] = NULL;

	// ------------------------------------  
	//  partition the data in each thread by particles type
	//

       #pragma omp single
	printf("partitioning particles by type in each thread..\n");
        
	{
	  memset( type_positions, 0, sizeof(ull_t)*NTYPES );
	  ull_t start = 0;
	  ull_t stop = myNp;

	  telapsed = 0;
	  for( int t = 0; t < NTYPES; t++ ) {
	    type_positions[t] = partition_P_by_type( start, stop, t );
	   #if defined(DEBUG)
	    {
	      ull_t check = check_partition(start, stop, type_positions[t], (PID_t)t, 1);
	      if ( check )
		printf("[err] %d partitioning by type is broken at %llu\n",
		       me, check);
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

	ull_t ret = 0;
	ull_t start = 0;    
	for( int t = 0; t < NTYPES; t++ ) {
	  ull_t stop = type_positions[t];
	  ret += check_sorting( start, stop, 0);
	  start = stop; }	
	if( ret )
	  printf("[err] thread %d :: sorting of particles by IDs is broken (%llu)\n",
		 me, ret);
       #endif

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
      fprintf( timings, "%35s%6.2g s\n", "sorting", telapsed );

      
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

	  catalog_header_t catalog_header;
	  catalog_header.id_size = (int)sizeof(PID_t);
	  catalog_header.nfiles = 1;
	  catalog_header.particle_t_size = (int)sizeof(particle_t);
	  catalog_header.Nparts = Np_all[t];                               // particles in this file
									   // writing of multiple files has yet to be
									   // implemented
	  catalog_header.Nparts_total = Np_all[t];                         // total number of particles

	  fwrite( &catalog_header, sizeof(catalog_header_t), 1, catalog_files[t] );        

	}

      for ( int type = 0; type < NTYPES; type ++ )
	{
	  FILE *file = catalog_files[type];
      
	 #pragma omp parallel for ordered
	  for( int th = 0; th < Nthreads; th++ )
	   #pragma omp ordered                  // force a pid-ordered writing
	    write_particles( file, type );

	  fclose( catalog_files[type] );
	}

      telapsed = CPU_TIME - tstart;
      fprintf( timings, "%35s%6.2g s\n", "writing catalogs",
	       telapsed );
      
      free(catalog_files);
  
  
      if( !( action | SEARCH_PARTICLES) ){
	for( int t = 0; t < NTYPES; t++ )
	  free( P_all[t] );
       #pragma omp parallel
	{
	  free(P[0]);
	}}
    }

  

  if ( action | SEARCH_PARTICLES )                            /* ------------------------ *
                                                               *                          *
                                                               *    search in catalogs    *
                                                               *                          *
							       * ------------------------ */

    {
      if( !(action | CREATE_CATALOGS ) )
	// In case we just created catalogs, we assume to reuse
	//   P, Nparts and myNp
	// Otherwise, we read in the required catalog
	//
	get_catalog_data( catalog_name, input_types );

      IDdomains[0] = (ull_t*)malloc( NTYPES * Nthreads * sizeof(ull_t) );
      for ( int i = 1; i < NTYPES; i++ )
	IDdomains[i] = IDdomains[i-1] + Nthreads;
      
     #pragma omp parallel
      {
	for( int t = 0; t < NTYPES; t++ )
	  if( Nparts[t] > 0 )
	    IDdomains[t][me] = P[t][Nparts[t]-1].pid;
	  else
	    IDdomains[t][me] = 0;
      }
      
      // get the ids from the list
      //
      int typed = get_list_ids( list_name );

      ull_t out_of_range = 0;
      ull_t fails = 0;
      
      for( ull_t j = 0; j < myNl; j++ )
	{
	  int type_start, type_end;
	  if ( typed )
	    type_start = type_end = List[j].type;
	  else
	    type_start = 0, type_end = NTYPES;
	    
	  for( int t = type_start; t < type_end; t++ )
	    {
	      int target_thread = 0;
	      while( (target_thread < Nthreads) && (List[j].pid > IDdomains[t][target_thread]) )
		target_thread++;
	  
	      if( target_thread < Nthreads )
		{
		  particle_t *res;
		 #if !defined(USE_LIBC_BSEARCH)
		  int err;
		  ull_t pos = mybsearch_in_P(P_all[t][target_thread], Nparts_all[t][target_thread], List[j].pid, &err);
		  res = ( err == 0? &P_all[t][target_thread][pos] : NULL );
		  
		 #else
		  res = bsearch( (void*)&(List[j].pid), P_all[t][target_thread], Nparts_all[t][target_thread],
				 sizeof(particle_t), compare_pid_with_particle_t );
		 #endif
		  
		  if( res != NULL ) {
		    if( ! typed)
		      List[j].type = res->type;
		    List[j].gen    = res->gen;
		    List[j].fofid  = res->fofid;
		    List[j].gid    = res->gid; }
		  else {
		    fails++;
		    List[j].type = -1; }
		}
	      else
		out_of_range++;
	    }
	}

      
    }



  telapsed = CPU_TIME - tbegin;
  fprintf( timings, "%35s%6.2g\n", "total time",
	   telapsed );

  fclose( details );
  fclose( timings );
  
 #if defined(USE_MPI)
  MPI_Finalize();
 #endif
  return 0;
}



