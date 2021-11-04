#include "cidx.h"


/*
 * Table of Contents
 * 
 * - routines for id masking
 *
 * - routines for subfind data input
 *
 * - routines for id and type input from snapshot
 *
 * - routines for catalog output
 *
 * - routines for catalog input
 *
 */


/*  -------------------------------------------------------------
 * 
 *   routines for id masking
 *
 *  ------------------------------------------------------------- */



PID_t get_stellargenerations_mask( int n_gen, int *bitshift )
/*
 * This function builds the bit mask necessary to recover
 * the original ID of a star particle when n_gen generations
 * of stars are spawned from a single gas particle.
 * The mask is returned with type PID_t, while the number of
 * bits it needs is returned in *bitshift
 */
{
 #if !defined(LONG_IDS)
 #define BASE 32
 #else
 #define BASE 64
 #endif
  
  PID_t mask = 0;
  int   log2g = 0;
  for( ; (1 << log2g) < n_gen; mask += (1<<log2g++));
  *bitshift = log2g;
  return ~(mask << (BASE-log2g));
}



int get_stellargenerations( PID_t id, PID_t mask, int bitshift )
/*
 * This function returns the stellar generation to which
 * a given particle with id id belongs
 */
{
  return (int)((id & mask) >> bitshift);  
}




/*  -------------------------------------------------------------
 * 
 *   routines for subfind data input
 *
 *  ------------------------------------------------------------- */


typedef struct { int tag1; char name[4]; int len; int tag2; } head_t;
typedef struct {
  int npart[6];			/*!< number of particles of each type in this file */
  double mass[6];		
  double time;
  double redshift;
  int flag_sfr;
  int flag_feedback;
  unsigned int npartTotal[6];	/*!< total number of particles of each type in this snapshot. This can be
				  different from npart if one is dealing with a multi-file snapshot. */
  int flag_cooling;
  int num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  int flag_stellarage;
  int flag_metals;
  unsigned int npartTotalHighWord[6];
  int flag_entropy_instead_u;
  int flag_doubleprecision;
  
  int flag_ic_info;
  float lpt_scalingfactor;
  
  char fill[18];
  char names[15][2];
} snapheader_t;


int seek_block ( FILE *file, char name[5] )
/*
 * This function get to the block specified by name
 * in a format-2 snapshot file pointed by *file
 *
 * RETURN VALUE:
 * 0 if the seeking is not successful, otherwiase
 * the 4-bytes long tag at the begin of the data block
 */
{
  head_t        head;
  size_t        ret;
  
  rewind( file );

  ret = fread( &head, sizeof(head_t), 1, file );
  while ( (!feof(file)) && (strncmp( head.name, name, 4 ) != 0) && (ret == 1))
    {
      fseeko( file, (off_t)head.len, SEEK_CUR );
      if ( ! feof(file) )
	ret = fread( &head, sizeof(head_t), 1, file );
    }

  if ( !feof(file) && (ret == 1))
    {
      int tag;
      ret = fread( &tag, sizeof(int), 1, file);
      return tag;
    }
  else  
    return 0;
}



int get_subfind_data(char *working_dir, char *subf_base)
/*
 * this function gets the subfind data, i.e. the
 * arrays of pairs (galaxy, starid )
 *
 * temporarly, for development purposes, we
 * randomly generate them
 */
{

  FILE * filein;
  size_t ret = 0;
  unsigned int multifile = 0, nfiles = 1;
  char   fname[ strlen(subf_base)+strlen(working_dir)+5 ];

  // as first, let's find out how many files there are
  // and how many particles do they have in substructures
  // so that we can distribute the file reading among threads
  //

  
  // check whether we are multi-file
  //
  sprintf( fname, "%s/%s", working_dir, subf_base );  
  if( (filein = fopen( fname, "r" )) == NULL )
    {
      multifile = 1;
      sprintf( fname, "%s/%s.0", working_dir, subf_base );
      if( (filein = fopen( fname, "r" )) == NULL )
	{
	  // unable to find the specified subfind files
	  printf("unable to find the subfind file both as %s and as %s.0\n",
		 subf_base, subf_base);
	  return -2;
	}
      
      do
	{
	  fclose(filein);
	  sprintf( fname, "%s/%s.%d", working_dir, subf_base, nfiles );
	  nfiles += ((filein = fopen( fname, "r" )) != NULL);
	}
      while( filein != NULL) ;  
    }


  // now reconstruct how many particles there are in total in
  // the haloes and in the subhaloes and assign equal bunches
  // of them to each thread
  //

  #define Fof   0
  #define SubH  1
  #define FofP  2
  #define SubHP 3
  ull_t HowMany[4] = {0};
  {
    if ( multifile )
      sprintf( fname, "%s/%s.0", working_dir, subf_base );
    else
      sprintf( fname, "%s/%s", working_dir, subf_base );
    filein = fopen( fname, "r" );
    ret = seek_block( filein, "HEAD");
    if ( ret == 0)
      {
	printf("unable to find the HEAD block in subfind files\n");
	return -3;
      }
    // skip some bytes in the header and get to the total number
    // of particles (which actually is:
    //   type 0: Fof haloes
    //   type 1: Sub Haloes
    //   type 2: how many particles in Fof haloes
    fseeko(filein, (off_t)(sizeof(int)*8+sizeof(double)*8), SEEK_CUR);

    // load the low part of particles' number
    unsigned int low[3];
    unsigned int high[3];
    ret = fread( &low[0], sizeof(int), 3, filein );

    // skip to the number of file the snapshot is splitted into
    fseeko(filein, (off_t)(sizeof(int)*4), SEEK_CUR);
    unsigned int numfiles;
    ret = fread( &numfiles, sizeof(int), 1, filein);

    if( numfiles != nfiles )
      {
	fclose(filein);
	printf("the number of file specified in the header (%d) "
	       "is different than the number of files found (%d)\n",
	       numfiles, nfiles );
	return -3;
      }

    // load the high part of particles' number
    fseeko(filein, (off_t)(sizeof(double)*4+sizeof(int)*2), SEEK_CUR);
    ret = fread( &high[0], sizeof(int), 3, filein );

    // reconstruct the total number of particles
    for( int i = 0; i < 3; i++ )
      HowMany[i] = (high[i]<<31) + low[i];
    
    fclose(filein);
  }  
  
  typedef struct { FILE *ptr; ull_t npart, npart_inc; omp_lock_t lock;} record_file_t;
  typedef struct { ull_t npart, offset; int parent, index; } record_halo_t;
  typedef struct { ull_t npart, offset, nsubh; } record_fof_t;
  
  record_file_t *files  = (record_file_t*)calloc( nfiles, sizeof(record_file_t) );
  record_halo_t *haloes = (record_halo_t*)calloc( HowMany[SubH]+1, sizeof(record_halo_t) );
  record_fof_t  *fofs   = (record_fof_t*)calloc( HowMany[Fof]+1, sizeof(record_fof_t) );

  // now we need to read all the haloes to discover how many particles in total live
  // therein, because that useful number is not present in the header
  
  PID_t Nhaloes            = 0;
  PID_t Nfofs              = 0;
  unsigned int hidx        = 0;  
  unsigned int prev_parent = 0;
  unsigned int go          = 0;
  do
    {
      #define input files[go].ptr

      unsigned int nfof;
      unsigned int nsubh;
      unsigned int nfofP;
      
      if ( multifile )
	sprintf( fname, "%s/%s.%d", working_dir, subf_base, go );
      else
	sprintf( fname, "%s/%s", working_dir, subf_base );
      
      DPRINT(2,-1, "opening file %d: %s\n", go, fname);
      
      files[go].ptr  = fopen( fname, "r" );
      omp_init_lock(&files[go].lock);
		
      ret = seek_block( input, "HEAD");
      ret = fread( &nfof   , sizeof(int), 1, input );
      ret = fread( &nsubh  , sizeof(int), 1, input );
      ret = fread( &nfofP , sizeof(int), 1, input );
      
      if( nfof > 0 )
	// if no haloes or subhaloes are present,
	// we don't have anything to do here
	//
	{	  
	  unsigned int *buffer = (unsigned int*)calloc( (nsubh>nfof?nsubh:nfof), sizeof(int) );

	  seek_block( input, "GLEN" );
	  ret = fread( buffer, sizeof(int), nfof, input );
	  for( unsigned int i = 0; i < nfof; i++ )
	    fofs[ Nfofs + i ].npart       = buffer[i];
	  
	  seek_block( input, "GOFF" );
	  ret = fread( buffer, sizeof(int), nfof, input );
	  for( unsigned int i = 0; i < nfof; i++ )
	    fofs[ Nfofs + i ].offset = buffer[i];
	  
	  seek_block( input, "NSUB" );
	  ret = fread( buffer, sizeof(int), nfof, input );
	  for( unsigned int i = 0; i < nfof; i++ )
	    fofs[ Nfofs + i ].nsubh = buffer[i];
	  
	  seek_block( input, "SLEN" );
	  ret = fread( buffer, sizeof(int), nsubh, input );
	  for( unsigned int i = 0; i < nsubh; i++ ) {
	    HowMany[SubHP]                  += buffer[i];
	    haloes[ Nhaloes + i ].npart      = buffer[i]; }
	  
	  seek_block( input, "SOFF" );
	  ret = fread( buffer, sizeof(int), nsubh, input );
	  for( unsigned int i = 0; i < nsubh; i++ )
	    haloes[ Nhaloes + i ].offset = buffer[i];

	  seek_block( input, "GRNR" );
	  ret = fread( buffer, sizeof(int), nsubh, input );
	  for( unsigned int i = 0; i < nsubh; i++ )
	    {
	      haloes[ Nhaloes + i ].parent = buffer[i];
	      hidx = ( buffer[i] == prev_parent ? hidx+1 : 0);
	      prev_parent = buffer[i];
	      haloes[ Nhaloes + i ].index  = hidx;
	    }

	  free(buffer);
	  Nfofs   += nfof;
	  Nhaloes += nsubh;
	}

      files[go].npart     = nfofP;
      files[go].npart_inc = files[go].npart + (go > 0 ? files[go-1].npart_inc : 0 );
      
      go++;      
    }
  while( go < nfiles);

  fofs[ HowMany[Fof] ].npart  = 0;
  fofs[ HowMany[Fof] ].offset = HowMany[FofP];

  haloes[ HowMany[SubH] ].npart  = 0;
  haloes[ HowMany[SubH] ].offset = HowMany[FofP];
  
  // determines which is the number of
  // particles to be distributed among threads
  //
  Np = HowMany[FofP];
  
  printf("%llu total particles belong to %llu sub-haloes\n", HowMany[SubHP], HowMany[SubH]);

  int failures = 0;
  #pragma omp parallel
  {

    int rem  = (int)(Np % (ull_t)Nthreads);
    myNp     = Np / Nthreads + (me < rem);
    PPP      = (particle_t*)calloc( myNp, sizeof(particle_t));
    
    PID_t avg_Np  = Np / Nthreads;
    ll_t  fofn    = 0;
    ll_t  hn      = 0;
    ul_t  fn      = 0;
    ull_t myoff   = avg_Np * me;
    myoff  += ( me > rem ? rem : me );

    ull_t in_fof_off  = 0;
    ull_t halo_read   = 0;
    ull_t fof_read    = 0;
    ull_t read        = 0;

    int  _failures_ = 0;
    ull_t particles_off;
    
    while( (fofn < Nfofs) && (myoff < fofs[fofn+1].offset) )   // note: usinf fofn+1 is safe because
      fofn++;						   // fofs has Nfofs+1 elements
    // find the offset inside the fof
    in_fof_off = myoff - fofs[fofn].offset;
    fof_read   = in_fof_off;
    read       = in_fof_off;
    
    while( (hn < Nhaloes) && (myoff < haloes[hn+1].offset) )   // note: same here as for hn+1
      hn++;
    hn += ( myoff > haloes[hn].offset + haloes[hn].npart );    // in case myoff falls at the end of a fof
                                                               // where its haloes are finished and those
                                                               // of the next fof did not start yet.
    if( hn < Nhaloes )
      {
	ul_t in_halo_off = 0;
	// find the offset inside the halo
	if ( myoff > haloes[hn].offset )
	  in_halo_off = myoff - haloes[hn].offset;
	else
	  hn = -hn;
	halo_read = in_halo_off;
      }
    
    particles_off = myoff;
    
    while( (read < myNp) && (hn < Nhaloes) && !_failures_ )
      {

	particles_off = myoff + read;
	  
	while( (fn < nfiles) && ( particles_off >= files[fn].npart_inc) )
	  fn++;
	
	if ( fn < nfiles )
	  {

	    ull_t in_file_off = particles_off - (files[fn].npart_inc - files[fn].npart); // that is zero when this is not the first file we read
	    
	    ull_t howmany_toread = files[fn].npart - in_file_off;
	    
	    PID_t *buffer = (PID_t*)calloc( howmany_toread, sizeof(PID_t) );
	    
	    omp_set_lock( &files[fn].lock );
	    
	    /* DPRINT( 2, -1, "[ I/O ] [ th %d ] is opening file %u to read %llu particles starting from %llu from halo %u with halo_read %llu and in_halo_offset %llu: total read particles are %llu -> %llu / %llu\n",		     */
	    /* 	    me, fn, howmany_toread, in_file_off, hn, halo_read, in_halo_off, read, read+howmany_toread, (ull_t)myNp); */
	    /* DPRINT(2, -1, "[ I/O ] [ th %d ] ---- file %d --- %llu %llu %llu %llu\n", me, fn, */
	    /* 	   haloes[hn].offset, in_halo_off, files[fn].npart_inc, files[fn].npart ); */
	    
	    ret = seek_block( files[fn].ptr, "PID " );
	    if ( ret == 0 ) {                                                             // check that the block exists
	      printf("unable to find the PID block in subfind files\n");
	      free(buffer); free(P[0]); _failures_++; }
	    
	    else {
	      ret = (size_t)((ull_t)ret / files[fn].npart);
	      if ( ret != sizeof(PID_t) ) {                                               // check that the ID type lenght matches
		printf("thread %d: ID block in subfind file %d "
		       "has %lu-bytes long data while my ID type has %lu bytes\n",
		       me, fn, ret, sizeof(PID_t) ); ret = 0; _failures_++; }

	      if ( in_file_off > 0 )                                                      // skip initial data if needed
		{
		  off_t before = ftello( files[fn].ptr);                                  
		  fseeko( files[fn].ptr, (off_t)(in_file_off*sizeof(PID_t)), SEEK_CUR);   
		  off_t after = ftello( files[fn].ptr);
		  
		  if( (ull_t)(after-before) != (in_file_off*sizeof(PID_t) )) {            // check that we read the right amount of data
		    printf("%d has got problem, in file %d it is at offset "
			   "%llu instead of %llu due to %llu \n",
			   me, fn, (ull_t)after, (ull_t)before+(in_file_off*sizeof(PID_t)), in_file_off);
		    _failures_++; }
		}

	      ret = fread( buffer, sizeof(PID_t), howmany_toread, files[fn].ptr);         // actually read the data
	      if ( ret != howmany_toread ) {                                              // check that everything went fine
		printf("%d has read %llu instead of %llu from file %d\n",
		       me, (ull_t)read, (ull_t)howmany_toread, fn);
		_failures_++; }
	    }

	    if ( _failures_ )                                                             // if we've got some failure
	      {                                                                           // signal it
	       #pragma omp atomic update
		failures++;
		omp_unset_lock( &files[fn].lock);	
		continue;
	      }
	    
	    omp_unset_lock( &files[fn].lock);
	    
	    for( ull_t i = 0; (i < howmany_toread) && (read < myNp); i++, read++ )
	      {

		// fill particle's info
		//
		PPP[read].pid   = buffer[i];
		PPP[read].fofid = fofn;
		PPP[read].gid   = -1;

		// keep track with the proceeding of
		// this fof's reading
		//
		fof_read++;
		if( fof_read == fofs[fofn].npart )
		  {
		    fofn++;
		    if( fofn < Nfofs )
		      in_fof_off = 0;
		  }

		if ( hn < 0 )
		  // check whether we entered in a subhalo if
		  // we were outside of them
		  if( particles_off+i >= haloes[-hn].offset )
		    hn = -hn;
		  
		if ( hn >= 0 )
		  {
		    // we are in a subhalo,
		    // keep track of it
		    halo_read ++;
		    PPP[read].gid = haloes[hn].index;
		    
		    if( halo_read == haloes[hn].npart ) 
		      {
			hn++;
			halo_read = 0;
			if( haloes[hn].offset > particles_off + i )
			  // if we are in fof-onyl region
			  // disable the subhalo
			  hn = -hn;
		      }
		  }

	      }

	    free(buffer);	    
	    if( fn < nfiles-1 ) {
	      fn++;
	      in_file_off = 0; }
	  }
	else
	  {
	   #pragma omp atomic update
	    failures++;
	    
	    printf("[ E ] [ th %d ] we've got a problem (failures=%d): "
		   "the sub files seem not to be enough, got to inspect nr %u "
		   "which clearly is unconsistent. Let me abort and signal..\n", me, _failures_,fn);
	  }

       #pragma omp atomic read
	_failures_ = failures;

      }

    fprintf( details, "thread %d got %llu particles from subfind data\n", me, myNp);
  }

  fprintf( details, "\n");
  
  for( unsigned int i = 0; i < nfiles; i++ )
    fclose( files[i].ptr);
  free(haloes);  
  free(files);
  
  return failures;
}





int get_id_data(char *working_dir, char *name)
/*
 * this function gets the subfind data, i.e. the
 * arrays of pairs (galaxy, starid )
 *
 * temporarly, for development purposes, we
 * randomly generate them
 */
{

  FILE *filein;
  size_t ret = 0;
  unsigned int   multifile = 0, nfiles = 1;
  char  fname[ strlen(subf_base)+strlen(working_dir)+5 ];
  ull_t AllN = 0;

  // as first, let's find out how many files there are
  // and how many particles do they have in substructures
  // so that we can distribute the file reading among threads
  //

  
  // check whether we are multi-file
  //  
  sprintf( fname, "%s/%s", working_dir, name );  
  if( (filein = fopen( fname, "r" )) == NULL )
    {
      multifile = 1;
      sprintf( fname, "%s/%s.0", working_dir, name );
      if( (filein = fopen( fname, "r" )) == NULL )
	{
	  // unable to find the specified subfind files
	  printf("unable to find the snapshot file both as %s and as %s.0 under dir %s\n",
		 name, name, working_dir);
	  return -2;
	}
      
      do
	{
	  fclose(filein);
	  sprintf( fname, "%s/%s.%d", working_dir, name, nfiles );
	  nfiles += ((filein = fopen( fname, "r" )) != NULL);
	}
      while( filein != NULL) ;  
    }
  else
    fclose(filein);
  
  {
    if ( multifile )
      sprintf( fname, "%s/%s.0", working_dir, name );
    else
      sprintf( fname, "%s/%s", working_dir, name );
    filein = fopen( fname, "r" );
    ret = seek_block( filein, "HEAD");
    if ( ret == 0)
      {
	printf("unable to find the HEAD block in snapshot files\n");
	return -3;
      }
    fseeko(filein, (off_t)(sizeof(int)*8+sizeof(double)*8), SEEK_CUR);

    unsigned int low[NTYPES];
    unsigned int high[NTYPES];
    ret = fread( &low[0], sizeof(int), NTYPES, filein );

    fseeko(filein, (off_t)(sizeof(int)*1), SEEK_CUR);
    unsigned int numfiles;
    ret = fread( &numfiles, sizeof(int), 1, filein);

    if( numfiles != nfiles )
      {
	fclose(filein);
	printf("the number of file specified in the header (%d) of the snapshot file #0 is different than the number of files found (%d)\n",
	       numfiles, nfiles );
	return -3;
      }
    
    fseeko(filein, (off_t)(sizeof(double)*4+sizeof(int)*2), SEEK_CUR);
    ret = fread( &high[0], sizeof(int), NTYPES, filein );

    for( int i = 0; i < NTYPES; i++ )
      AllN += (high[i]<<31) + low[i];
    
    fclose(filein);
  }  

 #pragma omp parallel
  {
    int rem  = (int)(AllN % (ull_t)Nthreads);
    myNID    = AllN / Nthreads + (me < rem);
    IDs      = (pidtype_t*)calloc( myNID, sizeof(pidtype_t));
  }
  
  typedef struct { FILE *ptr; ul_t npart[NTYPES]; ull_t npart_all, npart_all_inc; omp_lock_t lock;} record_file_t;  
  record_file_t *files  = (record_file_t*)calloc( nfiles, sizeof(record_file_t) );

  unsigned int go = 0;
  do
    {
     #define input files[go].ptr
      
      if ( multifile )
	sprintf( fname, "%s/%s.%d", working_dir, name, go );
      else
	sprintf( fname, "%s/%s", working_dir, name );
      
      files[go].ptr  = fopen( fname, "r" );
      omp_init_lock(&files[go].lock);
		
      seek_block( input, "HEAD");

      ret = fread( &files[go].npart[0], sizeof(int), NTYPES, input );

      files[go].npart_all_inc = (go > 0 ? files[go-1].npart_all_inc : 0);
      files[go].npart_all = files[go].npart[0];      
      for( int i = 1; i < NTYPES; i++ )
	files[go].npart_all += files[go].npart[i];
      files[go].npart_all_inc += files[go].npart_all;
      
      go++;
     #undef input
    }
  while( go < nfiles);

  int failures = 0;
  
 #pragma omp parallel
  {
    int rem       = (int)(AllN % (ull_t)Nthreads);
    PID_t avg_NID = AllN / Nthreads;
    ull_t myoff   = avg_NID * me;
    myoff  += ( me > rem ? rem : me );

    ul_t  file_start = 0;
    ull_t file_start_offset = 0;

    while( (file_start < nfiles) && (myoff >= files[file_start].npart_all_inc) )
      file_start++;
    file_start_offset = myoff - (files[file_start].npart_all_inc - files[file_start].npart_all);

    if( file_start < nfiles )
      {
	unsigned int _failures_ = 0;
	unsigned int fn = file_start;
	PID_t read = 0;

	while( (read < myNID) && (fn < nfiles) && !_failures_ )
	  {
	    ull_t npart_local_inc[NTYPES] = {0};
	    ull_t howmany_toread = 0;
	    int   type = 0;
		
	    npart_local_inc[0] = files[fn].npart[0];
	    for( int i = 1; i < NTYPES; i++ )
	      npart_local_inc[i] = files[fn].npart[i] + npart_local_inc[i-1];
	
	    while( file_start_offset > npart_local_inc[type] )
	      type++;

	    howmany_toread = myNID - read;
	    howmany_toread = (howmany_toread > files[fn].npart_all - file_start_offset ?
			      files[fn].npart_all - file_start_offset : howmany_toread);
	    
	    PID_t *buffer = (PID_t*)calloc( howmany_toread, sizeof(PID_t) );

	    omp_set_lock( &files[fn].lock );
	
	    ret = seek_block(files[fn].ptr, "ID  ");
	    if ( ret )
	      {
		ret = (size_t)((ull_t)ret / files[fn].npart_all);
		if ( ret != sizeof(PID_t) ) {
		  printf("thread %d: ID block in snap file %d has %lu-bytes long data while my ID type has %lu bytes\n",
			 me, fn, ret, sizeof(PID_t) ); ret = 0; }
	      }
	    if ( ret )
	      {
		ret = fseeko(files[fn].ptr, (off_t)(sizeof(PID_t)*file_start_offset), SEEK_CUR);
		fread(buffer, sizeof(PID_t), howmany_toread, files[fn].ptr);

		for ( ull_t i = 0; i < howmany_toread; read++, i++ )
		  {
		    type += (file_start_offset+i > npart_local_inc[type]);
		    IDs[read].pid  = buffer[i];
		    IDs[read].type = type;
		  }
		
		omp_unset_lock( &files[fn].lock );
		free(buffer);

		fn += (read < myNID);
		file_start_offset = 0;

		if( fn == nfiles ) {
		  _failures_++;
		 #pragma omp atomic update
		  failures++; }
		else {
		 #pragma omp atomic read
		  _failures_ = failures; }
	      }
	    else {
	      _failures_++;
	      free(buffer);
	     #pragma omp atomic update
	      failures++; }	    
	  }
      }
    else
      {
       #pragma omp atomic update
	failures++;
      }

    fprintf( details, "thread %d got %llu particles from snapshot data\n", me, myNID);
  }

  fprintf( details, "\n");
  
  for( unsigned int n = 0; n < nfiles; n++ )
    fclose(files[n].ptr);
  
  free( files );

  return failures;
  
}


/*  -------------------------------------------------------------
 * 
 *   routines for data output
 *
 *  ------------------------------------------------------------- */


int write_particles( FILE *file, int type )
{
  ull_t start = (type > 0? type_positions[type-1] : 0);
  ull_t stop  = type_positions[type];

  if( stop > start )
    fwrite( (void*)&PPP[start], sizeof(particle_t), (size_t)(stop-start), file);
  
  return 0;
}



/*  -------------------------------------------------------------
 * 
 *   routines for catalog inpunt
 *
 *  ------------------------------------------------------------- */

int get_catalog_numparticles( char *, catalog_header_t * );
ull_t get_catalog_data_from_file(char *, int, ull_t, ull_t, particle_t *);

int get_catalog_data(char *name, int *types)
{


  // find how many particles are in files for each type
  // of particles
  //

  catalog_header_t header;
  char fname[ CATALOG_NAME_SIZE + 5];

  Nparts_all[0] = (ull_t*)malloc( NTYPES*Nthreads*sizeof(ull_t));
  for ( int i = 1; i < NTYPES; i++ )
    Nparts_all[i] = Nparts_all[i-1] + Nthreads;

  Np = 0;
  for ( int t = 0; t < NTYPES; t++ )
    if( types[t] )
      {
	sprintf( fname, "%s%d", name, t);
	int nfiles = get_catalog_numparticles( fname, &header);
	if ( nfiles == 0 )
	  {
	    printf("unable to find the subfind file both as %s and as %s.0\n",
		   fname, fname);
	    return -1;
	  }
	if ( nfiles != header.nfiles )
	  {
	    printf("There a problem with catalogs for type %d: the expected number of files is %d but I've found %d files\n",
		   t, header.nfiles, nfiles);
	    return -2;
	  }
	Np_all[t] = header.Nparts_total;
	Np += Np_all[t];
       #pragma omp parallel
	{
	  int rem            = (int)(Np_all[t] % (ull_t)Nthreads);
	  Nparts[t]          = Np_all[t] / Nthreads + (me < rem);
	  Nparts_all[t][me]  = Nparts[t];
	  myNp              += Nparts[t];
	}
      }

  // allocate memory to hold all the particles
  // and set-up pointers for each particles type
  //
  int failures = 0;
 #pragma omp parallel
  {
    int t0 = 0;
    while ( Nparts[t0] == 0 ) t0++;
    if( t0 < NTYPES )
      {
	Pbase = (particle_t*)calloc( myNp, sizeof(particle_t));   // we use a single allocation
								  // to hold all the particles
	
	P[t0] = Pbase;                                            // then we set-up pointers to
	P_all[t0][me] = P[t0];                                    // each single type
	ull_t sum = Nparts[t0];
	for ( int t = t0+1 ; t < NTYPES; t++ ) {
	  P[t] = Pbase + sum;
	  P_all[t][me] = P[t]; }
      }
    else
     #pragma omp atomic update
      failures++;
  }

  if ( failures )
    {
      #pragma omp parallel
      {
	if( Pbase != NULL )
	  free(Pbase);
      }
      printf("I've got a problem in memory allocation\n");
      return -3;
    }

  for ( int t = 0; t < NTYPES; t++ )
    if( types[t] )
      {
	sprintf( fname, "%s%d", name, t);
	get_catalog_data_from_file( fname, 1, Nparts[t], Np_all[t], P[t]);
      }

  
  
  return 0;
}


int get_catalog_numparticles( char *name, catalog_header_t *header )
{
  FILE * filein;
  size_t ret;
  char   fname[ strlen(subf_base)+strlen(working_dir)+5 ];
  int    nfiles = 1;
  
  // check whether we are multi-file
  //
  sprintf( fname, "%s", name );  
  if( (filein = fopen( fname, "r" )) == NULL )
    {
      sprintf( fname, "%s.0", name );
      if( (filein = fopen( fname, "r" )) == NULL )
	// unable to find the specified subfind files
	return 0;
      
      do
	{
	  fclose(filein);
	  sprintf( fname, "%s.%d", name, nfiles );
	  nfiles += ((filein = fopen( fname, "r" )) != NULL);
	}
      while( filein != NULL) ;

      sprintf( fname, "%s.0", name );
      filein = fopen( fname, "r" );
    }


  ret = fread( header, sizeof(catalog_header_t), 1, filein );
  fclose(filein);

  if( ret != sizeof(catalog_header_t) )
    return -1;

  if ( nfiles != header->nfiles )
    return -2;

  return nfiles;

}

ull_t get_catalog_data_from_file(char *name, int nfiles, ull_t AllN, ull_t myN, particle_t *target)
{
  size_t ret;
  typedef struct { FILE *ptr; ull_t npart, npart_inc; omp_lock_t lock;} record_file_t;  
  record_file_t *files  = (record_file_t*)calloc( nfiles, sizeof(record_file_t) );
  
  int go = 0;
  do
    {
      catalog_header_t header;
      char fname[strlen(name)+5];
     #define input files[go].ptr
      
      if ( nfiles > 0 )
	sprintf( fname, "%s.%d", name, go );
      else
	sprintf( fname, "%s", name );
      
      files[go].ptr  = fopen( fname, "r" );
      omp_init_lock(&files[go].lock);
		
      ret = fread( &header, sizeof(header), 1, input );
      if ( ret != sizeof header )
	break;

      files[go].npart_inc = files[go].npart;
      files[go].npart_inc += (go > 0 ? files[go-1].npart_inc : 0);
      
      go++;
     #undef input
    }
  while( go < nfiles);

  if ( go <nfiles ) {
    printf("there was an I/O problem in reading the catalog files\n");
    return -1; }

  unsigned int failures = 0;
 #pragma omp parallel
  {
    int   rem     = (int)(AllN % (ull_t)Nthreads);
    PID_t avg_NID = AllN / Nthreads;
    ull_t myoff   = avg_NID * me;
    myoff  += ( me > rem ? rem : me );

    int   file_start        = 0;
    ull_t file_start_offset = 0;

    while( (file_start < nfiles) && (myoff >= files[file_start].npart_inc) )
      file_start++;
    file_start_offset = myoff - (files[file_start].npart_inc - files[file_start].npart);

    unsigned int _failures_ = 0;
    if( file_start < nfiles )
      {
        int fn = file_start;
	PID_t read = 0;
	
	while( (read < myNID) && (fn < nfiles) && !_failures_ )
	  {
	    ull_t howmany_toread = myN - read;
	    howmany_toread = (howmany_toread > files[fn].npart - file_start_offset ?
			      files[fn].npart - file_start_offset : howmany_toread);

	    omp_set_lock( &files[fn].lock );
	    fread( target+read, sizeof(particle_t), howmany_toread, files[fn].ptr);
	    omp_unset_lock( &files[fn].lock );

	    fn += (read < myNID);
	    file_start_offset = 0;
	    
	    if( fn == nfiles ) {
	      _failures_++;
	     #pragma omp atomic update
	      failures++; }
	    else {
	     #pragma omp atomic read
	      _failures_ = failures; }	    
	  }
      }
    else {
     #pragma omp atomic update
      failures++; }
  }

  
  free(files);
  
  return 0;
}




/*  -------------------------------------------------------------
 * 
 *   routines for list inpunt
 *
 *  ------------------------------------------------------------- */


int get_list_numparticles( char *, list_header_t * );
ull_t get_list_data_from_file(char *, int, ull_t, ull_t, list_t *);


int get_list_ids ( char *name )
{

  // find how many particles are in files 
  //

  list_header_t header;

  Nl = 0;
  int nfiles = get_list_numparticles( name, &header);
  if ( nfiles == 0 )
    {
      printf("unable to find the list file both as %s and as %s.0\n",
	     name, name);
	    return -1;
    }
  
  if ( nfiles != header.nfiles )
    {
      printf("There a problem with list: the expected number "
	     "of files is %d but I've found %d files\n",
	     header.nfiles, nfiles);
      return -2;
    }

  if ( header.id_size != sizeof(PID_t) )
    {
      printf("the id size in provided list is incompatible "
	     "with the one used here: %d bytes instead of %lu\n",
	     header.id_size, sizeof(PID_t) );
      return -3;
    }
  
  Nl = header.Nparts_total;

 #pragma omp parallel
  {
    int rem = (int)(Nl % (ull_t)Nthreads);
    myNl    = Nl / Nthreads + (me < rem);
  }


  // allocate memory to hold all the particles
  // and set-up pointers for each particles type
  //
  int failures = 0;
 #pragma omp parallel
  {
    List = (list_t*)calloc( myNl, sizeof(list_t));
   #pragma omp atomic update
    failures += (List == NULL);
  }

  if ( failures )
    {
      #pragma omp parallel
      {
	if( List != NULL )
	  free(List);
      }
      printf("I've got a problem in memory allocation\n");
      return -3;
    }

  int typed = get_list_data_from_file( name, 1, myNl, Nl, List);
  
  return typed;
}


int get_list_numparticles( char *name, list_header_t *header )
{
  FILE * filein;
  size_t ret;
  char   fname[ strlen(subf_base)+strlen(working_dir)+5 ];
  int    nfiles = 1;
  
  // check whether we are multi-file
  //
  sprintf( fname, "%s", name );  
  if( (filein = fopen( fname, "r" )) == NULL )
    {
      sprintf( fname, "%s.0", name );
      if( (filein = fopen( fname, "r" )) == NULL )
	// unable to find the specified subfind files
	return 0;
      
      do
	{
	  fclose(filein);
	  sprintf( fname, "%s.%d", name, nfiles );
	  nfiles += ((filein = fopen( fname, "r" )) != NULL);
	}
      while( filein != NULL) ;

      sprintf( fname, "%s.0", name );
      filein = fopen( fname, "r" );
    }


  ret = fread( header, sizeof(list_header_t), 1, filein );
  fclose(filein);

  if( ret != sizeof(list_header_t) )
    return -1;

  if ( nfiles != header->nfiles )
    return -2;

  return nfiles;

}


ull_t get_list_data_from_file(char *name, int nfiles, ull_t myN, ull_t AllN, list_t *target)
{
  size_t ret;
  list_header_t header;
  typedef struct { PID_t pid; int type; } list_pidtype_t;
  typedef struct { PID_t pid; } list_pid_t;
  
  typedef struct { FILE *ptr; ull_t npart, npart_inc; omp_lock_t lock;} record_file_t;
  record_file_t *files  = (record_file_t*)calloc( nfiles, sizeof(record_file_t) );
  
  int go = 0;
  do
    {      
      char fname[strlen(name)+5];
     #define input files[go].ptr
      
      if ( nfiles > 0 )
	sprintf( fname, "%s.%d", name, go );
      else
	sprintf( fname, "%s", name );
      
      files[go].ptr  = fopen( fname, "r" );
      omp_init_lock(&files[go].lock);
		
      ret = fread( &header, sizeof(header), 1, input );
      if ( ret != sizeof header )
	break;

      files[go].npart_inc = files[go].npart;
      files[go].npart_inc += (go > 0 ? files[go-1].npart_inc : 0);
      
      go++;
     #undef input
    }
  while( go < nfiles);

  if ( go <nfiles ) {
    printf("there was an I/O problem in reading the catalog files\n");
    return -1; }

  unsigned int failures = 0;
 #pragma omp parallel
  {
    int   rem     = (int)(AllN % (ull_t)Nthreads);
    PID_t avg_NID = AllN / Nthreads;
    ull_t myoff   = avg_NID * me;
    myoff  += ( me > rem ? rem : me );

    int   file_start        = 0;
    ull_t file_start_offset = 0;

    while( (file_start < nfiles) && (myoff >= files[file_start].npart_inc) )
      file_start++;
    file_start_offset = myoff - (files[file_start].npart_inc - files[file_start].npart);

    unsigned int _failures_ = 0;
    if( file_start < nfiles )
      {
        int fn = file_start;
	PID_t read = 0;
	
	while( (read < myNID) && (fn < nfiles) && !_failures_ )
	  {
	    ull_t howmany_toread = myN - read;
	    howmany_toread = (howmany_toread > files[fn].npart - file_start_offset ?
			      files[fn].npart - file_start_offset : howmany_toread);

	    int typesize = ( header.type_is_present ?
			     sizeof(list_pid_t) :
			     sizeof(list_pidtype_t) );
	    char *buffer = (char*)malloc( howmany_toread * typesize);
	    
	    omp_set_lock( &files[fn].lock );
	    ret = fread( buffer, typesize, howmany_toread, files[fn].ptr);
	    omp_unset_lock( &files[fn].lock );

	    list_t *_target_ = target+read;
	    switch(header.type_is_present )
	      {
	      case 0: {
		for( ull_t i = 0; i < howmany_toread; i++, _target_++ )
		  _target_->pid = ((list_pid_t*)buffer)[i].pid, _target_->type = -1; } break;
	      default: {
		for( ull_t i = 0; i < howmany_toread; i++, _target_++ )
		  _target_->pid = ((list_pidtype_t*)buffer)[i].pid, _target_->type = ((list_pidtype_t*)buffer)[i].type; } break;		
	      }
	    free( buffer );
	    
	    fn += (read < myNID);
	    file_start_offset = 0;
	    
	    if( fn == nfiles ) {
	      _failures_++;
	     #pragma omp atomic update
	      failures++; }
	    else {
	     #pragma omp atomic read
	      _failures_ = failures; }	    
	  }
      }
    else {
     #pragma omp atomic update
      failures++; }
  }

  
  free(files);
  
  return (!header.type_is_present);
}
