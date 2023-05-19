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


// these 2 defines are used to invoke get_particles_num()
//
#define FILE_PARTICLES 0
#define ALL_PARTICLES  1


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


PID_t get_unmasked_pid( PID_t masked_pid, int n_gen, int bitshift )
{
 #define SHIFT_BASE (BASE-1)
  
  PID_t mask = n_gen << (SHIFT_BASE-bitshift);
  return masked_pid | mask;
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

int check_number_of_files ( char *wdir, char *basename )
{
  
  int   nfiles    = 1;
  FILE *filein = NULL;

  char *fname = (char*)malloc(strlen(wdir)+strlen(basename)+10);
  // check whether we are multi-file
  //
  sprintf( fname, "%s/%s", wdir, basename );  
  if( (filein = fopen( fname, "r" )) == NULL )
    {
      sprintf( fname, "%s/%s.0", wdir, basename );
      if( (filein = fopen( fname, "r" )) == NULL )
	{
	  // unable to find the specified subfind files
	  printf("[io] unable to find the file under %s/ "
		 "both as %s and as %s.0\n",		 
		 wdir, basename, basename);
	  free(fname);
	  return -2;
	}
      
      do
	{
	  fclose(filein);
	  sprintf( fname, "%s/%s.%d", wdir, basename, nfiles );
	  nfiles += ((filein = fopen( fname, "r" )) != NULL);
	}
      while( filein != NULL) ;  
    }
  else
    fclose(filein);

  free(fname);
  return nfiles;
}



		      
num_t seek_block ( FILE *file, char name[5] )
/*
 * This function get to the block specified by name
 * in a format-2 file pointed by *file.
 *
 * RETURN VALUE:
 * 0 if the seeking is not successful, otherwise
 * the 4-bytes long tag at the begin of the data block
 * (i.e. the block's length in bytes).
 * The file position is set to the begin of the data 
 * block, right after the length tag.
 *
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
      int   tag1, tag2;
      ret = fread( &tag1, sizeof(int), 1, file);
      off_t pos = ftello( file );
      ret = fseeko( file, tag1, SEEK_CUR );
      ret = fread( &tag2, sizeof(int), 1, file);
      ret = fseeko ( file, pos, SEEK_SET);
      if( (tag1 != tag2) ) {
	printf("[SEEKBLOCK] error in tags of block %s: %d vs %d\n",
	       name, tag1, tag2 );
	return 0; }
	
      if( tag1 == 0 )
	printf("[SEEKBLOCK] error 0-valued tags for block %s\n", name);
	
      return (num_t)tag1;
    }
  else  
    return 0;
}


int get_snap_detail( FILE *file, int detail, void *result )
{
  int ret = seek_block( file, "HEAD");
  
  if ( ret == 0)
    {
      printf("unable to find the HEAD block in snapshot files\n");
      return -1;
    }

  snapheader_t header;
  
  ret = fread( &header, sizeof(snapheader_t), 1, file);

  switch(detail)
    {
    case EXP_FACTOR: *(double*)result = header.time; break;
    }

  return 0;
}


int get_particles_num( FILE *file, nparts_t *nparts, int mode )
{
  
  int ret = seek_block( file, "HEAD");
  
  if ( ret == 0)
    {
      printf("unable to find the HEAD block in snapshot files\n");
      return -1;
    }

  snapheader_t header;
  
  ret = fread( &header, sizeof(snapheader_t), 1, file);
  
  (*nparts)[ALL] = 0;

  switch( mode )
    {
    case  ALL_PARTICLES: {       // get the total number of particles
      for( int i = 0; i < NTYPES; i++ )
	(*nparts)[ALL] += ((*nparts)[i] =
			   (header.npartTotalHighWord[i]<<31) + header.npartTotal[i]); } break;
    case FILE_PARTICLES: {      // get the snapshot's number of particles    
    for( int i = 0; i < NTYPES; i++ )
      (*nparts)[ALL] += ((*nparts)[i] = header.npart[i]); } break;
    }

  return header.num_files;
}



int get_subfind_data(char *working_dir, char *subf_base, num_t HowMany[4] )
/*
 * this function gets the subfind data, i.e. the
 * arrays of pairs (galaxy, starid )
 *
 */
{
  
  FILE * filein;
  size_t ret = 0;
  int    multifile = 0, nfiles = 1;
  char   fname[ strlen(subf_base)+strlen(working_dir)+5 ];

  // as first, let's find out how many files there are
  // and how many particles do they have in substructures
  // so that we can distribute the file reading among threads
  //

  
  // check whether we are multi-file
  //
  multifile = ( (nfiles = check_number_of_files(working_dir, subf_base)) > 1);

  // now reconstruct how many particles there are in total in
  // the haloes and in the subhaloes and assign equal bunches
  // of them to each thread
  //

  #define Fof   0
  #define SubH  1
  #define FofP  2
  #define SubHP 3
  memset( HowMany, 0, 4*sizeof(num_t) );
  
  {
    nparts_t nparts;
    
    if ( multifile )
      sprintf( fname, "%s/%s.0", working_dir, subf_base );
    else
      sprintf( fname, "%s/%s", working_dir, subf_base );
    filein = fopen( fname, "r" );

    // skip some bytes in the header and get to the total number
    // of particles (which actually is:
    //   type 0: Fof haloes
    //   type 1: Sub Haloes
    //   type 2: how many particles in Fof haloes
    
    int numfiles = get_particles_num( filein, &nparts, ALL_PARTICLES );
    for( int i = 0; i < 3; i++ )
      HowMany[i] = nparts[i];
    
    if( numfiles != nfiles )
      {
	fclose(filein);
	printf("\t[sf data] the number of files specified in the header (%d) "
	       "is different than the number of files found (%d)\n",
	       numfiles, nfiles );
	return -3;
      }
    
    fclose(filein);
  }  
  
  typedef struct { FILE *ptr; num_t npart, npart_inc; omp_lock_t lock;} record_file_t;
  typedef struct { num_t npart, offset; int parent, index; } record_halo_t;
  typedef struct { num_t npart, offset, nsubh; } record_fof_t;
  
  record_file_t *files  = (record_file_t*)calloc( nfiles, sizeof(record_file_t) );
  record_halo_t *haloes = (record_halo_t*)calloc( HowMany[SubH]+1, sizeof(record_halo_t) );
  record_fof_t  *fofs   = (record_fof_t*)calloc( HowMany[Fof]+1, sizeof(record_fof_t) );

  // now we need to read all the haloes to discover how many particles in total live
  // therein, because that useful number is not present in the files' header

  // >>> NOTE: for the sake of simplicity, we load up all the fofs and haloes even
  //           if at command line a fof/halo masking has been specified
  
  PID_t Nhaloes            = 0;  
  PID_t Nfofs              = 0;
  int   go                 = 0;
  unsigned int hidx        = 0;  
  unsigned int prev_parent = 0;
  
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
      
      DPRINT(2,-1, "\t[sf data] opening file %d: %s\n", go, fname);
      
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
	      hidx = ( buffer[i] == prev_parent ? hidx : 0); // if the parent has changed, sets hidx to zero
	      haloes[ Nhaloes + i ].index  = hidx;           // hidx is the ordinal number of this halo in the parent fof
	      hidx++;
	      prev_parent = buffer[i];
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


  assert( Nfofs == HowMany[Fof] );

  assert( Nhaloes == HowMany[SubH] );
  
  fofs[ HowMany[Fof] ].npart  = 0;
  fofs[ HowMany[Fof] ].offset = HowMany[FofP];

  haloes[ HowMany[SubH] ].npart  = 0;
  haloes[ HowMany[SubH] ].offset = HowMany[FofP];
  
  // determines which is the number of
  // particles to be distributed among threads
  //
  Np = HowMany[FofP];
  
  printf("\t[sf data] %llu total particles belong to %llu sub-haloes\n", HowMany[SubHP], HowMany[SubH]);
  
  int failures = 0;
  #pragma omp parallel
  {
    // the avg # of particles per thread
    //
    PID_t avg_Np  = Np / Nthreads;

    // the actual number accounts for the reminder
    //
    int rem       = (int)(Np % (num_t)Nthreads);
    myNp          = avg_Np + (me < rem);
    PPP           = (particle_t*)calloc( myNp, sizeof(particle_t));    

    // myoff is the starting particle of my bunch
    // of particles
    //
    num_t myoff   = avg_Np * me;
    myoff  += ( me > rem ? rem : me );

    
    ll_t  fofn    = 0;
    ll_t  hn      = 0;
    int   fn      = 0;

    num_t in_fof_off  = 0;
    num_t halo_read   = 0;
    num_t fof_read    = 0;
    num_t read        = 0;

    int  _failures_ = 0;
    num_t particles_off;

    // find the fof to which my first particle belongs
    //    
    while( (fofn < (ll_t)Nfofs) && (myoff > fofs[fofn+1].offset) )   // note: using fofn+1 is safe because
      fofn++;						             // fofs has Nfofs+1 elements
    
    // find the offset inside that fof
    in_fof_off = myoff - fofs[fofn].offset;
    fof_read   = in_fof_off;

    // find the halo to which my first particle belongs
    //
    while( (hn < (ll_t)Nhaloes) && (myoff > haloes[hn+1].offset) )   // note: same here as for hn+1
      hn++;
    hn += ( myoff > haloes[hn].offset + haloes[hn].npart );    // in case myoff falls at the end of a fof
                                                               // where its haloes are finished and those
                                                               // of the next fof did not start yet.

    // find the offset inside that halo
    // hn is set negative if my first particle falls among
    // the unbound particles of a fof
    //
    if( hn < (ll_t)Nhaloes )
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

    dprint(2,me, "\t[sf data] thread %d is getting %llu particles, "
	   "starting from halo %lld off %llu , shalo %lld off %llu\n",
	   me, myNp, fofn, in_fof_off, hn, halo_read);


    // now read my particles from the appropriate files
    //
    while( (read < myNp) && (fofn < (ll_t)Nfofs) && (hn < (ll_t)Nhaloes) && !_failures_ )
      {

	// particles_off is the current position of reading 
	// i.e. the "coordinate" of the last particles that I've read
	// from the begininning of the particles (i.e. 0)
	//
	particles_off = myoff + read;

	// find which is the file that I must read from
	//
	while( (fn < nfiles) && ( particles_off >= files[fn].npart_inc) )
	  fn++;
	
	if ( fn < nfiles )
	  {

	    // find the offset inside the file
	    // remind: npart_inc is the incremental number of particles,
	    //         accounting also for this file's particles;
	    //         npart is the num of particles in this file
	    //
	    num_t in_file_off = particles_off - (files[fn].npart_inc - files[fn].npart); // that is zero when this is not the first file we read

	    // determine how many particles I must read from this file
	    //
	    num_t howmany_toread = files[fn].npart - in_file_off;

	    // allocate enough room for data
	    //
	    PID_t *buffer = (PID_t*)calloc( howmany_toread, sizeof(PID_t) );
	    
	    omp_set_lock( &files[fn].lock );
	    
	    /* DPRINT( 2, -1, "[ I/O ] [ th %d ] is opening file %u to read %llu particles starting " */
            /*         "from %llu from halo %u with halo_read %llu and in_halo_offset %llu: "" */
	    /*         "total read particles are %llu -> %llu / %llu\n",		     */
	    /* 	    me, fn, howmany_toread, in_file_off, hn, halo_read, in_halo_off, read, read+howmany_toread, (num_t)myNp); */
	    /* DPRINT(2, -1, "[ I/O ] [ th %d ] ---- file %d --- %llu %llu %llu %llu\n", me, fn, */
	    /* 	   haloes[hn].offset, in_halo_off, files[fn].npart_inc, files[fn].npart ); */
	    
	    ret = seek_block( files[fn].ptr, "PID " );
	    if ( ret == 0 ) {                                                             // check that the block exists
	      printf("unable to find the PID block in subfind files\n");
	      free(buffer); free(P[0]); _failures_++; }
	    
	    else {
	      int save_ret = ret;
	      ret = (size_t)((num_t)ret / files[fn].npart);
	      if ( ret != sizeof(PID_t) ) {                                               // check that the ID type lenght matches
		printf("thread %d: ID block in subfind file %d "
		       "has %lu-bytes long data while my ID type has %lu bytes (%d %llu)\n",
		       me, fn, ret, sizeof(PID_t), save_ret, files[fn].npart ); ret = 0; _failures_++; }

	      if ( in_file_off > 0 )                                                      // skip initial data if needed
		{
		  // this returns the current position, i.e. the begin of the block
		  off_t before = ftello( files[fn].ptr );
		  // jump ahead in_file_off data
		  fseeko( files[fn].ptr, (off_t)(in_file_off*sizeof(PID_t)), SEEK_CUR );
		  // get the position again
		  off_t after = ftello( files[fn].ptr );
		  
		  if( (num_t)(after-before) != (num_t)(in_file_off*sizeof(PID_t) )) {     // check that we skipped the right amount of data
		    printf("%d has got problem, in file %d it is at offset "
			   "%llu instead of %llu due to %llu \n",
			   me, fn, (num_t)after, (num_t)before+(in_file_off*sizeof(PID_t)), in_file_off);
		    _failures_++; }
		}

	      ret = fread( buffer, sizeof(PID_t), howmany_toread, files[fn].ptr);         // actually read the data
	      if ( ret != (ull_t)howmany_toread ) {                                       // check that everything went fine
		printf("%d has read %llu instead of %llu from file %d\n",
		       me, (num_t)read, (num_t)howmany_toread, fn);
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
	    
	    for( num_t i = 0; (i < howmany_toread) && (read < myNp); i++, read++ )
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
		    DPRINT(1, -1, "-- chf thread %d has exhausted fof %lld (%llu) with %llu parts\n",
			   me, fofn, fofs[fofn].npart, read );
		    fofn++;
		    //if( fofn < Nfofs )
		    in_fof_off = 0;
		    fof_read = 0;
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
			  // if we are in fof-only region
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

   /* #if defined(DEBUG) */
   /*  fof_table_t *fof_table = calloc( Nfofs, sizeof(fof_table_t)); */
   /*  make_fof_table( fof_table, 0); */
   /*  { */
   /*    char fname[200]; */
   /*    sprintf( fname, "check_fof_table.%d", me); */
   /*    FILE *file = fopen( fname, "w" ); */
      
   /*    for( num_t i = PPP[0].fofid; i <= PPP[myNp-1].fofid; i++ ) { */
   /* 	fprintf(file, "%10llu\t%llu", i, fof_table[i].TotN); */
   /* 	for( int j = 0; j < NTYPES; j++ ) */
   /* 	  fprintf(file, "\t%10llu", fof_table[i].Nparts[j] ); */
   /* 	fprintf(file, "\n"); }	     */
      
   /*    fclose(file);	   */
   /*  } */
   /* #endif */
    
    
    fprintf( details, "thread %d got %llu particles from subfind data\n", me, myNp);
  }

  fprintf( details, "\n");

  
  
  for( int i = 0; i < nfiles; i++ )
    fclose( files[i].ptr);
  free(fofs);
  free(haloes);  
  free(files);
  
  return failures;
}





int get_id_data(char *working_dir, char *name)
/*
 * this function gets the ids from the snapshot
 *
 */
{

  FILE   *filein;
  size_t  ret = 0;
  int     multifile = 0, nfiles = 1;
  char    fname[ strlen(subf_base)+strlen(working_dir)+5 ];
  num_t   AllN = 0;

  // as first, let's find out how many files there are
  // and how many particles do they have in substructures
  // so that we can distribute the file reading among threads
  //

  
  // check whether we are multi-file
  //
  multifile = ( (nfiles = check_number_of_files(working_dir, name)) > 1);
  
  {
    nparts_t nparts;
    
    if ( multifile )
      sprintf( fname, "%s/%s.0", working_dir, name );
    else
      sprintf( fname, "%s/%s", working_dir, name );
    filein = fopen( fname, "r" );

    int numfiles = get_particles_num( filein, &nparts, ALL_PARTICLES );
    AllN = nparts[ALL];
    
    if( numfiles != nfiles )
      {
	fclose(filein);
	printf("the number of file specified in the header (%d) of the snapshot file #0 is different than the number of files found (%d)\n",
	       numfiles, nfiles );
	return -3;
      }
    
    fclose(filein);
  }  

 #pragma omp parallel
  {
    int rem  = (int)(AllN % (num_t)Nthreads);
    myNID    = AllN / Nthreads + (me < rem);
    IDs      = (pidtype_t*)calloc( myNID, sizeof(pidtype_t));
  }
  
  typedef struct { FILE *ptr; ul_t npart[NTYPES]; num_t npart_all, npart_all_inc; omp_lock_t lock;} record_file_t;  
  record_file_t *files  = (record_file_t*)calloc( nfiles, sizeof(record_file_t) );

  int go = 0;
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
    // in this parallel region each thread determines which section
    // of ID data it should load from files
    //


    // calculate the limits for this thread
    //
    int rem       = (int)(AllN % (num_t)Nthreads);
    PID_t avg_NID = AllN / Nthreads;
    num_t myoff   = avg_NID * me;
    myoff  += ( me > rem ? rem : me );

    int   file_start        = 0;
    num_t file_start_offset = 0;

    // determine the first file of interest for this thread
    //
    while( (file_start < nfiles) && (myoff >= files[file_start].npart_all_inc) )
      file_start++;
    // find the offset in this file
    file_start_offset = myoff - (files[file_start].npart_all_inc - files[file_start].npart_all);

    if( file_start < nfiles )
      {
	unsigned int _failures_ = 0;
	int          fn         = file_start;
	PID_t        read       = 0;

	while( (read < myNID) && (fn < nfiles) && !_failures_ )
	  // continue reading until the desired number of
	  // data have not been loaded
	  // 
	  {
	    num_t npart_local_inc[NTYPES] = {0};  // this array stores the incremental nr. of
						  // particles in this file, by their type
	    num_t howmany_toread = 0;
	    int   type = 0;

	                                          // calculate the incremental number of
	                                          // particles of each type
	    npart_local_inc[0] = files[fn].npart[0];
	    for( int i = 1; i < NTYPES; i++ )
	      npart_local_inc[i] = files[fn].npart[i] + npart_local_inc[i-1];

	                                          // find the type of the first particle
	    while( file_start_offset >= npart_local_inc[type] )
	      type++;
	                                          // find how many particles must be
	                                          // read
	    
	    howmany_toread = myNID - read;        // these are all the particles that remains

						  // check whether there are less particles 
						  // available in the current file	    
	    howmany_toread = (howmany_toread > files[fn].npart_all - file_start_offset ?
			      files[fn].npart_all - file_start_offset : howmany_toread);

	                                          // allocate the memory needed 
	    PID_t *buffer = (PID_t*)calloc( howmany_toread, sizeof(PID_t) );

	                                          // protect the file access 
	    omp_set_lock( &files[fn].lock );

	                                          // find the entry point of the ID block
	    size_t  myret = seek_block(files[fn].ptr, "ID  ");
	    if ( myret )
	      {
		myret = (size_t)((num_t)myret / files[fn].npart_all);
		if ( myret != sizeof(PID_t) ) {
		  printf("thread %d: ID block in snap file %d has %lu-bytes long data while my ID type has %lu bytes\n",
			 me, fn, ret, sizeof(PID_t) ); myret = 0; }
	      }
	    if ( myret )
	      {
		// move to the first position for this thread
		myret = fseeko(files[fn].ptr, (off_t)(sizeof(PID_t)*file_start_offset), SEEK_CUR);
		// load the data
		fread(buffer, sizeof(PID_t), howmany_toread, files[fn].ptr);

		// assign ids to our data structure
		for ( num_t i = 0; i < howmany_toread; read++, i++ )
		  {
		    type += (file_start_offset+i >= npart_local_inc[type]);  // check whether the type has changed
		    IDs[read].pid  = buffer[i];
		    IDs[read].type = type;
		    IDs[read].file = fn;
		    IDs[read].pos  = file_start_offset + i;
		  }

		// release the file access
		omp_unset_lock( &files[fn].lock );
		// release memory               
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
  
  for( int n = 0; n < nfiles; n++ )
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
  num_t start = (type > 0? type_positions[type-1] : 0);
  num_t stop  = type_positions[type];

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
num_t get_catalog_data_from_file(char *, int, int, num_t );

int get_catalog_data(char *name, int *types)
{


  // find how many particles are in files for each type
  // of particles
  //

  catalog_header_t header;
  char fname[ CATALOG_NAME_SIZE + 5];

  Nparts_all[0] = (num_t*)malloc( NTYPES*Nthreads*sizeof(num_t));
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
	    printf("There was a problem with catalogs for type %d: "
		   "the expected number of files is %d but I've found %d files\n",
		   t, header.nfiles, nfiles);
	    return -2;
	  }
	dprint(1, 0, "%llu particles of type %d from catalog %s\n",
	       header.Nparts_total, t, fname );
	Np_all[t] = header.Nparts_total;
	Np += Np_all[t];
       #pragma omp parallel
	{
	  int rem            = (int)(Np_all[t] % (num_t)Nthreads);
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
	Pbase = (particle_t*)calloc( myNp + n_stars_generations,  // we use a single allocation
				     sizeof(particle_t));         // to hold all the particles
	
	P[t0] = Pbase;                                            // then we set-up pointers to
	P_all[t0][me] = P[t0];                                    // each single type
	num_t sum = Nparts[t0];
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
	get_catalog_data_from_file( fname, 1, t, Np_all[t]);	
      }

  
  
  return 0;
}


int get_catalog_numparticles( char *name, catalog_header_t *header )
{
 #define NSIZE (DIR_SIZE+NAME_SIZE+NUM_SIZE+5)
  FILE * filein = NULL;
  size_t ret;
  char   fname[ NSIZE ];
  int    nfiles = 1;
  
  // check whether we are multi-file
  //
  snprintf( fname, NSIZE, "%s", name );  
  if( (filein = fopen( fname, "r" )) == NULL )
    {
      snprintf( fname, NSIZE, "%s.0", name );
      if( (filein = fopen( fname, "r" )) == NULL )
	// unable to find the specified subfind files
	return 0;
      
      do
	{
	  fclose(filein);
	  snprintf( fname, NSIZE, "%s.%d", name, nfiles );
	  nfiles += ((filein = fopen( fname, "r" )) != NULL);
	}
      while( filein != NULL) ;

      snprintf( fname, NSIZE, "%s.0", name );
      filein = fopen( fname, "r" );
    }


  ret = fread( header, sizeof(catalog_header_t), 1, filein );
  fclose(filein);

  if( ret != 1 )
    return -1;

  if ( nfiles != header->nfiles )
    return -2;

  return nfiles;
 #undef NSIZE
}


num_t get_catalog_data_from_file(char *name, int nfiles, int type, num_t AllN )
{
  size_t ret;
  typedef struct { FILE *ptr; num_t npart, npart_inc; omp_lock_t lock;} record_file_t;  
  record_file_t *files  = (record_file_t*)calloc( nfiles, sizeof(record_file_t) );
  
  int go = 0;
  do
    {
      catalog_header_t header;
      char fname[strlen(name)+5];
      
      if ( nfiles > 1 )
	sprintf( fname, "%s.%d", name, go );
      else
	sprintf( fname, "%s", name );
      
      files[go].ptr  = fopen( fname, "r" );
      if( files[go].ptr == NULL ) break;
      omp_init_lock(&files[go].lock);
		
      ret = fread( &header, sizeof(header), 1, files[go].ptr );
      if ( ret != 1 )
	break;

      files[go].npart     = header.Nparts;
      files[go].npart_inc = files[go].npart + (go > 0 ? files[go-1].npart_inc : 0);
      
      go++;
    }
  while( go < nfiles);

  if ( go < nfiles ) {
    printf("[read catalog] there was an I/O problem in opening the catalog files\n");
    for ( int f = 0; f < go; f++ )
      if( files[f].ptr != NULL )
	fclose(files[f].ptr);
    if( files != NULL )
      free(files);
    return -1; }

  int          signal   = 0;
  unsigned int failures = 0;
 #pragma omp parallel
  {
    int   rem          = (int)(AllN % (num_t)Nthreads);
    PID_t avg_NP       = AllN / Nthreads;
    num_t myoff        = avg_NP * me;   
    myoff  += ( me > rem ? rem : me );

    int   file_start        = 0;
    num_t file_start_offset = 0;

    while( (file_start < nfiles) && (myoff > files[file_start].npart_inc) )
      file_start++;
    file_start_offset = myoff - (files[file_start].npart_inc - files[file_start].npart);

    particle_t *_target_ = P[type];
    
    unsigned int _failures_ = 0;
    if( file_start < nfiles )
      {
        int fn = file_start;
	PID_t read = 0;
	
	while( (read < Nparts[type]) && (fn < nfiles) && !_failures_ )
	  {
	    num_t howmany_toread = Nparts[type] - read;

	    howmany_toread = (howmany_toread > files[fn].npart - file_start_offset ?
			      files[fn].npart - file_start_offset : howmany_toread);

	    dprint( 2, me, "[read catalog] TH %d - read %llu data from position %llu of file %d/%d\n",
		    me, (num_t)howmany_toread, sizeof(catalog_header_t)+sizeof(particle_t)*file_start_offset,
		    fn, nfiles );

	    omp_set_lock( &files[fn].lock );
	    {
	      int seek = fseeko(files[fn].ptr, (off_t)(sizeof(catalog_header_t) +
						       sizeof(particle_t)*file_start_offset), SEEK_SET);
	      size_t v = fread(_target_+read, sizeof(particle_t),
				howmany_toread, files[fn].ptr);
	      
	      if( (seek != 0) || (v != (size_t)howmany_toread) ) {
	       #pragma omp atomic update
		failures++; }
	      else
		read += howmany_toread;
	    }
	    omp_unset_lock( &files[fn].lock );


	    fn += (read < myNID);
	    file_start_offset = 0;
	    
	    if( fn == nfiles ) {
	     #pragma omp atomic update
	      failures++; }
	    
	   #pragma omp atomic read
	    _failures_ = failures;
	  }
      }
    else {
     #pragma omp atomic update
      failures++; }

    // in the following block we make sure that the group of
    // particle_t entries that refer to the same masked_id
    // belong to the same task, in order to ease the search
    // operations afterwards.
    //    
   #pragma omp barrier
   #pragma omp for ordered schedule(static)
    for( int th = 0; th < Nthreads; th++ )
      {
	#pragma omp ordered
	{
	  if( signal ) {
	    // some of my initial entries
	    // have been moved by the previous
	    // me-1 thread;
	    // the number of entries is stored
	    // in the variable signal.
	    // Hence, here we discard the first
	    // [signal] entries
	    //
	    myNp            -= signal;
	    Nparts[type]    -= signal;
	    P[type]         += signal;
	    P_all[type][me] += signal; }
	  
	  if( th < Nthreads-1 ) {
	    // not the least thread,
	    // check whether we must unify
	    // the last particle and the
	    // first ones of the next thread
	    //
	    num_t n    = 0;
	    num_t np   = Nparts[type];
	    while( P[type][np-1].pid == P_all[type][th+1][n].pid )
	      P[type][np++] = P_all[type][th+1][n++];
	    Nparts[type]          = np;
	    Nparts_all[type][me]  = np;
	    myNp                 += n;
	    signal                = n; }
	}
      }
    
   #if defined(MASKED_ID_DBG)
    // track a given masked id for debugging purposes
    for( num_t i = 0; i < myNp; i++ )
      if( P[type][i].pid == MASKED_ID_DBG )
	dprint(0, me, "[ID DBG][S][th %d] found particle with masked id %llu, "
	       "type %d, gen %d at pos %llu\n", me,
	       (ull_t)P[type][i].pid, P[type][i].type, P[type][i].gen, i);	
   #endif

  }

  for( int f = 0; f < nfiles; f++ )
    if( files[f].ptr != NULL )
      fclose(files[f].ptr);
  free(files);
  
  return 0;
}




/*  -------------------------------------------------------------
 * 
 *   routines for list inpunt
 *
 *  ------------------------------------------------------------- */


int get_list_numparticles( char *, list_header_t *, int );
num_t get_list_data_from_file( char *, int, int, num_t );


int get_list_ids ( char *name, int file_type )
{

  // find how many particles are in files 
  //

  list_header_t header;

  Nl = 0;
  int nfiles = get_list_numparticles( name, &header, file_type);
  if ( nfiles == 0 )
    {
      printf("unable to find the list file both as %s and as %s.0\n",
	     name, name);
	    return -1;
    }
  
  if ( nfiles != header.nfiles )
    {
      printf("There was a problem with list: the expected number "
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
  dprint(1, 0, "[read list] %llu particles from list file %s\n",
	 header.Nparts_total, name );

  // allocate memory to hold all the particles
  // and set-up pointers for each particles type
  //
  int failures = 0;
  
 #pragma omp parallel
  {
    int rem = (int)(Nl % (num_t)Nthreads);
    myNl    = Nl / Nthreads + (me < rem);

    List = (list_t*)calloc( myNl, sizeof(list_t));
    if( file_type == -1 )
      Types = (int*)calloc( myNl, sizeof(int));
   #pragma omp atomic update
    failures += (List == NULL);
    dprint(2, me, "[read list] thread %d allocated %llu bytes at %p\n",
	   me, myNl * sizeof(list_t), List);
  }

  if ( failures )
    {
      #pragma omp parallel
      {
	if( List != NULL )
	  free(List);
	if( Types != NULL )
	  free(Types);
      }
      printf("I've got a problem in memory allocation\n");
      return -3;
    }

  int ret = get_list_data_from_file( name, 1, file_type, Nl);
  
  return ret;
}


int get_list_numparticles( char *name, list_header_t *header, int type )
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
	// unable to find the files
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

  if( type == -1 )
    ret = fread( header, sizeof(list_header_t), 1, filein );
  else {
    ret = fread ( &(header->id_size), sizeof(int), 1, filein);
    ret = fread ( &(header->Nparts), sizeof(num_t), 1, filein);
    header->Nparts_total = header->Nparts;                     // multiple list files not
    header->nfiles = 1; }                                      // yet implemented
    
  fclose(filein);

  if( (type == -1 ) && (ret != sizeof(list_header_t)) )
    return -1;

  if ( nfiles != header->nfiles )
    return -2;

  return nfiles;

}


num_t get_list_data_from_file(char *name, int nfiles, int file_type, num_t AllN)
{
  size_t ret;
  list_header_t header;
  typedef struct __attribute__((packed)) { PID_t pid; int type; } list_pidtype_t;
  
  typedef struct { FILE *ptr; num_t npart, npart_inc; omp_lock_t lock;} record_file_t;
  record_file_t *files  = (record_file_t*)calloc( nfiles, sizeof(record_file_t) );
  
  int go          = 0;
  int header_size = 0;
  do
    {      
      char fname[strlen(name)+5];
      
      if ( nfiles > 1 )
	sprintf( fname, "%s.%d", name, go );
      else
	sprintf( fname, "%s", name );
      
      files[go].ptr  = fopen( fname, "r" );
      if( files[go].ptr == NULL ) break;
      omp_init_lock(&files[go].lock);

      if( file_type == -1 ) {
	header_size = sizeof(header);
	ret = fread( &header, sizeof(header), 1, files[go].ptr );
	if ( ret != sizeof header )
	  break; }
      else {
	header_size = sizeof(int)+sizeof(num_t);
	header.type_is_present = 0;
	ret = fread ( &header.id_size, sizeof(int), 1, files[go].ptr );
	ret = fread ( &header.Nparts, sizeof(num_t), 1, files[go].ptr );
	header.Nparts_total = header.Nparts; }

      files[go].npart      = header.Nparts;
      files[go].npart_inc  = files[go].npart;
      files[go].npart_inc += (go > 0 ? files[go-1].npart_inc : 0);
      
      go++;
    }
  while( go < nfiles);

  if ( go < nfiles ) {
    printf("[read list ] there was an I/O problem in opening the list files\n");
    for ( int f = 0; f < go; f++ )
      if( files[f].ptr != NULL )
	fclose(files[f].ptr);
    return -1; }

  int typesize = ( header.type_is_present ?
		   sizeof(list_pidtype_t) :
		   sizeof(PID_t) );

  unsigned int failures = 0;
 #pragma omp parallel
  {
    int   rem     = (int)(AllN % (num_t)Nthreads);
    PID_t avg_NL  = AllN / Nthreads;
    num_t myoff   = avg_NL * me;
    myoff  += ( me > rem ? rem : me );

    int   file_start        = 0;
    num_t file_start_offset = 0;

    while( (file_start < nfiles) && (myoff >= files[file_start].npart_inc) )
      file_start++;
    file_start_offset = myoff - (files[file_start].npart_inc - files[file_start].npart);    
    
    unsigned int _failures_ = 0;
    if( file_start < nfiles )
      {
        int fn = file_start;
	PID_t read = 0;
	
	while( (read < myNl) && (fn < nfiles) && !_failures_ )
	  {
	    num_t howmany_toread = myNl - read;
	    howmany_toread = (howmany_toread > files[fn].npart - file_start_offset ?
			      files[fn].npart - file_start_offset : howmany_toread);

	    dprint( 2, me, "[read list] TH %d - read %llu data from position %llu of file %d/%d\n",
		    me, (num_t)howmany_toread, (header_size + typesize*file_start_offset), fn, nfiles );
	    
	    char *buffer = (char*)malloc( howmany_toread * typesize);

	    if( buffer != NULL )
	      {
		omp_set_lock( &files[fn].lock );
		
		int seek = fseeko(files[fn].ptr,
				  (off_t)(header_size + typesize*file_start_offset), SEEK_SET);
		ret = fread( buffer, typesize, howmany_toread, files[fn].ptr);
		
		if( (seek > 0) || ((num_t)ret != howmany_toread) ) {
		  _failures_++;
		 #pragma omp atomic update
		  failures++; }
		
		omp_unset_lock( &files[fn].lock );

		if( !_failures_)
		  {
		    list_t *_target_ = List+read;
		    int    *_types_ = NULL;
		    if( Types != NULL )
		      _types_ = Types+read;
		    switch(header.type_is_present )
		      {
		      case 0: {
			for( num_t i = 0; i < howmany_toread; i++, _target_++ )
			  _target_->pid = ((PID_t*)buffer)[i]; } break;
		      default: {
			for( num_t i = 0; i < howmany_toread; i++, _target_++, _types_++ )
			  _target_->pid = ((list_pidtype_t*)buffer)[i].pid, *_types_ = ((list_pidtype_t*)buffer)[i].type; } break;		
		      }
		    free( buffer );
		
		    read += howmany_toread;
		    fn += (read < myNID);
		    file_start_offset = 0;
		
		    if( fn == nfiles ) {
		      _failures_++;
		     #pragma omp atomic update
		      failures++; }
		    else {       		      
		     #pragma omp atomic read            // update the failures value
		      _failures_ = failures; }          // from other threads
		  }  // close < if( !_failures_ ) >
	      }  // close < if(buffer != NULL) >
	    else
	      {
		// there was a problem allocating buffer
		dprint( 0, me, "[read list] TH %d - unable to allocate %llu bytes "
			"to read data from %d/%d\n",
			me, (num_t)(howmany_toread*typesize), fn, nfiles );
		_failures_++;
	       #pragma omp atomic update
		failures++;
	      }
	  }  // close while loop on files
      }
    else {
     #pragma omp atomic update
      failures++; }
  }


  for( int f = 0; f < nfiles; f++ )
    if( files[f].ptr != NULL )
      fclose(files[f].ptr);
  free(files);
  
  return 0;
}



/*  -------------------------------------------------------------
 * 
 *   routines for generating snapshot-like files of fofs and
 *   galaxies
 *
 *  ------------------------------------------------------------- */



typedef struct { FILE *ptr; num_t nparts; } file_t;
typedef struct { int type, fofnum, ngals; } fofheader_t; 


int write_block_header( FILE *file, head_t *header )
// writes a block header into a file
// note: the label and the size must be set
//
{
  if( file == NULL )
    return 1;
  header -> tag1 = 8;
  header -> tag2 = 8;
  fwrite( header, sizeof(head_t), 1, file );
  return 0;
}



int write_block( FILE *file, head_t *header, char *block, int size )
// writes a block into a file
// note: the label in hte header must be set
//       size is the data size
{
  if( (file == NULL) ||
      (block == NULL ) ||
      (header == NULL) ||
      (size == 0) )
    return 1;
  
  header->len = size + 8;    // add the size of the 2 tags in the block
  write_block_header( file, header );
  fwrite( &size, sizeof(int), 1, file );
  fwrite( block, 1, size, file );
  fwrite( &size, sizeof(int), 1, file );
  return 0;
}



num_t get_block_from_file( char field[5], FILE *ptr, void *array )
{

  size_t n   = seek_block(ptr, field );
  size_t ret = fread( array, 1, n, ptr );

  return (num_t)ret;
}


int write_fofgal_ids( char *wdir, char *snapname, int fofid )
{
  UNUSED(wdir);
  FILE *file;
  {
    int   len = strlen(snapname);
    char *filename = (char*)malloc( len + 100 );
    char  string[99] = {0};
    sprintf( string, "_ids_fof_%04d", fofid );
    sprintf( filename, "%s%s%s", snapname, "_masked", string );
    file = fopen( filename, "w" );
  }

  
 #pragma omp parallel for ordered
  for( int t = 0; t < Nthreads; t++ )
    {
      for( num_t i = 0; i < myNp; i++ ) {
	switch(PPP[i].fofid == fofid) {
	case 1: fwrite(&PPP[i].pid, sizeof(PID_t), 1, file ); break; }
      }
    }

  fclose(file);

  return 0;
}

int write_fofgal_snapshots( char *wdir, char *snapname,
			    mask_galaxies_in_fof_t *masks, int nmasks )

{

  // [ 1 ]
  // create all the new snap-like files
  //
  FILE **FofGalFiles = (FILE**)malloc( nmasks * sizeof(FILE*) );
  int    failures    = 0;
  {
    int     len = strlen(snapname);
    char   *filename = (char*)malloc( len + 100 );
    for( int n = 0; n < nmasks; n++ )
      {
	char string[99] = {0};
	sprintf( string, "_fof_%04d", masks[n].fof_num );	
	if( masks[n].ngal >= 0 )
	  sprintf( &string[strlen(string)], "_gal_%05d", masks[n].ngal);
	sprintf( filename, "%s%s%s", snapname, "_masked", string );

	FofGalFiles[n] = fopen( filename, "w" );
	failures += ( FofGalFiles[n] == NULL );
      }
  }

  if( failures > 0 )
    {
      printf ( "%d mask files failed to create\n", failures );
      free( FofGalFiles );
      return -1;
    }

  int nfiles = check_number_of_files(wdir, snapname);
  
  char  *p_is_masked[Nthreads];
  num_t *file_positions[Nthreads];
  num_t *snapout_sizes = calloc(nmasks, sizeof(num_t));

  nparts_t masked_parts[Nthreads];
  num_t    masked_in_file[Nthreads][nfiles];
  memset( masked_parts, 0, sizeof(nparts_t)*Nthreads );
  memset( masked_in_file, 0, sizeof(num_t)*Nthreads*nfiles );
  
 #define p_idx(i) ((i)/8)
 #define p_bytemask(i,s) (char)((s) << (i%8))
  
 #pragma omp parallel reduction(+:snapout_sizes[0:nmasks])
  {
    // [ 2 ]
    //
    // partition particles by the file they belong to
    //

    num_t *fpositions = (num_t*)calloc( nfiles, sizeof(num_t) );
    memset(fpositions, 0, sizeof(num_t)*nfiles);
    file_positions[me] = fpositions;

    num_t franges[nfiles];
    for( int n = 0; n < nfiles; n++ )
      franges[n] = n;
    
    k_way_partition(0, myNp, 0, nfiles-2, franges, fpositions, 3);

   #if defined(DEBUG)
    {
      num_t fails = 0;
      for( num_t i = 1; i < myNp; i++ )
	fails += (PPP[i].file < PPP[i-1].file);

      if( fails > 0 )
	printf("\t[thread %d] has %llu errors in by-file partitioning\n",
	       me, fails );
    }
   #endif


    // [ 3 ]
    // allocate room for particle structures to be written
    // in the files at each thread
    //
    // allocate an array to check whether a given particle
    // is masked or not        

    char *_pmasked_ = (char*)calloc( (myNp/8+1), 1);
    p_is_masked[me] = _pmasked_;
    
    nparts_t _masked_parts_ = {0};
    num_t    _masked_in_file_[nfiles];
    memset( _masked_in_file_, 0, sizeof(num_t)*nfiles );
    
    for ( num_t i = 0; i < myNp; i++ )
      {
	int first = 1;
	for( int m = 0; m < nmasks; m++ )
	  {
	    int s = ( (PPP[i].fofid == masks[m].fof_num) &&
		      ((PPP[i].gid == masks[m].ngal) ||
		       (masks[m].ngal == -1)) );
	    snapout_sizes[m]          += s;
	    
	    if( s && first ) {
	      _masked_parts_[PPP[i].type] += s;
	      _masked_parts_[ALL]         += s;	    	  
	      _pmasked_[p_idx(i)]         |= p_bytemask(i,s);
	      first = 0; }
	  }
      }
    
    outsnap_data = (outsnap_t*)malloc( _masked_parts_[ALL] * sizeof(outsnap_t) );
    for ( int i = 0; i <= ALL; i++ )
      masked_parts[me][i] = _masked_parts_[i];

    // [ 4 ]
    //
    // fill some info for every particle
    
    num_t p = 0;
    for ( num_t i = 0; i < myNp; i++ )
      {
	int is_masked = ((_pmasked_[p_idx(i)] & (p_bytemask(i,1))) != 0);
	switch(is_masked) {
	case 1: {
	  PID_t pid = PPP[i].pid;
	  _masked_in_file_[PPP[i].file]++;
	  switch( PPP[i].gen > 0 ) {
	  case 1: pid = get_unmasked_pid( pid, PPP[i].gen, id_bitshift ); break; }
	  outsnap_data[p].idx.pid   = pid;
	  outsnap_data[p].idx.fofid = PPP[i].fofid;
	  outsnap_data[p].idx.gid   = PPP[i].gid;
	  outsnap_data[p].type      = PPP[i].type; p++; } break; }
      }

    for ( int i = 0; i < nfiles; i++ )
      masked_in_file[me][i] = _masked_in_file_[i];

   #if defined(DEBUG)
    if ( p != _masked_parts_[ALL] )
      printf(">>>>>> W AR N I N G <<<<<<< thread %d : mismatch in masked parts: "
	     "found %llu instead of %llu (line %d of file %s)",
	     me, p, _masked_parts_[ALL], __LINE__, __FILE__ );
   #endif
      
    dprint(1, -1, "\t[thread %d] %llu particles masked\n", me, (unsigned long long)p);  // just for debug

    
    }
  
  // [ 5 ]
  // allocate buffer for storing blocks from snapshot files
  //
  // here we should parse the required blocks
  // to individuate the largest one
  // ...
  // for the moment, we skip this and assume
  // that the largest is the pos / vel which
  // contain 3 entries per particle
  //
  
                                                                     // note: in the following,
						                     // - sizeof_out_data is either 4 or 8 bytes
						                     // depending on the desired precision
						                     // - max_data_multiplicity is the largest
						                     // amount of data per particles per block
						                     // (i.e. 3 for positions or velocities, 1
						                     // for most quantities, NMet for metallicities
						                     // and so on)
  
  int max_data_multiplicity = 3;  
  int max_size_per_particle = max_data_multiplicity*sizeof_out_data; 
  
  // [ 6 ]
  // loops over snapshot files and get particles
    
  double  exp_factor;
  char   *fname = (char*)malloc( strlen(wdir) + strlen(snapname) + 10 );  

  num_t   progress[Nthreads];
  memset( progress, 0, sizeof(num_t)*Nthreads );


  for ( int f = 0; f < nfiles; f++ )
    {
      FILE  *filein;      
      
      if ( nfiles > 1 )
	sprintf( fname, "%s/%s.%d", wdir, snapname, f );
      else
	sprintf( fname, "%s/%s", wdir, snapname );
      filein = fopen( fname, "r" );

      if ( f == 0 )
	get_snap_detail( filein, EXP_FACTOR, &exp_factor);
      
      nparts_t nparts;
      get_particles_num( filein, &nparts, FILE_PARTICLES ); 

      char *block = (char*)malloc(max_size_per_particle * nparts[ALL]);

      snapheader_t snapheader;
      get_block_from_file( "HEAD", filein, (void*)&snapheader );

      get_block_from_file( "POS ", filein, (void*)block );      
      
     #pragma omp parallel
      {
	if ( masked_in_file[me][f] > 0 )
	  {	
	    char  *_pmasked_ = p_is_masked[me];
	    num_t  p         = progress[me];
	    num_t  i         = (f == 0 ? 0 : file_positions[me][f-1]);
	    
	    float_out halfbox = (float_out)snapheader.BoxSize / 2.0;
	    
	    for( ; (i < myNp) && (PPP[i].file==f); i++ )
	      {
		int is_masked = ((_pmasked_[p_idx(i)] & (p_bytemask(i,1))) != 0);
		switch( is_masked )
		  {
		  case 1: { switch( sizeof_in_data == sizeof(float) ) {
		      case 0: { for( int j = 0; j < 3 ; j++ )
			    outsnap_data[p].pos[j] = (float_out)*((double*)block + PPP[i].pos*3 + j); } break;
		      case 1: { for( int j = 0; j < 3 ; j++ )
			    outsnap_data[p].pos[j] = (float_out)*((float*)block + PPP[i].pos*3 + j); } break; }
		      
		      switch( shiftbox ) {
		      case 1: { for( int j = 0; j < 3 ; j++ ) outsnap_data[p].pos[j] -= halfbox; } break; }
		      
		      p++; } break;
		  default: break;
		  }
	      }
	    dprint(2,-1, "\t\tth %d has got %llu masked particles from file %d\n", me, p-progress[me], f);
	  }	
      }

      get_block_from_file( "VEL ", filein, (void*)block );
      
     #pragma omp parallel
      {
	if ( masked_in_file[me][f] > 0 )
	  {		    
	    char  *_pmasked_ = p_is_masked[me];
	    num_t  i         = (f == 0 ? 0 : file_positions[me][f-1]);
	    num_t  p         = progress[me];
	    for( ; (i < myNp) && (PPP[i].file==f); i++ )
	      {
		int is_masked = ((_pmasked_[p_idx(i)] & (p_bytemask(i,1))) != 0);
		switch( is_masked )
		  {
		  case 1: { double v2 = 0;
		      if( sizeof_in_data == sizeof(float) ) {
			for( int j = 0; j < 3 ; j++ )
			  v2 += *((float*)block + PPP[i].pos*3 + j)* *((float*)block + PPP[i].pos*3 + j) * exp_factor; }
		      else { for( int j = 0; j < 3 ; j++ )
			  v2 += *((double*)block + PPP[i].pos*3 + j)* *((double*)block + PPP[i].pos*3 + j) * exp_factor; }
		      outsnap_data[p].vel2 = v2; p++; }
		  default: break;
		  }
	      }
	  }
      }


      get_block_from_file( "MASS", filein, (void*)block );
      
     #pragma omp parallel
      {
	if ( masked_in_file[me][f] > 0 )
	  {		    
	    char *_pmasked_ = p_is_masked[me];
	    double masstable[NTYPES];
	    memcpy( masstable, snapheader.mass, sizeof(double)*NTYPES );

	    num_t offset[NTYPES] = {0};
	    for( int t = 1; t < NTYPES; t++ )
	      {
		for ( int j = 0; j < t; j++ )
		  offset[t] += nparts[j]*(masstable[j]>0);
	      }
	
	    num_t i = (f == 0 ? 0 : file_positions[me][f-1]);
	    num_t p = progress[me];
	    for( ; (i < myNp) && (PPP[i].file==f); i++ )
	      {
		int is_masked = ((_pmasked_[p_idx(i)] & (p_bytemask(i,1))) != 0);
		double mass = masstable[PPP[i].type];
		switch( is_masked )
		  {
		  case 1: {
		    switch( mass == 0 )
		      {
		      case 1: {		    
			num_t n = PPP[i].pos - offset[PPP[i].type];
			if( sizeof_in_data == sizeof(float) )
			  mass = *((float*)block + n);
			else 
			  mass = *((double*)block + n);
		      } break;
		      default: break;
		      }
		    outsnap_data[p].vel2 *= mass; p++;
		  }
		  default: break;
		  }
	      }
	  }
      }

      get_block_from_file( "POT ", filein, (void*)block );
      
     #pragma omp parallel
      {
	if ( masked_in_file[me][f] > 0 )
	  {		    	    
	    char  *_pmasked_ = p_is_masked[me];
	    num_t  i         = (f == 0 ? 0 : file_positions[me][f-1]);
	    num_t  p         = progress[me];
	    for( ; (i < myNp) && (PPP[i].file==f); i++ )
	      {
		int is_masked = ((_pmasked_[p_idx(i)] & (p_bytemask(i,1))) != 0);
		switch( is_masked )
		  {
		  case 1: { double pot = 0;
		      if( sizeof_in_data == sizeof(float) )
			pot = (double)*((float*)block + PPP[i].pos);
		      else 
			pot = *((double*)block + PPP[i].pos);
		      outsnap_data[p].pot = pot; p++; }
		  default: break;
		  }
	      }

	    /* -----------------------
	     *
	     * This is the last block
	     * ==>> UPDATE progress[me]
	     * ----------------------- */
	    progress[me] += (p-progress[me]);
	  }
      }

      free(block);
      fclose(filein);


    }

  for( int f = 0; f < nmasks; f++ )
    {
      int size;
      printf("\twriting %llu particles in mask file %d, with %lu bytes per particle\n",
	     (unsigned long long)snapout_sizes[f], f, sizeof(outsnap_t));

      
      // write how many particles will be in the file
      fwrite( &snapout_sizes[f], sizeof(num_t), 1, FofGalFiles[f] );      

      // write the size of IDs
      size = sizeof(PID_t);
      fwrite( &size, sizeof(int), 1, FofGalFiles[f] );

      // write the size of floats
      size = sizeof_out_data;
      fwrite( &size, sizeof(int), 1, FofGalFiles[f] );
    }

  
  #pragma omp parallel
  {
    char *_pmasked_ = p_is_masked[me];
    
   #pragma omp for ordered
    for( int t = 0; t < Nthreads; t++ )
      {
       #pragma omp ordered
	{
	  num_t p = 0;
	  num_t w = 0;
	  for( num_t i = 0; i < myNp; i++ ) {
	    int is_masked = ((_pmasked_[p_idx(i)] & (p_bytemask(i,1))) != 0);	  
	    switch( is_masked ) {
	    case 1: {
	      for( int m = 0; m < nmasks; m++ )
		{
		  int s = ( (PPP[i].fofid == masks[m].fof_num) &&
			    ((PPP[i].gid == masks[m].ngal) ||
			     (masks[m].ngal == -1)) );
		switch( s != 0 ) {
		case 1: {fwrite( &outsnap_data[p], sizeof(outsnap_t), 1, FofGalFiles[m] ); w++;} break; }
		} p++; } break; }
	  }	
	  dprint(2,-1, "\t\t thread %d has found %llu masked particles and written %llu\n", me, p, w );
	}
      }
    
    free( p_is_masked[me] );
    free( outsnap_data );
  }
  
  for( int f = 0; f < nmasks; f++ )
    fclose( FofGalFiles[f] );
  
  free( fname );
  free( snapout_sizes );

  
  return 0;
}



