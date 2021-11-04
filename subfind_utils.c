
#include "cidx.h"


#define XCHG_P(P1, P2) {particle_t temp1 = *P1; particle_t temp2 = *P2; *P1 = temp2; *P2 = temp1;}
#define XCHG_ID(P1, P2) {pidtype_t temp1 = *P1; pidtype_t temp2 = *P2; *P1 = temp2; *P2 = temp1;}


int mask_ids_and_find_ranges           ( const PID_t, const int, ull_t *, int );




int mask_ids_and_find_ranges( const PID_t id_mask, const int bitshift, ull_t *ranges, int mode )
/*
 * this routine is getting through the particles masking their
 * ID by the stellar generation mask, so that the id of each particle
 * is brought back to the original ID and its generation is stored
 * in the .gen field.
 * While doing that, the min and max id are found so that it is possible
 * to span the ids among the threads
 *
 * the routine is called by the main process from a serial region and
 * opens a multithread region.
 *
 * This routine is called 
 * (i)  to run over particles found in subfind files,
 *      from the routine distribute_particles()
 * (ii) to run over the ids read from the snapshot,
 *      from the routine distribute_ids()
 *
 * id_mask  : that is the stellar-generation-mask to be used
 * bitshift : that is the needed bitshift to recover the generation
 * ranges   : the 2-entries array where the min and max id are placed
 * mode     : 0 means "run over subfind particles", 
 *            1 means "run over ids from snap"
 *
 * NOTE: for mode 0 the loop is over PPP, which stands for P[0]. At this stage,
 *       subfind data have not yet been sorted by type, so everything is still in
 *       P[0]
 */
{
  PID_t min_bound = (PID_t)-1;
  PID_t max_bound = 0;
  int id_size_alert   = 0;
  
  DPRINT(1,0, "id_mask is %llu\n", (ull_t)id_mask);
 #pragma omp parallel reduction(min:min_bound) reduction(max:max_bound)
  {
    int   bs = sizeof(PID_t)*8 - bitshift;
    ull_t N;

    switch(mode)
      {
      case 0: N = myNp; break;
      case 2: N = myNID; break;
      default: N = 0;
      }
    
    for ( ull_t i = 0; i < N; i++ )
      {
	PID_t ID, masked_id;
	switch( mode )
	  {
	  case 0: ID = PPP[i].pid; break;
	  case 2: ID = IDs[i].pid; break;
	  default: ID = 0;
	  }
	//id_size_alert += ( ID > id_mask );
	masked_id = ID & id_mask;
	int  gen       = ID >> bs;
	min_bound      = ( min_bound > masked_id ? masked_id : min_bound );
	max_bound      = ( max_bound < masked_id ? masked_id : max_bound );
	switch( mode )
	  {
	  case 0: PPP[i].pid = masked_id; PPP[i].gen = gen; break;
	  case 2: IDs[i].pid = masked_id; IDs[i].gen = gen; break;
	  }
      }
  }

  /* if ( id_size_alert ) */
  /*   printf("* ---> [ err ] it seems that you have less than %d stellar generations\n", n_stars_generations); */
  DPRINT(1, 0, "min and max are %llu %llu\n", (ull_t)min_bound, (ull_t)max_bound);
  PID_t delta = (max_bound - min_bound) / Nthreads;
  
  for ( int i = 1; i < Nthreads; i++ )
    ranges[i-1] = min_bound + delta*i;
  ranges[Nthreads-1] = max_bound;
  
  return id_size_alert;
}




inline ull_t partition_P_by_pid( const ull_t start, const ull_t stop, const PID_t V)
/*
 * This function separates the elements of an array which are smaller and larger than V.
 * ALl the elements that are smaller or equal than V are moved in the first part of the
 * array, and all the elements that are strictly larger than V are moved in the second
 * part of the array.
 * The function returns the index of the first element strictly larger than V.
 *
 * INPUT:
 * start is the first element to be considered in the array
 * stop is the second element to be considered in the array
 * V is the value to be compared with
 *
 * OUTPUT:
 * the index of the first element strictly larger than V
 * note: if the returned value is larger than stop, it means that all the elements
 * resulted to be smaller or equal than V
 *
 * NOTE: the run is over PPP, i.e. P[0], since at this point the subfind particles
 * hae not yet been sorted by type
 *
 */
{
  if( start >= stop )
    return start + (start < myNp ? PPP[start].pid < V : 0);

  ull_t j = start;
  register PID_t myV = V;

  while( (j < stop) && (PPP[j].pid <= myV) )
    j++;

  if ( j == stop )
    return stop;

  ull_t      stop_1 = stop-1;
  PID_t      A;
  register particle_t *ptr;
  
  A   = PPP[j].pid;
  ptr = &PPP[j];
  
  do
    {      
      if(A > myV)
	{
	  ++j;
	  __builtin_prefetch(&PPP[j+1].pid, 0, 3);
	  A = PPP[j].pid;	  
	}
      else
	{
	  __builtin_prefetch(&PPP[j+1].pid, 0, 3);
	  XCHG_P(ptr, &PPP[j]);	  
	  
	  ptr++;
	  A = PPP[++j].pid;
	}
    }
  while( j < stop_1 );

  if( PPP[stop_1].pid <= myV ) {
    XCHG_P( ptr, &PPP[stop_1]); ptr++; }

  return (ull_t)(ptr - PPP);
}


inline ull_t partition_IDs_by_pid( const ull_t start, const ull_t stop, const PID_t V)
/*
 * This function separates the elements of an array which are smaller and larger than V.
 * ALl the elements that are smaller or equal than V are moved in the first part of the
 * array, and all the elements that are strictly larger than V are moved in the second
 * part of the array.
 * The function returns the index of the first element strictly larger than V.
 *
 * INPUT:
 * start is the first element to be considered in the array
 * stop is the second element to be considered in the array
 * V is the value to be compared with
 *
 * OUTPUT:
 * the index of the first element strictly larger than V
 * note: if the returned value is larger than stop, it means that all the elements
 * resulted to be smaller or equal than V
 *
 */
{
  if( start >= stop )
    return start + (start < myNID ? IDs[start].pid < V : 0);
  
  ull_t j = start;
  register PID_t myV = V;
  
  while( (j < stop) && (IDs[j].pid <= myV) )
    j++;

  if ( j == stop )
    return stop;
  
  ull_t      stop_1 = stop-1;
  PID_t      A;
  register pidtype_t *ptr;

  A   = IDs[j].pid;
  ptr = &IDs[j];
  
  do
    {      
      if(A > myV)
	{
	  ++j;
	  __builtin_prefetch(&IDs[j+1].pid, 0, 3);
	  A = IDs[j].pid;
	}
      else
	{
	  __builtin_prefetch(&IDs[j+1].pid, 0, 3);
	  XCHG_ID(ptr, &IDs[j]);
	  
	  ptr++;
	  A = IDs[++j].pid;
	}
    }
  while( j < stop_1 );

  if( IDs[stop_1].pid <= myV ) {
    XCHG_ID( ptr, &IDs[stop_1]); ptr++; }


  return (ull_t)(ptr - IDs);
}


inline ull_t partition_P_by_type( const ull_t start, const ull_t stop, const unsigned int V)
/*
 * This routine does the partitioning of the particles of each thread by a given type.
 * At the end, the particles in the PPP (i.e. P[0]) array are divided in 2 bunches: the first
 * one contains all the particles of a given type while the second one contains all the
 * other particles.
 *
 * start : the starting point in the PP array
 * end   : the end point in the PPP array
 * v     : the value to partition by
 */
{
  if( start >= stop )
    return start + (start < myNp ? PPP[start].type < V : 0);
  
  ull_t j = start;
  register unsigned int myV = V;

  while( (j < stop) && (PPP[j].type <= myV) )
    j++;
  
  if ( j == stop )
    return stop;

  ull_t        stop_1 = stop - 1;
  unsigned int A;
  register particle_t *ptr;

  A   = PPP[j].type;
  ptr = &PPP[j];
  
  do
    {      
      if(A > myV)
	{
	  ++j;
	  __builtin_prefetch(&PPP[j+1].type, 0, 3);
	  A = PPP[j].type;	  
	}
      else
	{
	  __builtin_prefetch(&PPP[j+1].type, 0, 3);
	  XCHG_P(ptr, &PPP[j]);
	  
	  ptr++;
	  A = PPP[++j].type;
	}
    }
  while( j < stop_1 );

  if( PPP[stop_1].type <= myV ) {
    XCHG_P( ptr, &PPP[stop_1]); ptr++; }


  return (ull_t)(ptr - PPP);
}


int k_way_partition( const ull_t start, const ull_t stop, const int start_range, const int stop_range,
		     const ull_t *ranges, ull_t *restrict positions, int mode)
/*
 * This routine does the partitioning of the particles by multiple values.
 * At the end, the particles in the PPP (i.e. P[0]) array are divided in bunches of
 * homogeneous value.
 *
 * the indexes at which each partition begins are stored in *positions
 *
 * start        : the PPP position at which the multipartitioning should start
 * end          : the PPP last position to be considered
 * start_range  : the first partitioning value to consider (in *ranges)
 * stop_range   : the last partitioning value to consider (in *ranges)
 * ranges       : the partitioning values
 * positions    : an array as long as ranges that stores the beginning position
 *                of each partition
 * mode         : 0 means partition PPP by id
 *                1 means partition PPP by type
 *                2 measn partitions ids by id
 *
 * note: this function is recursive
 */
{
  int Nhalf = start_range + (stop_range - start_range)/2;
  
  PAPI_START;
  switch(mode)
    {
    case 0: positions[Nhalf] = partition_P_by_pid(start, stop, ranges[Nhalf]); break;
    case 1: positions[Nhalf] = partition_P_by_type(start, stop, ranges[Nhalf]); break;
    case 2: positions[Nhalf] = partition_IDs_by_pid(start, stop, ranges[Nhalf]); break;
    }
  PAPI_STOP;

 #if defined(DEBUG)
  ull_t fails = check_partition(start, stop, positions[Nhalf], ranges[Nhalf], mode);
  if( fails )
    printf("--------- thread %d has %llu failures with mode %d in partition from %llu to %llu with pivot %llu\n",
	   me, fails, mode, start, stop, positions[Nhalf]);
 #endif

  if (start_range < Nhalf) {
    if (positions[Nhalf] > start) 
      k_way_partition(start, positions[Nhalf], start_range, Nhalf-1, ranges, positions, mode);
    else for( int j = Nhalf-1; j>=start_range;  j-- )
	   positions[j] = start; }
    
  if (stop_range > Nhalf) {
    if (positions[Nhalf] < stop )
      k_way_partition(positions[Nhalf], stop, Nhalf+1, stop_range, ranges, positions, mode);
    else for( int j = Nhalf+1; j <= stop_range; j++ )
	   positions[j] = stop; }

  return 0;
}


#if defined(DEBUG)
ull_t check_partition(const ull_t start, const ull_t stop, const ull_t pivot, const PID_t V, const int mode)
/*
 * when DEBUG is set, ths routine checks that a partitioning has been done correctly
 *
 */
{
  ull_t ii;
  ull_t tmp = 0;
  
  switch( mode )
    {
    case 0:
      {	
	for(ii = start; ii < pivot; ii++)     // check the first part
	  tmp += (PPP[ii].pid > V);

	if( ii < stop - 1)                // check the second part
	  for(; ii < stop; ii++)
	    tmp += (PPP[ii].pid < V);
      } break;
    case 1:
      {
	unsigned int iV = (int)V;
	for(ii = start; ii < pivot; ii++){     // check the first part
	  ull_t f = (PPP[ii].type > iV);
	  if(f)
	    tmp ++; }
	
	if( ii < stop - 1)                // check the second part
	  for(; ii < stop; ii++) {
	    ull_t f = (PPP[ii].type < iV);
	    if(f)
	      tmp++; }
      } break;      
    case 2:
      {
	for(ii = start; ii < pivot; ii++)     // check the first part
	  tmp += (IDs[ii].pid > V);
	
	if( ii < stop - 1)                // check the second part
	  for(; ii < stop; ii++)
	    tmp += (IDs[ii].pid < V);
      } break;
    }

  return tmp;
}



ull_t check_sorting(const ull_t start, const ull_t stop, const int mode)
/*
 * when DEBUG is set, ths routine checks that a sorting has been done correctly
 *
 */
{
  ull_t ii;
  ull_t tmp = 0;
  
  switch( mode )
    {
    case 0:
      {	
	for(ii = start+1; ii < stop; ii++) {
	  int eq = (PPP[ii].pid == PPP[ii-1].pid);
	  tmp += (PPP[ii].pid < PPP[ii-1].pid) + (eq ? (PPP[ii].gen < PPP[ii-1].gen) :0); }	  
      } break;
    case 2:
      {
	for(ii = start+1; ii < stop; ii++)
	  tmp += (IDs[ii].pid < IDs[ii-1].pid);	
      } break;
    }

  return tmp;
}

#endif



int distribute_particles( void )
{
  int         fails = 0;
  ull_t       matrix[Nthreads][Nthreads];
  ull_t      *Pranges;
    
  memset( matrix, 0, sizeof(ull_t) * Nthreads * Nthreads );
  
  /* ---------------------------------------------------------
   * 
      DATA FOR IN-PLACE REDISTRIBUTION ALGORITHM

  #define NBUFFERS    3
  #define BUFFER_SIZE 0
  #define BUFFER_HEAD 1
  #define SIGNALS     2
  int flags[NBUFFERS][Nthreads];
  memset( flags, 0, NBUFFERS * sizeof(int) * Nthreads );
  
      END of DATA FOR IN-PLACE DISTRIBUTION AlGORITHM
  * -----------------------------------------------------------
  */

  // mask the ids and find the ranges
  //
  Pranges = calloc( Nthreads, sizeof(ull_t) );
  mask_ids_and_find_ranges( id_mask, id_bitshift, Pranges, 0 );

  if( Nthreads == 1 )
    {
      P_all[0][0] = P[0];
      return 0;
    }

  ull_t min_p = myNp, max_p = 0;
 #pragma omp parallel reduction(min:min_p) reduction(max:max_p)
  {
    ull_t positions[Nthreads];
    
    // partition the particles in each thread
    // so that at the end they are divided in bunches
    // ready to be copied into other threads
    //
    k_way_partition(0, myNp, 0, Nthreads-2, Pranges, positions, 0);
    positions[Nthreads-1] = myNp;
    
    matrix[0][me] = positions[0];
    for ( int i = 1; i < Nthreads; i++ )
      matrix[i][me] = positions[i] - positions[i-1];
    //matrix[Nthreads-1][me] = myNp - positions[Nthreads-2];

   #if defined(DEBUG)
    ull_t _sum_ = 0;
    for ( int i = 0; i < Nthreads; i++ ) {
      DPRINT(2,-1, "[ particles ] th %d will send %llu particles to %d\n", me, (ull_t)matrix[i][me], i); _sum_ += matrix[i][me];}
    DPRINT(2,-1, "[ particles ] th %d : total p is %llu\n", me, _sum_);
   #endif
    
   #pragma omp barrier  // ensures that matrix is complete

    // check that I do have enough memory
    //
    ull_t NN   = 0;
    for( int i = 0; i < Nthreads; i++ )
      {
	ull_t amount = matrix[me][i];
	matrix[me][i] = NN;
	NN += amount;
      }
    DPRINT(2,-1, "[ particles ] th %d : old total is %llu, new total is %llu\n", me, myNp, NN);

    // allocate space for the final amount of particles
    //
    //particle_t *my_P = (particle_t*)aligned_alloc( ALIGN, (NN+1)*sizeof(particle_t) );
    particle_t *my_P = (particle_t*)calloc( (NN+1), sizeof(particle_t) );
    P_all[0][me] = my_P;


   #pragma omp barrier  // ensures that everybody did allocate the memory region
    
    if( my_P != NULL )
      {

	// each thread here is copying its initial particles
	// that belong to another thread into the destination thread
	//
	particle_t *source = P[0];

	for ( int i = 0; i < Nthreads; i++ )
	  {
	    
	    
	    particle_t *target = P_all[0][i] + matrix[i][me];
	    ull_t amount;
	    if ( i == 0 )
	      amount = positions[0];
	    else
	      amount = positions[i] - positions[i-1];

	    DPRINT(2,-1, "[ D4 ] th %d -> th %d %llu\n", me, i, amount );
	    for( ull_t j = 0; j < amount; j++ )
	      {
		*target++ = *source++;
	      }

	  }

	// update relevant variables
	//
	myNp = NN;
	free(P[0]);
	P[0] = my_P;
      }
    else
      #pragma omp atomic update
      fails++;

    fprintf( details, "thread %d got %llu subfind particles after re-distribution\n", me, myNp);

    min_p = (min_p > myNp ? myNp : min_p);
    max_p = (max_p < myNp ? myNp : max_p);
  }

  fprintf( details, "---- maximum imbalance in subfind particles distribution is %3.1f%%\n\n",
	   (double)(max_p - min_p)*100/max_p);
    
    /* ---------------------------------------------------------
     * 
       IN-PLACE REDISTRIBUTION ALGORITHM

    // find the maximum bunch that I will receive
    // and allocate a buffer memory region
    //
    unsigned int largest_bunch = 0;
    for( int i = 0; i < Nthreads; i++ )
      largest_bunch = ( largest_bunch < matrix[me][i] ? matrix[me][i] : largest_bunch );
    data_t *buffer = calloc( largest_bunch, sizeof(data_t) );
    buffer_size[me] = largest_bunch;
    
    // exchange data with other threads
    //
    

    typedef struct {
      char          *ptr;
      long long int  size;
    } memslot_t;
    memslot_t memslots[Nthreads];
    memset( memslots, 0, sizeo(memslot_t)*Nthreads);

    int Ptasks;                                                                                                                                                                                                  
    for( Ptasks = 0; (1 << Ptasks) < Ntasks; Ptasks++ )                                                                                                                                                          
      ;  
    Ptasks = (1 << Ptasks);

    

    for ( int ngrp = 1; ngrp < Ptasks; ngrp++ )
      {
	int target = me ^ ngrp;
	if( target < Ntasks )
	  {
	    int copied_data = 0;
	    int available_in_target_buffer = largest_bunch[target] - buffer_head[target];
	    
	    
	  }
      }
    
    free(buffer);

   }
       END OF IN-PLACE REDISTRIBUTION ALGORITHM
    * ----------------------------------------------------------
    */
    

  return fails;
}



int distribute_ids( void )
{

  int   fails = 0;
  ull_t matrix[Nthreads][Nthreads];  
  memset( matrix, 0, sizeof(ull_t) * Nthreads * Nthreads );

  all_IDs = (pidtype_t**)calloc( Nthreads, sizeof(pidtype_t*) );
  all_NID = (ull_t*)calloc( Nthreads, sizeof(ull_t) );
  
  // mask the ids and find the ranges
  //
  IDranges = (ull_t*)calloc( Nthreads, sizeof(ull_t) );  
  mask_ids_and_find_ranges( id_mask, id_bitshift, IDranges, 2 );

  if ( Nthreads == 1 )
    {
      all_IDs[0] = IDs;
      all_NID[0] = myNID;
      return 0;
    }

  ull_t min_p = myNID, max_p = 0;
 #pragma omp parallel reduction(min:min_p) reduction(max:max_p)
  {
    ull_t NN = 0;
    ull_t positions[Nthreads];
    memset(positions, 0, sizeof(ull_t)*Nthreads);
    
    // partition the particles in each thread
    // so that at the end they are divided in bunches
    // ready to be coped into other threads
    //
    k_way_partition(0, myNID, 0, Nthreads-2, IDranges, positions, 2);
    positions[Nthreads-1] = myNID;

    matrix[0][me] = positions[0];
    for ( int i = 1; i < Nthreads; i++ )
      matrix[i][me] = positions[i] - positions[i-1];


   #if defined(DEBUG)
   #pragma omp barrier
   #pragma omp master
    {
      fflush(stdout);
      dprint(2, 0, "---------------\nto-be-sent particles\n\n");
      dprint(2, 0, "%7s", "thread");
      for( int i = 0; i < Nthreads; i++ )
	dprint(2, 0, "    %2d   ", i);
      dprint(2, 0, "\n");
      for( int i = 0; i < Nthreads; i++ ) {
	ull_t _sum_ = 0;
	dprint(2, 0, "   %2d  ", i);
	for( int j = 0; j< Nthreads; j++ ) {
	  _sum_ += matrix[i][j]; dprint(2, 0, "%9llu", matrix[i][j]);}
	dprint(2, 0, "\t%12llu\n", _sum_);}
    }
    ull_t _sum_ = 0;
    for ( int i = 0; i < Nthreads; i++ ) {
      _sum_ += matrix[i][me];}
   #endif


    // calculate how many particles I will have
    //
   #pragma omp barrier
    for( int i = 0; i < Nthreads; i++ )
      {
	ull_t amount = matrix[me][i];
	matrix[me][i] = NN;
	NN += amount;
      }

   #if defined(DEBUG)
   #pragma omp barrier
   #pragma omp master
    {
      fflush(stdout);
      dprint(2, 0, "---------------\nIDss offset\n\n");
      dprint(2, 0, "%7s", "thread");
      for( int i = 0; i < Nthreads; i++ )
	dprint(2, 0, "    %2d   ", i);
      dprint(2, 0, "\n");
      for( int i = 0; i < Nthreads; i++ ) {
	dprint(2, 0, "   %2d  ", i);
	for( int j = 0; j< Nthreads; j++ )
	  dprint(2, 0, "%9llu", matrix[i][j]);
	dprint(2, 0, "\n");}
    }
   #endif

    
    // allocate space for the final amount of particles
    //
    pidtype_t *my_ID = (pidtype_t*)calloc( NN, sizeof(pidtype_t) );
    all_IDs[me] = my_ID;
    
    if( my_ID != NULL )
      {

	
       #pragma omp barrier

	// each thread here is copying its initial particles
	// that belong to another thread into the destination thread
	//
	pidtype_t *source = IDs;
	for ( int i = 0; i < Nthreads; i++ )
	  {
	    pidtype_t *target = all_IDs[i] + matrix[i][me];
	    ull_t amount;
	    if ( i == 0 )
	      amount = positions[0];
	    else
	      amount = positions[i] - positions[i-1];

	    DPRINT(2,-1, "[ IDs ] th %d -> th %d %llu ** %llu %llu\n", me, i, amount, (i>0?positions[i-1]:0), positions[i] );
	    for ( ull_t j = 0; j < amount; j++ )
	      *target++ = *source++;
	  }

       #pragma omp barrier
	
	// update relevant variables
	//
	myNID = NN;
	all_NID[me] = NN;
	free(IDs);
	IDs = my_ID;
      }
    else
      {
       #pragma omp atomic update
	fails++;
      }

    fprintf( details, "thread %d got %llu snapshot particles after re-distribution\n", me, myNID);

    min_p = (min_p > myNID ? myNID : min_p);
    max_p = (max_p < myNID ? myNID : max_p);

  }

  fprintf( details, "---- maximum imbalance in snapshot particles distribution is %3.1f%%\n\n",
	   (double)(max_p - min_p)*100/max_p);

  return fails;
}


int compare_idsgen_in_particle_t( const void *A, const void *B )
{
  PID_t idA = ((particle_t*)A)->pid;
  PID_t idB = ((particle_t*)B)->pid;
  int   genA = ((particle_t*)A)->gen;
  int   genB = ((particle_t*)B)->gen;
  
  int  res = ( idA > idB ) - ( idA < idB );

  res += (res == 0 ? ( genA - genB ) : 0);
  
  return res;
}

int compare_ids_in_pidtype_t( const void *A, const void *B )
{
  PID_t idA = ((pidtype_t*)A)->pid;
  PID_t idB = ((pidtype_t*)B)->pid;
  int  res = ( idA > idB ) - ( idA < idB );

  return res;
}

int compare_pid_with_particle_t( const void *A, const void *B )
{
  PID_t idA = *(PID_t*)A;
  PID_t idB = ((particle_t*)B)->pid;
  int  res = ( idA > idB ) - ( idA < idB );

  return res;
}


int sort_thread_particles( ull_t *positions )
{

  ull_t start = 0;
  for ( int t = 0; t < NTYPES; t++ )
    {
      ull_t stop = positions[t];
     #if !defined(USE_LIBC_QSORT)
      inline_qsort_particle_id_gen( (void*)(PPP+start), (ull_t)(stop-start));
     #else  
      qsort( (void*)(PPP+start), (size_t)(stop-start), sizeof(particle_t), compare_idsgen_in_particle_t);
     #endif
      start = stop;
    }
  
  return 0;
}

int sort_thread_idtype( void )
{

 #if !defined(USE_LIBC_QSORT)
  inline_qsort_idtype_id( (void*)IDs, myNID);
 #else
  qsort( IDs, (size_t)myNID, sizeof(pidtype_t), compare_ids_in_pidtype_t);
 #endif

  return 0;
}


int assign_type_to_subfind_particles( ull_t *o_of_r, ull_t *failures)
{
  /* FILE *FILE; */
  /* char name[100]; */
  /* sprintf(name, "fails.%d", me); */
  /* FILE = fopen( name, "w"); */
  
  ull_t out_of_range = 0;
  ull_t fails = 0;
  
  for( ull_t j = 0; j < myNp; j++ )
    {
      
      int target_thread = 0;

      while( (target_thread < Nthreads) && (PPP[j].pid > IDranges[target_thread]) )
	target_thread++;

      if( target_thread < Nthreads )
	{
	  pidtype_t *res;
	  
	 #if !defined(USE_LIBC_BSEARCH)
	  int err;
	  ull_t pos = mybsearch_in_ids(all_IDs[target_thread], all_NID[target_thread], PPP[j].pid, &err);
	  res = ( err == 0? &all_IDs[target_thread][pos] : NULL );
	  
	 #else
	  res = bsearch( (void*)&(PPP[j].pid), all_IDs[target_thread], all_NID[target_thread],
			 sizeof(pidtype_t), compare_pid_with_particle_t );
	 #endif
	  
	  if( res != NULL )
	    PPP[j].type = res->type;
	  else {
	    fails++;
	    //DPRINT(1, -1, "* particle %llu pid %llu gen %d has no counterpart in thread %d\n", j, (ull_t)PPP[j].pid, PPP[j].gen, target_thread );
	  }
	}
      else
	out_of_range++;
    }

  *o_of_r = out_of_range;
  *failures = fails;
  //fclose(FILE);
  return (out_of_range > 0) + (fails > 0);
}
