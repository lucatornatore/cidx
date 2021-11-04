#define _FILE_OFFSET_BITS 64

#if defined(__STDC__)
#  if (__STDC_VERSION__ >= 199901L)
#     define _XOPEN_SOURCE 700
#  endif
#endif


#ifndef INLINE
# if __GNUC__ && !__GNUC_STDC_INLINE__
#  define INLINE extern 
# else
#  define INLINE inline
# endif
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <omp.h>
#if defined(USE_MPI)
#include <mpi.h>
#endif

#define ALIGN 32

/* --------------------------
 *   define PAPI calls
 * -------------------------- */
#ifdef PAPI_MEASURE
// define PAPI calls
//
#define PAPI_NEVENTS 4
#define PAPI_START {int retval = PAPI_start_counters(papi_events, PAPI_NEVENTS); PCHECK(retval);}
#define PAPI_STOP  {int retval = PAPI_stop_counters(papi_values, PAPI_NEVENTS); PCHECK(retval);}
#define PCHECK(e) \
  if (e!=PAPI_OK)							\
    {printf("Problem in papi call, line %d\n",__LINE__); return 1;}

#else
//define empty calls
//

#define PAPI_START
#define PAPI_STOP 

#endif

typedef unsigned int       ul_t;
typedef unsigned long long ull_t;
typedef long long int      ll_t;

#if defined(LONG_IDS)
typedef unsigned long long PID_t;
#else
typedef unsigned int       PID_t;
#endif

typedef struct {PID_t pid; unsigned int type, gen, fofid, gid;} particle_t;     // pid   is the gadget's ID of the particle
                                                                                // fofid is the cardinal number of the fof the particle belongs to, sequentially read from the sub_ files 
								                // gid   is the halo number the particle belong to, sequentially read from the sub_ files
								                // type  is the gadget's type of the particle
								                // gen   is the stellar generation, for type 4 particles

typedef struct {PID_t pid; int type; unsigned int gen, fofid, gid;} list_t;

typedef struct { PID_t pid; int gen, type; } pidtype_t;

typedef struct { int id_size, particle_t_size, nfiles; ull_t Nparts, Nparts_total; } catalog_header_t;

typedef struct { int id_size, nfiles, type_is_present; ull_t Nparts, Nparts_total; } list_header_t;

#define NTYPES 6

extern int     n_stars_generations, n_stars_generations_def;
extern int     id_bitshift;
extern int     verbose;
extern FILE   *details;

extern int     Nthreads, me;   // the number of omp threads, the thread id

//  these are used to create catalogs
extern ull_t     *IDranges;
extern pidtype_t *IDs, **all_IDs;
extern ull_t      type_positions[NTYPES];
extern ull_t      myNID;
extern ull_t      *all_NID;

// thise is used when we search particles in the catalogs
extern ull_t *IDdomains[NTYPES];

extern particle_t *P[NTYPES], **P_all[NTYPES], *Pbase;
extern list_t     *List;

extern ull_t    Nparts[NTYPES];
extern ull_t    *Nparts_all[NTYPES];
extern ull_t    Np_all[NTYPES];
extern ull_t    Np, myNp;
extern ull_t    Nl, myNl;

extern PID_t    id_mask;

#pragma omp threadprivate(me, P, Pbase, myNp, IDs, myNID, type_positions, Nparts, myNl, List)


int   distribute_particles        ( void );
int   distribute_ids              ( void );
int   sort_thread_particles       ( ull_t * );
int   sort_thread_idtype          ( void );
int   assign_type_to_subfind_particles( ull_t *, ull_t * );
PID_t get_stellargenerations_mask ( int, int * );
int   get_stellargenerations      ( PID_t, PID_t, int );
int   get_subfind_data            ( char *, char *);
int   get_id_data                 ( char *, char *);
int   get_catalog_data            ( char *, int *);
int   get_list_ids                ( char * );

ull_t partition_P_by_pid   ( const ull_t, const ull_t, const PID_t);
ull_t partition_IDs_by_pid ( const ull_t, const ull_t, const PID_t);
ull_t partition_P_by_type  ( const ull_t, const ull_t, const unsigned int  );

int k_way_partition        ( const ull_t, const ull_t, const int, const int,
			     const ull_t *, ull_t *restrict, int );


int write_particles( FILE *, int );

#if !defined(USE_LIBC_QSORT)
int inline_qsort_particle_id_gen( void*, ull_t);
int inline_qsort_idtype_id( void*, ull_t);
#endif

#if !defined(USE_LIBC_BSEARCH)
extern ull_t mybsearch_in_ids     (const pidtype_t *, const ull_t, const PID_t, int *);
extern ull_t mybsearch_in_P       (const particle_t *, const ull_t, const PID_t, int *);
#endif

#if defined(DEBUG)
ull_t check_partition                  (const ull_t, const ull_t, const ull_t, const PID_t, const int);
ull_t check_sorting                    (const ull_t, const ull_t, const int);
#endif



#define dprint( L, T, ... ) if( ((L)<=verbose) && ((T)==-1 || (T)==me) ) printf(__VA_ARGS__);
#if defined(DEBUG)
#define DPRINT( L, T, ... ) dprint(L, T, __VA_ARGS__)
#else
#define DPRINT( L, T, ... ) 
#endif



extern char working_dir_def[];
extern char working_dir_subf_base_def[];
extern char working_dir_snap_base_def[];
extern char snap_base_def[];
extern char subf_base_def[];
extern char catalog_base_def[];

#define DIR_SIZE 99
#define NAME_SIZE 49
#define NUM_SIZE 6
#define CATALOG_NAME_SIZE (DIR_SIZE+NAME_SIZE+2*NUM_SIZE+20)

extern char working_dir[DIR_SIZE];
extern char working_dir_snap_base[DIR_SIZE];
extern char working_dir_subf_base[DIR_SIZE];
extern char working_dir_snap[DIR_SIZE+NUM_SIZE];
extern char working_dir_subf[DIR_SIZE+NUM_SIZE];
extern char snap_base[NAME_SIZE];
extern char subf_base[NAME_SIZE];
extern char catalog_base[NAME_SIZE];
extern char snap_name[NAME_SIZE+NUM_SIZE];
extern char subf_name[NAME_SIZE+NUM_SIZE];
extern char catalog_name[CATALOG_NAME_SIZE];
extern char list_name[CATALOG_NAME_SIZE];
extern char snapnum[NUM_SIZE];


#define PPP   P[0]
