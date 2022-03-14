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
#include <float.h>
#include <omp.h>
#if defined(USE_MPI)
#include <mpi.h>
#endif

#define CPU_RTIME ({_Alignas(32) struct timespec ts; double v;		\
      v=(clock_gettime( CLOCK_REALTIME, &ts ),				\
	 (double)ts.tv_sec +						\
	 (double)ts.tv_nsec * 1e-9); v;})

#define CPU_TIME ({_Alignas(32) struct timespec ts; double v;		\
      v=(clock_gettime( CLOCK_THREAD_CPUTIME_ID, &ts ),			\
	 (double)ts.tv_sec +						\
	 (double)ts.tv_nsec * 1e-9); v;})

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
typedef unsigned long long int      num_t;

#if defined(LONG_IDS)
typedef unsigned long long PID_t;
#else
typedef unsigned int       PID_t;
#endif

#define NTYPES 6
#define ALL NTYPES

typedef struct {PID_t pid; unsigned int type, gen, fofid, gid;} particle_t;     // pid   is the gadget's ID of the particle
                                                                                // fofid is the cardinal number of the fof the particle belongs to, sequentially read from the sub_ files 
								                // gid   is the halo number the particle belong to, sequentially read from the sub_ files
								                // type  is the gadget's type of the particle
								                // gen   is the stellar generation, for type 4 particles

typedef struct {PID_t pid; unsigned int fofid, gid;} list_t;

typedef struct { PID_t pid; int gen, type; } pidtype_t;

typedef struct { int id_size, particle_t_size, nfiles, dummy; num_t Nparts, Nparts_total; } catalog_header_t;

typedef struct { int id_size, nfiles, type_is_present; num_t Nparts, Nparts_total; } list_header_t;

typedef struct { num_t fof_id; num_t TotN, Nparts[NTYPES]; } fof_table_t;


extern int     n_stars_generations, n_stars_generations_def;
extern int     id_bitshift;
extern int     verbose;
extern FILE   *details;
//extern FILE   *pdetails;

extern int     Nthreads, me;   // the number of omp threads, the thread id

//  these are used to create catalogs
extern num_t       *IDranges;
extern pidtype_t   *IDs, **all_IDs;
extern num_t        type_positions[NTYPES];
extern num_t        myNID;
extern num_t       *all_NID;
extern fof_table_t *fof_table;

// thise is used when we search particles in the catalogs
extern num_t *IDdomains[NTYPES];

extern particle_t *P[NTYPES], **P_all[NTYPES], *Pbase;
extern list_t     *List;
extern int        *Types;

extern num_t    Nparts[NTYPES+1];
extern num_t   *Nparts_all[NTYPES+1];
extern num_t    Np_all[NTYPES+1];
extern num_t    Nl, myNl;

#define  Np Np_all[ALL]
#define  myNp Nparts[ALL]


extern PID_t    id_mask;

#pragma omp threadprivate(me, P, Pbase, IDs, myNID, type_positions, Nparts, myNl, List, Types)


int   distribute_particles        ( void );
int   distribute_ids              ( void );
int   sort_thread_particles       ( num_t * );
int   sort_thread_idtype          ( void );
int   assign_type_to_subfind_particles( num_t *, num_t *, num_t *);
PID_t get_stellargenerations_mask ( int, int * );
int   get_stellargenerations      ( PID_t, PID_t, int );
int   get_subfind_data            ( char *, char *, num_t [4]);
int   get_id_data                 ( char *, char *);
int   get_catalog_data            ( char *, int *);
int   get_list_ids                ( char *, int );
int   make_fof_table              ( fof_table_t*, int );

num_t partition_P_by_pid   ( const num_t, const num_t, const PID_t);
num_t partition_IDs_by_pid ( const num_t, const num_t, const PID_t);
num_t partition_P_by_type  ( const num_t, const num_t, const unsigned int  );

int k_way_partition        ( const num_t, const num_t, const int, const int,
			     const num_t *, num_t *restrict, int );


int write_particles( FILE *, int );

#if !defined(USE_LIBC_QSORT)
int inline_qsort_particle_id_gen( void*, num_t);
int inline_qsort_idtype_id( void*, num_t);
#endif

#if !defined(USE_LIBC_BSEARCH)
extern num_t mybsearch_in_ids     (const pidtype_t *, const num_t, const PID_t, int *);
extern void* mybsearch_in_P       (const particle_t *, const num_t, const PID_t );
#endif

#if defined(DEBUG)
num_t check_partition             (const num_t, const num_t, const num_t, const PID_t, const int);
num_t check_sorting               (const num_t, const num_t, const int);
#define MASKED_ID_DBG             0
#endif




#define dprint( L, T, ... ) {if( ((L)<=verbose) && ((T)==-1 || (T)==me) ) printf(__VA_ARGS__);}
#if defined(DEBUG)
#define DPRINT( L, T, ... ) dprint(L, T, __VA_ARGS__)
#else
#define DPRINT( L, T, ... ) 
#endif

#define PRINT_TIMINGS( S, s, T ) fprintf( timings, "%40s %6.2g %s\n", S, (T), s );

#define UNUSED(x) (void)(x)

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
extern char catalog_table_name[CATALOG_NAME_SIZE];
extern char snapnum[NUM_SIZE];
#define MAX_N_LISTS 10
extern char *list_names[MAX_N_LISTS];
extern int   list_types[MAX_N_LISTS];
extern int   Nlists;

#define PPP   P[0]

#define ARG_SIZE      20
#define ARG_HELP_SIZE 120
#define ARG_DFLT      30
#define ARG_EX_SIZE   120

#define UNDEFINED           0
#define CREATE_CATALOGS     1
#define SEARCH_PARTICLES    2
#define BUILD_FOF_TABLE     4


typedef int (*action_t)(char **, int, int);

typedef struct {
  char arguments_family[ARG_HELP_SIZE*2];
  char example[ARG_EX_SIZE]; } arg_family_t;
  
typedef struct {
  char       argument[ARG_SIZE];
  char       dflt[ARG_DFLT];
  char       help[ARG_HELP_SIZE];
  char       example[ARG_EX_SIZE];
  int        family;
  action_t   action; } arg_t;

extern int   action;
extern int   help_given;
extern int   n_args_families;
extern int   n_args;
extern arg_t args[];
extern arg_family_t args_f[];
