#define INLINE
#include "cidx.h"

int     n_stars_generations, n_stars_generations_def = 4;
int     id_bitshift;
int     verbose = 0;
FILE   *details = NULL;

int     Nthreads, me;   // the number of omp threads, the thread id

ull_t  *IDranges;
_Alignas(32) pidtype_t *IDs, **all_IDs;
ull_t   type_positions[NTYPES];
ull_t   myNID;
ull_t  *all_NID;

ull_t  *IDdomains[NTYPES];

_Alignas(32) particle_t *P[NTYPES], **P_all[NTYPES], *Pbase;
_Alignas(32) list_t *List;

ull_t   Nparts[NTYPES] = {0};
ull_t   *Nparts_all[NTYPES];
ull_t   Np_all[NTYPES] = {0};
ull_t   Np = 0, myNp = 0;
ull_t   Nl, myNl;

PID_t   id_mask = 0;

char working_dir_def[] = "./";
char working_dir_subf_base_def[] = "groups_";
char working_dir_snap_base_def[] = "snapdir_";
char snap_base_def[] = "snap_";
char subf_base_def[] = "sub_";

char working_dir[DIR_SIZE];
char working_dir_snap_base[DIR_SIZE];
char working_dir_subf_base[DIR_SIZE];
char working_dir_snap[DIR_SIZE+NUM_SIZE];
char working_dir_subf[DIR_SIZE+NUM_SIZE];
char snap_base[NAME_SIZE];
char subf_base[NAME_SIZE];
char catalog_base[NAME_SIZE] = "catalog";
char snap_name[NAME_SIZE+NUM_SIZE];
char subf_name[NAME_SIZE+NUM_SIZE];
char catalog_name[CATALOG_NAME_SIZE];
char list_name[CATALOG_NAME_SIZE];
char snapnum[NUM_SIZE];
