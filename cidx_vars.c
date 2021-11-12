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

_Alignas(32) particle_t *P[NTYPES]={NULL}, **P_all[NTYPES]={NULL}, *Pbase = NULL;
_Alignas(32) list_t *List = NULL;
_Alignas(32) int    *Types = NULL;

ull_t   Nparts[NTYPES] = {0};
ull_t   *Nparts_all[NTYPES] = {NULL};
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

char **list_names = NULL;
int   *list_types = NULL;
int    Nlists = 0;
int    action = UNDEFINED;
int    help_given = 0;

// note: in the following functions, "mode" should be passed
// as argc from the main, so that is is > 0 and also conveys
// infos on the number of arguments which is needed for
// options that scroll through next arguments
//

int action_create_catalogs( char **argv, int n, int mode )
{
  UNUSED(argv);
  UNUSED(n);
  UNUSED(mode);
  if( mode > 0 )
    action |= CREATE_CATALOGS;
  return n;
}

int action_search_particles( char **argv, int n, int mode )
{
  UNUSED(argv);
  UNUSED(n);
  UNUSED(mode);
  if( mode > 0 )
    action |= SEARCH_PARTICLES;
  return n;
}

int set_snapf( char **argv, int n, int mode)
{
  if( mode > 0 )
    snprintf( snap_base, NAME_SIZE, "%s", *(argv + ++n) );
  else
    sprintf( args[n].dflt, "%s+snapnum", snap_base_def );
  return n;
}

int set_subf( char **argv, int n, int mode )
{
  if( mode > 0 )
    snprintf( subf_base, NAME_SIZE, "%s", *(argv + ++n) );
  else
    sprintf( args[n].dflt, "%s+snapnum", subf_base_def );
  return n;
}

int set_stellargenerations( char **argv, int n, int mode)
{
  if( mode > 0 )
    n_stars_generations = atoi( *(argv + ++n) );
  else
    sprintf( args[n].dflt, "%d", n_stars_generations_def);
  return n;
}

int set_wdir( char **argv, int n, int mode )
{
  if( mode > 0 )
    snprintf( working_dir, DIR_SIZE, "%s", *(argv + ++n) );
  else
    sprintf( args[n].dflt, "%s", working_dir_def );
  return n;
}

int set_snapdir( char **argv, int n, int mode )
{
  if( mode > 0 )
    snprintf( working_dir_snap_base, DIR_SIZE, "%s", *(argv + ++n) );
  else
    sprintf( args[n].dflt, "%s+snapnum", working_dir_snap_base_def);
  return n;
}

int set_subfdir( char **argv, int n, int mode )
{
  if( mode > 0 )
    snprintf( working_dir_subf_base, DIR_SIZE, "%s", *(argv + ++n) );
  else
    sprintf( args[n].dflt, "%s+snapnum", working_dir_subf_base_def);
  return n;
}

int set_num( char **argv, int n, int mode )
{
  if( mode > 0 )
    snprintf( snapnum, NUM_SIZE, "%s", *(argv + ++n) );
  else
    args[n].dflt[0] = '\0';
  return n;
}

int set_list( char **argv, int n, int mode )
{
  if( mode > 0 )
    {
      Nlists++;
      list_names = (char **)realloc( list_names, sizeof(char *)*Nlists );
      list_types = (int *)realloc( list_types, sizeof(int)*Nlists );
      int len = strlen(*(argv + ++n)) + 1;
      list_names[Nlists-1] = (char*)malloc( len );
      snprintf( list_names[Nlists-1], len, "%s", *(argv + n ));
      if ( (n+1 < mode ) && (**(argv+n+1) != '-') )
	list_types[Nlists-1] = atoi( *(argv + ++n) );
      else
	list_types[Nlists-1] = -1;
    }
  else if ( mode < 0 )
    {
      free( list_types );
      for ( int l = Nlists-1; l >= 0; l-- )
	free( list_names[l] );
      free( list_names );
    }
  
  return n;
}

int set_verbose( char **argv, int n, int mode)  
{
  if( mode > 0 )
    verbose = atoi( *(argv + ++n) );
  else
    sprintf( args[n].dflt, "%d", verbose );
  
  return n;
}

int ask( char **argv, int n, int mode)
{
  UNUSED(argv);
  UNUSED(n);
  UNUSED(mode);

  if( !mode )
    return 0;
  
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

  return n;
}

int help( char **, int, int);

arg_family_t args_f[] = { {"read subfind files and create the inverse catalogs", ""},
			  {"search a list of IDs from snapshot in inverse catalogs", ""},
			  {"get infos on me", ""}};

arg_t args[] = {
		{"-c", "", "", "", 0, action_create_catalogs },
		{"-cc", "", "", "", 0, action_create_catalogs },
		{"-create", "", "", "", 0, action_create_catalogs },
		{"-c", "", "", "", 0, action_create_catalogs },
		{"-catalogs", "", "", "", 0, action_create_catalogs },

		{"-s", "", "", "", 1, action_search_particles },
		{"-sp", "", "", "", 1, action_search_particles },
		{"-search", "", "", "", 1, action_search_particles },

		{"-snapf", "%%s+snapnum", "specify the basename of snapshot files", "[snap name]", -1, set_snapf },

		{"-subf", "%%s+snapnum", "specify the basename of subfind files", "[subf_name]", -1, set_subf },

		{"-wdir", "%%s+snapnum", "specify the working dir", "[path]", -1, set_subf },

		{"-snapdir", "%%s+snapnum", "specify the snapshot folder", "[path]", -1, set_snapdir },

		{"-subfdir", "%%s+snapnum", "specify the subfind files folder", "[path]", -1, set_subfdir },

		{"-num", "", "specify the snap/subfind number", "[num]", -1, set_num },

		{"-g", "4", "specify the number of stellar generation to consider", "[num]", -1, set_stellargenerations },

		{"-list", "", "add list file", "[path/to/file_name] [type]\n\t"
		 "where type is in <-2,-1,[0..5]>\n\t"
		 "with -2 = structured file, -1 = mixed, [0..5] typed", -1, set_list }, 

		{"-v", "0", "set verbosity level", "[num]", -1, set_verbose },

		{"-ask", "", "ask info", "", 2, ask },

		{"-h", "", "ask for help", "", 2, help },
		{"-help", "", "ask for help", "", 2, help } };
	       

int n_args = (int)((&args)[1] - args);
int help( char **argv, int n, int mode)
{  
  UNUSED(argv);
  UNUSED(n);
  UNUSED(mode);

  if( !mode )
    return 0;
  
  printf( "Here it goes a little help from my side:\n\n");
  
  for ( int i = 0; i < n_args; i++ ) {    
    args[i].action( NULL, i, 0);                                // set the default

    if( args[i].family < 0 )
      printf("* %s\n\t%s ",                                     // print help
	     args[i].help, args[i].argument);
    else {
      int f = args[i].family;
      printf("* %s\n\t< ", args_f[args[i].family].arguments_family);
      while( (i < n_args) && (args[i].family == f) ) {
	char c = ((i+1 < n_args) && (args[i+1].family == f))? '|' : '>';
	printf("%s %c ", args[i].argument, c); i++;}
      i--;
    }
    if( args[i].example[0] != '\0' )
      printf("%s ", args[i].example );
    if( args[i].dflt[0] != '\0' )
      printf(" --> default is %s", args[i].dflt);
    printf("\n\n");}
  printf("\n");

  printf("typical command lines may look like the following:\n\n"
	 
	 "[*] ./cidx -c -num 090 -snapdir ./\n"
	 "    create catalogs for snapshot num 090 looking for snapshot files in the current dir\n"
	 "    subfind files are searched under ./grups_090/sub_090*\n\n"
	 
	 "[*] ./cidx -c -num 101 -snapf mysnap.le_ -snapdir sdir_ -subf subhaloes_ -subfdir data_\n"
	 "    create catalogs for snapshot files sdir_101/mysnap.le_101* and subfind outputs are\n"
	 "    searched as data_101/subhaloes_101*\n\n"
	 
	 "[*] ./cidx -s -num 045 -list list_of_stars 4 -list mixed_list\n"
	 "    it will search the particles whose ids are in file \"list_of_stars\" assuming that they all are stars\n"
	 "    it will then search the particles whose ids are in file \"mixed_list\" assuming that they type is whatsoever\n"
	 "    the catalogs loaded are supposed to be named \"catalog_045.?\" with ? ranging in [0..5]\n"
	 );
	 
  printf("\n");
  help_given = 1;
  return ++n;
}

#undef UNUSED
