#define INLINE
#include "cidx.h"

int     n_stars_generations, n_stars_generations_def = 4;
int     id_bitshift;
int     shiftbox = 0;
int     verbose = 0;
FILE   *details = NULL;
//FILE   *pdetails = NULL;

int     Nthreads, me;   // the number of omp threads, the thread id

num_t  *IDranges;
num_t   type_positions[NTYPES];
num_t   myNID;
num_t  *all_NID;
fof_table_t *fof_table;

num_t  *IDdomains[NTYPES];

_Alignas(32) pidtype_t  *IDs;
_Alignas(32) pidtype_t **all_IDs;
_Alignas(32) particle_t  *P[NTYPES]     ={NULL};
_Alignas(32) particle_t **P_all[NTYPES] ={NULL};
_Alignas(32) particle_t   *Pbase        = NULL;
_Alignas(32) list_t *List  = NULL;
_Alignas(32) int    *Types = NULL;

num_t    Nparts[NTYPES+1]     = {0};
num_t   *Nparts_all[NTYPES+1] = {NULL};
num_t    Np_all[NTYPES+1]     = {0};


num_t  Nl, myNl;

mask_galaxies_in_fof_t *mask_fofgal = NULL;
int    mask_fofgal_n = 0;

int    sizeof_in_data = 4;

#if defined(DOUBLE_OUT)
int    sizeof_out_data = sizeof(double);
#else
int    sizeof_out_data = sizeof(float);
#endif

outsnap_t *outsnap_data = NULL;

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
char catalog_table_name[CATALOG_NAME_SIZE];
char list_name[CATALOG_NAME_SIZE];
char snapnum[NUM_SIZE];

num_t diagnostic[2];

char *list_names[MAX_N_LISTS] = {NULL};
int   list_types[MAX_N_LISTS] = {0};
int   Nlists = 0;
int   action = UNDEFINED;
int   help_given = 0;

// note: in the following functions, "mode" should be passed
// as argc from the main, so that it is > 0 and also conveys
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

int action_make_foftable( char **argv, int n, int mode )
{
  UNUSED(argv);
  UNUSED(n);
  UNUSED(mode);
  if( mode > 0 )
    action |= BUILD_FOF_TABLE;
  return n;
}

int action_create_snapshots( char **argv, int n, int mode )
{
  UNUSED(argv);
  UNUSED(n);
  UNUSED(mode);
  if( mode > 0 )
    action |= WRITE_MASK_FILES;
  return n;
}

int set_mask_fofgal( char **argv, int n, int mode )
{
  if( mode > 0 )
    {
      int fof_id = -1;
      int first = 1;
      n++;
      while( (n<mode) && (**(argv + n) != '-' ) )
	{	  
	  if ( isdigit( **(argv+n) ) )
	    {
	      if( (!first ) || (fof_id < 0) ) {
		mask_fofgal = (mask_galaxies_in_fof_t*)realloc( mask_fofgal, sizeof(mask_galaxies_in_fof_t)*(++mask_fofgal_n));
		memset( &mask_fofgal[mask_fofgal_n-1], 0, sizeof(mask_galaxies_in_fof_t) ); }
	      
	      int gal_id = -1;
	      if( fof_id < 0 ) fof_id = atoi( *(argv + n++) );
	      else { gal_id = atoi( *(argv + n++) ); first = 0; }
	      
	      mask_fofgal[mask_fofgal_n-1].fof_num = fof_id;
	      mask_fofgal[mask_fofgal_n-1].ngal    = gal_id;
	      
	    }
	  else if ( **(argv+n) != '-' )
	    printf("arg %d [\" %s \" ]is meaningless\n", n, *(argv+n) );
	}

      n--;
    }
  else
    {
      if( mask_fofgal != NULL ) {
	free( mask_fofgal );
	mask_fofgal = NULL; }	
    }

  return n;
}


int set_shiftbox( char **argv, int n, int mode )
{
  UNUSED(argv);
  shiftbox = ( mode > 0 );
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
      if( Nlists == MAX_N_LISTS-1 )
	{
	  printf("[reading args][set list] too many list files: "
		 "you can specify %d lists at maximum\n"
		 "%s file will be ignored\n",
		 MAX_N_LISTS,
		 *(argv + ++n));
	  return 0;
	}
      Nlists++;
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
      for ( int l = Nlists-1; l >= 0; l-- )
	free( list_names[l] );
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
		
		{"-t", "not active", "write a table with fof properties", "", -1, action_make_foftable },

		{"-w", "not active", "create snapshots with particles from the masked fof/galaxies", "", -1, action_create_snapshots },

		{"-m", "none", "set the masked fof/galaxies, you may have as many as you want", "-m FOFN [gal num1] [gal num2] [gal num 3] ...", -1, set_mask_fofgal },

		{"-snapf", "%%s+snapnum", "specify the basename of snapshot files", "[snap name]", -1, set_snapf },

		{"-subf", "%%s+snapnum", "specify the basename of subfind files", "[subf_name]", -1, set_subf },

		{"-wdir", "%%s+snapnum", "specify the working dir", "[path]", -1, set_subf },

		{"-snapdir", "%%s+snapnum", "specify the snapshot folder", "[path]", -1, set_snapdir },

		{"-subfdir", "%%s+snapnum", "specify the subfind files folder", "[path]", -1, set_subfdir },

		{"-num", "", "specify the snap/subfind number", "[num]", -1, set_num },

		{"-g", "4", "specify the number of stellar generation to consider", "[num]", -1, set_stellargenerations },

		{"-shiftbox", "0", "add -(boxsize/2) to positions", "[0|1]", -1, set_shiftbox },

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
