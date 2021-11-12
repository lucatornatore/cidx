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

#if defined(DEBUG)
#define dprintf(...) printf(__VA_ARGS__)
#else
#define dprintf(...)
#endif

#define NTYPES   6
#define NAMESIZE 100

int                get_infos( char *, int *, unsigned long long int *, unsigned long long int *, int * );
unsigned long long get_ids  ( char *, int *, int, int, char *, double );
int                dump_file( char *, int , unsigned long long, int *, char * );

int main( int argc, char **argv )
{
  
  if ( argc < 4 ) {
    printf( "please give me the needed cmd-line arguments\n" );
    exit(1); }
  
  // initialize
  //
  
  char name[NAMESIZE];
  int type;
  int types[NTYPES] = {0};
  int id_size;
  int nfiles;
  double selection;
  unsigned long long int Nall = 0, Nparts[NTYPES] = {0}, nids = 0;

  {
    snprintf( name, NAMESIZE, "%s", *(argv+1));
    selection = atof(*(argv+2));
    int j = 3;
    while( j < argc )
      {
	int type = atoi(*(argv+j));
	
	if( type < 0 ) { j = argc;
	  for( int i = 0; i < NTYPES; i++ )
	    types[i] = 1; }
	else
	  types[type] = 1;
	j++;
      }
  }
  
  dprintf("getting infos.. "); fflush(stdout);
  
  nfiles = get_infos( name, types, &Nall, Nparts, &id_size );
  nids   = Nall*selection;
  
 #if defined(DEBUG)
  printf("got infos:\nID's size is %d bytes\n", id_size);
  for ( int t = 0; t < NTYPES; t++ )
    printf("%c type %d: %llu\n", (types[t]?'>':' '), t, Nparts[t]);
  printf("%llu particles of selected types\n", Nall );
  if( selection > 0 )
    printf("about %llu will be selected\n", (unsigned long long)(Nall*selection) );  
 #endif

  // get ids
  //

  dprintf("getting the ids.. "); fflush(stdout);
  
  char *ids = (char*)malloc( id_size * nids );
  nids = get_ids( name, types, nfiles, id_size, ids, selection );
  
  dprintf("done, got %llu particles.\nwriting the file.. ", nids); fflush(stdout);
  
  dump_file(ids, id_size, nids, types, name );
  
  dprintf("done\n");
  
  free(ids);
  return 0;
}


typedef struct { int tag1; char name[4]; int len; int tag2; } head_t;

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


int get_infos( char *namebase, int *types, unsigned long long int *Nall, unsigned long long int *Nparts, int *id_size )
{
  FILE *file;
  int   nfiles = 1;
  char  name[NAMESIZE+5];


  // find in how many files the snapshot is splitted in
  //
  snprintf( name, NAMESIZE+5, "%s", namebase );
  file = fopen( name, "r" );
  if ( file == NULL )
    {
      snprintf( name, NAMESIZE+5, "%s.0", namebase );
      file = fopen( name, "r" );
      if ( file == NULL )
	return -1;

      do
	{
	  fclose(file);
	  snprintf( name, NAMESIZE+5, "%s.%d", namebase, nfiles );
	  nfiles += ((file = fopen( name, "r" )) != NULL);
	}
      while( file != NULL);
    }

  if( nfiles > 1 )
    snprintf( name, NAMESIZE+5, "%s.0", namebase );
  else
    snprintf( name, NAMESIZE+5, "%s", namebase );
  file = fopen( name, "r" );

  // get how many particles are present
  // in all the snapshot and in the first file
  //
  
  unsigned long long int nparts_local_all = 0;
  int                    nparts_local[NTYPES];
  int                    ret;

  ret = seek_block( file, "HEAD");
  if ( ret == 0 )
    {
      printf("unable to find the HEAD block in snapshot files\n");
      return -2;
    }

  // local particles
  //
  ret = fread( nparts_local, sizeof(int), NTYPES, file );
  for( int i = 0; i < NTYPES; i++ )
    nparts_local_all += nparts_local[i];

  // get to total particles
  //
  unsigned int low[NTYPES];
  unsigned int high[NTYPES];
  
  fseeko(file, (off_t)(sizeof(int)*2+sizeof(double)*8), SEEK_CUR);
  ret = fread( &low[0], sizeof(int), NTYPES, file );

  // get to high words of total particles
  //
  fseeko(file, (off_t)(sizeof(int)*4+sizeof(double)*4), SEEK_CUR);
  ret = fread( &high[0], sizeof(int), NTYPES, file );

  // sum-up
  //
  for( int i = 0; i < NTYPES; i++ )
    if( types[i] )
      *Nall += (Nparts[i] = ((unsigned long long)high[i]<<31) + low[i]);

  // get the ID size
  //
  ret = seek_block( file, "ID  ");
  if ( ret )
    *id_size = (size_t)((unsigned long long)ret / nparts_local_all);
  
  fclose( file );

  return nfiles;
}


unsigned long long get_ids(char *namebase, int *types, int nfiles, int id_size, char *ids, double selection )
{

  unsigned long long amount = 0;
  char *data = ids;

  if ( selection > 0 )
    srand48( time(NULL) );
  
  for( int ff = 0; ff < nfiles; ff++ )
    {
      char name[NAMESIZE+5];
      // open the current file
      //
      if( nfiles > 1 )
	sprintf( name, "%s.%d", namebase, ff );
      else
	sprintf( name, "%s", namebase );
      FILE *file = fopen( name, "r" );

      // get the local number of particles
      //
      int ret = seek_block( file, "HEAD");
      int nparts[NTYPES];
      ret = fread( nparts, sizeof(int), NTYPES, file );

      // move to IDS
      //
      ret = seek_block( file, "ID  ");

      for( int t = 0; t < NTYPES; t++ )
	{
	  if( !types[t] )
	    // this type won't be loaded
	    fseeko( file, (off_t)(id_size * nparts[t]), SEEK_CUR );
	  else
	    {
	      if( selection == 0 ) {
		// load all
		amount += nparts[t];
		ret = fread( data, id_size, nparts[t], file );
		data += (off_t)(id_size * nparts[t]); }
	      
	      else {
		char *buffer = (char *)malloc( nparts[t] * id_size );		
		ret = fread( buffer, id_size, nparts[t], file );
		int max = selection*(nparts[t]+0.5);
		
		switch( id_size )
		  {
		  case 4: {
		    int *ibuffer  = (int*)buffer;
		    int *idata    = (int*)data;
		    int *max_data = idata+max;
		    
		    for( int k = 0; k < nparts[t] && idata < max_data; k++ ) {		      
		      int get = (drand48() < selection);
		      amount += get;
		      *idata = ( get ? ibuffer[k] : 0);
		      idata += ( get ? 1 : 0); }}
		    break;
		  case 8: {
		    unsigned long long *lbuffer  = (unsigned long long*)buffer;
		    unsigned long long *ldata    = (unsigned long long*)data;
		    unsigned long long *max_data = ldata+max;
		    for( int k = 0; k < nparts[t] && ldata < max_data; k++ ) {		      
		      int get = (drand48() < selection);
		      amount += get;
		      *ldata = ( get ? lbuffer[k] : 0);
		      ldata += ( get ? 1 : 0); }}
		    break;
		  }
		free(buffer);		
	      }
	    }
	}
      
      
      fclose(file);
    }

  
  return amount;
}


int dump_file( char *ids, int id_size, unsigned long long nids, int *types, char *namebase )
{
  FILE *file;
  char name[NAMESIZE+30];

  snprintf( name, NAMESIZE+10, "%s.ids", namebase );
  for( int t = 0; t < NTYPES; t++ )
    if( types[t] )
      snprintf( name+strlen(name), NAMESIZE+10-strlen(name)-1, ".%d", t );

  file = fopen( name, "w" );

  if( file == NULL ) {
    printf("unable to create the file %s\n", name );
    return -1; }

  fwrite( &id_size, sizeof(int), 1, file );
  fwrite( &nids, sizeof(long long), 1, file );
    
  fwrite( ids, id_size, nids, file );

  fclose(file);

  return 0;
}
