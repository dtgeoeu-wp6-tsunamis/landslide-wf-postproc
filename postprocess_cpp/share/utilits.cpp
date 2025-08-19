#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#include "utilits.h"

#define iscomment(ch) (((ch) == ';' || (ch) == '#') ? 1 : 0)
//int iscomment( int ch ) { if( ch == ';' || ch == '#' ) return 1; else return 0; }
#define EOS '\0'

// ==================================================================
// =======================  Class Option List  ======================
cOptions::cOptions(const int num_opt, const char *lst_label, const char *lst_about)
{
  int k;
  char *cp0,*cp;
  char *buf = new char [256];

  nopt = num_opt;

  label = new char* [nopt];
  for( cp0=(char *)lst_label, k=0; k<nopt-1; k++ ) {
    cp = strchr( cp0, '|' ); if(!cp) {nopt=0; return;}
    sscanf( cp0, " %s", buf );
    label[k] = strdup(buf);
    cp0 = cp+1;
  }
  sscanf( cp0, " %s", buf );
  label[nopt-1] = strdup(buf);

  about = new char* [nopt];
  for( cp0=(char *)lst_about, k=0; k<nopt-1; k++ ) {
    cp = strchr( cp0, '|' ); if(!cp) {nopt=0; return;}
    memset(buf,0,256);
    strncpy(buf,cp0,cp-cp0);
    about[k] = strdup(buf);
    cp0 = cp+1;
  }
  about[nopt-1] = strdup(cp0);

}

cOptions::~cOptions()
{
}

int cOptions::getI( const char *given_label )
{
  int iopt = -1;

  for( int k=0; k<nopt; k++ ) {
    if( !stricmp( given_label, label[k] ) ) {
      iopt = k;
      break;
    }
  }

  return iopt;
}

int cOptions::chkAvailability( const char *given_label )
{
  int iopt = 0;

  for( int k=0; k<nopt; k++ ) {
    if( !stricmp( given_label, label[k] ) ) {
      iopt = 1;
      break;
    }
  }

  return iopt;
}

char *cOptions::listAll()
{
  char buf[nopt*16];
  memset( buf,0,nopt*16 );

  for( int k=0; k<nopt; k++ ) {
    strcat( buf, " " );
    strcat( buf, label[k] );
  }

  return strdup(buf);
}

// ==================================================================


#define ERRORFILE "error.msg"
cMsg Msg;
cErrMsg Err;
cLog Log;

// ==================================================================
// ==========================  Class Message  =======================
cMsg::cMsg()
{
  enabled = 1;
  timestamp_enabled = 0;
  setchannel(MSG_OUTSCRN);
  setfilename( "default.msg" );
}

cMsg::~cMsg()
{

}

void cMsg::enable()
{
  enabled = 1;
}

void cMsg::disable()
{
  enabled = 0;
}

void cMsg::timestamp_enable()
{
  timestamp_enabled = 1;
}

void cMsg::timestamp_disable()
{
  timestamp_enabled = 0;
}

void cMsg::setfilename( const char* newfilename )
{
  sprintf( fname, "%s", newfilename );
}

void cMsg::setchannel( int newchannel )
{
  channel = newchannel;
}

int cMsg::print( const char* fmt, ... )
{
  if(!enabled) return 0;

  if( channel & MSG_OUTSCRN ) {
    va_list arglst;
    va_start(arglst, fmt);
    vprintf( fmt, arglst );
    printf( "\n" );
    va_end(arglst);
  }
  if( channel & MSG_OUTFILE ) {
    va_list arglst;
    va_start(arglst, fmt);
    FILE *fp = fopen( fname, "at" );
    if( timestamp_enabled )fprintf( fp, "%s --> ", utlCurrentTime() );
    vfprintf( fp, fmt, arglst );
    fprintf( fp, "\n" );
    fclose( fp );
    va_end(arglst);
  }

  return 0;
}
// =================== END Class Message =================


// ==================================================================
// =======================  Class Error Message  ====================
cErrMsg::cErrMsg()
{
  enabled = 1;
  timestamp_enabled = 0;
  setchannel(MSG_OUTSCRN|MSG_OUTFILE);
  setfilename( "error.msg" );
  flag = 0;
}


cErrMsg::~cErrMsg()
{

}


void cErrMsg::set()
{
  flag = 1;
}


void cErrMsg::unset()
{
  flag = 0;
}


int cErrMsg::exist()
{
  if( flag ) return 1; else return 0;
}


int cErrMsg::post( const char* fmt, ... )
{
  if(!enabled) return 0;

  /* Simple re-direction to the method of the parent class does not work
  va_list arglst;
  va_start(arglst, fmt);
  print( fmt, arglst );
  va_end(arglst);
  */

  if( channel & MSG_OUTSCRN ) {
    va_list arglst;
    va_start(arglst, fmt);
    vprintf( fmt, arglst );
    printf( "\n" );
    va_end(arglst);
  }
  if( channel & MSG_OUTFILE ) {
    va_list arglst;
    va_start(arglst, fmt);
    FILE *fp = fopen( fname, "at" );
    if( timestamp_enabled )fprintf( fp, "%s --> ", utlCurrentTime() );
    vfprintf( fp, fmt, arglst );
    fprintf( fp, "\n" );
    fclose( fp );
    va_end(arglst);
  }

  return 1;
}


// We need this routine to overcome the stdlib 'exit' function overloading by the current class
void my_exit()
{
  exit(EXIT_FAILURE);
}

void cErrMsg::exit( const char* fmt, ... )
{
  if(!enabled) return;

  /* Simple re-direction to the method of the parent class does not work
  va_list arglst;
  va_start(arglst, fmt);
  print( fmt, arglst );
  va_end(arglst);
  */

  if( channel & MSG_OUTSCRN ) {
    va_list arglst;
    va_start(arglst, fmt);
    printf( "\nCritical error! Aborting execution...\n" );
    vprintf( fmt, arglst );
    printf( "\n" );
    va_end(arglst);
  }
  if( channel & MSG_OUTFILE ) {
    va_list arglst;
    va_start(arglst, fmt);
    FILE *fp = fopen( fname, "at" );
    if( timestamp_enabled )fprintf( fp, "%s --> ", utlCurrentTime() );
    fprintf( fp, "Critical error! Aborting execution...\n" );
    vfprintf( fp, fmt, arglst );
    fprintf( fp, "\n" );
    fclose( fp );
    va_end(arglst);
  }

  my_exit();
}

char* cErrMsg::msgAllocateMem( void )
{
  char msg[256];

  sprintf( msg, "Error allocating memory" );

  return strdup(msg);
}

char* cErrMsg::msgOpenFile( const char *fnam )
{
  char msg[256];

  sprintf( msg, "Cannot open file %s", fnam );

  return strdup(msg);
}

char* cErrMsg::msgReadFile( const char *fnam, int line, const char *errmsg )
{
  char msg[256];

  if( line )
    sprintf( msg, "Error reading file %s, line number %-d : %s", fnam, line, errmsg );
  else
    sprintf( msg, "Error reading file %s : %s", fnam, errmsg );

  return strdup(msg);
}
// =================== END Class Message =================


// ==================================================================
// ==========================  Class Log  ===========================
cLog::cLog()
{
  enabled = 0;
  timestamp_enabled = 1;
  setchannel(MSG_OUTFILE);
  setfilename( "default.log" );
  fpMsg = NULL;
}

cLog::~cLog()
{
  if( fpMsg ) fclose(fpMsg);
}

void cLog::start( const char* filename )
{
  int ifenabled=0;

  enabled = 1;
  setfilename( filename );

  if( timestamp_enabled ) {
    ifenabled = 1;
    timestamp_disable();
  }
  print("\n\n =============  Starting new log at %s", utlCurrentTime() );
  if( ifenabled ) timestamp_enable();
}


void cLog::hold()
{
  fpMsg = fopen( fname, "at" );
}


void cLog::release()
{
  if(fpMsg) { fclose(fpMsg); fpMsg = NULL; }
}


int cLog::print( const char* fmt, ... )
{
  if(!enabled) return 0;

  if( channel & MSG_OUTSCRN ) {
    va_list arglst;
    va_start(arglst, fmt);
    vprintf( fmt, arglst );
    printf( "\n" );
    va_end(arglst);
  }
  if( channel & MSG_OUTFILE ) {
    va_list arglst;
    va_start(arglst, fmt);
    if( fpMsg ) {
      if( timestamp_enabled )fprintf( fpMsg, "%s --> ", utlCurrentTime() );
      vfprintf( fpMsg, fmt, arglst );
      fprintf( fpMsg, "\n" );
    }
    else {
      FILE *fp = fopen( fname, "at" );
      if( timestamp_enabled ) fprintf( fp, "%s --> ", utlCurrentTime() );
      vfprintf( fp, fmt, arglst );
      fprintf( fp, "\n" );
      fclose( fp );
    }
    va_end(arglst);
  }

  return 0;
}

// =================== END Class Log =================



// ==================================================================
// Command-line processing
int utlCheckCommandLineOption( int argc, char **argv, const char *option, int letters_to_compare )
{
  int k;

  for( k=1; k<argc; k++ ) {
    if( argv[k][0] != '-' && argv[k][0] != '/' ) continue;
    if( !strncmp( argv[k]+1, option, letters_to_compare ) ) break;
  }

  if( k == argc )
    return 0;
  else
    return k;
}



// ******************************************************************
// ***                                                            ***
// ***            S T R I N G   H A N D L I N G                   ***
// ***                                                            ***
// ******************************************************************

// ==================================================================
//  Reads specified option
//  Option: [option-name]=contents
//  Returns 0 by success, otherwise -1
int utlStringReadOption( char *record, char *option_name, char *contents )
{
  int found, length;
  char *cp, *cpe, buf[64];

  cp = record;
  found = 0;
  contents[0] = EOS;

  while( *cp != EOS ) {

    cp = strchr( cp, '[' );
    if( cp == NULL ) break;

    cpe = strchr( cp+1, ']' );
    if( cpe == NULL ) break;

    length = cpe - cp - 1;
    strncpy( buf, cp+1, length );
    buf[length] = EOS;

    if( !strcmp( buf, option_name ) ) {   // found option name

      // skip to option value
      cp = strchr( cpe+1, '=' );
      if( cp == NULL ) break;
      while( isspace( *(++cp) ) ) ;

      if( *cp == EOS )  // no value assigned
        ;
      else if( *cp == '[' )  // no value assigned
        ;
      else if( *cp == '"' ) {  // quoted string
        cpe = strchr( cp+1, '"' );
        if( cpe == NULL ) break;
        length = cpe - cp - 1;
        strncpy( contents, cp+1, length );
        contents[length] = EOS;
      }
      else // string without quotation
        sscanf( cp, "%s", contents );

      found = 1;
      break;

    }  // found option name

    cp++;
  }

  if( found )
    return 0;
  else
    return -1;
}


// ==================================================================
// Reads sequence of integers from a string. Returns number of read
int utlReadSequenceOfInt( char *line, int *value )
{
  int N=0;
  int itmp;
  char *cp=line;

  while( 1 ) {
    while( isspace(*cp) ) cp++;
    if( *cp == EOS || iscomment(*cp) ) return N;
    if( sscanf( cp, "%d", &itmp ) )
      value[N++] = itmp;
    while( !isspace(*cp) && *cp != EOS ) cp++;
  }
}


// ==================================================================
// Reads sequence of doubles from a string. Returns number of read
int utlReadSequenceOfDouble( char *line, double *value )
{
  int N=0;
  double dtmp;
  char *cp=line;

  while( 1 ) {
    while( isspace(*cp) ) cp++;
    if( *cp == EOS || iscomment(*cp) ) return N;
    if( sscanf( cp, "%lf", &dtmp ) )
      value[N++] = dtmp;
    while( !isspace(*cp) && *cp != EOS ) cp++;
  }
}


// ==================================================================
// Picks up n-th word from string (starting from 1). By error (number words < n) returns -1
int utlPickUpWordFromString( const char *str, int n, char *word )
{
  const char *cp;
  int i;

  for( cp=str, i=1; 1; i++ ) {
    while( isspace(*cp) ) cp++;
    if( *cp == EOS ) {
      word[0] = EOS;
      return -1;
    }
    if( i == n ) break;
    while( !isspace(*cp) && *cp != EOS ) cp++;
  }

  sscanf( cp, "%s", word );
  return 0;
}


// ==================================================================
//   Picks up n-th word from string. By error (number words < n) returns NULL
char *utlPickUpWordFromString( const char *str, int n )
{
  const char *cp;
  char buf[256];
  int i;

  for( cp=str, i=1; 1; i++ ) {
    while( isspace(*cp) ) cp++;
    if( *cp == EOS ) return NULL;
    if( i == n ) break;
    while( !isspace(*cp) && *cp != EOS ) cp++;
  }

  sscanf( cp, "%s", buf );
  return strdup(buf);
}


// ==================================================================
// Picks up a word from string starting from position pos.
// Return value: pointer to next char after the word or NULL if failed
char *utlPickUpWordFromString( const char *str, char *pos, char *word )
{
  char *cp;

  if( pos < str ) { word[0] = EOS; return NULL; }
  if( pos > (str+strlen(str)) ) { word[0] = EOS; return NULL; }

  for( cp = pos; isspace(*cp); cp++ ) ;
  if( *cp == EOS ) { word[0] = EOS; return NULL; }
  sscanf( cp, "%s", word );
  cp += strlen(word);

  return cp;
}


// ==================================================================
// Count words in a string
int utlCountWordsInString( const char *line )

{
  int nwords=0;
  const char *cp=line;

  while( 1 ) {
    while( isspace(*cp) ) cp++;
    if( *cp == EOS || iscomment(*cp) ) return nwords;
    nwords++;
    while( !isspace(*cp) && *cp != EOS ) cp++;
  }
}


// ==================================================================
// Write a sequence of words into a string
char *utlWords2String( int nwords, char **word )
{
  char *buf;
  int k,lenstr;

  for( lenstr=0,k=0; k < nwords; k++ )
    lenstr += strlen( word[k] );

  lenstr += (nwords + 1);  // add separators plus final null character

  buf = new char[lenstr];

  memset( buf, 0, lenstr );

  for( k=0; k < nwords; k++ ) {
    if( k>0 ) strcat( buf, " " );
    strcat( buf, word[k] );
  }

  return buf;
}


// ==================================================================
// reads substring from p1 to p2 into a buffer substr
int utlSubString( char *str, int p1, int p2, char *substr )
{
  char *cp;
  int k;

  if( p1<0 || p2<0 || p1>p2 || (unsigned)p2 > (strlen(str)-1) ) return -1;

  for( k=0,cp = &str[p1]; cp <= &str[p2]; cp++ )
    substr[k++] = *cp;
  substr[k] = EOS;

  return 0;
}


// ==================================================================
// Prints string to another starting from given position
void utlPrintToString( char *str, int position, char *insstr )
{
  char *cp;
  int i;

  for( i=0; i<position; i++ )
    if( str[i] == EOS ) str[i] = ' ';
  for( cp=insstr, i=0; *cp != EOS; i++, cp++ )
    str[position+i] = *cp;
}


// ==================================================================
// Case-insensitive string comparison
int stricmp(const char *s1, const char *s2)
{
   if(s1 == NULL) return s2 == NULL ? 0 : -(*s2);
   if(s2 == NULL) return *s1;

   int c1, c2;
   while( (c1 = tolower(*s1)) == (c2 = tolower (*s2)) ) {
     if (*s1 == '\0') break;
     ++s1; ++s2;
   }

   return (c1 - c2);
}


// ==================================================================
// Case-insensitive string comparison of the first n-chars
int strincmp(const char *s1, const char *s2, int nchars)
{
   if(s1 == NULL) return s2 == NULL ? 0 : -(*s2);
   if(s2 == NULL) return *s1;

   int c1, c2, n=1;
   while( (c1 = tolower(*s1)) == (c2 = tolower (*s2)) ) {
     if(*s1 == '\0' || n == nchars) break;
     ++s1; ++s2; ++n;
   }

   return (c1 - c2);
}


// ******************************************************************
// ***                                                            ***
// ***                F I L E   H A N D L I N G                   ***
// ***                                                            ***
// ******************************************************************

// ==================================================================
//  Search for next input field in a file skipping blanks, tabs, empty lines, comments (from ';' or '#' to EOL)
//  Returns number of lines scanned until input field found
int utlFindNextInputField( FILE *fp )
{
  int ch = 0;

L1:
  while( isspace( (ch=fgetc(fp)) ) ) ;

  if( iscomment(ch) ) {
    while( (ch=fgetc(fp)) != '\n' && ch != EOF ) ;
    if( ch == '\n' )
      goto L1;
    else
      return EOF;
  }
  else if( ch == EOF )
    return EOF;

  ungetc( ch, fp );

  return 0;
}


// ==================================================================
// Reads from the next record according to given format string
int fscanfl( FILE *fp, const char *fmt, ... )
{
  char ch;
  int nread(0);

  if( utlFindNextInputField( fp ) == EOF ) return EOF;

  va_list args;
  va_start(args, fmt);
  nread = vfscanf(fp, fmt, args);
  va_end(args);

  while( (ch=fgetc(fp)) != '\n' && ch != EOF ) ;

  return nread;
}


// ==================================================================
// Reads next non-commented and non-empty record from a file
int utlReadNextRecord( FILE *fp, char *record )
{
  int found = 0;
  char *cp;

  while( !found && fgets(record, MaxFileRecordLength, fp) ) {
    for( cp=record; isspace(*cp); cp++ ) ;
    if( !iscomment(*cp) && *cp != EOS ) found = 1;
  }

  if( !found )
    return EOF;
  else {
    while( *cp != '\n' && *cp != '\r' && *cp != '\0' ) cp++;
    *cp = '\0';
    return 0;
  }
}


// ==================================================================
// Reads next non-commented and non-empty record from a file. Counts number of lines
int utlReadNextRecord( FILE *fp, char *record, int &line )
{
  int found = 0;
  char *cp;

  while( !found && fgets(record, MaxFileRecordLength, fp) ) {
    line++;
    for( cp=record; isspace(*cp); cp++ ) ;
    if( !iscomment(*cp) && *cp != EOS ) found = 1;
  }

  if( !found )
    return EOF;
  else {
    while( *cp != '\n' && *cp != '\r' && *cp != '\0' ) cp++;
    *cp = '\0';
    return 0;
  }
}


// ==================================================================
// Reads next non-commented and non-empty record starting with numbers from a file
int utlReadNextDataRecord( FILE *fp, char *record )
{
  int found=0;
  float ftmp;

  while( !found && utlReadNextRecord( fp, record ) != EOF ) {
    if( sscanf(record, "%f", &ftmp) == 1 ) found = 1;
  }

  if( !found )
    return EOF;
  else
    return 0;
}


// ==================================================================
// Reads next non-commented and non-empty record starting with numbers from a file. Counts number of lines
int utlReadNextDataRecord( FILE *fp, char *record, int &line )
{
  int found=0;
  float ftmp;

  while( !found && utlReadNextRecord( fp, record, line ) != EOF ) {
    if( sscanf(record, "%f", &ftmp) == 1 ) found = 1;
  }

  if( !found )
    return EOF;
  else
    return 0;
}


// ==================================================================
// Find record containing given pattern in a file
// Read the whole record into the passed buffer 'record' and return 1 or 0
int utlFileFindRecord( const char *fname, const char *pattern, int skip_comments, char *record )
{
  int found = 0;
  char *cp=NULL;

  FILE *fp = fopen(fname, "rt");  if(!fp) { Err.post(Err.msgOpenFile(fname)); return 0; }

  while( !found && fgets(record, MaxFileRecordLength, fp) ) {
    if( (cp = strstr(record, pattern) ) == NULL ) continue;
    if( skip_comments ) {
      for( cp=record; isspace(*cp); cp++ ) ;
      if( iscomment(*cp) ) found = 0;
    }
    else
      found = 1;
  }

  fclose(fp);

  if( !found )
    return 0;
  else {
    for( cp=record; *cp != '\n' && *cp != '\r' && *cp != '\0'; cp++ ) ;
    *cp = '\0';
    return 1;
  }
}


// ==================================================================
// Parce file from current position to a given pattern
// Returns 0 if pattern found. File pointer remains at the first position after the pattern!
int utlFileParceToString( FILE *fp, const char *pattern )
{
  int ch,pos=0;


  while( utlFindNextInputField(fp) != EOF ) {
    ch = fgetc(fp);
    if( ch != pattern[0] ) continue;

    while( (ch=fgetc(fp)) != EOF ) {
      pos++;
      if( ch != pattern[pos] ) break;
      if( pattern[pos+1] == '\0' ) return 0;
    }

    pos = 0;
  }

  return EOF;
}


// ==================================================================
int utlGetNumberOfRecords( const char *fname )
{
  FILE *fp = fopen(fname, "rt");  if(!fp) { Err.post(Err.msgOpenFile(fname)); return 0; }
  char *record = (char *)malloc(1024*8);

  int nrec = 0;
  while( utlReadNextRecord( fp, record ) != EOF )
    nrec++;

  fclose( fp );
  free(record);

  return nrec;
}


// ==================================================================
int utlGetNumberOfDataRecords( const char *fname )
{
  FILE *fp = fopen(fname, "rt");  if(!fp) { Err.post(Err.msgOpenFile(fname)); return 0; }
  char *record = new char[1024*32];

  int nrec = 0;
  while( utlReadNextDataRecord( fp, record ) != EOF )
    nrec++;

  delete [] record;
  fclose( fp );

  return nrec;
}


// ==================================================================
//  Option-   [option-name]=contents  (can be in quotation marks)
//  Returns 0 by success, otherwise -1
int utlFileReadOption( const char *fname, const char *option_name, char *contents )
{
  FILE *fp;
  int ierr,ch;
  char pattern[128];

  fp = fopen( fname, "rt" );
  if( !fp ) return Err.post(Err.msgOpenFile(fname));

  sprintf( pattern, "[%s]", option_name );
  ierr = utlFileParceToString( fp, pattern ); if(ierr) return -1;

  ierr = utlFindNextInputField(fp); if(ierr==EOF) return -1;
  ch = fgetc(fp); if(ch != '=') return -1;

  ierr = utlFindNextInputField(fp); if(ierr==EOF) return -1;
  ch = fgetc(fp);

  if(ch == '"') { // quoted string
    int k=0;
    while( (ch=fgetc(fp)) != '"' && ch != EOF ) {
      contents[k++] = ch;
    }
    if( ch == EOF ) return -1;
    contents[k] = EOS;
  }
  else { // option value without quotation
    ungetc(ch,fp);
    ierr = fscanf( fp, "%s", contents ); if(ierr!=1) return -1;
    }

  fclose( fp );

  return 0;

  /*
  char *cp,pattern[128];
  char *record = new char[1024*8];

  sprintf( pattern, "[%s]", option_name );
  cp = utlFileFindRecord( fname, pattern, record ); if( cp == NULL ) return -1;

  for( cp=strchr(cp,']')+1; isspace(*cp); cp++ ) ;
  if(*cp != '=') return -1;

  for( cp=cp+1; isspace(*cp); cp++ ) ;
  if(*cp == '"') { // quoted string
    int k=0;
    contents[0] = '\0';
    for( cp=cp+1; *cp != '"' && *cp != EOS; cp++,k++ )
      contents[k] = *cp;
    if( *cp == EOS ) return -1;
    contents[k] = '\0';
  }
  else { // option value without quotation
    if( sscanf(cp, "%s", contents) !=1 ) return -1;
  }

  delete [] record;

  return 0;
  */
}

/*
// ==================================================================
//  Read INTEGER option from a file
//  On success returns option value, otherwise returns 0 and sets Err flag
int utlFileReadOption( const char *fname, const char *option_name )
{
  char buf[32]; // by integers short string is enough
  int ipar;

  if( utlFileReadOption(fname, option_name, buf) ) {Err.set(); return 0;}
  if( sscanf(buf, "%d", &ipar) != 1 ) {Err.set(); return 0;}

  return ipar;
}


// ==================================================================
//  Read DOUBLE option from a file
//  On success returns option value, otherwise returns 0 and sets Err flag
double utlFileReadOption( const char *fname, const char *option_name )
{
  char buf[32]; // by integers short string is enough
  double dpar;

  if( utlFileReadOption(fname, option_name, buf) ) {Err.set(); return 0;}
  if( sscanf(buf, "%lf", &dpar) != 1 ) {Err.set(); return 0;}

  return dpar;
}
*/

// ==================================================================
//   Change file extension. Return new file name
char *utlFileChangeExtension( const char *fname, const char *newext )
{
  char buf[256], *cp;
  int len,extlen;

  len = strlen(fname);
  extlen = strlen(newext);
  memset( buf, 0, len+extlen+1 );
  strcpy( buf, fname );

  for( cp=&buf[len-1]; *cp!='.' && cp>=&buf[0]; cp-- ) ;

  if( *cp=='.' )
    sprintf( cp+1, "%s", newext );
  else
    sprintf( &buf[len], ".%s", newext );

  return strdup( buf );
}


// ==================================================================
//   Add new extension. Return new file name
char *utlFileAddExtension( const char *fname, const char *addext )
{
  char buf[256];
  int len,extlen;

  len = strlen(fname);
  extlen = strlen(addext);
  memset( buf, 0, len+extlen+1 );
  sprintf( buf, "%s.%s", fname, addext );

  return strdup( buf );
}


// ==================================================================
//   Remove file extension. Return file name without extension
char *utlFileRemoveExtension( const char *fname )
{
  char buf[256], *cp;

  strcpy( buf, fname );
  for( cp=&buf[strlen(fname)-1]; *cp!='.' && cp>=&buf[0]; cp-- ) ;
  if( *cp=='.' ) *cp = '\0';

  return strdup( buf );
}


// ==================================================================
//   Split file name into path and name
void utlFileSplitName( char *fullname, char *fpath, char *fname )
{
  char *cp;

  for( cp = &fullname[strlen(fullname)-1]; cp > &fullname[0]; cp-- ) {
    if( *cp == '/' || *cp == '\\' ) break;
  }

  if( cp == &fullname[0] ) {
    strcpy( fname, fullname );
    if( fpath ) strcpy( fpath, "./" );
  }
  else {
    strcpy( fname, cp+1 );
    if( fpath ) strncpy( fpath, fullname, cp-fullname+1 );
  }

}


// ==================================================================
// Reads XY ASCII data file. Returns number of data points read
int utlReadXYdata( const char *fname, double **x, double **y )
{
  char record[256];
  int nRec,line;
  double xv,yv;

  FILE *fp = fopen(fname, "rt");  if(!fp) Err.post(Err.msgOpenFile(fname));

  // Calculate number of data lines
  for( nRec=0,line=0; utlReadNextDataRecord( fp, record, line ) != EOF; ) {
    if( sscanf( record, "%lf %lf", &xv, &yv ) != 2 ) return Err.post(Err.msgReadFile(fname, line, "X Y"));
    nRec++;
  }

  // Allocate memory
  *x = new double [nRec];
  *y = new double [nRec];

  // Read data
  rewind( fp );
  line=0;
  for( int k=0; k<nRec; k++ ) {
    utlReadNextDataRecord( fp, record, line );
    sscanf( record, "%lf %lf", &(*x)[k], &(*y)[k] );
  }

  fclose(fp);

  return nRec;
}



// ****************************************************************
// ***                                                          ***
// ***                        MATH & GEOMETRY                   ***
// ***                                                          ***
// ****************************************************************

// ==================================================================
// UNIFORM RANDOM NUMBER GENERATOR (MIXED CONGRUENTIAL METHOD)
double utlRandom( int& seed )
{
  seed = seed * 65539;
  if( seed < 0 ) seed = (seed + 2147483647) + 1;

  return( (double)seed * 0.4656613e-9 );
}


// ==================================================================
// returns random double value between v_from and v_to
// be sure to call srand() before calling this function
double utlRandom( double v_from, double v_to )
{
  return( (double)rand()/RAND_MAX*(v_to-v_from) + v_from );
}


// ==================================================================
// returns random integer value between i_from and i_to
// be sure to call srand() before calling this function
int utlRandom( int i_from, int i_to )
{
  return( rand() % (i_to-i_from+1) + i_from );
}


//========================================================
// NORMAL RANDOM NUMBER GENERATOR (BOX-MULLER METHOD)
double utlNormal( double avg, double stdev )
{
  double r1, r2, wnorm;

  r1 = utlRandom( 0., 1. );
  r2 = utlRandom( 0., 1. );
  wnorm = sqrt(-2*log(r1)) * sin(2*My_PI*r2);

  return( avg + stdev*wnorm );
}


//========================================================
int utlPointInPolygon( int nvrtx, double *xvrtx, double *yvrtx, double x, double y )
// Check if point lies inside a polygon
{
  int l,ln,nsec;
  double yside;

  for( nsec=0, l=0; l<nvrtx; l++ ) {

    ln = l + 1;
    if( ln == nvrtx) ln = 0;

    if( (xvrtx[l]-x)*(xvrtx[ln]-x) > 0 ) continue;

    if( (xvrtx[l]-x)*(xvrtx[ln]-x) == 0 )
      return utlPointInPolygon( nvrtx, xvrtx, yvrtx, (x*(1+1.e-8)), y );

    yside = yvrtx[l] + (yvrtx[ln]-yvrtx[l])/(xvrtx[ln]-xvrtx[l])*(x-xvrtx[l]);

    if( yside > y ) nsec++;
  }

  if( (nsec/2)*2 != nsec )
    return 1;
  else
    return 0;

}


//========================================================
double utlDistPoint2Line( double x0, double y0, double x1, double y1, double x2, double y2 )
// Distance from point (x0,y0) to the line defined by the two points (x1,y) and (x2,y2)
{

  return( fabs(-(y2-y1)*x0 + (x2-x1)*y0 + x1*y2 - x2*y1) / sqrt((y2-y1)*(y2-y1)+(x2-x1)*(x2-x1)) );

}


//========================================================
// Haversine formula for distance between any two points on the Earth surface
double utlSphereDistance( const double lon1, const double lat1, const double lon2, const double lat2 )
{
  double a = pow(sindeg((lat2-lat1)/2), 2.) + cosdeg(lat1)*cosdeg(lat2)*pow(sindeg((lon2-lon1)/2), 2.);

  return( REARTH * 2 * asin(sqrt(a)) );
}


//========================================================
double utlSphereStrike( double lon1, double lat1, double lon2, double lat2 )
{
  double strike;

  if( (lat1 == lat2) && (lon1 == lon2) ) {
    strike = 0.;
  }
  else if( lon1 == lon2 ) {
    if( lat1 > lat2 )
      strike = 180.;
    else
      strike = 0.;
  }
  else {
    double distRad = utlSphereDistance( lon1,lat1,lon2,lat2 ) / REARTH;
    strike = r2g( asin( cosdeg(lat2)*sindeg(lon2-lon1)/sin(distRad) ) );

    if( (lat2 > lat1) && (lon2 > lon1) ) {
    }
    else if( (lat2 < lat1) && (lon2 < lon1) ) {
      strike = 180.0 - strike;
    }
    else if( (lat2 < lat1) && (lon2 > lon1) ) {
      strike = 180.0 - strike;
    }
    else if( (lat2 > lat1) && (lon2 < lon1) ) {
      strike += 360.0;
    }
  }

//  if( strike > 180.0 ) strike -= 180.0;
  return strike;
}


// ****************************************************************
// ***                                                          ***
// ***            T I M E  &  T I M E  L A B E L S              ***
// ***                                                          ***
// ****************************************************************

// ==================================================================
// Returns current time as a string
char *utlCurrentTime()
{
  time_t timer;
  char *cp;

  timer = time(NULL);
  cp = asctime(localtime(&timer));
  cp[strlen(cp)-1] = '\0';

  return cp;
}


// ==================================================================
// Splits elapsed time into hours:minutes:seconds.milliseconds
// ctime is time in seconds with milliseconds as fraction
void utlSplitElapsedTime( double ctime, int& nHour, int& nMin, int& nSec, int& nMsec )
{
  int fullSec;

  fullSec = (int)ctime;

  nMsec = (int)((ctime - fullSec)*1000 + 0.1);
  nHour = fullSec/3600;
  nMin = (fullSec - nHour*3600)/60;
  nSec = fullSec - nHour*3600 - nMin*60;

  return;
}


// ==================================================================
// Splits elapsed time into hours:minutes:seconds.milliseconds and returns formatted string
// ctime is time in seconds with milliseconds as fraction
char *utlFormatElapsedTime( double ctime )
{
  int nHour,nMin,nSec,nMsec;
  static char timebuf[32];

  utlSplitElapsedTime( ctime, nHour, nMin, nSec, nMsec );

  if( nMsec > 0 )
    sprintf( timebuf, "%2.2d:%2.2d:%2.2d.%3.3d", nHour,nMin,nSec,nMsec );
  else
    sprintf( timebuf, "%2.2d:%2.2d:%2.2d", nHour,nMin,nSec );

  return timebuf;
}


// ==================================================================
char *utlFormatElapsedTime( int ctime )
{
  int nHour,nMin,nSec,nMsec;
  static char timebuf[32];

  utlSplitElapsedTime( (double)ctime, nHour, nMin, nSec, nMsec );

  sprintf( timebuf, "%2.2d:%2.2d:%2.2d", nHour,nMin,nSec );

  return timebuf;
}


// ==================================================================
char *utlTimeStamp( const char *fmt )
// 'fmt' argument contains formatting rules valid for the C-library strftime() function
{
  time_t rawtime;
  struct tm * timeinfo;
  char buffer[128];

  time (&rawtime);
  timeinfo = localtime (&rawtime);

  strftime( buffer, 128, fmt, timeinfo );

  return strdup(buffer);
}
