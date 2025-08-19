// last modified 11.11.2014
#ifndef UTIL_H
#define UTIL_H

#include <math.h>

#define MaxFileRecordLength 16384
#define RealMax 1.e+30
#define RealMin 1.e-30
#define My_max(a,b)  (((a) > (b)) ? (a) : (b))
#define My_min(a,b)  (((a) < (b)) ? (a) : (b))

#define My_PI  3.14159265358979
#define g2r(x)  (((double)(x))*My_PI/180)
#define r2g(x)  (((double)(x))/My_PI*180)
#define cosdeg(x) cos(g2r(x))
#define sindeg(x) sin(g2r(x))
#define tandeg(x) tan(g2r(x))

#define REARTH 6371.e+3


// Clast list of options
class cOptions {
public:
  int nopt;
  char **label;
  char **about;
  cOptions( int num_opt, const char *lst_label, const char *lst_about );
  ~cOptions();
  int getI( const char *given_label );
  int chkAvailability( const char *given_label );
  char *listAll();
};


// Class Message
#define MSG_OUTSCRN 1
#define MSG_OUTFILE 2
class cMsg
{
protected:
  int enabled;
  int timestamp_enabled;
  int channel;
  char fname[64];

public:
  cMsg();
  ~cMsg();
  void enable();
  void disable();
  void timestamp_enable();
  void timestamp_disable();
  void setchannel( int newchannel );
  void setfilename( const char* newfilename );
  int print( const char* fmt, ... );
  };
extern cMsg Msg;


// Class Error message
class cErrMsg : public cMsg
{
protected:
  int flag;

public:
  cErrMsg();
  ~cErrMsg();
  void set();
  void unset();
  int exist();
  int post( const char* fmt, ... );
  void exit( const char* fmt, ... );
  char* msgAllocateMem();
  char* msgOpenFile( const char *fname );
  char* msgReadFile( const char *fname, int line, const char *expected );
  };
extern cErrMsg Err;


// Class Log
class cLog : public cMsg
{
protected:
  FILE *fpMsg;

public:
  cLog();
  ~cLog();
  void start( const char* filename );
  void hold();
  void release();
  int print( const char* fmt, ... );
  };
extern cLog Log;


// *** Command-line processing ***
int utlCheckCommandLineOption( int argc, char **argv, const char *option, int letters_to_compare );


// ***  String handling  ***
int utlStringReadOption( char *record, char *option_name, char *contents );
int utlReadSequenceOfInt( char *line, int *value );
int utlReadSequenceOfDouble( char *line, double *value );
int utlPickUpWordFromString( const char *string, int n, char *word );
char *utlPickUpWordFromString( const char *string, int n );
char *utlPickUpWordFromString( const char *string, char *pos, char *word );
int utlCountWordsInString( const char *record );
char *utlWords2String( int nwords, char **word );
int utlSubString( char *logstr, int p1, int p2, char *substr );
void utlPrintToString( char *string, int position, char *insert );
int stricmp(const char *s1, const char *s2);
int strincmp(const char *s1, const char *s2, int nchars);


// ***  File handling  ***
int utlFindNextInputField( FILE *fp );
int fscanfl( FILE *fp, const char *fmt, ... );
int utlReadNextRecord( FILE *fp, char *record );
int utlReadNextRecord( FILE *fp, char *record, int &line );
int utlReadNextDataRecord( FILE *fp, char *record );
int utlReadNextDataRecord( FILE *fp, char *record, int &line );
int utlFileFindRecord( const char *fname, const char *pattern, int skip_comments, char *record );
int utlFileParceToString( FILE *fp, const char *pattern );
int utlGetNumberOfRecords( const char *fname );
int utlGetNumberOfDataRecords( const char *fname );
int utlFileReadOption( const char *fname, const char *option_name, char *contents );
//int utlFileReadOption( const char *fname, const char *option_name );
//double utlFileReadOption( const char *fname, const char *option_name );
char *utlFileChangeExtension( const char *fname, const char *newext );
char *utlFileRemoveExtension( const char *fname );
char *utlFileAddExtension( const char *fname, const char *addext );
void utlFileSplitName( char *fullname, char *fpath, char *fname );
int utlReadXYdata( const char *fname, double **x, double **y );


// ***   MATH & GEOMETRY   ***
double utlRandom( int& seed );
double utlRandom( double v_from, double v_to );
int utlRandom( int i_from, int i_to );
double utlNormal( double avg, double stdev );
int utlPointInPolygon( int nvrtx, double *xvrtx, double *yvrtx, double x, double y );
double utlDistPoint2Line( double x0, double y0, double x1, double y1, double x2, double y2 );
double utlSphereDistance( double lon1, double lat1, double lon2, double lat2 );
double utlSphereStrike( double lon1, double lat1, double lon2, double lat2 );


// ***  Time & Time labels   ***
char *utlCurrentTime();
void utlSplitElapsedTime( double ctime, int& nHour, int& nMin, int& nSec, int& nMsec );
char *utlFormatElapsedTime( double ctime );
char *utlFormatElapsedTime( int ctime );
char *utlTimeStamp( const char *fmt );

#endif  // UTIL_H
