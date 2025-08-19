#include <stdio.h>
#include <string.h>
#include <math.h>
#include "utilits.h"
#include "cXY.h"



//=========================================================================
// Zero Constructor
cXY::cXY( )
{
  nx = 0;
  x = y = NULL;
}


//=========================================================================
// Destructor
cXY::~cXY()
{
  if( x != NULL ) delete [] x;
  if( y != NULL ) delete [] y;
}


//=========================================================================
// copy constructor similar as operator =
cXY::cXY( cXY& xy )
{
  nx = xy.nx;
  x = new double[nx]; memcpy( x, xy.x, nx*sizeof(double) );
  y = new double[nx]; memcpy( y, xy.y, nx*sizeof(double) );
}


//=========================================================================
// Constructor with initialization
cXY::cXY( int nx0 )
{
  nx = nx0;
  x = new double[nx];
  y = new double[nx];

  for( int k=0; k<nx; k++ ) {
    x[k] = y[k] = nan("");
  }
}


//=========================================================================
// Another constructor with initialization by x-vector
cXY::cXY( int nx0, double* x0 )
{
  nx = nx0;
  x = new double[nx]; memcpy( x, x0, nx*sizeof(double) );
  y = new double[nx];
  for( int k=0; k<nx; k++ ) y[k] = nan("");
}


//=========================================================================
// Another constructor with initialization by data arrays
cXY::cXY( int nx0, double* x0, double* y0 )
{
  nx = nx0;
  x = new double[nx]; memcpy( x, x0, nx*sizeof(double) );
  y = new double[nx]; memcpy( y, y0, nx*sizeof(double) );
}


//=========================================================================
// Another constructor with initialization of x-vector
cXY::cXY( double xmin0, double xmax0, int nx0 )
{
  nx = nx0;
  x = new double[nx];
  y = new double[nx];

  double dx = (xmax0-xmin0)/(nx0-1);
  for( int k=0; k<nx; k++ ) {
    x[k] = xmin0 + k*dx;
    y[k] = nan("");
    }
}


//=========================================================================
// Initialization from ycol-column of XY-file
int cXY::read( char *fname, int ycol )
{
  FILE *fp;
  char *record,*fmt;
  int k,line=0;
  double xs,ys;


  record = new char [MaxFileRecordLength];
  fmt = new char [MaxFileRecordLength];

  // Form read format string which depends on ycol
  strcpy(fmt, "%lf" );
  for( k=1; k<ycol; k++ )
    strcat( fmt, " %*s" );
  strcat( fmt, " %lf" );

  // Open plain XY-file
  if( (fp = fopen( fname, "rt" )) == NULL ) return Err.post( Err.msgOpenFile(fname) );

  // Check xy-file format, define number of points
  rewind( fp );
  for( line=0,nx=0; utlReadNextDataRecord( fp, record, line ) != EOF; ) {
    if( sscanf( record, fmt, &xs, &ys ) != 2 ) return Err.post( Err.msgReadFile( fname, line, fmt ) );
    nx++;
  }

  // Allocate memory
  if( x != NULL ) {delete [] x;} x = new double[nx];
  if( y != NULL ) {delete [] y;} y = new double[nx];

  // Read xy-values
  rewind( fp );
  line = 0;
  for( k=0; k<nx; k++ ) {
    utlReadNextDataRecord( fp, record, line );
    sscanf( record, fmt, &x[k], &y[k] );
  }

  fclose( fp );

  delete [] fmt;
  delete [] record;

  return 0;
}


//=========================================================================
// Initialization from two-column XY-file
int cXY::read( char *fname )
{
  return read( fname, 1 );
}


//=========================================================================
cXY& cXY::operator= ( const cXY& xy )
{
  nx = xy.nx;
  if( x != NULL ) {delete [] x;} x = new double[nx]; memcpy( x, xy.x, nx*sizeof(double) );
  if( y != NULL ) {delete [] y;} y = new double[nx]; memcpy( y, xy.y, nx*sizeof(double) );

  return *this;
}


//=========================================================================
cXY& cXY::operator= ( double val )
{
  for( int k=0; k<nx; k++ )
    y[k] = val;

  return *this;
}


//=========================================================================
cXY& cXY::operator= ( double* val )
{
  for( int k=0; k<nx; k++ )
    y[k] = val[k];

  return *this;
}


//=========================================================================
cXY& cXY::operator= ( float* val )
{
  for( int k=0; k<nx; k++ )
    y[k] = (double)val[k];

  return *this;
}


//=========================================================================
cXY& cXY::operator*= ( double coeff )
{
  for( int k=0; k<nx; k++ )
    y[k] *= coeff;

  return *this;
}


//=========================================================================
cXY& cXY::operator+= ( cXY& xy )
{
  int i;

  for( i=0; i<nx; i++ )
    y[i] += xy.getY( x[i] );

  return *this;
}


//=========================================================================
cXY& cXY::operator+= ( double shift )
{
  int i;

  for( i=0; i<nx; i++ )
    y[i] += shift;

  return *this;
}


//=========================================================================
double& cXY::operator() ( int i )
{
  return( y[i] );
}


//=========================================================================
// Reset Y-values
void cXY::resetY()
{
  for( int k=0; k<nx; k++)
    y[k] = nan("");
}


//=========================================================================
// Reset complete timeseries
void cXY::reset( int n, double *x0, double *y0 )
{
  nx = n;
  if( x != NULL ) {delete [] x;} x = new double[nx]; memcpy( x, x0, nx*sizeof(double) );
  if( y != NULL ) {delete [] y;} y = new double[nx]; memcpy( y, y0, nx*sizeof(double) );
}


//=========================================================================
void cXY::shift( double xshift, double yshift )
{
  int i;

  for( i=0; i<nx; i++ ) {
    x[i] += xshift;
    y[i] += yshift;
  }

  return;
}


//=========================================================================
void cXY::scale( double xscale, double yscale )
{

  for( int k=0; k<nx; k++ ) {
    x[k] *= xscale;
    y[k] *= yscale;
  }

  return;
}



//=========================================================================
// Get minimal x-value
double cXY::minX()
{
  return x[0];
}


//=========================================================================
// Get maximal x-value
double cXY::maxX()
{
  return x[nx-1];
}


//=========================================================================
// Get minimal y-value
double cXY::min()
{
  double ymin(RealMax);
  for( int k=0; k<nx; k++ )
    if( y[k] < ymin )
      ymin = y[k];

  return ymin;
}


//=========================================================================
// Get minimal y-value within given interval
double cXY::min( int k1, int k2 )
{
  if( k1 > k2 ) return nan("");

  int kbeg = k1; if( kbeg < 0 ) kbeg = 0;
  int kend = k2; if( kend > nx-1 ) kend = nx-1;

  double ymin(RealMax);
  for( int k=kbeg; k<=kend; k++ )
    if( y[k] < ymin )
      ymin = y[k];

  return ymin;
}


//=========================================================================
// Get position of minimal y-value
int cXY::posMin()
{
  int kmin(0); double ymin(RealMax);
  for( int k=0; k<nx; k++ ) {
    if( y[k] < ymin ) {
      ymin = y[k];
      kmin = k;
    }
  }

  return kmin;
}


//=========================================================================
// Get position of minimal y-value within given interval
int cXY::posMin( int k1, int k2 )
{
  int kbeg = k1; if( kbeg < 0 ) kbeg = 0;
  int kend = k2; if( kend > nx-1 ) kend = nx-1;

  int kmin(0); double ymin(RealMax);
  for( int k=kbeg; k<=kend; k++ ) {
    if( y[k] < ymin ) {
      ymin = y[k];
      kmin = k;
    }
  }

  return kmin;
}


//=========================================================================
// Next local minimum
double cXY::minNext( int ks, int& kmin )
{
  int k;

  if( ks < 0 || ks > nx-1 ) { kmin = -1; return nan(""); }

  // Going right
  for( k=ks; k<nx-1; k++ )
    if( y[k-1] > y[k] && y[k] < y[k+1] )
      break;

  if( k == nx-1 ) { kmin = -1; return nan(""); }

  kmin = k;
  return y[kmin];
}


//=========================================================================
// Previous local minimum
double cXY::minPrev( int ks, int& kmin )
{
  int k;

  if( ks < 0 || ks > nx-1 ) { kmin = -1; return nan(""); }

  // Going left
  for( k=ks; k>0; k-- )
    if( y[k-1] > y[k] && y[k] < y[k+1] )
      break;

  if( k == 0 ) { kmin = -1; return nan(""); }

  kmin = k;
  return y[kmin];
}


//=========================================================================
// Get maximal y-value
double cXY::max()
{
  double ymax(-RealMax);
  for( int k=0; k<nx; k++ )
    if( y[k] > ymax )
      ymax = y[k];

  return ymax;
}


//=========================================================================
// Get maximal y-value within given interval
double cXY::max( int k1, int k2 )
{
  if( k1 > k2 ) return nan("");

  int kbeg = k1; if( kbeg < 0 ) kbeg = 0;
  int kend = k2; if( kend > nx-1 ) kend = nx-1;

  double ymax(-RealMax);
  for( int k=kbeg; k<=kend; k++ )
    if( y[k] > ymax )
      ymax = y[k];

  return ymax;
}


//=========================================================================
// Get maximal y-value within given interval
double cXY::max( int k1, int k2, int& kmax )
{
  if( k1 > k2 ) { kmax = -1; return nan(""); }

  int kbeg = k1; if( kbeg < 0 ) kbeg = 0;
  int kend = k2; if( kend > nx-1 ) kend = nx-1;

  double ymax(-RealMax); kmax = -1;
  for( int k=kbeg; k<=kend; k++ ) {
    if( y[k] > ymax ) {
      ymax = y[k];
      kmax = k;
    }
  }

  return ymax;
}


//=========================================================================
// Get maximal absolute y-value
double cXY::maxAbs()
{
  double ymax(-RealMax);
  for( int k=0; k<nx; k++ )
    if( fabs(y[k]) > ymax )
      ymax = fabs(y[k]);

  return ymax;
}


//=========================================================================
// Get nearest maximum
double cXY::maxNear( double xs, int& kmax )
{
  int k,k0,k_l,k_r;

  if( xs < x[0] || xs > x[nx-1] ) return nan("");
  for( k0=0; xs > x[k0+1]; k0++ ) ;

  // Going left
  for( k=k0; k>0; k-- )
    if( y[k-1] < y[k] && y[k] > y[k+1] )
      break;
  k_l = k;

  // Going right
  for( k=k0; k<nx-1; k++ )
    if( y[k-1] < y[k] && y[k] > y[k+1] )
      break;
  k_r = k;

  if( k_l == 0 && k_r == nx-1 ) {
    kmax = -1;
    return nan("");
  }
  else if( k_l == 0 ) {
    kmax = k_r;
    return y[kmax];
  }
  else if( k_r == nx-1 ) {
    kmax = k_l;
    return y[kmax];
  }
  else {
    if( (x[k_r]-x[k0]) > (x[k0]-x[k_l]) )
      kmax = k_l;
    else
      kmax = k_r;
    return y[kmax];
  }

}


//=========================================================================
// Get position of maximal y-value
int cXY::posMax()
{
  int kmax(0); double ymax(-RealMax);
  for( int k=0; k<nx; k++ ) {
    if( y[k] > ymax ) {
      ymax = y[k];
      kmax = k;
    }
  }

  return kmax;
}


//=========================================================================
// Get position of maximal y-value within given interval
int cXY::posMax( int k1, int k2 )
{
  int kbeg = k1; if( kbeg < 0 ) kbeg = 0;
  int kend = k2; if( kend > nx-1 ) kend = nx-1;

  int kmax(0); double ymax(-RealMax);
  for( int k=kbeg; k<=kend; k++ ) {
    if( y[k] > ymax ) {
      ymax = y[k];
      kmax = k;
    }
  }

  return kmax;
}


//=========================================================================
// Analyse time series for extrema
int cXY::whichExtrema( int kExtr )
{
  if( kExtr == 0 || kExtr == nx-1 ) return 0;  // end node - not an extremum
  if( (y[kExtr]-y[kExtr-1])*(y[kExtr+1]-y[kExtr]) >= 0 ) return 0; // not an extremum
  if( y[kExtr] > y[kExtr-1] ) return 1;
  else return -1;
}


//=========================================================================
// Get y-value with interpolation (linear)
double cXY::getY( double xs )
{
  int i;
  double result;

  if( xs < x[0] || xs > x[nx-1] ) return nan("");

  for( i=0; xs > x[i+1]; i++ ) ;

  result = y[i] + (y[i+1]-y[i])/(x[i+1]-x[i])*(xs-x[i]);

  return( result );
}


//=========================================================================
// Resample (interpolate) y-values onto new arbitrary grid
void cXY::resample( int new_nx, double* new_x )
{
  double *new_y = new double [new_nx];

  for( int i=0; i<new_nx; i++ )
    new_y[i] = getY( new_x[i] );

  reset( new_nx, new_x, new_y );

  delete [] new_y;
}


//=========================================================================
// Resample (interpolate) y-values onto new equidistant grid
void cXY::resample( int new_nx )
{
  double dx = ( maxX()-minX() ) / (new_nx-1);
  double *new_x = new double [new_nx];
  for( int k=0; k<new_nx; k++ )
    new_x[k] = minX() + k*dx;

  resample( new_nx, new_x );

  delete [] new_x;
}


//=========================================================================
// Resample (interpolate) y-values onto new arbitrary grid
void cXY::resample( cXY& new_xy )
{
  resample( new_xy.nx, new_xy.x );
}


//=========================================================================
// write into simple XY-file
int cXY::write( const char *xyfile )
{
  FILE *fp = fopen( xyfile, "wt" );

  for( int k=0; k<nx; k++ ) {
    fprintf( fp, "%g %g\n", x[k], y[k] );
    //fprintf( fp, "%d %g\n", k, y[k] );
  }

  fclose( fp );

  return 0;
}


//=========================================================================
// write part of graph into simple XY-file
int cXY::write( int k1, int k2, const char *xyfile )
{
  FILE *fp = fopen( xyfile, "wt" );

  if( k1 >= k2 ) return -1;
  if( k1<0 ) k1 = 0;
  if( k2 > nx-1 ) k2 = nx-1;

  for( int k=k1; k<=k2; k++ ) {
    fprintf( fp, "%g %g\n", x[k], y[k] );
  }

  fclose( fp );

  return 0;
}


//=========================================================================
// write into simple XY-file
int cXY::write( const char *xyfile, const char *header )
{
  FILE *fp = fopen( xyfile, "wt" );

  fprintf( fp, "%s\n", header );

  for( int k=0; k<nx; k++ ) {
    fprintf( fp, "%g %g\n", x[k], y[k] );
  }

  fclose( fp );

  return 0;
}


//=========================================================================
// get first arrival
int cXY::get1st( double thresh )
{
  int pos1stArrival(-1);

  for( int i=0; i<nx; i++ )
    if( pos1stArrival < 0 && fabs(y[i]) > thresh ) {
      pos1stArrival = i;
      break;
    }

  return pos1stArrival;
}



//=========================================================================
// get average of non-zero values
double cXY::getAverage()
{
  double avg(0.); int navg(0);
  for( int k=0; k<nx; k++ ) {
    if( fabs(y[k]) > 1.e-10 ) {
      navg++;
      avg += y[k];
    }
  }
  avg /= navg;

  return avg;
}


//=========================================================================
// Count zero-crosses
int cXY::cntZero( int k1, int k2 )
{
  int kbeg = k1; if( kbeg < 0 ) kbeg = 0;
  int kend = k2; if( kend > nx-1 ) kend = nx-1;

  int cnt(0);
  for( int k=kbeg; k<kend; k++ )
    if( y[k]*y[k+1] < 0 ) cnt++;

  return cnt;
}


//=========================================================================
// Basic analysis for minima, maxima and zero-crossings

int cXY::analyse( int k1, int k2, double thresh, int& k1st, int& nmins, int *kmins, int& kmin, int& nmaxs, int *kmaxs, int& kmax, int& nzeros, int *kzeros )
{
  if( k1 > k2 ) return -1;

  int kbeg = k1; if( kbeg < 0 ) kbeg = 0;
  int kend = k2; if( kend > nx-1 ) kend = nx-1;
  double ymax(-RealMax), ymin(RealMax);

  nmins = nmaxs = nzeros = 0; k1st = kmin = kmax = -1;

  // Left edge
  if( fabs(y[kbeg]) > thresh ) {
    k1st = kbeg;
    if( y[kbeg+1] > y[kbeg] ) kmins[nmins++] = kbeg;
    if( y[kbeg+1] < y[kbeg] ) kmaxs[nmaxs++] = kbeg;
    if( y[kbeg]*y[kbeg+1] < 0 ) kzeros[nzeros++] = kbeg;
    kmin = kmax = kbeg; ymin = ymax = y[kbeg];
  }

  // Internal points
  for( int k=kbeg+1; k<=kend-1; k++ ) {
    if( k1st < 0 ) {
      if( fabs(y[k]) < thresh )
        continue;
      else
        k1st = kzeros[nzeros++] = k;
    }
    else
      if( y[k]*y[k+1] < 0 ) kzeros[nzeros++] = k;

    if( y[k-1] > y[k] && y[k] < y[k+1] ) kmins[nmins++] = k;
    if( y[k-1] < y[k] && y[k] > y[k+1] ) kmaxs[nmaxs++] = k;
    if( y[k] < ymin ) { ymin = y[k]; kmin = k; }
    if( y[k] > ymax ) { ymax = y[k]; kmax = k; }
  }

  // Right edge
  if( k1st < 0 && fabs(y[kend]) > thresh )
    k1st = kend;

  if( k1st >= 0 ) {
    if( y[kend] < y[kend-1] ) kmins[nmins++] = kend;
    if( y[kend] > y[kend-1] ) kmaxs[nmaxs++] = kend;
    if( y[kend] < ymin ) { ymin = y[kend]; kmin = kend; }
    if( y[kend] > ymax ) { ymax = y[kend]; kmax = kend; }
  }

  return( (k1st >= 0) ? 0 : -1 );
}

