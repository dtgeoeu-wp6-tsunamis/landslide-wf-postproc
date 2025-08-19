#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <cMrg.h>


//=========================================================================
// Zero Constructor
cMrg::cMrg() : cXY()
{
  debug_level = 0;
  id = k1st = kmin = kmax = nminabs = nmaxabs = -1;
  polarity = nmins = nmaxs = nzeros = 0;
  kmins = kmaxs = kzeros = NULL;
  period = 0.;
  par_arr_threshold = DEFAULT_ARR_THRESHOLD; par_precision = DEFAULT_PRECISION; par_lowess = DEFAULT_LOWESS;
}


//=========================================================================
// Another constructor with memory allocation
cMrg::cMrg( int nx0 ) : cXY( nx0 )
{
  debug_level = 0;
  id = k1st = kmin = kmax = -1;
  polarity = nmins = nmaxs = nzeros = 0;
  int nspec=nx0;
  kmins = new int[nspec]; kmaxs = new int[nspec]; kzeros = new int[nspec];
  for( int k=0; k<nspec; k++ ) {
    kmins[k] = kmaxs[k] = kzeros[k] = -1;
  }
  period = 0.;
  par_arr_threshold = DEFAULT_ARR_THRESHOLD; par_precision = DEFAULT_PRECISION; par_lowess = DEFAULT_LOWESS;
}


//=========================================================================
// Another constructor with initialization by data arrays
cMrg::cMrg( int nx0, double* x0, double* y0 ) : cXY( nx0, x0, y0 )
{
  debug_level = 0;
  id = k1st = kmin = kmax = -1;
  polarity = nmins = nmaxs = nzeros = 0;
  int nspec=nx0;
  kmins = new int[nspec]; kmaxs = new int[nspec]; kzeros = new int[nspec];
  for( int k=0; k<nspec; k++ ) {
    kmins[k] = kmaxs[k] = kzeros[k] = -1;
  }
  period = 0.;
  par_arr_threshold = DEFAULT_ARR_THRESHOLD; par_precision = DEFAULT_PRECISION; par_lowess = DEFAULT_LOWESS;
}


//=========================================================================
// Destructor
cMrg::~cMrg()
{
  if( kmins != NULL ) delete [] kmins;
  if( kmaxs != NULL ) delete [] kmaxs;
  if( kzeros != NULL ) delete [] kzeros;
}


//=========================================================================
// Set debug level
void cMrg::set_debug_level( int new_debug_level )
{
  debug_level = new_debug_level;
}


//=========================================================================
// Set arrival threshold
void cMrg::set_arr_threshold( double new_threshold )
{
  par_arr_threshold = new_threshold;
}

//=========================================================================
// Set parameter of the Lowess filter (this is a dynamic scaling parameter!)
void cMrg::set_lowess( double new_lowess )
{
  par_lowess = new_lowess;
}

//=========================================================================
// Set parameter of the precison (to search for minima and maxima)
void cMrg::set_precision( double new_precision )
{
  par_precision = new_precision;
}


//=========================================================================
// Get first arrival
int cMrg::get_1st_arrival()
{
  k1st = -1;
  for( int k=0; k<nx; k++ )
    if( k1st < 0 && fabs(y[k]) > par_arr_threshold ) {
      k1st = k;
      break;
    }

  if( k1st > -1 )
    return 0;
  else
    return -1;
}


//=========================================================================
int cMrg::get_stats()
{
  return get_stats( 0, nx-1 );
}

//========================================================
int compare2dbl(const void * a, const void * b)
{
  if ( *(double *)a <  *(double *)b ) return -1;
  if ( *(double *)a == *(double *)b ) return 0;
  if ( *(double *)a >  *(double *)b ) return 1;
}

//=========================================================================
// Basic analysis for mareogram statistics: first arrival, minima, maxima and zero-crossings
// Returns 0 if mareogram seems to be OK, otherwise -1
int cMrg::get_stats( int k1, int k2 )
{

  if( k1 > k2 ) return -1;

  int kbeg = k1; if( kbeg < 0 ) kbeg = 0;
  int kend = k2; if( kend > nx-1 ) kend = nx-1;
  double ssh_max(-1.e+30), ssh_min(1.e+30);

  nmins = nmaxs = nzeros = 0; k1st = kmin = kmax = -1; nminabs = nmaxabs = -1;
  period = 0.; polarity = 0;

  // Left edge
  if( fabs(y[kbeg]) > par_arr_threshold ) {
    k1st = kbeg;
    if( y[kbeg]*y[kbeg+1] < 0 ) kzeros[nzeros++] = kbeg;
    kmin = kmax = kbeg; ssh_min = ssh_max = y[kbeg];
  }

  // Internal points
  int kp_p(-1), kp_n(-1);
  for( int k=kbeg+1; k<=kend-1; k++ ) {

    if( k1st < 0 ) {
      if( fabs(y[k]) < par_arr_threshold )
        continue;
      else
        k1st = kzeros[nzeros++] = k;  // NOTE: k1st is also always considered as a 'first zero crossing'
    }

    if( y[k]*y[k+1] < 0 ) kzeros[nzeros++] = k;

    if( kp_p > 0 ) {
      if( fabs(y[k+1]-y[kp_p]) > par_precision ) {
        if( y[k+1] < y[kp_p] ) {
          kmaxs[nmaxs++] = (kp_p + k)/2;
        }
        kp_p = -1;
      }
    }
    else if( kp_n > 0 ) {
      if( fabs(y[k+1]-y[kp_n]) > par_precision ) {
        if( y[k+1] > y[kp_n] ) {
          kmins[nmins++] = (kp_n + k)/2;
        }
        kp_n = -1;
      }
    }
    else {
      if( (y[k]-y[k-1]) > par_precision && (y[k]-y[k+1]) > par_precision )
        kmaxs[nmaxs++] = k;
      else if( (y[k-1]-y[k]) > par_precision && (y[k+1]-y[k]) > par_precision )
        kmins[nmins++] = k;
      else if( (y[k]-y[k-1]) > par_precision && fabs(y[k]-y[k+1]) < par_precision )
        kp_p = k;
      else if( (y[k-1]-y[k]) > par_precision && fabs(y[k+1]-y[k]) < par_precision )
        kp_n = k;
    }

    if( y[k] < ssh_min ) { ssh_min = y[k]; kmin = k; }
    if( y[k] > ssh_max ) { ssh_max = y[k]; kmax = k; }
  }

  // Right edge
  if( k1st < 0 && fabs(y[kend]) > par_arr_threshold )
    k1st = kend;

  if( k1st > -1 ) {
    if( y[kend] < ssh_min ) { ssh_min = y[kend]; kmin = kend; }
    if( y[kend] > ssh_max ) { ssh_max = y[kend]; kmax = kend; }
  }

  if( k1st < 0 || (nmins < 2 && nmaxs < 2) ) return -1;

  double *per = new double [k2-k1+1]; memset( per, 0, sizeof(double)*(k2-k1+1) ); // [nmaxs>nmins ? nmaxs*2 : nmins*2];
  int l;

  // Now analyse maximums: find the global maximum and wave period
  double periodMax(0.);

  if( debug_level > 1 ) { printf( "Maximums: k, y, x, per, nmaxabs\n" ); printf( "%d: %d\t%g\t%g\n", 0, kmaxs[0], y[kmaxs[0]], x[kmaxs[0]] ); }

  if( nmaxs >= 2 ) {
    nmaxabs = 0;
    for( l=1; l<nmaxs; l++ ) {
      if( y[kmaxs[l]] > y[kmaxs[nmaxabs]] )
        nmaxabs = l;
      per[l-1] = x[kmaxs[l]] - x[kmaxs[l-1]];

      if( debug_level > 1 ) { printf( "%d: %d\t%.3f\t%.0f\t%.0f\t%d\n", l, kmaxs[l], y[kmaxs[l]], x[kmaxs[l]], per[l-1], nmaxabs ); }
    }
    // Take median value of periods
    qsort( (void *)per, nmaxs, sizeof(double), compare2dbl );
    if( nmaxs % 2 != 0 )
      periodMax = per[nmaxs/2];
    else
      periodMax = (per[(nmaxs - 1)/2] + per[nmaxs/2]) / 2.;
  }
  else
    nmaxabs = 0;

  if( periodMax > 0 ) period = periodMax;

  // Same for minimums
  double periodMin(0.);

  if( debug_level > 1 ) { printf( "Minimums: k, y, x, per, nminabs\n" ); printf( "%d: %d\t%g\t%g\n", 0, kmins[0], y[kmins[0]], x[kmins[0]] ); }

  if( nmins >= 2 ) {
    nminabs = 0;
    for( l=1; l<nmins; l++ ) {
      if( y[kmins[l]] < y[kmins[nminabs]] )
        nminabs = l;
      per[l-1] = x[kmins[l]] - x[kmins[l-1]];

      if( debug_level > 1 ) { printf( "%d: %d\t%.3f\t%.0f\t%.0f\t%d\n", l, kmins[l], y[kmins[l]], x[kmins[l]], per[l-1], nminabs ); }
    }
    qsort( (void *)per, nmins, sizeof(double), compare2dbl );
    if( nmins % 2 != 0 )
      periodMin = per[nmins/2];
    else
      periodMin = (per[(nmins - 1)/2] + per[nmins/2]) / 2.;
  }
  else
    nminabs = 0;

  if( periodMin > 0 ) period = (period + periodMin)/2;

  // Polarity as sign of minimum preceeding the absolute maximum
  for( l=nmins-1; l>= 0 && kmins[l] > kmaxs[nmaxabs]; l-- ) ;
  if( l >= 0 )
    polarity = (y[kmins[l]] < -0.25*y[kmaxs[nmaxabs]]) ? -1 : 1;
  else
    // if there is no minimum before the first maximum
    polarity = 1;

  delete [] per;

  return 0;
}


//=========================================================================
// Retrieving dominant mareogram parameters: period and polarity for later use with local amplification factors
void lowess( double *x, double *y, int n, double f, double *ys );

int cMrg::get_dominant()
{
  int ierr;

  int iret = 1; // return code. Signals for how successful was the mareogram analysis. Larger is better
    // 1: analysis failed
    // 2: parameters retrieved from the raw mareogram
    // 3: parameters retrieved from the filtered mareogram
    // not implemented in this version! 4: period estimated locally as difference between the two minima around the abs maximum at the filtered mareogram

  if( debug_level > 1 ) { write( "mrg-raw.dat" ); }

  // Original maregram: make analysis for minima, maxima, zero-crosses
  ierr = get_stats(); if( ierr ) return iret;  // mareogram not plausible

  if( debug_level > 0 ) { printf( "Raw: nmins: %d, nmaxs: %d, nzeros: %d, period: %0.f, polarity: %d\n", nmins, nmaxs, nzeros, period, polarity ); }

  iret = 2;  // analysis completed at the raw mareogram

  // Next stage: produce smoothed mareogram image with lowess filter

  // Define filtering window: two wave periods but at least 2 hours
  double wind = 2*period; if( wind < 7200 ) wind = 7200.;
  double tstep = (x[nx-1]-x[0])/(nx-1);
  int kw1 = kmax - (int)(wind/tstep); if( kw1 < k1st ) kw1 = k1st;
  int kw2 = kmax + (int)(wind/tstep); if( kw2 > nx-1 ) kw2 = nx-1;

  // Temporary storage for a filtered mareogram
  cMrg fltr( kw2-kw1+1, &x[kw1], &y[kw1] );

  // Run filter
  lowess( &x[kw1], &y[kw1], kw2-kw1+1, par_lowess*period/2/(x[kw2]-x[kw1]), &fltr.y[0] );

  // Remove eventual offset
  double avg = fltr.getAverage();
  fltr.shift( 0., -avg );

  if( debug_level > 1 ) { fltr.write( "mrg-flt.dat" ); }

  // Now run analysis on filtered maregram: make analysis for minima, maxima, zero-crosses
  ierr = fltr.get_stats();  if( ierr ) return iret;  // filtered mareogram not plausible

  if( debug_level > 0 ) { printf( "Fltr: nmins: %d, nmaxs: %d, nzeros: %d, period: %0.f, polarity: %d\n",fltr.nmins, fltr.nmaxs, fltr.nzeros, fltr.period, fltr.polarity ); }

  period = fltr.period;
  polarity = fltr.polarity;
  iret = 3;  // analysis completed on filtered mareogram but with 'bulk' period

  return iret;
}

/*
int get_dominant_save()
{
  int i, ierr;

  int iret = 1; // return code. Signals for how successful was the mareogram analysis. Larger is better
    // 1: analysis failed
    // 2: parameters retrieved from the raw mareogram only (failed to filter); period estimated from zero-crossings
    // 3: parameters retrieved from the filtered mareogram, but period estimated from zero-crossings only
    // 4: period estimated locally as difference between the two minima around the abs maximum at the filtered mareogram

  polarity = 0; period = 0.;

  // Original maregram: make analysis for minima, maxima, zero-crosses
  ierr = get_stats(); if( ierr ) return iret;  // 1st arrival not found for given threshold

  // If no zero-crossings, mareogram can have vertical offset. Try to remove vertical shift and repeat analysis
  if( nzeros < 2 ) {
    double avg = getAverage();
    shift( 0., -avg );
    ierr = get_stats();
    if( ierr || nzeros < 2 ) return iret;  // did not help
  }

  // Roughly evaluate period and polarity from the raw (unfiltered) mareogram
  // Period as distance between first and last zero-crosses divided by number of intervals and multiplied by 2
  period = (x[kzeros[nzeros-1]] - x[kzeros[0]]) / (nzeros-1) * 2;

  // Polarity as sign of minimum preceeding the absolute maximum
  for( i=nmins-1; i>= 0 && kmins[i] > kmax; i-- ) ;
  if( i >= 0 )
    polarity = (y[kmins[i]] < -0.25*ssh_max) ? -1 : 1;
  else
    // if there is no minimum before the first maximum
    polarity = 1;

  iret = 2;  // analysis completed at the raw mareogram

  // Next stage: produce smoothed mareogram image with lowess filter

  // Define filtering window: two wave periods but at least 2 hours
  double wind = 2*period; if( wind < 7200 ) wind = 7200.;
  double tstep = (x[nx-1]-x[0])/(nx-1);
  int kw1 = kmax - (int)(wind/tstep); if( kw1 < k1st ) kw1 = k1st;
  int kw2 = kmax + (int)(wind/tstep); if( kw2 > nx-1 ) kw2 = nx-1;

  // Temporary storage for a filtered mareogram
  cMrg fltr( kw2-kw1+1, &x[kw1], &y[kw1] );

  // Run filter
  lowess( &x[kw1], &y[kw1], kw2-kw1+1, par_lowess*period/2/(x[kw2]-x[kw1]), &fltr.y[0] );

  // Remove eventual offset
  double avg = fltr.getAverage();
  fltr.shift( 0., -avg );

  #ifdef DEBUG2
    fltr.write( "mrg-flt.dat" );
  #endif

  // Now run analysis on filtered maregram: make analysis for minima, maxima, zero-crosses
  ierr = fltr.get_stats();

  // If nor zero-crossings: exit
  if( fltr.nzeros < 2 ) return iret;

  // Another rough estimation of wave period as a distance between first and last zero-crosses divided by number of intervals, and finally multiplied by 2
  period = (fltr.x[fltr.kzeros[fltr.nzeros-1]] - fltr.x[fltr.kzeros[0]]) / (fltr.nzeros-1) * 2;

  iret = 3;  // analysis completed on filtered mareogram but with 'bulk' period

  // Find filtered maximum closest to the absolute maximum at the raw mareogram
  int knear(-1), kdist, kdistMin(nx);
  for( i=0; i<fltr.nmaxs; i++ ) {
    if( (kdist = abs(fltr.kmaxs[i]-(kmax-kw1))) < kdistMin ) {
      kdistMin = kdist;
      knear = fltr.kmaxs[i];
    }
  }
  // in case nearest filtered maximum seems not be the filtered image of the original raw maximum - exit
  if( fabs(x[kmax] - fltr.x[knear]) > period/4 ) return iret;

  // Polarity as sign of minimum preceeding this maximum
  for( i=fltr.nmins-1; i>= 0 && fltr.kmins[i] > knear; i-- ) ;
  if( i >= 0 )
    polarity = (fltr(fltr.kmins[i]) < -0.25*fltr(knear)) ? -1 : 1;
  else
    // if there is no minimum before the first maximum
    polarity = 1;

  // Finally try better period estimation as a distance between the two minima around the knear-maxima
  for( i=0; i<fltr.nmins-1; i++ ) {
    if( fltr.kmins[i] < knear && fltr.kmins[i+1] > knear ) {
      period = (fltr.kmins[i+1]-fltr.kmins[i])*tstep;
      iret = 4;  // full analysis
      break;
    }
  }

  return iret;
}
*/

//=========================================================================
//=========================================================================
//=========================================================================
// This is realization of the LOWESS (LOcally WEighted Scatterplot Smoothing) filter
/* Adopted and simplified by A.B. based on:
 *
 *  c++ implementation of Lowess weighted regression by
 *  Peter Glaus http://www.cs.man.ac.uk/~glausp/
 *
 *
 *  Based on fortran code by Cleveland downloaded from:
 *  http://netlib.org/go/lowess.f
 *  original author:
* wsc@research.bell-labs.com Mon Dec 30 16:55 EST 1985
* W. S. Cleveland
* Bell Laboratories
* Murray Hill NJ 07974
 *
 */

static inline double fcube(double x)
{
  return x * x * x;
}

static inline double fmax2(double x, double y)
{
  return (x < y) ? y : x;
}

static void plowest(double *x, double *y, int n, double xs, double& ys, int nleft, int nright, double *w, int *ok)
{
  int nrt, j;
  double a, b, c, h, h1, h9, r, range;

  range = x[n-1]-x[0];
  h = fmax2( xs-x[nleft], x[nright]-xs );
  h9 = 0.999*h; h1 = 0.001*h;

  // sum of weights
  a = 0.;
  j = nleft;
  for( j = nleft; j < n; j++ ) {
    // compute weights (pick up all ties on right)
    w[j] = 0.;
    r = fabs(x[j]-xs);
    if( r <= h9 ){
      if( r > h1 )
        w[j] = fcube(1.-fcube(r/h));
      else
        w[j] = 1.;
      a += w[j];
    }
    else if( x[j] > xs )
      break;  // get out at first zero wt on right
  }

  // rightmost pt (may be greater than nright because of ties)
  nrt = j-1;
  if( a <= 0. )
    *ok = 0;
  else {
    // weighted least squares
    *ok = 1;
    // normalize weights
    for( j = nleft; j <= nrt; j++ )
      w[j] /= a;

    if( h > 0. ) {
      // use linear fit
      a = 0.;
      for( j = nleft; j <= nrt; j++ )
        a += w[j] * x[j]; // weighted center of x values
      b = xs - a;
      c = 0.;
      for( j = nleft; j <= nrt; j++ )
        c += w[j]*(x[j]-a)*(x[j]-a);
      if( sqrt(c) > 0.001*range ) {         // points are spread out enough to compute slope
        b /= c;
        for( j = nleft; j <= nrt; j++ )
          w[j] *= (b*(x[j]-a) + 1.);
      }
    }

    ys = 0.;
    for( j = nleft; j <= nrt; j++ )
      ys += w[j] * y[j];
  }
}

void lowess( double *x, double *y, int n, double f, double *ys )
{
  int nleft, nright, ns, ok;
  double *dtmp = new double [n];

  // effective window: at least two, at most n points
  ns = (int)(f*n); if( ns > n ) ns = n; if( ns < 2 ) ns = 2;

  nleft = 0; nright = ns-1;

  for( int i = 0; i < n; i++ ) {

    while( nright < n-1 ) {
      // move nleft, nright to right while radius decreases
      if( (x[i]-x[nleft]) <= (x[nright+1]-x[i]) ) break;
      nleft++; nright++;
    }

    // fitted value at x[i]
    plowest( x, y, n, x[i], ys[i], nleft, nright, dtmp, &ok );
    if( !ok ) ys[i] = y[i];
  }

  delete [] dtmp;
  return;
}
