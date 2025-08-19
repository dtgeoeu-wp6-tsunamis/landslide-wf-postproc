#ifndef MRG_H
#define MRG_H

#include <cXY.h>

#define DEFAULT_ARR_THRESHOLD 1.e-5  // in meters
#define DEFAULT_LOWESS 3.5
#define DEFAULT_PRECISION 1.e-5  // in meters

class cMrg : public cXY
{

public:
  int debug_level;   // level of debugging. Currently 0,1, or 2.
  int id;   // eventually ID of the time series (e.g., POI number)
  int k1st;  // position of first arrival
  int polarity; // sign of the first significant arrival (analysed by applying lowess filering)
  int nmins; // number of local minima
  int *kmins; // positions of local minima
  int nminabs; // which min is the deepest
  int kmin; // position of absolute minimum
  int nmaxs; // number of local maxima
  int *kmaxs; // positions of local maxima
  int nmaxabs; // which max is the highest
  int kmax; // position of absolute maximum
  int nzeros; // number of zero-crossings
  int *kzeros; // positions of zero-crossings
  double period; // dominant period (analysed by applying lowess filering)
  double par_arr_threshold;
  double par_precision;
  double par_lowess;

  cMrg();
  cMrg( int nx0 );
  cMrg( int nx0, double* x0, double* y0 );
  ~cMrg();

  void set_debug_level( int new_debug_level );
  void set_arr_threshold( double new_threshold );
  void set_lowess( double new_lowess );
  void set_precision( double new_precision );
  int get_1st_arrival();
  int get_stats( int k1, int k2 );
  int get_stats();
  int get_dominant();

};


#endif  // MRG_H
