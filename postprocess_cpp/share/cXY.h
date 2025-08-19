#ifndef XY_H
#define XY_H


class cXY
{

public:
  int nx;
  double *x;
  double *y;

  cXY();
  cXY( cXY& );
  cXY( int nx0 );
  cXY( int nx0, double* x0 );
  cXY( int nx0, double* x0, double* y0 );
  cXY( double xmin0, double xmax0, int nx0 );
  ~cXY();

  int read( char *xyfile, int ycol );
  int read( char *xyfile );
  cXY& operator= ( const cXY& xy );
  cXY& operator= ( double val );
  cXY& operator= ( double* val );
  cXY& operator= ( float* val );
  cXY& operator*= ( double coeff );
  cXY& operator+= ( cXY& xy );
  cXY& operator+= ( double shift );
  double& operator() ( int i );
  void reset( int n, double *x, double *y );
  void resetY();
  double minX();
  double maxX();
  double min();
  double min( int i1, int i2 );
  double minNext( int ks, int& kmin );
  double minPrev( int ks, int& kmin );
  int posMin();
  int posMin( int i1, int i2 );
  double max();
  double max( int i1, int i2 );
  double max( int i1, int i2, int& kmax );
  double maxNear( double x0, int& kmax );
  int posMax();
  int posMax( int i1, int i2 );
  double maxAbs();
  int whichExtrema( int kExtr );
  double getY( double xs );
  void shift( double xshift, double yshift );
  void scale( double xscale, double yscale );
  void resample( int new_nx, double* new_x );
  void resample( int newnx );
  void resample( cXY& newxy );
  int get1st( double thresh );
  double getAverage();
  int cntZero( int k1, int k2 );
  int analyse( int k1, int k2, double thresh, int& k1st, int& nmins, int *kmins, int& kmin, int& nmaxs, int *kmaxs, int& kmax, int& nzeros, int *kzeros );
  int write( const char *fname );
  int write( int k1, int k2, const char *fname );
  int write( const char *xyfile, const char *header );

};

#endif  // XY_H
