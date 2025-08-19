/* Class to read HySea ver 3.8.2 netCDF output with ssh-timeseries encoded as short integers
and to convert it into: (1) compact output format aka-Tolkova and (2) calculate and write netCDF with dominant wave  parameters like period and polarity for later hazard calculations.
*/


#ifndef NCHYSEA_H
#define NCHYSEA_H

#define VER_HYSEA_381   "3.8.1"      // eta encoded as floats
#define VER_HYSEA_382   "3.8.2"      // eta encoded as short integers
#define VER_HYSEA_383   "3.8.3"      // added by roberto t.
#define VER_HYSEA_384   "3.8.4"      // added by roberto t.
#define VER_HYSEA_390   "3.9.0"      // added by roberto t.
#define VER_HYSEA_400MC "4.0.0 MC"   // added by sylfest g.
#define VER381 381
#define VER382 382

class cNcHySea
{

public:
  int version;      // HySea output version: only "3.8.1" (=1) or "3.8.2" (=2) allowed
  char *fname;      // name of the (input) file; tracked to create output file names
  int ntime;        // number of timesteps (length of 'time' dimension in HySea output netCDF)
  int npoi;         // number of POIs (length of 'grid_npoints' dimension in HySea output netCDF)
  float ssh_scale;  // ->
  float ssh_offset; // To convert 2-byte short ssh into floats: ssh_f = scale*ssh_short + offset
  short *sshs;      // Time series encoded by short integers: matrix of ntime rows X npoi columns  (HySea Version 3.8.2)
  float *sshf;      // Time series encoded by floats: matrix of ntime rows X npoi columns  (HySea Version 3.8.1)
  float *depth;     // POI depth in [m], dim=npois
  float *tstamp;    // Time vector [sec], dim=ntime
  float *sshmin;    // Minimal wave height [m], dim=npois
  float *sshmax;    // Maximal wave height [m], dim=npois
  double *lon;      // Longitude of POI [deg], dim=npois
  double *lat;      // Latitude of POI [deg], dim=npois

  cNcHySea();
  ~cNcHySea();
  int read_hysea( const char *ncfname );           // Read HySea netCDF output ver=VER_HYSEA and store data in class variables
  int write_hysea( const char *ncfname );          // Write class data into netCDF according to the HySEA ver=VER_HYSEA format
  int get_mareogram( int poi, float *sshPoi );   // Pick a single mareogram from the full 2D storage (ssh) for a given POI
  int get_mareogram( int poi, double *sshPoi );
    int get_snapshot( float t0, float *sshTime );  // Pick a single time snapshot from the full 2D storage (ssh)

  int write_tolkova();        // Compact timeseries (ssh) using Tolkova-aka algorithm and write them to a binary file
  int read_tolkova( const char *fname );         // Read compacted timeseries from binary file and store them into the class variables

  int write_offshore( const char *outFmt );         // Analyse mareograms for their dominant characteristics (period, polarity) and write into the specific netCDF format
  int analyse_offshore( int poi );  // offshore analysis for one poi
};

#endif  // NCHYSEA_H

