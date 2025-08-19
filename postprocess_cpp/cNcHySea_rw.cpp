/* Class to read HySea version 3.8.2 netCDF output with ssh-timeseries encoded as short integers
 * Imposed file structure is listed below:
netcdf grid30_ts {
dimensions:
	grid_npoints = ;
	time = UNLIMITED ; // (961 currently)
variables:
	float deformed_bathy(grid_npoints) ;
	float time(time) ;
	double longitude(grid_npoints) ;
	double latitude(grid_npoints) ;
	float min_height(grid_npoints) ;
	float max_height(grid_npoints) ;
	short eta(time, grid_npoints) ;
		eta:scale_factor = 0.004577707f ;
		eta:add_offset = 0.002288853f ;
// global attributes:
*/

#include <stdio.h>
#include <stdlib.h>  //added by roberto t.
#include <string.h>
#include <netcdf.h>
#include <utilits.h>
#include "cNcHySea.h"
#define ERRSYS(err) {Err.post("%s", nc_strerror(err)); return(-1);}
#define ERRFMT() {Err.post("Unexpected netCDF HySea file format. Should be 3.8.1 or 3.8.2"); return(-1);}

//=========================================================================
// Constructor
cNcHySea::cNcHySea()
{
  version = 0;
  fname = NULL;
  ntime = npoi = 0;
  ssh_scale = ssh_offset = 0.;
  sshf = NULL; sshs = NULL;
  depth = tstamp = sshmin = sshmax = NULL;
  lon = lat = NULL;
}


//=========================================================================
// Destructor
cNcHySea::~cNcHySea()
{
  if( fname != NULL ) delete [] fname;
  if( sshf != NULL ) delete [] sshf;
  if( sshs != NULL ) delete [] sshs;
  if( tstamp != NULL ) delete [] tstamp;
  if( depth != NULL ) delete [] depth;
  if( sshmin != NULL ) delete [] sshmin;
  if( sshmax != NULL ) delete [] sshmax;
  if( lon != NULL ) delete [] lon;
  if( lat != NULL ) delete [] lat;
}


//=========================================================================
// Open NetCDF file, check it and read arrays
int cNcHySea::read_hysea( const char *ncfname )
{
  char buf[NC_MAX_NAME];
  int ierr;
  size_t length;

  // GENERAL INFO
  int ncid,ncDims,ncVars;

  fname = new char [strlen(ncfname)+1];
  strcpy( fname, ncfname );

  ierr = nc_open( fname, NC_NOWRITE, &ncid ); if(ierr) ERRSYS(ierr);

  // First of all check the HySea version by reading the corresponding global attribute
  // Should be 3.8.1 or 3.8.2
  memset( buf, 0, NC_MAX_NAME );
  ierr = nc_get_att_text( ncid, NC_GLOBAL, "Tsunami-HySEA_version", buf ); if(ierr) ERRSYS(ierr);
  if( !strcmp( buf, VER_HYSEA_381 ) )
    version = VER381;
  else if( !strcmp( buf, VER_HYSEA_382 ) )
    version = VER382;
  else if( !strcmp( buf, VER_HYSEA_383 ) )    // added by roberto t.
    version = VER382;                         // added by roberto t.
  else if( !strcmp( buf, VER_HYSEA_384 ) )    // added by roberto t.
    version = VER381;                         // added by roberto t.
  else if( !strcmp( buf, VER_HYSEA_390 ) )    // added by roberto t.
    version = VER381;                         // added by roberto t.
  else if( !strcmp( buf, VER_HYSEA_400MC ) )  // added by sylfest g.
    version = VER381;                         // added by sylfest g.
  else
    {Err.post("Valid only for HySea versions %s or %s!", VER_HYSEA_381, VER_HYSEA_382); exit(-1);}

  ierr = nc_inq( ncid, &ncDims, &ncVars, NULL, NULL ); if(ierr) ERRSYS(ierr);
  //if( ncDims != 2 || ncVars != 7 ) ERRFMT();


  // DIMENSIONS
  int dimid_poi, dimid_time;

  // Dim 0: "grid_npoints"
  ierr = nc_inq_dimid( ncid, "grid_npoints", &dimid_poi ); if(ierr) ERRFMT();
  ierr = nc_inq_dimlen( ncid, dimid_poi, &length ); if(ierr) ERRSYS(ierr);
  npoi = length;

  // Dim 1: "time"
  ierr = nc_inq_dimid( ncid, "time", &dimid_time ); if(ierr) ERRFMT();
  ierr = nc_inq_dimlen( ncid, dimid_time, &length ); if(ierr) ERRSYS(ierr);
  ntime = length;


  // VARIABLES
  int typeVar,ndimVar,varid,dimidVar[2];
  nc_type att_type;

  // Var 0: "deformed_bathy", float #0
  ierr = nc_inq_varid( ncid, "deformed_bathy", &varid ); if(ierr) ERRFMT();
  ierr = nc_inq_var( ncid, varid, NULL, &typeVar, &ndimVar, dimidVar, NULL ); if(ierr) ERRSYS(ierr);
  if( typeVar != NC_FLOAT ) ERRFMT();
  if( ndimVar != 1 ) ERRFMT();
  if( dimidVar[0] != dimid_poi ) ERRFMT();
  depth = new float [npoi];
  ierr = nc_get_var_float( ncid, varid, depth ); if(ierr) ERRSYS(ierr);

  // Var 1: "time", float #1
  ierr = nc_inq_varid( ncid, "time", &varid ); if(ierr) ERRFMT();
  ierr = nc_inq_var( ncid, varid, NULL, &typeVar, &ndimVar, dimidVar, NULL ); if(ierr) ERRSYS(ierr);
  if( typeVar != NC_FLOAT ) ERRFMT();
  if( ndimVar != 1 ) ERRFMT();
  if( dimidVar[0] != dimid_time ) ERRFMT();
  tstamp = new float [ntime];
  ierr = nc_get_var_float( ncid, varid, tstamp ); if(ierr) ERRSYS(ierr);

  // Var 2: "longitude", double #0
  ierr = nc_inq_varid( ncid, "longitude", &varid ); if(ierr) ERRFMT();
  ierr = nc_inq_var( ncid, varid, NULL, &typeVar, &ndimVar, dimidVar, NULL ); if(ierr) ERRSYS(ierr);
  if( typeVar != NC_DOUBLE ) ERRFMT();
  if( ndimVar != 1 ) ERRFMT();
  if( dimidVar[0] != dimid_poi ) ERRFMT();
  lon = new double [npoi];
  ierr = nc_get_var_double( ncid, varid, lon ); if(ierr) ERRSYS(ierr);

  // Var 3: "latitude", double #0
  ierr = nc_inq_varid( ncid, "latitude", &varid ); if(ierr) ERRFMT();
  ierr = nc_inq_var( ncid, varid, NULL, &typeVar, &ndimVar, dimidVar, NULL ); if(ierr) ERRSYS(ierr);
  if( typeVar != NC_DOUBLE ) ERRFMT();
  if( ndimVar != 1 ) ERRFMT();
  if( dimidVar[0] != dimid_poi ) ERRFMT();
  lat = new double [npoi];
  ierr = nc_get_var_double( ncid, varid, lat ); if(ierr) ERRSYS(ierr);

  // Var 4: "min_height", float #0
  ierr = nc_inq_varid( ncid, "min_height", &varid ); if(ierr) ERRFMT();
  ierr = nc_inq_var( ncid, varid, NULL, &typeVar, &ndimVar, dimidVar, NULL ); if(ierr) ERRSYS(ierr);
  if( typeVar != NC_FLOAT ) ERRFMT();
  if( ndimVar != 1 ) ERRFMT();
  if( dimidVar[0] != dimid_poi ) ERRFMT();
  sshmin = new float [npoi];
  ierr = nc_get_var_float( ncid, varid, sshmin ); if(ierr) ERRSYS(ierr);

  // Var 5: "max_height", float #0
  ierr = nc_inq_varid( ncid, "max_height", &varid ); if(ierr) ERRFMT();
  ierr = nc_inq_var( ncid, varid, NULL, &typeVar, &ndimVar, dimidVar, NULL ); if(ierr) ERRSYS(ierr);
  if( typeVar != NC_FLOAT ) ERRFMT();
  if( ndimVar != 1 ) ERRFMT();
  if( dimidVar[0] != dimid_poi ) ERRFMT();
  sshmax = new float [npoi];
  ierr = nc_get_var_float( ncid, varid, sshmax ); if(ierr) ERRSYS(ierr);

  // Var 6: "eta", either float (ver 3.8.1), or short (ver 3.8.2) #1-#0
  // Read all time series at once into one array. Values are stored by POIs then by time, i.e., POIs vary faster
  ierr = nc_inq_varid( ncid, "eta", &varid ); if(ierr) ERRFMT();
  ierr = nc_inq_var( ncid, varid, NULL, &typeVar, &ndimVar, dimidVar, NULL ); if(ierr) ERRSYS(ierr);
  if( ndimVar != 2 ) ERRFMT();
  if( dimidVar[0] != dimid_time || dimidVar[1] != dimid_poi ) ERRFMT();
  if( version == VER381 ) { // as float
    if( typeVar != NC_FLOAT ) ERRFMT();
    sshf = new float [npoi*ntime];
    ierr = nc_get_var_float( ncid, varid, sshf ); if(ierr) ERRSYS(ierr);
  }
  else if( version == VER382 ) { // as short
    if( typeVar != NC_SHORT ) ERRFMT();
    sshs = new short [npoi*ntime];
    ierr = nc_get_var_short( ncid, varid, sshs ); if(ierr) ERRSYS(ierr);
    // also read the two attributes
    ierr = nc_inq_att( ncid, varid, "scale_factor", &att_type, &length ); if(ierr) ERRFMT();
    if( att_type != NC_FLOAT ) ERRFMT();
    if( length != 1 ) ERRFMT();
    ierr = nc_get_att_float( ncid, varid, "scale_factor", &ssh_scale ); if(ierr) ERRSYS(ierr);
    ierr = nc_inq_att( ncid, varid, "add_offset", &att_type, &length); if(ierr) ERRFMT();
    if( att_type != NC_FLOAT ) ERRFMT();
    if( length != 1 ) ERRFMT();
    ierr = nc_get_att_float( ncid, varid, "add_offset", &ssh_offset ); if(ierr) ERRSYS(ierr);
  }
  else
    ERRFMT();

  ierr = nc_close( ncid ); if(ierr) ERRSYS(ierr);

  return 0;
}


//=========================================================================
// Write output NetCDF file compatible with HySea output
int cNcHySea::write_hysea( const char *ncfname )
{
  char buf[NC_MAX_NAME];
  int ierr;


  // Create netCDF file
  int ncid;
  ierr = nc_create( ncfname, NC_CLOBBER|NC_NETCDF4, &ncid ); if(ierr) ERRSYS(ierr);


	// Define dimensions
	int dimid_poi,dimid_time,dimidVar[2];

	ierr = nc_def_dim( ncid, "grid_npoints", npoi, &dimid_poi ); if(ierr) ERRSYS(ierr);
	ierr = nc_def_dim( ncid, "time", ntime, &dimid_time); if(ierr) ERRSYS(ierr);


	// Define variables and attributes
	// Note: not all HySea 3.8.2 variables but only timeseries-related
	int varid_time, varid_lon, varid_lat, varid_eta;

  ierr = nc_def_var( ncid, "time", NC_FLOAT, 1, &dimid_time, &varid_time ); if(ierr) ERRSYS(ierr);
  ierr = nc_put_att_text( ncid, varid_time, "units", 23, "seconds simulation time" ); if(ierr) ERRSYS(ierr);
  ierr = nc_def_var( ncid, "longitude", NC_DOUBLE, 1, &dimid_poi, &varid_lon ); if(ierr) ERRSYS(ierr);
  ierr = nc_put_att_text( ncid, varid_lon, "units", 7, "degrees" ); if(ierr) ERRSYS(ierr);
  ierr = nc_def_var( ncid, "latitude", NC_DOUBLE, 1, &dimid_poi, &varid_lat ); if(ierr) ERRSYS(ierr);
  ierr = nc_put_att_text( ncid, varid_lat, "units", 7, "degrees" ); if(ierr) ERRSYS(ierr);
	dimidVar[0] = dimid_time;
	dimidVar[1] = dimid_poi;
	if( version == VER381 ) {
    ierr = nc_def_var( ncid, "eta", NC_FLOAT, 2, dimidVar, &varid_eta ); if(ierr) ERRSYS(ierr);
    ierr = nc_def_var_deflate( ncid, varid_eta, NC_SHUFFLE, 1, 5 ); if(ierr) ERRSYS(ierr);
    ierr = nc_put_att_text( ncid, varid_eta, "units", 6, "meters" ); if(ierr) ERRSYS(ierr);
	}
	else if( version == VER382 ) {
    ierr = nc_def_var( ncid, "eta", NC_SHORT, 2, dimidVar, &varid_eta ); if(ierr) ERRSYS(ierr);
    ierr = nc_def_var_deflate( ncid, varid_eta, NC_SHUFFLE, 1, 5 ); if(ierr) ERRSYS(ierr);
    ierr = nc_put_att_text( ncid, varid_eta, "units", 20, "meters after scaling" ); if(ierr) ERRSYS(ierr);
    ierr = nc_put_att_float( ncid, varid_eta, "scale_factor", NC_FLOAT, 1, &ssh_scale ); if(ierr) ERRSYS(ierr);
    ierr = nc_put_att_float( ncid, varid_eta, "add_offset", NC_FLOAT, 1, &ssh_offset ); if(ierr) ERRSYS(ierr);
  }

  // Global attributes
  if( version == VER381 ) {
    ierr = nc_put_att_text( ncid, NC_GLOBAL, "Compatible_HySEA_version", strlen(VER_HYSEA_381), VER_HYSEA_381 ); if(ierr) ERRSYS(ierr);
  }
  else if( version == VER382 ) {
    ierr = nc_put_att_text( ncid, NC_GLOBAL, "Compatible_HySEA_version", strlen(VER_HYSEA_382), VER_HYSEA_382 ); if(ierr) ERRSYS(ierr);
  }

  sprintf( buf, "selected variables only" );
  ierr = nc_put_att_text( ncid, NC_GLOBAL, "comments", strlen(buf), buf ); if(ierr) ERRSYS(ierr);

  // End definitions
  ierr = nc_enddef( ncid ); if(ierr) ERRSYS(ierr);

  // Write down variables
  ierr = nc_put_var_float( ncid, varid_time, tstamp ); if(ierr) ERRSYS(ierr);
  ierr = nc_put_var_double( ncid, varid_lon, lon ); if(ierr) ERRSYS(ierr);
  ierr = nc_put_var_double( ncid, varid_lat, lat ); if(ierr) ERRSYS(ierr);
  if( version == VER381 ) {
    ierr = nc_put_var_float( ncid, varid_eta, sshf ); if(ierr) ERRSYS(ierr);
  }
  else if( version == VER382 ) {
    ierr = nc_put_var_short( ncid, varid_eta, sshs ); if(ierr) ERRSYS(ierr);
  }

  ierr = nc_close( ncid ); if(ierr) ERRSYS(ierr);

  return 0;
}



//=========================================================================
// Extract time series for a given poi
int cNcHySea::get_mareogram( int poi, float *sshPoi )
{

  if( version == VER381 ) {
    for( int k=0; k<ntime; k++ )
      sshPoi[k] = sshf[npoi*k+poi];
  }
  else if( version == VER382 ) {
    for( int k=0; k<ntime; k++ )
      sshPoi[k] = sshs[npoi*k+poi]*ssh_scale + ssh_offset;
  }

  return 0;
}


//=========================================================================
// Extract time series for a given poi
int cNcHySea::get_mareogram( int poi, double *sshPoi )
{

  if( version == VER381 ) {
    for( int k=0; k<ntime; k++ )
      sshPoi[k] = (double)sshf[npoi*k+poi];
  }
  else if( version == VER382 ) {
    for( int k=0; k<ntime; k++ )
      sshPoi[k] = (double)( sshs[npoi*k+poi]*ssh_scale + ssh_offset );
  }

  return 0;
}


//=========================================================================
// Extract time slice for all pois
int cNcHySea::get_snapshot( float t0, float *sshTime )
{

  int k0=0; float tdiff_min=1.e15;
  for( int k=0; k<ntime; k++ ) {
    float tdiff = fabs(tstamp[k]-t0);
    if( tdiff < tdiff_min ) {tdiff_min = tdiff; k0 = k;}
  }

  if( version == VER381 ) {
    for( int n=0; n<npoi; n++ )
      sshTime[n] = sshf[npoi*k0+n];
  }
  else if( version == VER382 ) {
    for( int n=0; n<npoi; n++ )
      sshTime[n] = sshs[npoi*k0+n]*ssh_scale + ssh_offset;
  }

  return 0;
}
