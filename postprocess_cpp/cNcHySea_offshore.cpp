#include <stdio.h>
#include <string.h>
#include <netcdf.h>
#include <utilits.h>
#include <cMrg.h>
#include "cNcHySea.h"

#define ERRSYS(err) {Err.post("%s", nc_strerror(err)); return(-1);}

//=========================================================================
// Compute and write dominant wave information for hazard calculations:
// maximum, dominant period, polarity, offset, etc.
int cNcHySea::write_offshore( const char *outFmt )
{
  int ierr(0);
  char buf[MAX_NC_NAME];

  cMrg mrg( ntime );
  for( int i=0; i<ntime; i++ )
    mrg.x[i] = (double)tstamp[i];
  mrg.set_arr_threshold( ssh_scale );
  mrg.set_precision( 0.9*ssh_scale );

  // First make necessary calculations incluing mareogram analysis
  float *offset = new float [npoi];
  float *max_ssh = new float [npoi];
  float *min_ssh = new float [npoi];
  unsigned short *period = new unsigned short [npoi];
  signed char *polarity = new signed char [npoi];
  signed char *ret = new signed char [npoi];

  for( int n=0; n<npoi; n++ ) {

    ierr = get_mareogram( n, mrg.y ); if(ierr) return Err.post("Cannot retrieve mareogram for poi %d", n);

    if( mrg.debug_level == 1 ) { printf( "#### n= %d\n", n+1 ); }

    ierr = mrg.get_dominant();

    if( mrg.debug_level == 1 ) { printf( "Return code= %d\n", ierr ); }

    ret[n] = (signed char)ierr;
    offset[n] = mrg(0);
    max_ssh[n] = sshmax[n];
    min_ssh[n] = sshmin[n];
    period[n] = (unsigned short)mrg.period;
    polarity[n] = (signed char)mrg.polarity;
  }

  if( strstr( outFmt, "txt" ) ) {
    sprintf( buf, "%s.offshore.txt", fname );
    FILE *fp_lst = fopen( buf, "wt" );
    fprintf( fp_lst, "ID lon lat depth max_ssh min_ssh period polarity return_code\n" );
    for( int n=0; n<npoi; n++ ) {
      //fprintf( fp_lst, "%5d %+8.4f %+8.4f %3.0f %.3f %.3f %5d %2d %1d\n", n+1, lon[n], lat[n], depth[n], max_ssh[n], min_ssh[n], period[n], polarity[n], ret[n] );
      fprintf( fp_lst, "%5d\t%+8.4f\t%+8.4f\t%3.0f\t%.3f\t%.3f\t%5d\t%2d\t%1d\n", n+1, lon[n], lat[n], depth[n], max_ssh[n], min_ssh[n], period[n], polarity[n], ret[n] );
    }
    fclose( fp_lst );
  }


  if( strstr( outFmt, "netcdf" ) ) {

    // Create netCDF file
    int ncid(0);

    sprintf( buf, "%s.offshore.nc", fname );
    ierr = nc_create( buf, NC_CLOBBER|NC_NETCDF4, &ncid ); if(ierr) ERRSYS(ierr);

    // Define dimensions
    int dimid_poi;

    ierr = nc_def_dim( ncid, "pois", npoi, &dimid_poi ); if(ierr) ERRSYS(ierr);

    // Define variables and attributes
    int varid_lon,varid_lat,varid_depth,varid_offset,varid_max,varid_min,varid_period,varid_polarity,varid_ret;

    ierr = nc_def_var( ncid, "longitude", NC_DOUBLE, 1, &dimid_poi, &varid_lon ); if(ierr) ERRSYS(ierr);
    ierr = nc_put_att_text( ncid, varid_lon, "units", 7, "degrees" ); if(ierr) ERRSYS(ierr);
    ierr = nc_def_var( ncid, "latitude", NC_DOUBLE, 1, &dimid_poi, &varid_lat ); if(ierr) ERRSYS(ierr);
    ierr = nc_put_att_text( ncid, varid_lat, "units", 7, "degrees" ); if(ierr) ERRSYS(ierr);
    ierr = nc_def_var( ncid, "depth", NC_FLOAT, 1, &dimid_poi, &varid_depth ); if(ierr) ERRSYS(ierr);
    ierr = nc_put_att_text( ncid, varid_depth, "units", 6, "meters" ); if(ierr) ERRSYS(ierr);
    ierr = nc_def_var( ncid, "offset", NC_FLOAT, 1, &dimid_poi, &varid_offset ); if(ierr) ERRSYS(ierr);
    ierr = nc_put_att_text( ncid, varid_offset, "units", 6, "meters" ); if(ierr) ERRSYS(ierr);
    ierr = nc_def_var( ncid, "max_height", NC_FLOAT, 1, &dimid_poi, &varid_max ); if(ierr) ERRSYS(ierr);
    ierr = nc_put_att_text( ncid, varid_max, "units", 6, "meters" ); if(ierr) ERRSYS(ierr);
    ierr = nc_def_var( ncid, "min_height", NC_FLOAT, 1, &dimid_poi, &varid_min ); if(ierr) ERRSYS(ierr);
    ierr = nc_put_att_text( ncid, varid_min, "units", 6, "meters" ); if(ierr) ERRSYS(ierr);
    ierr = nc_def_var( ncid, "period", NC_USHORT, 1, &dimid_poi, &varid_period ); if(ierr) ERRSYS(ierr);
    ierr = nc_put_att_text( ncid, varid_period, "units", 7, "seconds" ); if(ierr) ERRSYS(ierr);
    ierr = nc_def_var( ncid, "polarity", NC_BYTE, 1, &dimid_poi, &varid_polarity ); if(ierr) ERRSYS(ierr);
    ierr = nc_put_att_text( ncid, varid_polarity, "values", 7, "-1 or 1" ); if(ierr) ERRSYS(ierr);
    ierr = nc_def_var( ncid, "return_code", NC_BYTE, 1, &dimid_poi, &varid_ret ); if(ierr) ERRSYS(ierr);
    ierr = nc_put_att_text( ncid, varid_ret, "values", 3, "1-4" ); if(ierr) ERRSYS(ierr);

    // Global attributes
    ierr = nc_put_att_text( ncid, NC_GLOBAL, "Original HySEA output", strlen(fname), fname ); if(ierr) ERRSYS(ierr);

    // End definitions
    ierr = nc_enddef( ncid ); if(ierr) ERRSYS(ierr);

    // Write down variables
    ierr = nc_put_var_double( ncid, varid_lon, lon ); if(ierr) ERRSYS(ierr);
    ierr = nc_put_var_double( ncid, varid_lat, lat ); if(ierr) ERRSYS(ierr);
    ierr = nc_put_var_float( ncid, varid_depth, depth ); if(ierr) ERRSYS(ierr);
    ierr = nc_put_var_float( ncid, varid_offset, offset ); if(ierr) ERRSYS(ierr);
    ierr = nc_put_var_float( ncid, varid_max, max_ssh ); if(ierr) ERRSYS(ierr);
    ierr = nc_put_var_float( ncid, varid_min, min_ssh ); if(ierr) ERRSYS(ierr);
    ierr = nc_put_var_ushort( ncid, varid_period, period ); if(ierr) ERRSYS(ierr);
    ierr = nc_put_var_schar( ncid, varid_polarity, polarity ); if(ierr) ERRSYS(ierr);
    ierr = nc_put_var_schar( ncid, varid_ret, ret ); if(ierr) ERRSYS(ierr);

    ierr = nc_close( ncid ); if(ierr) ERRSYS(ierr);
  }

  delete [] offset; delete [] max_ssh; delete [] min_ssh; delete [] period; delete [] polarity; delete [] ret;

  return 0;
}

//=========================================================================
// Compute dominant wave information for one POI
int cNcHySea::analyse_offshore( int poi )
{
  int ierr;

  cMrg mrg( ntime );
  for( int i; i<ntime; i++ )
    mrg.x[i] = (double)tstamp[i];
  mrg.set_arr_threshold( ssh_scale );
  mrg.set_precision( 0.9*ssh_scale );
  mrg.set_debug_level( 2 );

  ierr = get_mareogram( poi, mrg.y ); if(ierr) return Err.post("Cannot retrieve mareogram for poi %d", poi);

  if( mrg.debug_level > 0 ) { printf( "#### POI= %d\n", poi+1 ); }

  ierr = mrg.get_dominant();

  if( mrg.debug_level > 0 ) { printf( "Return code= %d\n", ierr ); }

  return 0;
}
