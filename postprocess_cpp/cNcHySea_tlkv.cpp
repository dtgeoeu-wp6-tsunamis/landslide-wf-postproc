/* Encode+write and read+decode mareograms in accord with 'Tolkova algorithm'

* Original data are to be presented as cNcHySea class variables (see 'cNcHySea.h' for description).
* At this stage timeseries (matrix 'ssh') are already presented by short integers after applying the netCDF data mapping with scale and offset.

* Compact data format is as follows:
*
* - 32 bytes: file header with version information
* - 64 bytes: reserved for scenario ID
* - 2 bytes | short: time step [s]
* - 2 bytes | short: number of time steps
* - 2 bytes | short: number of POIs (npoi)
* - 8*npoi bytes | double : array of poi longitudes
* - 8*npoi bytes | double : array of poi latitudes
* - 4 bytes | float : ssh scale_factor
* - 4 bytes | float : ssh add_offset
* Followed by loop by POIs
*   - 1 byte | char : flag of a next POI (=FLAG_NEWPOI)
*   - 2 bytes | short : ssh at time=0 (coseismic offset)
*   - 2 bytes | short : time step of the first arrival
* Optional:
*   - 1 byte | short : flag is next value (deviation from prediction) is coded with 2 bytes int (by default - with 1 byte) (=FLAG_2BYTE)
* If flag present:
*   - 2 bytes | short : ssh deviation as 2 byte int
* else
*   - 1 byte | char : ssh deviation as 1 byte int

*/


#include <stdio.h>
#include <string.h>
#include <netcdf.h>
#include <utilits.h>
#include "cNcHySea.h"

#define FLAG_NEWPOI -128
#define FLAG_2BYTE 127

#define ERRSYS(err) {Err.post("%s", nc_strerror(err)); return(-1);}
#define ERRCMPKT(err) {Err.post("Unexpected compact file format: %s", err); return(-1);}


//=========================================================================
// Compact time series with 'Tolkova algorithm'
// i.e., by storing deviation from extrapolated values as 1- or 2-bytes integers
// and write them together with additional data into the proprietary binary file
int cNcHySea::write_tolkova()
{
  typedef signed char schar;
  char buf[128];
  short sbuf;
  schar cbuf;
  const schar flagNewPoi(FLAG_NEWPOI);
  const schar flag2bytes(FLAG_2BYTE);

  if( version == VER381 ) return Err.post("Compressed output currently implemented for HySea ver 3.8.2 only");

  sprintf( buf, "%s.tkv", fname );
  FILE *fout = fopen( buf, "wb" );

  // 1: write file header. Version should correspond to that of HySea output -- 32 bytes
  memset( buf, 0, 64 );
  sprintf( buf, "TSCMPKT Ver:%s", VER_HYSEA_382 );
  fwrite( buf, sizeof(char), 32, fout );

  // 2: write scenario information -- 64 bytes
  memset( buf, 0, 64 );
  sprintf( buf, "64 bytes reserved for a scenario ID" );
  fwrite( buf, sizeof(char), 64, fout );

  // 3: write time step -- 2 bytes
  sbuf = 30;
  fwrite( &sbuf, sizeof(short), 1, fout );

  // 4: write number of time steps -- 2 bytes
  sbuf = ntime;
  fwrite( &sbuf, sizeof(short), 1, fout );

  // 5: write number of pois -- 2 bytes
  sbuf = npoi;
  fwrite( &sbuf, sizeof(short), 1, fout );

  // 6: write array of poi longitudes -- 8*npoi bytes
  fwrite( lon, sizeof(double), npoi, fout );
  // 7: write array of poi latitudes -- 8*npoi bytes
  fwrite( lat, sizeof(double), npoi, fout );

  // 8: write ssh scale_factor -- 4 bytes
  fwrite( &ssh_scale, sizeof(float), 1, fout );

  // 9: write ssh add_offset -- 4 bytes
  fwrite( &ssh_offset, sizeof(float), 1, fout );

  // 10: write flag for new poi -- 1 byte
  fwrite( &flagNewPoi, sizeof(schar), 1, fout );
  for( int n = 0; n<npoi; n++ ) {
    short u0,u1,u2,udev,k1st;

    u0 = sshs[npoi*0+n];
    // 11 x npoi: write first ssh -- 2 bytes
    fwrite( &u0, sizeof(short), 1, fout );

    k1st = 1;
    while( sshs[npoi*k1st+n] == u0 ) {
      k1st++;
      if( k1st == ntime-1 ) break;
    }
    // 12 x npoi: write time step of the first arrival -- 2 bytes
    fwrite( &k1st, sizeof(short), 1, fout );

    u2 = u1 = u0;
    for( int k=k1st; k<ntime; k++ ) {
      u0 = sshs[npoi*k+n];
      udev = u0 - (2*u1-u2);
      u2 = u1;
      u1 = u0;
      if( abs(udev) >= flag2bytes ) {
        // 13a x npoi: write flag for 2-bytes value which will follow -- 1 byte
        fwrite( &flag2bytes, sizeof(schar), 1, fout );
        // 13b x npoi: write deviation as -- 2 bytes
        fwrite( &udev, sizeof(short), 1, fout );
      }
      else {
        cbuf = (schar)udev;
        // 13 x npoi: write deviation as -- 1 byte
        fwrite( &cbuf, sizeof(schar), 1, fout );
      }
    }
  // 14 x npoi: write flag for new poi -- 1 byte
  fwrite( &flagNewPoi, sizeof(schar), 1, fout );
  }

  fclose(fout);

return 0;
}


//=========================================================================
// Read mareograms compacted with 'Tolkova algorithm',
// i.e., by storing deviation from extrapolated values as 1- or 2-bytes integers
// and fill the class members with this information
int cNcHySea::read_tolkova( const char *tkvfname )
{
  int ierr;
  typedef signed char schar;
  char buf[64], buf2[64];
  short sbuf;
  schar cbuf;
  const schar flagNewPoi(FLAG_NEWPOI);
  const schar flag2bytes(FLAG_2BYTE);

  FILE *fin = fopen( tkvfname, "rb" ); if( fin==NULL ) return Err.post(Err.msgOpenFile( tkvfname ));

  // 1: read file header. Version should correspond to that of HySea output -- 32 bytes
  memset( buf, 0, 64 ); memset( buf2, 0, 64 );
  sprintf( buf2, "TSCMPKT Ver:%s", VER_HYSEA_382 );
  fread( buf, sizeof(char), 32, fin );
  if( strcmp( buf, buf2 ) ) {Err.post("Valid for HySea version %s only!", VER_HYSEA_382); return(-1);}

  // 2: read scenario information -- 64 bytes
  memset( buf, 0, 64 );
  fread( buf, sizeof(char), 64, fin );

  // 3: read time step -- 2 bytes
  ierr = fread( &sbuf, sizeof(short), 1, fin ); if( ierr!=1 ) ERRCMPKT("time step");

  // 4: read number of time steps -- 2 bytes
  ierr = fread( &ntime, sizeof(short), 1, fin ); if( ierr!=1 ) ERRCMPKT("numer of time steps");

  // build tstamp vector using 3) and 4)
  tstamp = new float [ntime];
  for( int k=0; k<ntime; k++ )
    tstamp[k] = (float)(k*sbuf);

  // 5: read number of pois -- 2 bytes
  ierr = fread( &npoi, sizeof(short), 1, fin ); if( ierr!=1 ) ERRCMPKT("numer of pois");

  // 6: read array of poi longitudes -- 8*npoi bytes
  lon = new double [npoi];
  ierr = fread( lon, sizeof(double), npoi, fin ); if( ierr!=npoi ) ERRCMPKT("array of lon");
  // 7: read array of poi latitudes -- 8*npoi bytes
  lat = new double [npoi];
  ierr = fread( lat, sizeof(double), npoi, fin ); if( ierr!=npoi ) ERRCMPKT("array of lat");

  // 8: read ssh scale factor -- 4 bytes
  ierr = fread( &ssh_scale, sizeof(float), 1, fin ); if( ierr!=1 ) ERRCMPKT("scale_factor");

  // 9: read ssh add_offset -- 4 bytes
  ierr = fread( &ssh_offset, sizeof(float), 1, fin ); if( ierr!=1 ) ERRCMPKT("add_offset");

  sshs = new short [npoi*ntime];

  // 10: read flag for new poi -- 1 byte
  ierr = fread( &cbuf, sizeof(schar), 1, fin ); if( ierr!=1 || cbuf != flagNewPoi ) ERRCMPKT("start of poi data");

  for( int n = 0; n<npoi; n++ ) {
    short u0,u1,u2,udev,k1st;

    // 11 x npoi: read first ssh -- 2 bytes
    ierr = fread( &u0, sizeof(short), 1, fin ); if( ierr!=1 ) ERRCMPKT("u0");
    sshs[npoi*0+n] = u0;

    // 12 x npoi: read time step of the first arrival -- 2 bytes
    ierr = fread( &k1st, sizeof(short), 1, fin ); if( ierr!=1 ) ERRCMPKT("k1st");

    for( int k=1; k<k1st; k++ )
      sshs[npoi*k+n] = u0;

    u2 = u1 = u0;
    for( int k=k1st; k<ntime; k++ ) {

      // 12a x npoi: check if this is a flag for a 2-byte value to follow -- 1 byte
      ierr = fread( &cbuf, sizeof(schar), 1, fin ); if( ierr!=1 ) ERRCMPKT("try char");
      if( cbuf == flag2bytes ) {
        fseek( fin, SEEK_CUR, -1 );
        // 12b x npoi: read deviation as -- 2 bytes
        ierr = fread( &udev, sizeof(short), 1, fin ); if( ierr!=1 ) ERRCMPKT("2-bytes dev");
      }
      else {
        // 12b x npoi: deviation was already read as 1-byte
        udev = (short)cbuf;
      }

      sshs[npoi*k+n] = udev + (2*u1-u2);
      u2 = u1;
      u1 = sshs[npoi*k+n];
    }

    // 13 x npoi: read flag for new poi -- 1 byte
    ierr = fread( &cbuf, sizeof(schar), 1, fin ); if( ierr!=1 || cbuf != flagNewPoi ) ERRCMPKT("end of poi data");
  }

  fclose(fin);

return 0;
}
