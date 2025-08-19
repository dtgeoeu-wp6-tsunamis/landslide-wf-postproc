#include <stdio.h>
#include <utilits.h>
#include "cNcHySea.h"


int commandLineHelp(){
  printf("\n\nUsage: postprocess_HySea filename [-o fmtstring] [-poi poiID]\n");
  printf("\twhere 'filename' stands for HySea netCDF output file\n");
  printf("\t'fmtstring' can contain 'netcdf' or/and 'txt' to indicate output format, or both\n");
  printf("\t'poiID is POI number. If selected, analysis will be restricted to that POI only\n\n");
  return -1;
}


int main( int argc, char **argv )
{
  int ierr, argn, poi(-1);
  char fmtStr[32];

  if( argc == 1 ) return commandLineHelp();

  if( (argn = utlCheckCommandLineOption( argc, argv, "o", 1)) )
    sprintf( fmtStr, "%s", argv[argn+1] );
  else
    sprintf( fmtStr, "netcdf+txt" );

  if( (argn = utlCheckCommandLineOption( argc, argv, "poi", 3)) ) {
    if( sscanf( argv[argn+1], "%d", &poi ) != 1 ) return commandLineHelp();
    poi--;
  }

  cNcHySea hysea;

  ierr = hysea.read_hysea( argv[1] ); if(ierr) return ierr;

  if( poi > -1 )
    hysea.analyse_offshore( poi );
  else
    hysea.write_offshore( fmtStr );

  return 0;

}
