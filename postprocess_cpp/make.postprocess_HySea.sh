OTH=./share
NETCDF_INCLUDE=/usr/include
NETCDF_LIB=/usr/lib

g++-11 -I$OTH -I$NETCDF_INCLUDE -O2 -o postprocess_HySea postprocess_HySea.cpp cNcHySea_rw.cpp cNcHySea_tlkv.cpp cNcHySea_offshore.cpp $OTH/utilits.cpp $OTH/cXY.cpp $OTH/cMrg.cpp -L$NETCDF_LIB -lnetcdf

