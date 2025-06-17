#!/usr/bin/python3
import logging
import subprocess
import time

import numpy as np
# from rabbitmq_demo import report, exposure, config, waveform
# from rabbitmq_demo import az_blob_storage
# import netCDF4, shutil, glob
# import time as t
# import os
# import geopandas as gpd
# import pandas as pd
from scipy import interpolate
# import json

# from rabbitmq_demo.exposure import generate_shapefile_names
# from rabbitmq_demo.hysea_repository import HySEARepository, HttpResourceNotFoundError


# def tms_py(ncfile,af,ctr,HySEAid):
#     ampfact={}
#     ampfact["nonzero"]=0

#     file  = netCDF4.Dataset(ncfile, 'r')
#     latv  = file.variables["latitude"][:]
#     lonv  = file.variables["longitude"][:]
#     time  = file.variables["time"][:]
#     depthv= file.variables["deformed_bathy"][:]
#     eta   = file.variables["eta"][:]

#     tmscount=0
#     nonzero=0
    
#     for i in range(len(eta[0,:])):
#         e=eta[:,i]
#         lon=lonv[i]
#         lat=latv[i]
#         if lon>180.0:
#             lon-=360.0
#         keyll="%(lon).04f/%(lat).04f" %vars()   
#         try:
#             id=ctr[keyll]
#         except:
#             print("problem: keyll",keyll)
#             continue
#         depth=depthv[i]
#         period,max,polarity=waveform.waveform(time,e,i,id,HySEAid,treshold=0.05,plot=False,quiet=True) 



#         nonzero=0
        
#         if max > 0:
#             #get AF-factor
#             ampf=interp_ampfact(period,af["periods"],af[id][polarity])
#             #correct with Greens law
#             max *= pow(depth/50, 0.25)
#             #compute MIH
#             MIH = ampf*max
#             nonzero=1
#         elif abs(max)<0.05 and depth>0:
#             #this is points in sea with max below treshold (0.05m)
#             max=MIH=0.0
#             nonzero=1
#         else:
#             #rest (either points on land or strange wave signal)
#             MIH = -1
#             max = -1
#             period = -1
#             nonzero=0
#         ampfact[id]={}
#         ampfact[id]["count"]=tmscount
#         ampfact[id]["lon"]=lon
#         ampfact[id]["lat"]=lat
#         ampfact[id]["MIH"]=MIH
#         ampfact[id]["key"]=keyll
#         ampfact[id]["period"]=period
#         ampfact[id]["hmax"]=max
#         ampfact[id]["area"]=af[id]["area"]
#         ampfact["nonzero"]+=nonzero
#         tmscount+=1

#     return ampfact,tmscount

def tms_cpp(ncfile,ctrfile,af,ctr):

    ampfact={}
    ampfact["nonzero"]=0
    #output from cpp-program is a txt file:
    wavef=f"{ncfile}.offshore.txt"
    
    # read the ampfact into dictionary both the combination of output from cpp (waveforms) and controlfile (with correct id)
    # 1. read in control_dict: ctr_dict with key as the counting number of lines

    if ctrfile:
        count=0
        ctr_dict={}
        with open(ctrfile,'r') as f:
            for line in f:
                if count>0:
                    l=line.strip().split()
                    id,lon,lat=l[0],l[1],l[2]
                    ctr_dict[count]={}
                    ctr_dict[count]["id"]=id
                    ctr_dict[count]["lon"]=lon
                    ctr_dict[count]["lat"]=lat
                count+=1


    # 2. read in waveform dict, fill dict ampfact with calculated MIH from AF (key is global id, idXXXXX)
    tmscount=0
    nonzero=0
    with open(wavef,'r') as f:
        for line in f:
            l=line.strip().split()
            if tmscount>0:
                lon,lat,depth,max,min,period,polarity,return_code=float(l[1]),float(l[2]),float(l[3]),float(l[4]),float(l[5]),float(l[6]),int(l[7]),float(l[8])
                keyll="%(lon).04f/%(lat).04f" %vars()   
                if polarity==-1:
                    polarity="neg"
                elif polarity==1:
                    polarity="pos"
                if ctrfile:
                    id=ctr_dict[tmscount]["id"]
                else:
                    try:
                        id=ctr[keyll]
                    except:
                        print("problem: keyll",keyll)
                        continue
                nonzero=0
                if return_code > 1:
                    #get AF-factor
                    ampf=interp_ampfact(period,af["periods"],af[id][polarity])
                    #correct with Greens law
                    max *= pow(depth/50, 0.25)
                    #compute MIH
                    MIH = ampf*max
                    nonzero=1
                elif return_code == 1 and abs(max)<0.05 and depth>0:
                    #this is points in sea with max below treshold (0.05m)
                    max=MIH=0.0
                    nonzero=1
                else:
                    #rest (either points on land or strange wave signal)
                    MIH = -1
                    max = -1
                    period = -1
                    nonzero=0
                
                if isinstance(max,complex):
                    #bug in cpp, complex number due to negative depth, skip POI
                    #will be replaced with a neighboring POI in next step
                    MIH = -1
                    max = -1
                    period = -1
                    nonzero=0

                   
                ampfact[id]={}
                ampfact[id]["count"]=tmscount
                ampfact[id]["lon"]=lon
                ampfact[id]["lat"]=lat
                ampfact[id]["MIH"]=MIH
                ampfact[id]["key"]=keyll
                ampfact[id]["period"]=period
                ampfact[id]["hmax"]=max
                ampfact[id]["area"]=af[id]["area"]
                ampfact["nonzero"]+=nonzero
            tmscount+=1
    #fill points with MIH==-1 with values from closest point (new iteration over all items in results-dictionary)

    return ampfact,tmscount


def read_ampfacts(file):
    #fill a dictionary with  "lon/lat" as key - can substitute with idXXXXX later
    af={}
    ctr={}   #to be used for mapping id and lon/lat if file from HySEA mapping these are missing
    periods=np.array([120,200,300,600,1000,1800,3600])
    count=0
    af["periods"]=periods
    for l in open(file,'r').readlines()[1:]:

        l=l.strip().split()
        id,lon,lat=l[0],float(l[1]),float(l[2])
        neg,pos=np.array([float(i) for i in l[4:11]]),np.array([float(i) for i in l[12:19]])
        idll="%(lon).04f/%(lat).04f" %vars()

        #fill in into dictionary
        af[id]={}
        af[id]["id"]=id
        af[id]["lon"]=lon
        af[id]["lat"]=lat
        af[id]["neg"]=neg
        af[id]["pos"]=pos
        af[id]["area"]=str(l[20])

        ctr[idll]=id
        count+=1

    return af,ctr


def interp_ampfact(inper,per,af):
    #if period is shorter use T=120
    if inper<per[0]:
        AF=af[0]
    #if period is longer use T=3600
    elif inper>per[-1]:
        AF=af[-1]
    else:
        f = interpolate.interp1d(per,af)
        AF=f([inper])
        AF=AF[0]
    return AF


# def fill_empty_segments(ampfact):
#     """Segments without data (MIH=-1) will get copy of the values from the closest nonzero segment"""
#     from geopy.distance import geodesic as GD

#     def distance(ampfact,id1,id2):
#         """Compute great circle distance between two points"""
#         #check if segment id2 have MIH > -1
#         dist=False
#         if id2 in ampfact.keys():
#             if "MIH" in ampfact[id2].keys():
#                 if ampfact[id2]["MIH"]>0:
#                     lon1,lat1=ampfact[id1]["lon"],ampfact[id1]["lat"]
#                     lon2,lat2=ampfact[id2]["lon"],ampfact[id2]["lat"]
#                     A,B=(lat1,lon1),(lat2,lon2)
#                     dist=GD(A,B) #great circle distance
#         return dist


#     for k in ampfact.keys():
#         if  k != "nonzero" and ampfact[k]["MIH"]<0:
#             mindist=999999999
#             mindistid=""
#             #find closest segment with MIH>0, iterate over closest segments
#             idint=int(k[2:])
#             #iterate over closest segments (a maximum number) give a maximum allowed distance between segments
#             #running over ids +1 -1 +2 -2 +3 -3 .... (totally):
#             numbers = []
#             for i in range(1,10):
#                 numbers.extend([i,-i])
#             for i in numbers:
#                 id2="id{:05n}".format(idint+i)
#                 dist=distance(ampfact,k,id2)

#                 if dist and 0<dist<mindist:
#                     mindist=dist
#                     mindistid=id2
#                 if mindist< 20:
#                     #satisfied if the distance is less than 20 km
#                     break

#             if mindist<50:
#                 for ak in ampfact[k].keys():
#                     if ak in ["lon","lat"]:
#                         continue
#                     else:
#                         ampfact[k][ak]=ampfact[mindistid][ak]
                
                    
#                 #ampfact[k]["MIH"]       = ampfact[mindistid]["MIH"]
#                 #ampfact[k]["hmax"]      = ampfact[mindistid]["hmax"]
#                 #ampfact[k]["period"]    = ampfact[mindistid]["period"]
#                 #for documentation and tracking:
#                 ampfact[k]["mindistid"] = mindistid
#                 ampfact[k]["mindist"]   = mindist


           
#     return ampfact

# def NEAM_area(lon,lat):
#         #include Black Sea etc. located in Asia to NEAM region
#         if 23<lon<45 and 28<lat<50:
#             return True
#         else:
#             return False
        
# def add_MIH_percentiles(af, ampfact):
#     for k in ampfact.keys():
#         if k!="nonzero":
#             #print("cont",ampfact[k].keys())
#             #continent=ampfact[k]["continent"]
#             kll=ampfact[k]["key"]
#             lonp=af[k]["lon"]
#             latp=af[k]["lat"]
#             MIH=ampfact[k]["MIH"]
#             hmax=ampfact[k]["hmax"]
#             per=ampfact[k]["period"]
#             area=ampfact[k]["area"]
            

#             #correction factors for percenetiles and mean values (global)
#             mf5,mf16,mfm,mf84,mf95=0.1848,0.338,1.30,2.13,3.9

#             #NEAM correction factors
#             if area in ["NEA","MED","BLK"] or NEAM_area(lonp,latp):
#                 mf5,mf16,mfm,mf84,mf95=0.322,0.514,1.36,2.15,3.44

#             MIH5,MIH16,MIHm,MIH84,MIH95=mf5*MIH,mf16*MIH,mfm*MIH,mf84*MIH,mf95*MIH
#                         #update dictionary ampfact with percentiles and mean value
#             ampfact[k]["MIH5"]  = MIH5
#             ampfact[k]["MIH16"] = MIH16
#             ampfact[k]["MIHm"]  = MIHm
#             ampfact[k]["MIH84"] = MIH84
#             ampfact[k]["MIH95"] = MIH95 
#     return ampfact

def save_ampfacts_to_file(ampfact,outfile):

    #outfile=os.path.join(results_dir, f"results_{HySEAid}_MIH.txt")
    fil=open(outfile,'w')

    fil.write("#id lon lat MIH hmax period\n")
    for k in ampfact.keys():
        if k and k!="nonzero":        
            print("k",k)       
            #kll   = ampfact[k]["key"]
            lonp  = ampfact[k]["lon"]
            latp  = ampfact[k]["lat"]
            hmax  = ampfact[k]["hmax"]
            per   = ampfact[k]["period"]
            MIH   = ampfact[k]["MIH"]
            # MIH5  = ampfact[k]["MIH5"]  
            # MIH16 = ampfact[k]["MIH16"] 
            # MIHm  = ampfact[k]["MIHm"]   
            # MIH84 = ampfact[k]["MIH84"] 
            # MIH95 = ampfact[k]["MIH95"]
            # #print("id occup",k)





            fil.write("%(k)s %(lonp)f %(latp)f %(MIH).2f %(hmax).2f %(per).0f\n" %vars())

            #update dictionary ampfact with percentiles and mean value

    fil.close()
    return outfile



# def read_eqinfo(eqfile: str|bool) -> dict[str,list[str,str]]:
#     eqinfo={}
#     degr="\u00B0"

#     #fill default values
#     txt="<not available>"
#     eqinfo["lon"]       = txt
#     eqinfo["lat"]       = txt
#     eqinfo["depth"]     = txt
#     eqinfo["length"]    = txt
#     eqinfo["width"]     = txt
#     eqinfo["strike"]    = txt
#     eqinfo["dip"]       = txt
#     eqinfo["rake"]      = txt
#     eqinfo["slip"]      = txt
#     eqinfo["magnitude"] = txt
#     try:    
#         #read out data from file:
#         with open(eqfile, 'r') as file:
#             data = [float(x) for x in file.readlines()[1].strip().split()]
#         print("data",data)
#         eqinfo["lon"]       = f'{data[0]:.2f}{degr}'
#         eqinfo["lat"]       = f'{data[1]:.2f}{degr}' 
#         eqinfo["depth"]     = f'{data[2]:.0f} km'
#         eqinfo["length"]    = f'{data[3]:.0f} km'
#         eqinfo["width"]     = f'{data[4]:.0f} km'
#         eqinfo["strike"]    = f'{data[5]:.0f}{degr}' 
#         eqinfo["dip"]       = f'{data[6]:.0f}{degr}' 
#         eqinfo["rake"]      = f'{data[7]:.0f}{degr}' 
#         eqinfo["slip"]      = f'{data[8]:.1f} m'
#         eqinfo["magnitude"] = f'{data[9]:.1f} Mw'
#     except FileNotFoundError:
#         print(f"The file was not found {eqfile}.")
#     except PermissionError:
#         print(f"You don't have permission to read the file {eqfile}.")
#     except Exception as e:
#         print(f"An unexpected error occurred: {e}")
#     print("eqinfo",eqinfo)
#     return eqinfo

# def find_unique_values(obj, target_key, result_set=None):
#     if result_set is None:
#         result_set = set()

#     if isinstance(obj, dict):
#         for key, value in obj.items():
#             if key == target_key:
#                 result_set.add(value)
#             find_unique_values(value, target_key, result_set)

#     elif isinstance(obj, list):
#         for item in obj:
#             find_unique_values(item, target_key, result_set)

#     return sorted(result_set)

# def sum_values_for_key(obj, target_key):
#     total = 0

#     if isinstance(obj, dict):
#         for key, value in obj.items():
#             if key == target_key and isinstance(value, (int, float)):
#                 total += value
#             total += sum_values_for_key(value, target_key)

#     elif isinstance(obj, list):
#         for item in obj:
#             total += sum_values_for_key(item, target_key)

#     return total


# class InvalidHySEAJobError(Exception):
#     pass


# async def run(HySEAid: int, settings: config.Settings) -> list[str]:
#     logging.info(f"Connecting to storage account: {settings.account_url}, HySEA id {HySEAid}")
#     if settings.account_url:
#         storage_client = az_blob_storage.BlobStorageClient(
#             settings.account_url,
#             input_container=settings.input_blob_container,
#             output_container=settings.output_blob_container
#         )
#     else:
#         logging.warning("Using local storage only. This should only be used for local development.")
#         storage_client = None

#     if storage_client:
#         files_uploaded = []

#     try:
#         hysea_repository = HySEARepository(f"{settings.root_url}/{HySEAid}")
#     except HttpResourceNotFoundError as e:
#         raise InvalidHySEAJobError(str(e))

#     t0  = t.perf_counter()
#     ct0 = t.process_time()

#     ###################################################################33
#     #defaults and initiate
#     error=[]
#     results_dir=os.path.join(settings.results_dir,str(HySEAid))
    
#     if not os.path.exists(results_dir):
#         os.makedirs(results_dir)


#     ###################################################################33
#     #copy control-file for id and pos of timeseries
#     logging.info("Copy additional data from HySEA")
#     localctr = os.path.join(results_dir, f"ts_HySEA_ids_id{HySEAid}.txt")
#     extctr = f"ts_HySEA_ids_id{HySEAid}.txt"
#     try:
#         hysea_repository.extract_to_local_file(extctr, localctr)
#     except HttpResourceNotFoundError as e:
#         logging.info(f"Could not locate: {extctr}")  # TODO: what is done instead in this case?
#         # if controlfile (key between ids and timeseries from HySEA) is not present
#         localctr = False

#     #copy file with timeseries
#     tmslocalnc=os.path.join(results_dir, f"result_id{HySEAid}_ts.nc")
#     extnc=f"result_id{HySEAid}_ts.nc"
#     try:
#         hysea_repository.extract_to_local_file(extnc, tmslocalnc)
#     except HttpResourceNotFoundError as e:
#         raise InvalidHySEAJobError(str(e))


#     #file with exposure data for all segments
#     exposure_file=os.path.join(settings.data_dir,"World.txt")

#     ###################################################################33
#     #read in ampfactors
#     logging.info("Read ampfactors")
#     file=os.path.join(settings.data_dir, "global_ampf_v04.txt")

#     #ctr is a dictionary to be used if localctr is missing an id for a gauge
#     #is to be found by using lon/lat not id:
#     af,ctr=read_ampfacts(file)
    
#     ###################################################################
#     #
#     #   MAIN LOOP OVER ALL TIMESERIES, CALCULATE MIH FROM AF
#     #   Results are put into the results dictionary
#     #
#     ###################################################################
#     #extract waveforms from HySEA output (ncfile) - saved in file {ncfile}.offshore.txt
#     try:
#         logging.info("Run cpp for timeseries")
#         res = subprocess.run([settings.postprocess_script_path, tmslocalnc], check=True, capture_output=True)

#         logging.debug(f"res = {res}")
#         logging.info(f"stdout: {res.stdout!r}")
#         if res.stderr:
#             logging.error(f"stderr: {res.stderr!r}")
#         ampfact,tmscount=tms_cpp(tmslocalnc,localctr,af,ctr)  
#     except:
#         logging.info("Cpp failed, run python for timeseries") 
#         ampfact,tmscount=tms_py(tmslocalnc,af,ctr,HySEAid)   
            

#     #####################################################################

    
#     nonzero=ampfact["nonzero"]  
#     logging.info(f"Number timeseries processed: {tmscount}")
#     logging.info(f"Number of nonzero timeseries: {nonzero}")
#     logging.info(f"Netcdf: {extnc}")

#     #####################################################################
#     #
#     # Shapefile for coastal segments 
#     #
#     #####################################################################
    
#     # use file for segment line directly from data_dir (no modification needed) 
#     shape_seg=os.path.join(settings.data_dir,settings.shape_seg)

#     #####################################################################
#     #
#     # extract coordinates for domain and earthquake information
#     #
#     #####################################################################
#     # download file with corners and centrpoint of EQ
#     # extract corners of computational domain
#     logging.info("Download data from HySEA-simulations")
#     local = os.path.join(results_dir, f"corners_id{HySEAid}.txt")
#     ext = f"corners_id{HySEAid}.txt"
#     try:
#         hysea_repository.extract_to_local_file(ext, local)
#         logging.info("Downloaded coordinates for the domain for HySEA-simulations")
#         domain = [float(i) for i in open(local).readlines()[0].split()]
#     except HttpResourceNotFoundError as e:
#         #if corners-file not present:
#         #extract domain from timeseries file
#         logging.info(f"Create domain from timeseries location (file with corners is not present: {ext})")
#         import xarray as xr
#         ds = xr.open_dataset(tmslocalnc)

#         data=ds["longitude"]
#         df=data.to_dataframe().reset_index()
#         lon=df["longitude"].values
#         lon[lon>180]-=360
#         lonmin=lon.min()
#         lonmax=lon.max()

#         data=ds["latitude"]
#         df=data.to_dataframe().reset_index()
#         latmin=df["latitude"].min()
#         latmax=df["latitude"].max()
#         domain=[lonmin, lonmax, latmin, latmax]


#     #extract position of centerpoint of eq
#     local = os.path.join(results_dir, f"eq_position_id{HySEAid}.txt")
#     ext = f"eq_position_id{HySEAid}.txt"
#     try:
#         hysea_repository.extract_to_local_file(ext, local)
#         eqpos = open(local).readlines()[0].split()
#     except HttpResourceNotFoundError as e:
#         # if eq_position file not present, put eq to (0,0):
#         eqpos = [0, 0]
#         logging.info(f"Could not find eq position file: {ext}. Defaulting to {eqpos}")

#     #file with earthquake information
#     localeqinfo = os.path.join(results_dir, f"parameters_id{HySEAid}.txt")
#     exteqinfo = f"parameters_id{HySEAid}.txt"
#     try:
#         hysea_repository.extract_to_local_file(exteqinfo, localeqinfo)
#     except HttpResourceNotFoundError as e:
#         logging.info(f"File with parameters for the earthquake {exteqinfo} is not present")
#         localeqinfo = None

#     logging.info("Read earthquake information")
#     eqinfo=read_eqinfo(localeqinfo)

#     #####################################################################
#     #
#     #  add percentiles for MIH to ampfact dictonary
#     #
#     #####################################################################
#     logging.info("Adding percentiles")
#     ampfact=add_MIH_percentiles(af, ampfact)
    
#     #####################################################################
#     #
#     #  Replace segments without MIH (MIH=-1) with values from closest 
#     #  segment with a value
#     #
#     #####################################################################
#     logging.info("Fill empty segments")
#     ampfact = fill_empty_segments(ampfact)

#     #####################################################################
#     #
#     #  Add exposure data based on calculated MIH ()
#     #
#     #####################################################################
#     logging.info("Add exposure data into dictionary ampfact")
#     ampfact = exposure.add_exposure(ampfact,exposure_file,settings.mih_type)

#     #####################################################################
#     #
#     #  Add segment lengths into ampfact dictionary from shapefile
#     #     
#     #####################################################################
#     logging.info("Add segment length into dictionary ampfact")
#     ampfact = exposure.add_segment_lengths(shape_seg, ampfact)  

#     #####################################################################
#     #sum up exposed persons for entire domain and country for country
#     #####################################################################
    
#     total_occupants=0
#     exp_country={}

#     for id in ampfact.keys():
#         if "id" not in id:
#             continue
#         countr=ampfact[id]["country"]
#         occ=ampfact[id]["occupants"]
#         mih=ampfact[id][settings.mih_type]
#         if countr not in exp_country.keys():
#             exp_country[countr]={}
#             exp_country[countr]["exposed"]=0
#             exp_country[countr]["low_mih"]=999999
#             exp_country[countr]["high_mih"]=-999999
#         exp_country[countr]["exposed"]+=occ
#         low=exp_country[countr]["low_mih"]
#         if mih<low:
#             exp_country[countr]["low_mih"]=mih
#         high=exp_country[countr]["high_mih"]
#         if mih>high:
#             exp_country[countr]["high_mih"]=mih
#         total_occupants+=occ
#     logging.info(f"Exposed by country {exp_country}")

#     #####################################################################
#     #
#     #  extract inundation polygons and create kmz file
#     #  update exposure, countries etc to ampfact dict
#     #
#     #####################################################################
#     logging.info("Make kmz and calculate exposure")

    
#     if storage_client and settings.create_kmz:
#         polygon_treshold=settings.polygon_treshold
#         #shapefile_ids = [k for k in ampfact.keys() if k.startswith("id") and float(ampfact[k][settings.mih_type])>polygon_treshold and ampfact[k]["country"]=="Japan"]
#         shapefile_ids = [k for k in ampfact.keys() if k.startswith("id") and float(ampfact[k][settings.mih_type])>polygon_treshold]
#         shapefile_ext = [".sbn", ".sbx", ".cpg", ".dbf", ".prj", ".shp", ".shx"]
#         #shapefile_paths = ["shape/"+name for name in generate_shapefile_names(shapefile_ids, shapefile_ext)]
#         shapefile_names = generate_shapefile_names(shapefile_ids, shapefile_ext)


#         if not os.path.exists(settings.shapefile_dir_id):
#             os.makedirs(settings.shapefile_dir_id)

#         logging.info(f"Downloading {len(shapefile_names)} files...")
#         t0_download = time.perf_counter()

#         download_count=0
#         if len(shapefile_names)>0:
#             download_count = await storage_client.concurrent_download_multiple("shape", settings.shapefile_dir_id, shapefile_names)
#         # for p in shapefile_paths:
#         #     storage_client.download("shape/"+p, os.path.join(settings.shapefile_dir_id, p))
#         t1_download = time.perf_counter()
#         logging.info(f"Downloaded {download_count} files in {t1_download - t0_download} seconds")

#     if settings.create_kmz:
#         kmzfile_name_MIH = f"res_MIH_{HySEAid}.kmz"
#         kmzfile_local_path_MIH = os.path.join(results_dir, kmzfile_name_MIH)
#         exposure.run_kmz(results_dir, HySEAid, eqpos, domain, settings.shapefile_dir_id, shape_seg, ampfact, kmzfile_local_path_MIH, total_occupants, polygon_treshold=settings.polygon_treshold, colorthres=1,value_column=settings.mih_type, value_column_MIH=settings.mih_type)
#         kmzfile_name_exp = f"res_exp_{HySEAid}.kmz"
#         kmzfile_local_path_exp = os.path.join(results_dir, kmzfile_name_exp)
#         exposure.run_kmz(results_dir, HySEAid, eqpos, domain, settings.shapefile_dir_id, shape_seg, ampfact, kmzfile_local_path_exp, total_occupants, polygon_treshold=settings.polygon_treshold, colorthres=1,value_column="occupants", value_column_MIH=settings.mih_type)
#         if storage_client:
#             files_uploaded.append(storage_client.upload(f"job_{HySEAid}/{kmzfile_name_MIH}", kmzfile_local_path_MIH))
#             files_uploaded.append(storage_client.upload(f"job_{HySEAid}/{kmzfile_name_exp}", kmzfile_local_path_exp))

#     #####################################################################
#     #
#     #  Save ampfact to a text file (columnwise)
#     #
#     #####################################################################
#     logging.info(f"Save ampfact to file")
#     ampf_file=save_ampfacts_to_file(ampfact,results_dir,HySEAid,error)


#     # calculated total exposed:
#     exposed=0
#     exposed_ids=0
#     for k in ampfact.keys():
#         if "id" in k and "occupants" in ampfact[k].keys():
#             exposed+=ampfact[k]["occupants"]
#             exposed_ids+=1
    
#     logging.info(f"Totally {exposed} persons exposed by the tsunami, in totally {exposed_ids} segments.")



#     #####################################################################
#     #
#     #  Save exposed persons versus MIH levels and 
#     #
#     #####################################################################
#     logging.info(f"Save exposed persons versus MIH levels")   
#     # levels given as text string in config.py/.env file:
#     levels=np.array(settings.mih_levels.split(","))
#     # list must be numpy-array
#     levels=np.array(sorted([float(numb) for numb in levels]))

#     list_countries= list(find_unique_values(ampfact, "country"))
#     list_countries.append("all")

#     #initialize a dictionary
#     explev={}
#     for c in list_countries:
#         explev[c]={}
#         for i in levels:
#             explev[c][i]={}
#             explev[c][i]["km"]=0
#             explev[c][i]["segments"]=0
#             explev[c][i]["exposed"]=0
#             explev[c][i]["exposed16"]=0
#             explev[c][i]["exposed84"]=0
        
  
#     for k in ampfact.keys():
#         if "id" in k:
#             mih=float(ampfact[k][settings.mih_type])
#             c=ampfact[k]["country"]
#             #find correct mih interval to put results for this segment
#             filt_ind = np.where(levels < mih)[0]
#             if filt_ind.size > 0:
#                 closest_ind=filt_ind[np.argmax(levels[filt_ind])]
#                 closest_val=levels[closest_ind]
#                 explev[c][closest_val]["km"]            += ampfact[k]["seg_length"]
#                 explev[c][closest_val]["segments"]      += 1
#                 explev[c][closest_val]["exposed"]       += ampfact[k]["occupants"]

#                 explev["all"][closest_val]["segments"]  += 1
#                 explev["all"][closest_val]["km"]        += ampfact[k]["seg_length"]
#                 explev["all"][closest_val]["exposed"]   += ampfact[k]["occupants"]

            
#             # 16 percentile
#             mih=float(ampfact[k]["MIH16"])
#             c=ampfact[k]["country"]
#             #find correct mih interval to put results for this segment
#             filt_ind = np.where(levels < mih)[0]
#             if filt_ind.size > 0:
#                 closest_ind=filt_ind[np.argmax(levels[filt_ind])]
#                 closest_val=levels[closest_ind]
#                 explev[c][closest_val]["exposed16"]     += ampfact[k]["p16_occupants"]
#                 explev["all"][closest_val]["exposed16"] += ampfact[k]["p16_occupants"]

#             # 84 percentile
#             mih=float(ampfact[k]["MIH84"])
#             c=ampfact[k]["country"]
#             #find correct mih interval to put results for this segment
#             filt_ind = np.where(levels < mih)[0]
#             if filt_ind.size > 0:
#                 closest_ind=filt_ind[np.argmax(levels[filt_ind])]
#                 closest_val=levels[closest_ind]
#                 explev[c][closest_val]["exposed84"]     += ampfact[k]["p84_occupants"]
#                 explev["all"][closest_val]["exposed84"] += ampfact[k]["p84_occupants"]

 
 
#     #save sum of exposed persons by country to file (json)
#     outfile=os.path.join(results_dir, f"results_{HySEAid}_expcountries.json")
#     with open(outfile, "w", encoding="utf-8") as f:
#         json.dump(exp_country, f, indent=2, ensure_ascii=False)
#     #save explevels to json format
#     outfile=os.path.join(results_dir, f"results_{HySEAid}_explevels.json")
#     with open(outfile, "w", encoding="utf-8") as f:
#         json.dump(explev, f, indent=2, ensure_ascii=False)

#     #####################################################################
#     #
#     #  make plots
#     #
#     #####################################################################
#     logging.info("Make plots")
#     dx=5
#     xll,xur,yll,yur=float(eqpos[0])-dx,float(eqpos[0])+dx,float(eqpos[1])-dx,float(eqpos[1])+dx
#     zoom_domain=[xll,xur,yll,yur]
#     # make zoomed domain plot - output zoom.png
#     report.plot("zoom",results_dir,shape_seg,eqpos,zoom_domain,settings.boundary_file,ampfact,settings.mih_levels,value_column=settings.mih_type)
#     # make large domain plot - output large.png
#     report.plot("large",results_dir,shape_seg,eqpos,domain,settings.boundary_file,ampfact,settings.mih_levels,value_column=settings.mih_type,zoom=zoom_domain)

#     #potentially exposed:
#     # make zoomed domain plot - output zoom.png
#     report.plot("zoom",results_dir,shape_seg,eqpos,zoom_domain,settings.boundary_file,ampfact,False,value_column="occupants")
#     # make large domain plot - output large.png
#     report.plot("large",results_dir,shape_seg,eqpos,domain,settings.boundary_file,ampfact,False,value_column="occupants",zoom=zoom_domain)

#     #####################################################################
#     #
#     #  make report 
#     #
#     #####################################################################
#     # table with largest MIH etc.
#     rep=os.path.join(results_dir,f"report_{HySEAid}.pdf")
#     logging.info(f"Making report {rep}")
#     report.report(f"report_{HySEAid}",results_dir,HySEAid,ampfact,exp_country,eqinfo,explev, settings.mih_type)
    
#     ct1 = t.process_time()
#     t1  = t.perf_counter()
#     logging.info(f"Elapsed time: {t1 - t0} sec")
#     logging.info(f"Process time: {ct1 - ct0} sec")  

#     if storage_client:
#         blob_path = f"job_{HySEAid}/report_{HySEAid}.pdf"
#         files_uploaded.append(storage_client.upload(blob_path, rep))
#         return files_uploaded
#     else:
#         return [results_dir]


    
if __name__ == "__main__":
    # read in ampfactors
    af,ctrl=read_ampfacts("/home/vmg/NGI/P/2022/02/20220296/Calculations/testing-Messina-POIs-AmpFactor/landslide_tsunami_HySEA/factors/ampf_messina_new.txt")
    #af,ctrl=read_ampfacts("../factors/global_v04_messina.txt")
    
    # read in waveparam
    ampfact,tmscount=tms_cpp("/home/vmg/NGI/P/2022/02/20220296/Calculations/Messina_simulations_June2025/8391804/hysea_out/filterkajiura_res100_ts.nc",None,af,ctrl)
     
    # save ampfact-dict (with calculated MIH) to file                  
    save_ampfacts_to_file(ampfact,"/home/vmg/NGI/P/2022/02/20220296/Calculations/Messina_simulations_June2025/8391804/hysea_out/MIH.txt")
