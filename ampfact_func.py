#!/usr/bin/python3
import numpy as np
from scipy import interpolate

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

def save_ampfacts_to_file(ampfact,outfile):

    fil=open(outfile,'w')

    fil.write("#id lon lat MIH hmax period\n")
    for k in ampfact.keys():
        if k and k!="nonzero":        
            lonp  = ampfact[k]["lon"]
            latp  = ampfact[k]["lat"]
            hmax  = ampfact[k]["hmax"]
            per   = ampfact[k]["period"]
            MIH   = ampfact[k]["MIH"]

            fil.write("%(k)s %(lonp)f %(latp)f %(MIH).2f %(hmax).2f %(per).0f\n" %vars())

    fil.close()
    return outfile

    
def ampfact_func(ampfact_file, hyseaout_file, outfile_name):
    """
    Function to compute maximum inundation height (MIH) using amplification factors.
     - ampfact_file: name of the file with amplification factors for each POI 
     - hyseaout_file: name of the T-HySEA time series output file (without the .offshore.txt suffix)
     - outfile_name: name of the output file where the results will be saved
    """
    # read in ampfactors
    af,ctrl=read_ampfacts(ampfact_file)
    
    # read in waveparam
    ampfact,tmscount=tms_cpp(hyseaout_file,None,af,ctrl)
     
    # save ampfact-dict (with calculated MIH) to file                  
    save_ampfacts_to_file(ampfact,outfile_name)
