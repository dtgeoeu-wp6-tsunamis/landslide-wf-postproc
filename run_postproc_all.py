"""
Do postprocessing of all HySEA simulations:
- postprocess_HySea is C++ code that extract info from time series
- ampfact_func.py is a script that uses the amplification factors to compute the inundation
"""

import os
import sys
import pandas as pd
from ampfact_func import ampfact_func

parent_dir = "/home/vmg/NGI/P/2022/02/20220296/Calculations/Messina_simulations_June2025//"
volumes_dir = "/home/vmg/NGI/P/2022/02/20220296/Calculations/release-volume-sampler_20250602_kl1242/volumes/"
filename = os.path.join(volumes_dir,"Clusters_sorted_by_volume.csv")

# File with amplification factors
ampfact_file="/home/vmg/NGI/P/2022/02/20220296/Calculations/testing-Messina-POIs-AmpFactor/landslide_tsunami_HySEA/factors/ampf_messina_new.txt"

## Read the clusters CSV file to get the scenario IDs
## Use the ID to find the HySEA time series file and run the amp factor postprocessing
df = pd.read_csv(filename, usecols=["id"])
for i in range(0,2): #len(df)):
    scenario = str(df["id"][i])
    hysea_ts = os.path.join(parent_dir, scenario,'hysea_out','filterkajiura_res100_ts.nc')
    print(hysea_ts)
    # Run first step of the postprocessing that analyses the HySEA time series
    command = f"./postprocess_HySea {hysea_ts}"
    os.system(command)
    # Run the second step of the postprocessing that computes the inundation
    outfile_name = os.path.join(parent_dir, scenario, 'hysea_out', 'MIH.txt')
    ampfact_func(ampfact_file, hysea_ts, outfile_name)


