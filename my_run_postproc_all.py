"""
Script to run postprocessing of all HySEA simulations to compute maximum inundation height (MIH):
- postprocess_HySea is C++ code that extract info from time series
- ampfact_func.py is a script that uses the output of postprocess_HySea and the amplification factors to compute MIH
"""

import os
import sys
import pandas as pd
from ampfact_func import ampfact_func

parent_dir = "/home/vmg/NGI/P/2022/02/20220296/Calculations/Messina_simulations_August2025//"
volumes_dir = "/home/vmg/NGI/P/2022/02/20220296/Calculations/release-volume-sampler_20250806_kl1558/messina_20250806/volumes/"
filename = os.path.join(volumes_dir,"volume_representatives.csv")

# File with values of the amplification factors
ampfact_file="factors/ampf_messina.txt"

## Read the clusters CSV file to get the scenario IDs
## Use the ID to find the HySEA time series file and run the amp factor postprocessing
df = pd.read_csv(filename, usecols=["id"])
#idx_to_run = [7,8,10,11,19,24,25,31,38,46,48,50,52]

for i in range(0, len(df)):
    #for i in idx_to_run:
    scenario = str(int(df["id"][i]))
    hysea_ts = os.path.join(parent_dir, scenario,'hysea_out','filterkajiura_res100_ts.nc')
    print(hysea_ts)
    # Run the first step of the postprocessing that analyses the HySEA time series
    command = f"./postprocess_HySea {hysea_ts}"
    os.system(command)
    # Run the second step of the postprocessing that computes the inundation
    outfile_name = os.path.join(parent_dir, scenario, 'hysea_out', 'MIH.txt')
    ampfact_func(ampfact_file, hysea_ts, outfile_name)


