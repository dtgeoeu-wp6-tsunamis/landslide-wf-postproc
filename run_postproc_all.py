"""
Script to run postprocessing of all HySEA simulations to compute maximum inundation height (MIH):
- postprocess_HySea is C++ code that extract info from time series
- ampfact_func.py is a script that uses the output of postprocess_HySea and the amplification factors to compute MIH
Output: MIH_all.npy file with a 2D array with MIH values for all scenarios at all POIs
"""

import os
import pandas as pd
import numpy as np
from ampfact_func import ampfact_func

def read_hysea_results(sim_file):
    """
    Read simulation results after amplification factors postprocessing and return lists of id, longitude, latitude, and MIH values.
    - sim_file: name of file with MIH values at each POI for a specific scenario computed by the amplification factor postprocessing.
    """
    with open(sim_file, "r") as f:
        lines = f.readlines()[1:]  # skip first line

    id = []
    sim_lon = []
    sim_lat = []
    MIH = []

    for line in lines:
        parts = line.strip().split()
        if len(parts) == 6:
            id.append(parts[0])
            sim_lon.append(float(parts[1]))
            sim_lat.append(float(parts[2]))
            MIH.append(float(parts[3]))

    return id, sim_lon, sim_lat, MIH

def run_postproc_all():
    ### Define name and paths of input files and directories
    # Directory where the simulations are stored (e.g., parent_dir/hysea_out/filterkajiura_res100_ts.nc)
    # The structure of the simulations output folder is defined during the run of the bingclaw-to-hysea wrorkflow
    parent_dir = "path/to/simulations_output_directory/"

    # Directory where the volume representatives CSV file is located
    # This file is created during the run of the release-volume-sampler/preparational.py workflow
    volumes_dir = "path/to/volumes_directory/" # e.g., "messina_001/volumes/"
    filename = os.path.join(volumes_dir,"volume_representatives.csv")

    # File with values of the amplification factors
    ampfact_file="factors/ampf_messina.txt"

    ### Read the clusters CSV file to get the scenario IDs and
    ### Use the ID to find the HySEA time series file and run the amp factor postprocessing
    df = pd.read_csv(filename, usecols=["id"])

    ### Run the postprocessing for each scenario
    #idx_to_run = [7,8,10,11,19,24,25,31,38,46,48,50,52]
    #for i in idx_to_run:
    for i in range(0, len(df)):
        scenario = str(int(df["id"][i]))
        # Name of T-HySEA time series output file
        hysea_ts = os.path.join(parent_dir, scenario,'hysea_out','filterkajiura_res100_ts.nc')
        print(hysea_ts)
        # Run the first step of the postprocessing that analyses the HySEA time series
        command = f"./postprocess_HySea {hysea_ts}"
        os.system(command)
        # Run the second step of the postprocessing that computes the inundation
        outfile_name = os.path.join(parent_dir, scenario, 'hysea_out', 'MIH.txt')
        ampfact_func(ampfact_file, hysea_ts, outfile_name)


    ### Read the MIH results after amplification factors postprocessing for each scenarios,
    ### store them all in a 2D array with shape (n_scenarios,n_pois) and save it in a file
    n_scenarios = len(df)
    n_pois = 154  # Number of Points of Interest (POIs) (can be found in factors/ampf_messina.txt or hysea outputs)

    MIH_all = np.zeros((n_scenarios,n_pois))
    for i in range(0,n_scenarios): 
        scenario = str(int(df["id"][i]))
        print(f"Processing scenario {scenario}")

        # Read HySEA simulation results after amplification factors have been applied
        hysea_sim_file = os.path.join(parent_dir, scenario, 'hysea_out', 'MIH.txt')
        id, sim_lon, sim_lat, MIH = read_hysea_results(hysea_sim_file)
        MIH_all[i,:] = MIH

    np.save("MIH_all.npy",MIH_all)

if __name__ == '__main__':
    run_postproc_all()