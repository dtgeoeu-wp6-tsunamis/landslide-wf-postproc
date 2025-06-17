# **Postprocessing scripts: Amplification factors and plotting** 
 
This repo contains the scripts to execute the postprocessing of the HySEA simulation results using the amplification factors and to plot the runups to compare it with the Messina event.

The amplification factors scripts are from Sylfest Glimsdal.
Valentina Magni added scripts to run the workflow automatically

# Amplification factor postprocessing
## Description of the scripts
The script run_postproc_all.py is the main one that runs all the steps of the postprocessing. It reads the clusters CSV file to get the scenario IDs, then uses the ID to find the HySEA time series file and run the amp factor postprocessing. 
The amp factor postprocessing is made in two steps: 
- postprocess_HySea is C++ code that extract info from time series: this produces files with .offshore.tct (and .nc) extension that describe some key parameters of the wave.
- ampfact_func.py is a script that uses the amplification factors to compute the inundation. This produces a MIH.txt file with the inundation heights for each POI.

## Run
Open the run_postproc_all.py and write the right paths of the directories with the simulation results and amp factor files. Everything is on P, so in principle you just need to change the starting of the path with your home directory, the rest should be the same

```
poetry run python run_postproc_all.py
```

# Plotting run ups
## Description of the script
The script is run_postproc_all.py. It reads in the MIH.txt file of a specific scenario, the bingclaw outputs (time step 0 and 90), the bathymetry file and the Messina data. It plots the data, the MIH of the simulation, the difference between the data and the simulation and the landslide location at t0 and t90.

## Run
Change the name of the scenario and check the paths to the files are correct. Everything is on P, but I have ran this on windows, so maybe the paths do not work in other machines as they are now.
```
poetry run python run_postproc_all.py
```
Note that this script at the moment works for one single scenario. This will need to be automatized to run it for all scenarios. 