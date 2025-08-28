# **Postprocessing of T-HySEA output using amplification factors**
 
This repo contains the codes to execute the postprocessing of the T-HySEA simulation results and use the amplification factors to compute the Maximum Inundation Height (MIH) at specific Points of Interest (POIs).
      
It is the last step of the landslide-tsunami-workflow that produces the pre-computed database of landslide and tsunami simulations for the 1908 Messina event (but it can be applied to other areas). Therefore, it assumes that the results of the landslide and tsunami simulations are already available.


The steps to run the postprocessing are:
- The main script `run_postproc_all.py` reads the file with the ID of the scenarios to run the postprocessing of and launches the commands and run the functions to run the 2 steps of the postprocessing:
  - `postprocess_HySea` is C++ code that extracts info from time series (output of T-HySEA). It produces files with .offshore.txt (and .nc) extension that describe some key parameters of the wave. The code is in the directory `postprocess_cpp`. 
  - `ampfact_func.py` contains a python function that uses the results of the previous step and the amplification factors to compute the MIH. It produces a MIH.txt file with the inundation heights for each POI.   

  The output of the main script is a 2D array with all MIH values of all scenarios at each POI saved in the file `MIH_all.npy`. This is the main output of the precomputed database for the simulations part and it is used in the PTF landslide workflow.  
- (optional) Plot the MIH values, the results of the landslide simulation (output at the first and last time step), and the comparison with the historical runup data for each scenario with the script `plot_runups.py` 


## Requirements
### Python packages
The required python packages are in the `requirements.txt` file. The file `pyproject.toml` is also provided if Poetry is used instead of pip to install the python packages. Tested with python>=3.10.  
### Compile the C++ code
The executable `postprocess_HySea` created on the NGI cluster David is provided. However, each cluster has different compilers and settings, so unless you are running this on David@NGI, you need to recompile the C++ code to generate the executable from the `postprocess_cpp` directory. To do that:
```
cd postprocess_cpp
./make.postprocess_HySea.sh
```
Note that `make.postprocess_HySea.sh` assumes version 11 of the gcc compiler is installed. This should be changed accordingly to the gcc version used (or do not specify the version).

## Run
Open `run_postproc_all.py` and write the right paths of the directories with the simulation results and amplification factors file. This script reads the .csv file containing the list of scenarios produced by release-volume-sampler/preparational.py. 

```
python run_postproc_all.py
```
This will generate a `MIH.txt` file stored inside the directory where the T-HySEA outputs are for each scenario and a 2D array with MIH values of all scenarios at each POI saved in the file `MIH_all.npy`.  

Now, it is possible to plot the results (optional), by running the script `plot_runups.py`. First, check inside the script that the paths of the directories are correct, then run:
```
python plot_runups.py
```
This will produce a figure for each scenario stored where the T-HySEA outputs are with the following plots: 1908 Messina event historical runup data, MIH from the tsunami simulation at each POI, comparison between MIH and historical runup data, landslide thickness at the first and last time step of the landslide simulation. 

## Acknowledgments
This work is made in the context of the [DT-GEO](https://dtgeo.eu/) project.     

The amplification factors method is described in:   
Glimsdal, S., Løvholt, F., Harbitz, C. B., Romano, F., Lorito, S., Orefice, S., ... & Omira, R. (2019). A new approximate method for quantifying tsunami maximum inundation height probability. Pure and Applied Geophysics, 176(7), 3227-3246. https://doi.org/10.1007/s00024-019-02091-w  

The code that executes the postprocessing of the time series of T-HySEA is by Andrey Babeyko (GFZ). The amplification factors scripts are by Sylfest Glimsdal (NGI). Valentina Magni (NGI) added scripts to run the workflow automatically and plot the runups.

## License
This project is licensed under the Non-Profit Open Software License version 3.0 (NPOSL-3.0).
See the [LICENSE](./LICENSE.txt) file for details.
