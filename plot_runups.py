"""
Script to plot MIH values of the HySEA simulation, after being postprocessed with the amplification factor scripts.
This script reads the HySEA time series output, the landslide simulation results, the historical runup data and plots them on a map in one figure with multiple panels.
"""

import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import cartopy
from netCDF4 import Dataset
from math import radians, cos, sin, asin, sqrt


def read_bathymetry(bathyfile):
    """
    Read bathymetry data from a NetCDF file and return longitude, latitude, and depth arrays.
    - bathyfile: name of the NetCDF file containing bathymetry data.
    """
    with Dataset(bathyfile) as ds:
        bathylon = ds.variables["lon"][:]
        bathylat = ds.variables["lat"][:]
        bathyz = ds.variables["z"][:]
    
    return bathylon, bathylat, bathyz

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

def read_bingclaw_results(bingclaw_sim_file):
    """
    Read the output of one time step of the BingClaw landslide simulation
    - bingclaw_sim_file: name of the BingClaw output (one time step).
    """
    with open(bingclaw_sim_file) as f:
        header_tmp = [next(f).strip() for _ in range(8)]
        # Now read the rest as data
        bingclaw_res_tmp = np.loadtxt(f)
    header={}
    header['mx'] = int(header_tmp[2].split()[0])
    header['my'] = int(header_tmp[3].split()[0])
    header['xlow'] = float(header_tmp[4].split()[0])
    header['ylow'] = float(header_tmp[5].split()[0])
    header['dx'] = float(header_tmp[6].split()[0])
    header['dy'] = float(header_tmp[7].split()[0])  

    lon = header['xlow'] + np.arange(header['mx']) * header['dx']
    lat = header['ylow'] + np.arange(header['my']) * header['dy']
    X, Y = np.meshgrid(lon, lat)

    bingclaw_res = bingclaw_res_tmp[:,0].reshape(header['my'],header['mx'])

    return header, X, Y, bingclaw_res

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a))
    # Radius of earth in kilometers is 6371
    km = 6371* c
    return km

def find_closest_poi(data_lon, data_lat, sim_lon, sim_lat, data_mih, sim_mih):
    """
    Find the closest simulation points to the data points and compute the difference in MIH.
    - data_lon, data_lat: coordinates of the historical runup data points.
    - sim_lon, sim_lat: coordinates of the POIs used in the T-HySEA simulation.   
    - data_mih: array of MIH values (m) for the historical runup data points.
    - sim_mih: array of MIH values (m) from the simulations.
    
    Output: a dictionary with the coordinates of the data points, the coordinates of the POIs, and the difference in MIH (m).
    """

    data_sim_comp = {
        'data_lon': [],
        'data_lat': [],
        'sim_lon': [],
        'sim_lat': [],
        'delta_mih': []
    }
    i=0
    for dlon, dlat, dmih in (zip(data_lon, data_lat, data_mih)):
        # Calculate distances to all simulation points
        distances = [haversine(dlon, dlat, slon, slat) for slon, slat in zip(sim_lon, sim_lat)]
        # Find the index of the closest point
        closest_index = np.argmin(distances)
        # # Print the closest point's coordinates and distance
        # print(f"Closest point to data point ({dlon}, {dlat}) is at "
        #       f"({sim_lon[closest_index]}, {sim_lat[closest_index]}) with distance {distances[closest_index]:.2f} km")
        if distances[closest_index] < 5:  # threshold of 5 km
            data_sim_comp['data_lon'].append(dlon)
            data_sim_comp['data_lat'].append(dlat)
            data_sim_comp['sim_lon'].append(sim_lon[closest_index])
            data_sim_comp['sim_lat'].append(sim_lat[closest_index])
            data_sim_comp['delta_mih'].append(dmih - sim_mih[closest_index])

    return data_sim_comp

def plot_one_simulation(simulation, parent_dir, do_save_figures,
                        lon2d, lat2d, bathyz,
                        data_lon, data_lat, h):
    """
    Plot the results of one simulation in 4 subplots: 
    1. Messina historical runup data
    2. HySEA simulation results (MIH)
    3. Difference between data and simulation results
    4. Zoom in of the landslide location at t0 and t90

    - simulation: ID of the simulation
    - parent_dir: path to the parent directory where all simulation results are stored
    - do_save_figures: boolean to decide whether to save the figures or show them
    - lon2d, lat2d, bathyz: arrays of longitude, latitude, and bathymetry data
    - data_lon, data_lat, h: arrays of longitude, latitude, and runup values of historical runup data

    Output: one figure saved in the parent_dir/simulation/hysea_out folder.
    """

    # ======== INPUT ======== #
    # Read Bingclaw simulation results
    bingclaw_sim_file_0 = (parent_dir + simulation + r"\bingclaw_out\fort.q0000")
    bingclaw_sim_file_90 = (parent_dir + simulation + r"\bingclaw_out\fort.q0090")
    header, X, Y, bingclaw_res_0 = read_bingclaw_results(bingclaw_sim_file_0)
    _,_,_, bingclaw_res_90 = read_bingclaw_results(bingclaw_sim_file_90)

    # Read HySEA simulation results
    hysea_sim_file = (parent_dir + simulation + r"\hysea_out\MIH.txt")
    id, sim_lon, sim_lat, MIH = read_hysea_results(hysea_sim_file)

    # ======== PLOT ======== #
    proj = cartopy.crs.PlateCarree()  # Define the projection for the map
    fig = plt.figure(figsize=(16, 8))
    
    ### SUBPLOT 1: PLOT MESSINA DATA ###
    ax = fig.add_subplot(1, 4, 1, projection=cartopy.crs.Mercator())

    # Coast line from bathymetry #
    contour = ax.contour(
    lon2d, lat2d, bathyz, levels=[0], colors='k', linewidths=1,
    transform=cartopy.crs.PlateCarree())
    # Run up data #
    sc = ax.scatter(
        data_lon, data_lat, c=h, cmap='Blues', s=40, edgecolor='k',
        transform=cartopy.crs.PlateCarree(), zorder=10, vmax=10)
    
    # Plot settings #
    plt.colorbar(sc, ax=ax, orientation='horizontal', label='Runup height (m)')
    ax.set_extent([np.min(lon2d), np.max(lon2d), 36.6, np.max(lat2d)], crs=cartopy.crs.PlateCarree())
    gl = ax.gridlines(crs=proj, draw_labels=True, linewidth=1,
                    color="grey", alpha=0.4, linestyle='-')
    gl.top_labels = False
    gl.right_labels = False
    gl.bottom_labels = True
    gl.left_labels = True
    ax.set_title('Runup DATA', fontsize=14)


    ### SUBPLOT 2: PLOT HYSEA SIMULATION RESULTS ###
    ax2 = fig.add_subplot(1, 4, 2, projection=cartopy.crs.Mercator())

    # Coast line from bathymetry and bathymetry contours #
    contour = ax2.contour(lon2d, lat2d, bathyz, levels=[0], colors='k', linewidths=1,
                            transform=cartopy.crs.PlateCarree())
    con3 = ax2.contour(lon2d, lat2d, bathyz, levels=np.arange(-2500,-499,500), colors='k',linewidths=0.5, alpha=0.5, transform=cartopy.crs.PlateCarree())
    ax2.clabel(con3, fmt='%d', inline=True, fontsize=6, manual=False, rightside_up=True, use_clabeltext=True, 
           inline_spacing=1)
    # MIH tsunami simulation results #
    sc2 = ax2.scatter(
        sim_lon, sim_lat, c=MIH, cmap='Blues', s=40, edgecolor='k',
        transform=cartopy.crs.PlateCarree(), zorder=10, vmax=np.max(MIH))    
    # Bingclaw landslide location at t0 and t90 #
    masked_bingclaw_res_0 = np.ma.masked_where(bingclaw_res_0 == 0, bingclaw_res_0)
    masked_bingclaw_res_90 = np.ma.masked_where(bingclaw_res_90 == 0, bingclaw_res_90)
    con = ax2.contourf(X, Y, masked_bingclaw_res_0, alpha=1, cmap="Greys", 
                       transform=cartopy.crs.PlateCarree())
    con2 = ax2.contourf(X, Y, masked_bingclaw_res_90, alpha=1,cmap="Purples", 
                        transform=cartopy.crs.PlateCarree())

    # Find the bounding box of nonzero values in masked_bingclaw_res_90 for zoom in plot
    nonzero_indices = np.argwhere(masked_bingclaw_res_90.mask == False)
    if nonzero_indices.size > 0:
        y_min, x_min = nonzero_indices.min(axis=0)
        y_max, x_max = nonzero_indices.max(axis=0)
        x0, x1 = X[y_min, x_min], X[y_max, x_max]
        y0, y1 = Y[y_min, x_min], Y[y_max, x_max]
        # Add some padding
        pad_x = (x1 - x0) * 0.2
        pad_y = (y1 - y0) * 0.2
        x0, x1 = x0 - pad_x, x1 + pad_x
        y0, y1 = y0 - pad_y, y1 + pad_y
    else:
        # fallback to a default area if all masked
        x0, x1 = X.min(), X.max()
        y0, y1 = Y.min(), Y.max()
    
    # Draw rectangule of zoom in area
    rect = mpatches.Rectangle(
        (x0, y0), x1 - x0, y1 - y0,
        linewidth=1, edgecolor='orange', facecolor='none', zorder=20,
        transform=cartopy.crs.PlateCarree()
    )
    ax2.add_patch(rect)

    # Plot settings #
    plt.colorbar(sc2, ax=ax2, orientation='horizontal', label='Runup height (m)')
    ax2.set_extent([np.min(lon2d), np.max(lon2d), 36.6, np.max(lat2d)], crs=cartopy.crs.PlateCarree())
    gl = ax2.gridlines(crs=proj, draw_labels=True, linewidth=1,
                    color="grey", alpha=0.4, linestyle='-')
    gl.top_labels = False
    gl.right_labels = False
    gl.bottom_labels = True
    gl.left_labels = False
    ax2.set_title('Runup SIMULATION', fontsize=14)
    

    ### SUBPLOT 3: PLOT DIFFERENCE BETWEEN DATA and SIMULATION RESULTS ###
    # Compute the difference between data and simulation results
    # Find the closest simulation points to the data points and compute difference
    data_sim_comp = find_closest_poi(data_lon, data_lat, sim_lon, sim_lat, h, MIH)
    
    # PLOT DIFFERENCE ON THE MAP - SUBPLOT 1 #
    ax4 = fig.add_subplot(1, 4, 3, projection=cartopy.crs.Mercator())

    # Coast line from bathymetry #
    contour = ax4.contour(
    lon2d, lat2d, bathyz, levels=[0], colors='k', linewidths=1,
    transform=cartopy.crs.PlateCarree())

    # Run up data - simulation MIH #
    sc = ax4.scatter(data_sim_comp['data_lon'], data_sim_comp['data_lat'], 
                    c=data_sim_comp['delta_mih'], cmap='Reds', s=40, edgecolor='k',
                transform=cartopy.crs.PlateCarree(), zorder=10)

    # Plot settings #
    plt.colorbar(sc, ax=ax4, orientation='horizontal', label='Runup data - simulation (m)')
    ax4.set_extent([np.min(lon2d), np.max(lon2d), 36.6, np.max(lat2d)], crs=cartopy.crs.PlateCarree())
    gl = ax4.gridlines(crs=proj, draw_labels=True, linewidth=1,
                    color="grey", alpha=0.4, linestyle='-')
    gl.top_labels = False
    gl.right_labels = False
    gl.bottom_labels = True
    gl.left_labels = False
    ax4.set_title('Runup difference DATA-SIMULATION', fontsize=14)


    ### SUBPLOT 4: PLOT ZOOM IN of the LANDSLIDE ###
    ax3 = fig.add_subplot(1, 4, 4, projection=cartopy.crs.Mercator())

    # Plot Bingclaw landslide location at t0 and t90 #
    c1 = ax3.contourf(X, Y, masked_bingclaw_res_0, alpha=1, cmap="Greens", transform=cartopy.crs.PlateCarree())
    c2 = ax3.contourf(X, Y, masked_bingclaw_res_90, alpha=1, cmap="Purples", transform=cartopy.crs.PlateCarree())
    # Coast line from bathymetry and bathymetry contours#
    ax3.contour(lon2d, lat2d, bathyz, levels=[0], colors='k', linewidths=1, transform=cartopy.crs.PlateCarree())
    contour2 = ax3.contour(lon2d, lat2d, bathyz, levels=np.arange(-2400,-199,200), colors='k',linewidths=0.5, alpha=0.5, transform=cartopy.crs.PlateCarree())
    ax3.clabel(contour2, fmt='%d', inline=True, fontsize=6, manual=False, rightside_up=True, use_clabeltext=True, 
           inline_spacing=0.005)
    # MIH tsunami simulation results #
    ax3.scatter(sim_lon, sim_lat, c=MIH, cmap='Blues', s=10, edgecolor='k', transform=cartopy.crs.PlateCarree(), zorder=10, vmax=np.max(MIH))
    
    # Plot settings #
    gl = ax3.gridlines(crs=proj, draw_labels=True, linewidth=1,
                    color="grey", alpha=0.2, linestyle='-')
    gl.top_labels = False
    gl.right_labels = True
    gl.bottom_labels = True
    gl.left_labels = False
    ax3.set_extent([x0, x1, y0, y1], crs=cartopy.crs.PlateCarree())
    ax3.set_title('Landslide location \n at t0 and t 15min', fontsize=14)
    plt.colorbar(c1, ax=ax3, orientation='horizontal', label='H at t0 (m)')
    plt.colorbar(c2, ax=ax3, orientation='horizontal', label='H at t15min (m)')


    # ======== OUTPUT ======== #
    if do_save_figures:
        fig_name = parent_dir + simulation + r"\hysea_out\runups_" + simulation + ".png"
        fig.savefig(fig_name, bbox_inches='tight', dpi=300)
        print(f"Run up data and simulation figure saved as {fig_name}")
        plt.close()
    else:
        plt.show()


### MAIN FUNCTION
def plot_runups():
    ### Load all necessary files
    # Directory where the simulations are stored (e.g., parent_dir/hysea_out/filterkajiura_res100_ts.nc)
    # The structure of the simulations output folder is defined during the run of the bingclaw-to-hysea wrorkflow
    parent_dir = "path/to/simulations_output_directory/"
    do_save_figures = True

    # Directory where the volume representatives CSV file is located
    volumes_dir = "path/to/volumes_directory/" # e.g., "messina_001/volumes/"
    filename = os.path.join(volumes_dir,"volume_representatives.csv")
    # Load volumes IDs (these IDs correspond to the name of scenario folder where the landslide and tsunami simulation results are stored)
    df = pd.read_csv(filename, usecols=["id"])

    # Read bathymetry file #
    bathyfile = "data/GEBCO2024_100m_15_1650_36_3840.nc"
    bathylon, bathylat, bathyz = read_bathymetry(bathyfile)
    lon2d, lat2d = np.meshgrid(bathylon, bathylat)

    # Read Messina historical runups data #
    filename = "data/messina_runups_data.npz"
    runup_data = np.load(filename)

    ### Loop over each simulation to read the outputs and plot the results
    #idx_to_run = [7,8,10,11,19,24,25,31,38,46,48,50,52]
    #for i in idx_to_run:
    for i in range(0,len(df)):
        simulation = str(int(df["id"][i]))
        print(f"Plotting simulation {simulation}")
        plot_one_simulation(simulation, parent_dir, do_save_figures ,
                            lon2d, lat2d, bathyz,
                            runup_data["lon"], runup_data["lat"], runup_data["h"])


if __name__ == '__main__':
    plot_runups()