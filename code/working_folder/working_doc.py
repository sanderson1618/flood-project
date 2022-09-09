#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 08:43:12 2021

@author: SamAnderson
"""
from landlab.components import OverlandFlow, SinkFiller
from landlab.io import read_esri_ascii
from landlab import imshow_grid
import numpy as np
import copy
from matplotlib import pyplot as plt
#from landlab import multi_grain
import multi_grain

storm_flag = "1hr_1000yrStorm" # Different storm intensities from NOAA estimates
basin_flag = 'LC3' # 'Square' or Long'

# Set up LC1 and LC3
if basin_flag == 'LC1': #needs to be reclipped so that outlet is on a boundry
    watershed_dem = 'lc1_dem_test.txt'
    (rmg, z) = read_esri_ascii(watershed_dem, name='topographic__elevation')
    rmg.set_watershed_boundary_condition(z, nodata_value = -9999.)
    imshow_grid(rmg, 'topographic__elevation', vmin = 0, cmap ='twilight_r', color_for_closed = "black")
    plt.show()
elif basin_flag == 'LC3':
    watershed_dem = 'lc3_dem_test.txt'
    (rmg, z) = read_esri_ascii(watershed_dem, name='topographic__elevation')
    rmg.set_watershed_boundary_condition(z, nodata_value = -9999.)
    imshow_grid(rmg, 'topographic__elevation', vmin = 0, cmap ='twilight_r', color_for_closed = "black")
    plt.show()
 
# Precipitation inputs, 10min duration (600s)
# Used the precip_intensity sheet, which gives depth(in)/time(hr), 
# The conversion from (in/hr) to (m/s) is 7.05556 * 10 ** -6d
if storm_flag == "1hr_1000yrStorm":
    precip_inhr = 13.3
    precip_ms = precip_inhr * (7.05556 * 10 ** -6) 
    storm_duration = 600. # Seconds
elif storm_flag == '1hr_100yrStorm':
    precip_mmhr = 9.56
    precip_ms = precip_mmhr * (7.05556 * 10 ** -6)
    storm_duration = 600. # incorrect
elif storm_flag == '1hr_10yrStorm':
    precip_mmhr = 2.98
    precip_ms = precip_mmhr * (7.05556 * 10 ** -6)
    storm_duration = 600. # incorrect
    
# Fills sinks in case there are some
sf = SinkFiller(rmg, routing='D4', apply_slope=True, fill_slope=1.e-5)
sf.fill_pits()

# A copy of z for plotting purposes
z_initial = copy.deepcopy(z) 

# Pretty sure I need these lines...
rmg['link']['surface_water__discharge'] = np.zeros(rmg.number_of_links)
rmg['node']['water_surface__slope'] = np.zeros(rmg.number_of_nodes)
rmg['node']['surface_water__depth'] = np.zeros(rmg.number_of_nodes)
rmg['node']['velocity_shear_stress'] = np.zeros(rmg.number_of_nodes)
rmg['node']['depth_slope_shear_stress'] = np.zeros(rmg.number_of_nodes)
rmg['node']['vertical_velocity'] = np.zeros(rmg.number_of_nodes)
rmg['node']['cross_stream_velocity'] = np.zeros(rmg.number_of_nodes)

# Generate slope map of Last Chance canyon
rmg.at_node['topographic__slope'] = rmg.calc_slope_at_node(elevs='topographic__elevation')
imshow_grid(rmg,'topographic__slope');
plt.show()

# Initialization
of = OverlandFlow(rmg, mannings_n=0.03, steep_slopes=True, alpha = 0.3)
#mg = multi_grain(rmg, dt = dt)

# Creates time control settings
elapsed_time = 1.0
model_run_time = 700.0 # s

## Lists for saving data
discharge_at_outlet = []
hydrograph_time = []

## Running the overland flow component.
while elapsed_time < model_run_time:
    
    # Setting the adaptive time step
    of.dt = of.calc_time_step()
    
    ## The storm starts when the model starts. While the elapsed time is less
    ## than the storm duration, we add water to the system as rainfall.
    if elapsed_time < (storm_duration):

        of.rainfall_intensity =  precip_ms

    ## Then the elapsed time exceeds the storm duration, rainfall ceases.
    else:

        of.rainfall_intensity = 0.0

    ## Generating overland flow based on the deAlmeida solution.
    of.overland_flow()
    #mg.multi_grain()
    
    ## Append time and discharge to their lists to save data and for plotting.
    print(elapsed_time)
    
    ## Updating elapsed_time
    elapsed_time += of.dt
    
"""
imshow_grid(rmg, 'shear_stress', plot_name = 'Shear Stress at ...',
            var_name = 'Shear Stress', var_units = 'Pa', grid_units = ('m','m'),
            cmap = 'Blues')
"""
imshow_grid(rmg, 'tau', plot_name='Water depth at time = 2 hr', 
            var_name='Water Depth', var_units='m', grid_units=('m', 'm'), 
            cmap = 'Blues')
