#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 08:43:12 2021

@author: angelmonsalve
"""
from landlab.components import OverlandFlow, DetachmentLtdErosion, SinkFiller
from landlab.components.uniform_precip import PrecipitationDistribution
from landlab.io import read_esri_ascii
from landlab import imshow_grid
import numpy as np
import copy
from matplotlib import pyplot as plt

max_dt=10  # Using only the default time step calculator gives some unstable results. This low value ensures more stables results use it in seconds
tPlot = 500 # Just to get a plot every this seconds

max_dt = 0.5  
tPlot = 500 # seconds

basin_flag = 'sam_dem' # 'Square' or Long'
storm_flag = 'Base' # 'Base' or'HigherIntensity' or 'LongerDuration'

if basin_flag == 'Square':
    watershed_dem = '/Users/safiya/Documents/LocalPy/Square_TestBasin.asc'
    link_to_sample = 299
    node_to_sample = 300
elif basin_flag == 'Long':
    watershed_dem = 'Long_TestBasin.asc'
    link_to_sample = 149
    node_to_sample = 150
elif basin_flag == 'sam_dem':
    watershed_dem = 'steady_state_topography'
    link_to_sample = 1493
    node_to_sample = 1494
    
(rmg, z) = read_esri_ascii(watershed_dem, name='steady_state_topography')
imshow_grid(rmg, 'topographic__elevation', vmin=0, cmap='twilight_r')  # plots the DEM
plt.show()

rmg.set_watershed_boundary_condition_outlet_id((1493), z, -9999.)
z = np.where(z == 0, -9999.0, z)
sf = SinkFiller(rmg, routing='D4', apply_slope=True, fill_slope=1.e-5)
sf.fill_pits()
z_initial = copy.deepcopy(z)
rmg.set_watershed_boundary_condition(z)
plt.show()

# Fills sinks in case there are some
sf = SinkFiller(rmg, routing='D4', apply_slope=True, fill_slope=1.e-5)
sf.fill_pits()
z_initial = copy.deepcopy(z)                # A copy of z for plotting purposes

"""
Now we set the boundary conditions, instantiate the components, and set the appropriate storm parameters.
All NODATA nodes in the DEM are closed boundaries and the outlet is set to an open boundary.
This is all done in rmg.set_watershed_boundary_condition(z).
"""
rmg.set_watershed_boundary_condition(z)

"""
First, what components are requiered?
in OverlandFlow we need ('surface_water__depth', 'topographic__elevation')
in DetachmentLtdErosion we need ('surface_water__discharge', 'topographic__elevation')
This can be checked using the next commands
"""
OverlandFlow.input_var_names
DetachmentLtdErosion.input_var_names
rmg.add_zeros('surface_water__depth', at = 'node')      # All cells are dry at the beginning of the simulation
rmg.add_zeros('surface_water__discharge', at='link')    # All cells are dry at the beginning of the simulation - no discharge

# I don't know if these variables are used
rmg['node']['surface_water__discharge'] = np.zeros(rmg.number_of_nodes)
rmg['node']['water_surface__slope'] = np.zeros(rmg.number_of_nodes)

# Also, DetachmentLtdErosion requires the calculation of the 'topographic__slope'
rmg.at_node['topographic__slope'] = rmg.calc_slope_at_node(elevs='topographic__elevation')
imshow_grid(rmg,'topographic__slope');
plt.show()

# Initialization
of = OverlandFlow(rmg, mannings_n=0.03, steep_slopes=True, alpha = 0.3)   # instantiate OverlandFlow object
dle = DetachmentLtdErosion(rmg, K_sp = 1.259162261 * (10**-7)) #instantiate DetachmentLtdErosion object
"""
if storm_flag == 'Base':
    starting_precip_mmhr = 5.0
    starting_precip_ms = starting_precip_mmhr * (2.77778 * 10 ** -7)
    storm_duration = 7200.
elif storm_flag == 'HigherIntensity':
    starting_precip_mmhr = 20.0
    starting_precip_ms = starting_precip_mmhr * (2.77778 * 10 ** -7)
    storm_duration = 1200.
elif storm_flag == 'LongerDuration':
    starting_precip_mmhr = 5.0
    starting_precip_ms = starting_precip_mmhr * (2.77778 * 10 ** -7)
    storm_duration = 14400.
"""

uplift_rate = 3.170979 * (10**-10) # m/s

# Creates time control settings
elapsed_time = 1.0
model_run_time = 20000.0 # s

## Lists for saving data
discharge_at_outlet = []
hydrograph_time = []
incision_at_outlet = []



while elapsed_time < model_run_time:
    # Setting the adaptive time step
    of.dt = min(max_dt,of.calc_time_step())

    ## The storm starts when the model starts. While the elapsed time is less
    ## than the storm duration, we add water to the system as rainfall.
    
    np.random.seed(np.arange(10))
precip = PrecipitationDistribution(mean_storm_duration=1.5,
     mean_interstorm_duration=15.0, mean_storm_depth=0.5,
     total_t=100.0, delta_t=1.)
for (dt, rate) in precip.yield_storm_interstorm_duration_intensity():
    imshow_grid(
        rmg, 'rainfall__flux', cmap='gist_ncar', colorbar_label='Rainfall flux (m/h)'
    )

    plt.show()

    """
    if elapsed_time < (storm_duration):
        of.rainfall_intensity =  starting_precip_ms  
    else: # elapsed time exceeds the storm duration, rainfall ceases.
        of.rainfall_intensity = 0.0
"""
    of.overland_flow()                      # Generating overland flow based on the deAlmeida solution.
    of.dt = min(max_dt,of.calc_time_step())      # of.overland_flow changes the dt. This is to ensure that we are using the same in both components
    
    ## Mapping water discharge from links (m^2/s) to nodes (m^3/s) for use
    ## in the DetachmentLtdErosion component.
    #node_slope = of.discharge_mapper(rmg['link']['surface_water__discharge'])
   
    ## Mapping water discharge from links (m^2/s) to nodes (m^3/s) for use
    ## in the DetachmentLtdErosion component.
    node_slope = (of._water_surface_slope[rmg.links_at_node] * rmg.active_link_dirs_at_node)
    incision_Q = np.abs(of._q * rmg.dx)[rmg.links_at_node]
    rmg['node']['surface_water__discharge'] = (incision_Q[np.arange(len(node_slope)), np.argmax(node_slope, axis=1)])
   
    ## Calculating water surface slope from the OverlandFlow component.
    node_slope = node_slope.max(axis=1)
    rmg['node']['water_surface__slope'] = node_slope

    ## Eroding topographic__elevation using DetachmentLtdErosion component.
    #dle.erode(of.dt, slope='water_surface__slope')
    dle.run_one_step(of.dt)

    ## Updating topographic__elevation after surface was eroded in
    ## DetachmentLtdErosion component.
    rmg['node']['topographic__elevation'] += uplift_rate * of.dt

    ## Append time and discharge to their lists to save data and for plotting.
    hydrograph_time.append(elapsed_time)
    discharge_at_outlet.append(np.abs(of._q[link_to_sample]) * rmg.dx)
    incision_at_outlet.append(dle._I[node_to_sample])

    ## If one wanted to see the time evolve, uncomment below, but this
    ## produces A LOT of output, as time steps are very short.
    tPlot = tPlot-of.dt
    if tPlot <= 0:
        print('Elapsed time :',elapsed_time,' s. Current dt =',of.dt,' s')
        imshow_grid(rmg, 'surface_water__depth',cmap='Blues',vmin=0,vmax=1);
        plt.show()
       
        #imshow_grid(rmg,'topographic__elevation',vmin=0)  # plots the DEM
        #plt.show()
        tPlot = 500

    ## Updating elapsed_time
    elapsed_time += of.dt

## Plotting the hydrograph at the outlet
plt.figure(1)
plt.plot(hydrograph_time, discharge_at_outlet, color='mediumblue')
plt.ylabel('Discharge (cms)')
plt.xlabel('Time (seconds)')
plt.title('Hydrograph at Watershed Outlet')

## Below is if you want to save a figure.
# files=plt.savefig('HydrographAtBasinOutlet.png')

## Plotting the incision rate through time at the outlet
plt.figure(2)
plt.plot(hydrograph_time, incision_at_outlet, color='darkred')
plt.ylabel('incision rate (m/s)')
plt.xlabel('Time (seconds)')
plt.title('Incision at Watershed Outlet')

## Below is if you want to save a figure.
# files=plt.savefig('HydrographAtBasinOutlet.png')

z_diff = z - z_initial # the difference in elevation from start to finish
imshow_grid(rmg, z_diff, limits=(-0.00004, np.max(z_diff)))

### this is just fluvial incision, which will be negative
z_diff = z - z_initial - uplift_rate * model_run_time
imshow_grid(rmg, z_diff, limits=(-0.00004, 0.0))