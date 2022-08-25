# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 10:15:40 2021

@author: Sam
"""

from landlab.components.overland_flow import OverlandFlow
from landlab.plot.imshow import imshow_grid
from landlab.plot.colors import water_colormap
from landlab import RasterModelGrid
from landlab.io.esri_ascii import read_esri_ascii
from matplotlib.pyplot import figure
import numpy as np
from time import time

run_time = 300  # duration of run, (s)
h_init = 0.1  # initial thin layer of water (m)
n = 0.01  # roughness coefficient, (s/m^(1/3))
g = 9.8  # gravity (m/s^2)
alpha = 0.7  # time-step factor (nondimensional; from Bates et al., 2010)
u = 0.4  # constant velocity (m/s, de Almeida et al., 2012)
run_time_slices = (10, 50, 100, 300)

elapsed_time = 1.0

rmg, z = read_esri_ascii('Square_TestBasin.asc')
rmg.add_field('topographic__elevation', z, at='node')
rmg.set_closed_boundaries_at_grid_edges(True, True, True, True)

np.all(rmg.at_node['topographic__elevation'] == z)

my_outlet_node = 100  # This DEM was generated using Landlab and the outlet node ID was known
rmg.status_at_node[my_outlet_node] = 1  # 1 is the code for fixed value

rmg.add_zeros('surface_water__depth', at='node')  # water depth (m)

rmg.at_node['surface_water__depth'] += h_init

imshow_grid(rmg, 'topographic__elevation')

of = OverlandFlow(
    rmg, steep_slopes=True
)  #for stability in steeper environments, we set the steep_slopes flag to True

while elapsed_time < run_time:
    # First, we calculate our time step.
    dt = of.calc_time_step()
    # Now, we can generate overland flow.
    of.overland_flow()
    # Increased elapsed time
    print('Elapsed time: ', elapsed_time)
    elapsed_time += dt
    
imshow_grid(rmg, 'surface_water__depth', cmap='Blues')

elapsed_time = 1.
for t in run_time_slices:
    while elapsed_time < t:
        # First, we calculate our time step.
        dt = of.calc_time_step()
        # Now, we can generate overland flow.
        of.overland_flow()
        # Increased elapsed time
        elapsed_time += dt
    figure(t)
    imshow_grid(rmg, 'surface_water__depth', cmap='Blues')