# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 15:57:29 2021

@author: Sam
"""

from landlab.components.uniform_precip import PrecipitationDistribution
from landlab.components import SpatialPrecipitationDistribution
import numpy as np
import random
import matplotlib.pyplot as plt
from landlab.components import OverlandFlow
from landlab.plot.imshow import imshow_grid
from landlab.io import read_esri_ascii, write_esri_ascii
import os 

#initial topo and show it
(mg, z) = read_esri_ascii("steady_state_topography.txt", 
                          name="topographic__elevation")

plt.figure()
imshow_grid(mg, z, colorbar_label='Elevation (m)')
plt.show()

#storm in same locale and show it
rain = SpatialPrecipitationDistribution(mg)
"""
np.random.seed()
total_t_each_step = [
    (storm+interstorm) for (storm, interstorm) in rain.yield_storms()]
len(total_t_each_step)
np.isclose(sum(total_t_each_step)/24., 365.)
"""
#uniform rain
morerain = SpatialPrecipitationDistribution(mg)
np.random.seed(26)  # arbitrary to get a cool-looking storm out every time

# get the storm simulator to provide a storm
# There's only one storm generated here in the time series, so easy enough to do.
# first, check the directory we need for saving exists, and make it if not:
if not os.path.exists('./rainfall'):
    os.makedirs('./rainfall')
for (storm_t, interstorm_t) in rain.yield_storms(style='monsoonal'):  # storm lengths in hrs
    mg.at_node['rainfall__flux'] *= 0.001  # because the rainfall comes out in mm/h
    mg.at_node['rainfall__flux'] *= 10.  # to make the storm heavier and more interesting!
    plt.figure()
    # plot up this storm
    write_esri_ascii('./rainfall/rainfall.asc', mg, 'rainfall__flux', clobber=True)
 
"""
if not os.path.exists('./rainfall'):
    os.makedirs('./rainfall')
for (storm_t, interstorm_t) in rain.yield_storms(style='monsoonal'):  # storm lengths in hrs
    mg.at_node['rainfall__flux'] *= 0.001  # because the rainfall comes out in mm/h
    mg.at_node['rainfall__flux'] *= 10.  # to make the storm heavier and more interesting!
    plt.figure()
    # plot up this storm
    
"""


"""
mytotals = []
for yr in range(5):
     mytotals.append(rain.calc_annual_rainfall(style='whole_year'))
"""
#plot this storm
imshow_grid(
        mg, 'rainfall__flux', cmap='gist_ncar', colorbar_label='Rainfall flux (m/h)'
    )
plt.show()