"""
This component simulates sediment transport of differently sized grains using the 
numerical model from Parker, 1990. This compoent is meant to be coupled with the 
overland flow component from Adams et al, 2017.

.. codeauthor:: angel and sam 
"""

import numpy as np
import pandas as pd
from landlab.components import OverlandFlowVarRoughnessRainfall, BedloadTransport
from landlab.grid.mappers import map_mean_of_link_nodes_to_link
from landlab import RasterModelGrid

"""Create a grid on which to calculate sediment transport."""
grid = RasterModelGrid((9, 9))

"""
The grid will need some data to provide the bedload component. 
To check the names of the fields that provide input to
the detachment ltd transport component, use the *input_var_names* class
property.
"""

"""
Create fields of data for each of these input variables starting with topography.
"""
grid.at_node['topographic__elevation'] = np.array([
    1.5, 1.5, 1.08, 1.08, 1.08, 1.08, 1.08, 1.5, 1.5,
    1.5, 1.5, 1.07, 1.07, 1.07, 1.07, 1.07, 1.5, 1.5,
    1.5, 1.5, 1.06, 1.06, 1.06, 1.06, 1.06, 1.5, 1.5,
    1.5, 1.5, 1.05, 1.05, 1.05, 1.05, 1.05, 1.5, 1.5,
    1.5, 1.5, 1.04, 1.04, 1.04, 1.04, 1.04, 1.5, 1.5,
    1.5, 1.5, 1.03, 1.03, 1.03, 1.03, 1.03, 1.5, 1.5,
    1.5, 1.5, 1.02, 1.02, 1.02, 1.02, 1.02, 1.5, 1.5,
    1.5, 1.5, 1.01, 1.01, 1.01, 1.01, 1.01, 1.5, 1.5,
    1.5, 1.5, 1.01, 1.01, 1, 1.01, 1.01, 1.5, 1.5])

"""
From slopes at nodes calculates slopes at links as averaged over neighboring patches
"""
grid.at_node['topographic__slope'] = grid.calc_slope_at_node(elevs='topographic__elevation') 
# imshow_grid(rmg,'topographic__slope');plt.show() # Uncomment to display the slopes

"""
Now, add channel roughness (Manning's coefficient), grain size location (index 
for each grain size distribution), and rainfall (mm/hr) info to the grid.
"""

grid.at_node['bed_grain_size_distribution__location'] = np.array([
    0, 0, 1, 1, 1, 1, 1, 0, 0,
    0, 0, 1, 1, 1, 1, 1, 0, 0,
    0, 0, 1, 1, 1, 1, 1, 0, 0,
    0, 0, 1, 1, 1, 1, 1, 0, 0,
    0, 0, 1, 1, 1, 1, 1, 0, 0,
    0, 0, 1, 1, 1, 1, 1, 0, 0,
    0, 0, 1, 1, 1, 1, 1, 0, 0,
    0, 0, 1, 1, 1, 1, 1, 0, 0,
    0, 0, 1, 1, 1, 1, 1, 0, 0])


"""
Reads the grain size distribution from an xlsx file in the working directory
"""
gsd = pd.read_excel('GSD.xlsx',sheet_name='GSD',skiprows=0).values

"""
Configuring all topographic and flow variables related required by overland
flow, all are dry at the beginning of this example.
"""
grid.add_zeros('surface_water__depth', at = 'node')
grid.add_zeros('surface_water__discharge', at='link')
grid.add_zeros('surface_water__discharge', at='node')

"""
Instantiate the `bedload transport` component to work on this grid, and
run it. In this simple case, we need to pass it a time step ('dt') with
units in seconds. Here, we use a dt of 3 to avoid rounding errors, and a
max dt of 3 to ensure more stable results. This model will run for 60000 
seconds (model_run_time), and the storm will run for the duration of the 
model.
"""

dtPrecision = 3         
max_dt = 1             
model_run_time = 60000  
storm_starts = 0.0         
storm_ends = 60000       

""" Instantiation """
of = OverlandFlowVarRoughnessRainfall(grid, steep_slopes=True)
of.run_one_step()
qb = BedloadTransport(grid)
BedloadTransport.run_one_step(dt=dtPrecision)
""" Defines variables to store data """
elapsed_time = 0.0            

"""
The previous and current velocities are contained in each time step. With 
them, the du/dt term can be calculated in the shear stresses calculations
"""
u0 = np.vstack((np.zeros_like(of._q),np.zeros_like(of._q)))   

while elapsed_time < model_run_time:
    """ Calculates the rainfall """
    if (elapsed_time >= storm_starts) and (elapsed_time <= storm_ends):
        grid['node']['rainfall_intensity'] = 6.9e-6
    else: # if elapsed time exceeds the storm duration, rainfall ceases.
        grid['node']['rainfall_intensity'] = 0
   
    """ Runs overland flow for one time step """ 
    of.overland_flow()  # Generating overland flow based on the deAlmeida solution.
    qb._dt = of.dt      # Assigns overland flow dt to calculate time derivatives
    
    """ Updates u0 """   
    u0[0,:] = u0[1,:]               # What was the current discharge is now the previous
    u0[1,:] = of._q / of._h_links   # Now the velocity is updated    
    
    """ Runs bedload transport for one time step """ 
    qb.run_one_step(u0, gsd)

    ## Updating elapsed_time
    elapsed_time = round(elapsed_time + of.dt, dtPrecision)
    
"""
References
----------
**Required Software Citation(s) Specific to this Component**

None Listed

**Additional References**

G. Parker (1990): Surface-based bedload transport relation for gravel rivers,
Journal of Hydraulic Research, 28:4, 417-436,
http://dx.doi.org/10.1080/00221689009499058
"""