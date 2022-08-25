"""
NOTES:
date: 12/27
editor: Sam
objective: to start building "multi_grain" component and attept to find depth-slope product
while coupled to overland_flow component


INSERT DESCRIPTION

.. codeauthor:: ANGEL AND SAM

Examples
INSERT EXAMPLE

References
----------
**Required Software Citation(s) Specific to this Component**

INSERT REFERENCES 

ie Monsalve et al, Yager, Parker, probably some other science people etc.
"""

import numpy as np

from landlab import Component


class MultiGrain(Component):

    """Landlab component that simulates sediment transport

    This component calculates size and amount of sediment removed with different storms
    """

    _name = "MultiGrain"

    _unit_agnostic = True

    _info = {
        "surface_water__depth": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Volumetric discharge of surface water",
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
        "topographic__slope": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": "Gradient of the ground surface",
        },
    }

    def __init__(
        self,
        grid,
        fluid_density=1000,
        gravity=9.8,
        slope="topographic__slope",
    ):
        """Calculate sediment transport on nodes (eventually)...

        Parameters
        ----------
        grid : RasterModelGrid
            A landlab grid.
        fluid_density : float, optional
            Density of fluid, default set to water density of 1000 kg / m^3
        gravity : float, optional
            Acceleration due to gravity (m/s^2).
        slope : str
            Field name of an at-node field that contains the slope.
        """

"""
Later on...
# Calculate boundary shear stress (the weird velocity one)
boundary_shear_stress = fluid_density * drag_coefficient ('vertical_velocity' **2 * 'cross_stream_velocity' ** 2)
for now...
# Depth-slope shear stress
depth_slope_shear_stress = fluid_density * gravity * 'water_surface__slope' * 'surface_water__depth'
"""
        super(multi_grain, self).__init__(grid)

        assert slope in grid.at_node

        self._E = self._grid.zeros(at="node")
        self._g = g
        self._rho = fluid_density
        self._n = n_sp
        self._slope = slope

    def run_one_step(self, dt):
        """calculate shear stress (for now)

        For one time step, this calculates depth-slope product shear stress at nodes

        Parameters
        ----------
        dt : float
            Time step.
        """

        depth = self._grid.at_node["surface_water__depth"]
        S = self._grid.at_node[self._slope]

        self._tau = self._rho * self._g * depth * S

