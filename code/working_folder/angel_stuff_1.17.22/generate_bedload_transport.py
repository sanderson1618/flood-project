"""Calculates bedload transport using Parker 1990 equation and flow properties from OverlandFlow.

.. codeauthor:: Angel Monsalve

Examples
--------
>>> import numpy as np
>>> from landlab import RasterModelGrid
>>> from landlab.components import DepthSlopeProductErosion

Create an example here
"""
from landlab import Component
import numpy as np
import scipy.constants
from scipy.interpolate import interp1d

class BedloadTransport(Component):

    """Landlab component that predicts bedload transport using Parker 1990 surface based equation.
    This component predicts the bedload transport rate and fractional tranport at each link using the unsteady
    shear stress and Yager 2012 shear stress partitioning technique"""

    _name = "BedloadTransport"

    _unit_agnostic = True

    _info = {
        "bed_grain_size_distribution__location": {
        "dtype": float,
        "intent": "in",
        "optional": False,
        "units": "-",
        "mapping": "node",
        "doc": "Accounts for a grain size dsitribution that varies in space. Each value corresponds to a given GSD.",
        },
    }

    def __init__(
        self,
        grid,
        g=scipy.constants.g,    # Acceleration due to gravity (m/s^2).
        rho = 1000,             # Fluid density
        rho_s = 2650,           # Fluid density
        tauStar_rsgo = 0.0386,  # Reference dimensionless shear stress for the median size
        beta = 0.0951,          # Coefficient for the hiding function
        dt = 1e-3,              # time step
    ):
        
        """Calculate bedload transport and fractional rates on links using the discharge and 
        water depth predicted in overland flow. Then, shear stress equation as input, 
        but further partitioned using Yager 2012 shear stress partitioning technique.

        Landlab component that add more info here

        Parameters
        ----------
        grid :  RasterModelGrid
                A landlab grid. 
        g :     float, optional
                Acceleration due to gravity (m/s^2)
        rho :   float, optional
                Density of fluid, default set to water density of 1000 kg / m^3
        rho_s : float, optional
                Density of sediment, default set to sediment density of 2650 kg / m^3
        tauStar_rsgo : float, optional
                Reference dimensionless shear stress for the median size - default 0.0386 
        beta :  float, optional
                Coefficient for the hiding function - default 0.0951            
        dt :    float, set it equal to overland flow dt
                Time step for velocity derivatives in time - default 1e-3 
        """
        super().__init__(grid)

        self._g = g
        self._rho = rho
        self._rho_s = rho_s
        self._R = (rho_s-rho)/rho
        self._tauStar_rsgo = tauStar_rsgo
        self._beta = beta
        self._dt = dt
        
        self._z = self._grid.at_node["topographic__elevation"]
        self._h = self._grid.at_node["surface_water__depth"]
        self._h_links = self._grid.at_link["surface_water__depth"]
        self._q = self._grid.at_link["surface_water__discharge"]
        self._u = self._grid.zeros(at="link")
                
        self._dz_ds = self._grid.zeros(at='link')
        self._dh_ds = self._grid.zeros(at='link')
        self._du_ds = self._grid.zeros(at='link')
        self._du_dt = self._grid.zeros(at='link')
        
        self._Sf = self._grid.zeros(at="link")
        self._Rh = self._grid.zeros(at="link")
        self._tau = self._grid.zeros(at="link")
        
        self._gsd_location = self._grid.at_node["bed_grain_size_distribution__location"]       
        self._gsdInit_flag = True
        
        self._omega0 = self._grid.zeros(at='link')
        self._sigma0 = self._grid.zeros(at='link')
        self._omega_m = self._grid.zeros(at='link')
        
    def gsdProperties(self, gsd):
              
        grainSize = gsd[:,0] 
        D = np.zeros(grainSize.shape[0]-1)  # One for each equivalent grain size
        f = np.zeros([grainSize.shape[0]-1,gsd.shape[1]-1]) # One for each equivalent grain size at each different GSD
        geoMean_nodes = np.zeros_like(self._gsd_location)
        geoStd_nodes = np.zeros_like(self._gsd_location)
        
        for i in np.arange(0,gsd.shape[1]-1):

            grainPerc = gsd[:,i+1]/100  # Now working in frequency  
            
            """ The equivalent grain sizes and their frequencies are calculated """
            if i == 0:
                D = (grainSize[0:-1] * grainSize[1:])**0.5
            f[:,i] = grainPerc[0:-1] - grainPerc[1:]
            
            """ The geometric mean and geometric standard deviation are calculated """
            geoMean0 = 2**np.max(np.cumsum(f[:,i]*np.log2(D)))
            geoStd0 = 2**np.sqrt(np.max(np.cumsum((((np.log2(D)-(np.max(np.cumsum(f[:,i]*np.log2(D)))))**2)*f[:,i]))))
            
            """ finds all locations where i is equal to the specific gsd within the watershed """
            (gsdIndex,) = np.where(self._gsd_location == i)
            geoMean_nodes[gsdIndex] = geoMean0
            geoStd_nodes[gsdIndex] = geoStd0
        
        # At this point geoMean, geoStd, and the location of the different GSD 
        # are at nodes, we need them at links. To assing GSD properties the definition is: 
        # horizontal links are mapped from tail nodes
        # vertical links are mapped from head nodes
        
        geoMean = np.zeros_like(self._u)
        geoStd = np.zeros_like(self._u)
        gsd_loc_links = np.zeros_like(self._u)
        
        geoMean[self._grid.horizontal_links] = self._grid.map_link_tail_node_to_link(geoMean_nodes)[self._grid.horizontal_links]
        geoMean[self._grid.vertical_links] = self._grid.map_link_head_node_to_link(geoMean_nodes)[self._grid.vertical_links]

        geoStd[self._grid.horizontal_links] = self._grid.map_link_tail_node_to_link(geoStd_nodes)[self._grid.horizontal_links]
        geoStd[self._grid.vertical_links] = self._grid.map_link_head_node_to_link(geoStd_nodes)[self._grid.vertical_links]

        gsd_loc_links[self._grid.horizontal_links] = self._grid.map_link_tail_node_to_link(self._gsd_location)[self._grid.horizontal_links]
        gsd_loc_links[self._grid.vertical_links] = self._grid.map_link_head_node_to_link(self._gsd_location)[self._grid.vertical_links]
                 
        self._D = D
        self._f = f
        self._geoMean = geoMean
        self._geoStd = geoStd
        self._gsd_loc_links = gsd_loc_links
        self._gsdInit_flag = False
    
    def shearStress(self,u0):
        
        """ The shear stress is calculated at links according to
        tau = rho * g * Rh * Sf
        Where Sf is the unsteady friction slope which is calculated as
        Sf = S0 - dh/ds - U/g du/ds - 1/g du/dt

        The term ds indicates a certain direction, X and Y in this case. 
        All derivatives are first or second order approximations in each direction
        
        Most of the parameters within the grid are read again because they may 
        change in time in overland flow
        """
    
        """ First term in the friction slope - gradient of bed elevation S0 = - dz/ds """
        self._dz_ds = - self._grid.calc_grad_at_link(self._z)
        
        """ Second term in the friction slope - gradient of water depth """
        self._h = self._grid["node"]["surface_water__depth"]
        self._dh_ds = self._grid.calc_grad_at_link(self._h)
        
        self._q = self._grid["link"]["surface_water__discharge"]
        self._h_links = self._grid.at_link["surface_water__depth"]
        self._u = self._q / self._h_links
        
        """ Third term in the friction slope - gradient of flow velocity """
        # Velocity gradients are calculated in each direction
        self._du_ds = np.zeros_like(self._dz_ds)
        
        # In the X direction only horizontal gradients are valid. Here we use the 
        # map_mean_of_horizontal_links_to_node method. However, this tool is valid only for the X direction
        u_nodes_horizontal = self._grid.map_mean_of_horizontal_links_to_node(self._u)    # 2nd order interpolation
        self._du_ds[self._grid.horizontal_links] = self._grid.calc_grad_at_link(u_nodes_horizontal)[self._grid.horizontal_links]
        
        # In the Y direction we have to do it manually. A 2nd order interpolation is used     
        u_nodes_vertical = self._grid.map_mean_of_vertical_links_to_node(self._u)    # 2nd order interpolation
        u_nodes_vertical = np.flip(np.flip((np.reshape(u_nodes_vertical,(self._grid._shape[0],self._grid._shape[1])))),axis=1)
        
        du_dsV = np.zeros((u_nodes_vertical.shape[0]-1,u_nodes_vertical.shape[1]))
        
        for i in np.arange(u_nodes_vertical.shape[0]-2,-1,-1):
            du_dsV[i,:] = (u_nodes_vertical[i,:] - u_nodes_vertical[i+1,:]) / self._grid.dy
        
        # Three operations are in the right side. This is to recover landLab indexing order 
        self._du_ds[self._grid.vertical_links] = np.flip(du_dsV.T,axis=1).flatten(order='F')        
        
        """ Fourth term in the friction slope - rate of change of velocity in time"""
        self._du_dt = (u0[1,:] - u0[0,:]) / self._dt
        
        """ Friction slope calculation at links including unsteady effects """
        self._Sf = self._dz_ds - self._dh_ds - (u0[0,:]/self._g) * self._du_ds - 1/self._g * self._du_dt

        """ Now calculates the hydraulics ratio for the final shear stress calculation 
        Rh = wetted area / wetted perimeter 
        Different grid sizes (dx~=dy) are possible"""
        A = np.zeros_like(self._h_links)
        A[self._grid.horizontal_links] = self._h_links[self._grid.horizontal_links] * self._grid.dx
        A[self._grid.vertical_links] = self._h_links[self._grid.vertical_links] * self._grid.dy
        
        P = np.zeros_like(self._h_links)
        P[self._grid.horizontal_links] = self._grid.dx + 2 * self._h_links[self._grid.horizontal_links]
        P[self._grid.vertical_links] =  self._grid.dy + 2 * self._h_links[self._grid.vertical_links]     
        
        self._Rh = A / P
        
        """ Shear stress at links including unsteady effects 
        Equation is tau = Rho * g * Rh * Sf """
        self._tau = self._rho * self._g * self._Rh * self._Sf 

    def strainFunctions(self):
        po = np.array([0.6684,0.7639,0.8601,0.9096,0.9615,1,1.055,1.108,1.197,\
                       1.302,1.407,1.529,1.641,1.702,1.832,1.937,2.044,2.261,\
                       2.499,2.732,2.993,3.477,4.075,4.469,5.016,6.158,7.821,\
                       10.06,14.38,19.97,25.79,38.57,68.74,91.95,231.2,2320])
        oo = np.array([1.011,1.011,1.01,1.008,1.004,0.9997,0.9903,0.9789,0.9567,\
                       0.9273,0.8964,0.8604,0.8287,0.8123,0.7796,0.7554,0.7326,\
                       0.6928,0.6585,0.6345,0.615,0.5877,0.564,0.5523,0.5395,\
                       0.5209,0.5045,0.4917,0.479,0.4712,0.4668,0.462,0.4578,\
                       0.4564,0.4541,0.4527])
        so = np.array([0.8157,0.8157,0.8182,0.8233,0.8333,0.8439,0.8621,0.8825,\
                       0.9214,0.9723,1.025,1.083,1.13,1.153,1.196,1.225,1.25,1.287,\
                       1.313,1.333,1.352,1.38,1.403,1.414,1.426,1.444,1.458,1.469,\
                    1.48,1.486,1.49,1.493,1.497,1.498,1.499,1.5])
        
        omega0 = np.zeros_like(self._phi_sgom)
        sigma0 = np.zeros_like(self._phi_sgom)
            
        # There are three intervals where we can interpolate Omega and Sigma
        (I,) = np.where(self._phi_sgom <= 0.7639)
        if I.shape[0]>0:
            omega0[I] = oo[0]
            sigma0[I] = so[0]
        
        (I,) = np.where(self._phi_sgom > 231.2)
        if I.shape[0]>0:
            omega0[I] = np.interp(self._phi_sgom, po, oo)[I]
            sigma0[I] = np.interp(self._phi_sgom, po, oo)[I]
        
        (I,) = np.where( (self._phi_sgom > 0.7639) & (self._phi_sgom < 231.2))
        if I.shape[0]>0:        
            foo = interp1d(po, oo, kind='cubic')
            fso = interp1d(po, so, kind='cubic')
            omega0[I] = foo(self._phi_sgom[I])
            sigma0[I] = fso(self._phi_sgom[I])
        
        self._omega0 = omega0
        self._sigma0 = sigma0
        
    def bedloadParker1990(self):
        
        """ All predictions will be based on the absolute shear stress. Therefore,
        direction should be assigned later based on the velocity vector """
        
        tauT = np.abs(self._tau)
        self._R = (self._rho_s-self._rho)/self._rho
        
        tauStar_sgm = tauT/(self._rho*self._R*self._g*(self._geoMean/1000)) # geoMean is now in meters
        self._phi_sgom = tauStar_sgm/self._tauStar_rsgo
        
        self.strainFunctions()
        self._omega_m = 1+(np.log2(self._geoStd)/self._sigma0)*(self._omega0-1)
        
        # We need to obtain a phi_mi value for every equivalent grain size
        # Therefore, it is not possible to do it in a single step. We divided the
        # calculation and then everything is merged again
        
        phi_mi = np.zeros((self._omega_m.shape[0],self._D.shape[0]))
        G = np.zeros((self._omega_m.shape[0],self._D.shape[0]))
        Gf = np.zeros((self._omega_m.shape[0],self._D.shape[0]))
        p = np.zeros((self._omega_m.shape[0],self._D.shape[0]))
        
        for i in np.arange(0,self._D.shape[0]):
            phi_mi[:,i] = self._omega_m*self._phi_sgom*(self._D[i]/(self._geoMean))**-self._beta
            
            # There are three intervals where G is evaluated
            (I,) = np.where(phi_mi[:,i] > 1.59)
            if I.shape[0]>0:
                G[:,i] = 5474 * (1 - 0.853/phi_mi[:,i]) ** 4.5
            
            (I,) = np.where( (phi_mi[:,i] >= 1) & (phi_mi[:,i]<=1.59))
            if I.shape[0]>0:
                G[:,i] = np.exp(14.2*(phi_mi[:,i]-1)-9.28*(phi_mi[:,i]-1) ** 2)
            
            (I,) = np.where(phi_mi[:,i] < 1)
            if I.shape[0]>0:        
                G[:,i] = phi_mi[:,i] ** 14.2
            
        # Calculates Gf at each link - Links may have different f so it needs to identify each location        
        for i in np.arange(0,self._f.shape[1]):

            (I,) = np.where(self._gsd_loc_links == i)
            Gf[I,:] = self._f[:,i] * G[I,:]

        GfSum = np.sum(Gf,axis=1)
        Wstar_msi = 0.00218 * GfSum
        
        # Total bedload transport rate is calculated at each link
        self._bedloadRate = ((np.sqrt(tauT/self._rho))**3*Wstar_msi)/(self._R*self._g)
        
        # Frational bedload transport rate is calculated as a fraction of the
        # total bedload transport rate at each link
        
        for i in np.arange(0,self._D.shape[0]):
            (I,) = np.where(GfSum > 0)
            if I.shape[0]>0:
                p[I,i] = Gf[I,i] / GfSum[I]

        self._p =  p   
        
    def run_one_step(self, u0, gsd):
        
        """        
        Parameters
        ----------        
        u0 : float
            Flow velocity at links for the current and previous time steps
        gsd : array, float
            An array containing the grain size distribution for each different 
            GSD present at the watershed
        
        """
        # Unsteady shear stress is calculated every time step
        self.shearStress(u0)
        
        # Checks if the grain size distribution have been initialized
        # The very first time it should be empty, but after the first dt it does not change anymore
        if self._gsdInit_flag is True:
            self.gsdProperties(gsd)
        
        # Calculates bedload transport using Parker 1990 surface equation
        self.bedloadParker1990()        