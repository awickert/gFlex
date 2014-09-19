import numpy as np
import scipy
from scipy.io import loadmat
from scipy.sparse import linalg
import os
import sys

class iceflow(object):

  def __init__(self, boundary_condition='Dirichlet0', sliding=True, useGRASS=False, location=None, mapset='PERMANENT', gisbase='/usr/local/src/grass_trunk/dist.x86_64-unknown-linux-gnu', grassdata_dir='grassdata', C0=0.0006, T=273.15, bcap=1000.):
    """
    Instantiates the iceflow class. If useGRASS=True, other parameters are
    needed, and these connect to the GRASS GIS location from which the climate 
    and mass balance parameters are derived.
    
    BOUNDARY CONDITIONS
    The domain boundary should not intersect
    any major ice masses. These default boundary condition options are suited
    for dealing with small peripheral ice masses adjacent to the domain
    boundary.
    
    Dirichlet0 - prescribed head: boundary cells are prescribed as zero ice 
    thickness. "Eliminated" flux is recorded as a dynamic leackage term.
    
    Neumann0 - prescribed flux: boundary cells have zero flux leaving domain.
    Ice "piles up" where it flows against the domain boundary.
    
    BASAL SLIDING COEFFICIENT [m/a/Pa]
    Proportional to driving stress. Marshall et al. (2005) use 0.0006 m/a/Pa 
    for Vatnajokull Ice Cap.
    Colgan used 0.0003 m/a/Pa for the Front Range in Colorado.
    Defaults to 0.0006 m/a/Pa
    
    ICE TEMPERATURE [K]
    Defaults to T=273.15 K: isothermal and constant through time and space.
    Right now, is constant, and cannot be used with arrays
    So flow law parameter defined here alone
    """
    self.Option_BC = boundary_condition
    self.sliding = sliding # Whether or not sliding will happen here
    self.useGRASS = useGRASS
    if self.useGRASS:
      if location:
        print "Using GRASS location", location
      else:
        sys.exit("Must define a GRASS GIS location.")
      # Implicitly 
      os.environ['GISBASE'] = gisbase
      gisdbase = os.path.join(os.environ['HOME'], grassdata_dir)
      sys.path.append(os.path.join(os.environ['GISBASE'], "etc", "python"))
      import grass.script as grass
      import grass.script.setup as gsetup
      gsetup.init(gisbase, gisdbase, location, mapset)
      from grass.script import array as garray
    self.grass = grass
    self.garray = garray
    # define constants:
    self.secyr = 3600*24*365.24 # conversion factor 1 [s/a]
    self.g = 9.81 # gravitational acceleration [m/s2]
    self.rho = 917. # glacier density [kg/m3]
    self.n = 3 # flow law exponent [unitless]  
    self.C0 = C0/self.secyr # convert units of basal sliding coefficient [m/s/Pa]
    # Including those for effective ice viscosity
    R = 8.314 # ideal gas constant [J/mol/K]
    QcCold = 60000. # cold/warm ice creep activation energy [J/mol]
    QcWarm = 139000.
    AoCold = 1.14E-5 # cold/warm ice reference flow law parameter [/Pa3/a]
    AoWarm = 5.47E10
    # Combine cold and warm ice flow law parameters [/Pa3/a]
    self.T = T # for output
    self.A = AoCold*np.exp(-QcCold/R/T) * (T < 263.5) + AoWarm*np.exp(-QcWarm/R/T) * (T >= 263.5)
    self.bcap = bcap/self.secyr

  def initialize_elevation_grid_from_file(self, topo, X, Y):
    """
    # INPUT: Bedrock/earth surface topography.
    # XI_500m = local easting coordinates [m]
    # YI_500m = local northing coordinates [m]
    # Z_500m = elevation [m]
    topoin = loadmat('topography_500m.mat')
    x = X = topoin['XI_500m']
    y = Y = topoin['YI_500m']
    Z = topoin['Z_500m']
    """
    pass

  def initialize_modern_climate_grid_from_file(self, precip, temp):
    """
    # INPUT: Climatological mean contemporary precipitation and air temperature 
    # fields, gridded to save input domain as the topography.
    # Ta_500m = mean annual air temprature [K]
    # Pa_500m = mean annual precipitation [m/a]
    climatein = loadmat('climate_500m.mat')
    Tair = climatein['Ta_500m']
    Pair = climatein['Pa_500m']
    """
    pass
    
  def initialize_modern_climate_and_elevation_from_GRASS(self, temperature_summer, precipitation_annual, elevation, resolution, n, s, w, e, ice=None):
    self.mass_balance_type='PDD'
    """
    Get climate and elevation grids from GRASS
    Need to separate these if I will use different data for mass balance
    """
    # First, use the region extent from the topography, cropped (by hand) to remove nulls at borders
    self.grass.run_command('g.region', n=n, s=s, w=w, e=e)
    # Then apply the desired resolution. This won't be *quite* a square grid, but will be close.
    self.grass.run_command('g.region', res=resolution)
    # Next, import arrays of topography, climate, and (if applicable) ice
    self.Zb = self.garray.array()
    self.Zb.read(elevation)
    self.Zb_initial = self.Zb.copy() # in case isostasy is used
    self.Tair = self.garray.array()
    self.Tair.read(temperature_summer)
    self.Pair = self.garray.array()
    self.Pair.read(precipitation_annual)
    self.Pair /= 1000. # PRISM is in mm/yr, change to m/yr.
    # And then create the grids of the x and y values that go along with these positions
    self.dx = self.grass.region()['ewres']
    self.dy = self.grass.region()['nsres']
    self.grass.mapcalc('x = x()', overwrite=True, quiet=True)
    self.grass.mapcalc('y = y()', overwrite=True, quiet=True)
    #X = np.arange(0, Z.shape[1], self.dx)
    #Y = np.arange(0, Z.shape[0], self.dy)
    self.x = self.garray.array()
    self.y = self.garray.array()
    self.x.read('x')
    self.y.read('y')
    self.ny, self.nx = self.x.shape # number of easting and northing nodes [unitless]  
    self.H0 = self.garray.array() # Starts out with 0's at the right shape
    if ice:
      self.H0.read(ice)
    
  def initialize_modern_elevation_and_grids_from_GRASS(self, elevation, resolution, n, s, w, e, ice=None):
    """
    Get elevation grids from GRASS as well as the grid dimensions
    """
    # First, use the region extent from the topography, cropped (by hand) to remove nulls at borders
    self.grass.run_command('g.region', n=n, s=s, w=w, e=e)
    # Then apply the desired resolution. This won't be *quite* a square grid, but will be close.
    self.grass.run_command('g.region', res=resolution)
    # Next, import arrays of topography
    self.Zb = self.garray.array()
    self.Zb.read(elevation)
    self.Zb_initial = self.Zb.copy() # in case isostasy is used
    # And then create the grids of the x and y values that go along with these positions
    self.dx = self.grass.region()['ewres']
    self.dy = self.grass.region()['nsres']
    self.grass.mapcalc('x = x()', overwrite=True, quiet=True)
    self.grass.mapcalc('y = y()', overwrite=True, quiet=True)
    #X = np.arange(0, Z.shape[1], self.dx)
    #Y = np.arange(0, Z.shape[0], self.dy)
    self.x = self.garray.array()
    self.y = self.garray.array()
    self.x.read('x')
    self.y.read('y')
    self.ny, self.nx = self.x.shape # number of easting and northing nodes [unitless]  
    self.H0 = self.garray.array() # Starts out with 0's at the right shape
    # And ice array, if applicable
    if ice:
      self.H0.read(ice)

  def initialize_modern_climate_from_GRASS(self, temperature_summer, precipitation_annual):
    self.mass_balance_type='PDD'
    """
    Get grids from GRASS
    Run after elevation grids obtained (this needs to be fixed)
    """
    self.Tair = self.garray.array()
    self.Tair.read(temperature_summer)
    self.Pair = self.garray.array()
    self.Pair.read(precipitation_annual)
    self.Pair /= 1000. # PRISM is in mm/yr, change to m/yr.

  def climate_corrections(self, T_correction, P_factor, melt_season_length_days):
    self.T_correction = T_correction
    self.P_factor = P_factor
    self.melt_season_length_days = melt_season_length_days
    self.Ta = self.Tair + T_correction # Temperature average through melt season [degC]
    self.Pa = self.Pair * P_factor / self.secyr # Total annual precipitation (assume snow?) (doesn't have to be this) [m/s]
    
  def basic_mass_balance_with_GRASS(self, ELA, dbdz):
    """
    dbdx = change in mass balance w/ change in x
    dbdy = ditto for y
    ...
    """
    self.mass_balance_type='ELA'
    self.ela0 = ELA
    self.dbdz = dbdz # [mWE/a/m]
    self.dbdz /= self.secyr
    #self.bcap = bcap # accumulation cap for now

  def initialize_mass_balance_basic_PDD(self, melt_factor, precip_lapse, temp_lapse=-4.5E-3):
    """"
    Surface mass balance - An elevation-dependent surface mass 
    balance is parameterized from observed contemporary solid precipitation 
    and air temperature fields that are perturbed to glacial climate.
    """
    # Removed delete step -- having small arrays in memory not an issue
    self.precip_lapse = precip_lapse/self.secyr # convert units of precipitation lapse rate [m/s/m]
    self.mu = melt_factor/86400. # convert units of melt factor [m/s/K]
    self.temp_lapse=temp_lapse

  def setup(self, run_length, dt_years=0.5, record_frequency=None):
    """
    self.run_length: total simulation length [a]
    self.dt_years:   time step [a]    
    """
    # simulation length:
    self.run_length = run_length
    self.dt_years = dt_years
    tEnd = self.run_length*self.secyr # model run length [s] 
    self.dt = self.dt_years*self.secyr # time step [s]
    self.t = np.arange(0, tEnd+self.dt/2., self.dt) # time vector [s]

    # variables to record series and slices:
    # time interval between time steps selected for recording [a]
    if record_frequency:
      self.record_frequency = record_frequency
    else:
      self.record_frequency = np.min((self.run_length, 10.)) # time interval between time steps selected for recording [a]
      record_timesteps = np.arange(0, self.run_length+self.record_frequency/2., self.record_frequency) # time-steps to be record [a]
    self.record_index = 0 # index of recorded selected time steps [unitless]
    self.b_timestep = 0 # surface mass balance record counter [a]
    self.a_timestep = 0 # surface ablation record counter [a]
    self.c_timestep = 0 # surface accumulation record counter [a]
    self.outflow_timestep = 0 # domain outflow record counter [a]
    self.time_series = np.zeros((len(record_timesteps),5)) # record of selected variables at selected time step [variable]  
    self.H_record = [] # List to record ice thicknes at time slices [m]
    self.Zb_record = [] # List to record bed elevation at time slices [m] -- important for isostasy, maybe in future if erosion is included
    self.dz_record = [] # change in bed elev (isostasy, etc.)
    self.uS_record = [] # List to record sliding velocities at time slices [m/a]
    self.b_record = [] # List to record mass balance at time slices [m/a]
    self.uD_record = [] # List to record depth-averaged deformational velocities [m/a]
    self.t_record = [] # time step [yr]



  def update(self):
    self.build_sparse_array()
    self.solve_sparse_equation()

  def run(self, run_length, plot=False):
    self.setup(run_length=run_length)
    self.enforce_boundary_conditions()
    self.update()
    if plot:
      self.plot()

  """
  def initialize(self):
    pass

  def update(self):
    pass
  
  def finalize(self, plot=False, save=True):
    pass
  
  def run(self):
    self.initialize()
    self.run()
    self.finalize()
    
  def set_mass_balance_grid(self, bgrid):
    self.bgrid = bgrid

  def get_ice_elevation(self):
    return self.Z...
  
  def get_topo(self):
    return self.Z...
    
  def get_total_elevation(self):
    return self.Z... + ...
    
  def update_basal_elevation
  """    

  def enforce_boundary_conditions(self):
    # boundary conditions:  
    if self.Option_BC == 'Dirichlet0':
      # prescribed zero ice thickness
      # boundary condition mask [binary]
      self.BC = np.ones((self.ny,self.nx))
      self.BC[:1,:] = 0 
      self.BC[:,:1] = 0 
      self.BC[-1:,:] = 0 
      self.BC[:,-1:] = 0
      self.Zb_initial *= self.BC
    elif self.Option_BC == 'Neumann0':
      # type 2 (prescribed zero ice flux) 
      self.BC = np.ones((self.ny,self.nx)) # boundary condition mask [binary]
    
  def setup_sparse_array(self):

    self.R_term_yes = np.hstack((np.ones((self.ny,self.nx-1)), np.zeros((self.ny,1)))) # identify nodes where the RIGHT-hand matrix term is present [binary]
    self.L_term_yes = np.hstack((np.zeros((self.ny,1)), np.ones((self.ny,self.nx-1)))) # identify nodes where the LEFT-hand matrix term is present [binary]
    self.D_term_yes = np.vstack((np.ones((self.ny-1,self.nx)), np.zeros((1,self.nx)))) # identify nodes where the DOWN-ward matrix term is present [binary]  
    self.U_term_yes = np.vstack((np.zeros((1,self.nx)), np.ones((self.ny-1,self.nx)))) # identify nodes where the UP-ward matrix term is present [binary]  
    self.C_term_yes = np.ones((self.ny,self.nx)) # identify nodes where the CENTRE matrix term is present [binary]
       
    # bedrock/land surface field:
    #self.Zb = interp2(XI_500m,YI_500m,Z_500m,x,y);
    #  # resample bedrock elevation at nodes j,i [m] -- not doing here b/c should happen before this code -- Landlab integration
    self.Zb *= self.BC # constrain bedrock elevation [m]
    #
    # MUST MOVE if changing elevation (e.g., flexure, erosion, asteroid impacts, global thermonuclear war...)
    self.ZbP1 = np.hstack((self.Zb[:,1:], np.zeros((self.ny,1)))) # bedrock elevation at node j,i+1 [m]
    self.ZbM1 = np.hstack((np.zeros((self.ny,1)), self.Zb[:,:-1])) # bedrock elevation at node j,i-1 [m]
    self.Zb_jP1i = np.vstack((self.Zb[1:,:], np.zeros((1,self.nx)))) # bedrock elevation at node j+1,i [m]
    self.Zb_jM1i = np.vstack((np.zeros((1,self.nx)), self.Zb[:-1,:])) # bedrock elevation at node j-1,i [m]  
      
    self.A /= self.secyr # convert units of effective ice viscosity [/Pa3/s]  

    # initialize evolving ice geometry variables:
    self.H = self.H0.copy() # ice thickness at node j,i [m] 
    self.Zs = self.Zb + self.H # ice surface elevation at node j,i [m]
    
  def build_sparse_array(self, tt):
  
    # calculate ice form and flow variables:
    try:
      # If isostasy is involved, include projected future isostatic change in the ice flow calculations
      self.Zs = self.Zb + self.H# + self.dz # ice surface elevation at node j,i: projected to future w/ isostasy [m]
    except:
      self.Zs = self.Zb + self.H # ice surface elevation at node j,i [m]
    if self.mass_balance_type == 'PDD':
      #a = (self.Ta + (self.H*self.temp_lapse))*self.mu * self.melt_season_length_days/365.24 # surface ablation scaled to melt season [m/s]
      #c = self.Pa + (self.H*self.precip_lapse)/self.secyr # surface accumulation [m/s]
      self.dz_plus_H = self.H + self.Zb - self.Zb_initial
      a = (self.Ta + (self.dz_plus_H*self.temp_lapse))*self.mu * self.melt_season_length_days/365.24 # surface ablation scaled to melt season [m/s]
      c = self.Pa + (self.dz_plus_H*self.precip_lapse)/self.secyr # surface accumulation [m/s]
      self.b = (c - a) # surface mass balance [m/s]
    elif self.mass_balance_type == 'ELA':
      self.b = self.dbdz*(self.Zs-self.ela0)
    self.b[self.b > self.bcap] = self.bcap
      
    self.HPh = np.hstack(( (self.H[:,1:]+self.H[:,:-1])/2., np.zeros((self.ny,1)) )) # ice thickness at node j,i+1(?) [m]
    self.HMh = np.hstack(( np.zeros((self.ny,1)), (self.H[:,1:] + self.H[:,:-1])/2. )) # ice thickness at node j,i-1(?) [m]
    self.H_jPhi = np.vstack(( (self.H[1:,:]+self.H[:-1,:])/2., np.zeros((1,self.nx)) )) # ice thickness at node j+1(?),i [m]
    self.H_jMhi = np.vstack(( np.zeros((1,self.nx)), (self.H[1:,:]+self.H[:-1,:])/2. )) # ice thickness at node j-1(?),i [m]      
      
    self.alpha = (np.hstack(( np.zeros((self.ny,1)), (self.Zs[:,2:] - self.Zs[:,:-2]) / (self.x[:,2:] - self.x[:,:-2]), np.zeros((self.ny,1)) ))**2 + np.vstack(( np.zeros((1,self.nx)), (self.Zs[2:,:] - self.Zs[:-2,:])/(self.y[2:,:] - self.y[:-2,:]), np.zeros((1,self.nx)) ))**2)**0.5 # absolute ice surface slope at node i,j [m/m]
    self.alphaPh = np.hstack(( (self.alpha[:,1:]+self.alpha[:,:-1])/2., np.zeros((self.ny,1)) )) # absolute ice surface slope at node j,i+1(?) [m]
    self.alphaMh = np.hstack(( np.zeros((self.ny,1)), (self.alpha[:,1:]+self.alpha[:,:-1])/2. )) # absolute ice surface slope at node j,i-1(?) [m]
    self.alpha_jPhi = np.vstack(( (self.alpha[1:,:]+self.alpha[:-1,:])/2., np.zeros((1,self.nx)) )) # absolute ice surface slope at node j+1(?),i [m]
    self.alpha_jMhi = np.vstack(( np.zeros((1,self.nx)), (self.alpha[1:,:]+self.alpha[:-1,:])/2. )) # absolute ice surface slope at node j-1(?),i [m]      
        
    self.dZsdxPh = np.hstack(( (self.Zs[:,1:]-self.Zs[:,:-1]) / (self.x[:,1:]-self.x[:,:-1]), np.zeros((self.ny,1)) )) # directional ice surface slope at node j,i+1(?) [m/m]
    self.dZsdxMh = np.hstack(( np.zeros((self.ny,1)), (self.Zs[:,1:]-self.Zs[:,:-1]) / (self.x[:,1:]-self.x[:,:-1]) )) # directional ice surface slope at node j,i-1(?) [m/m]
    self.dZsdy_jPhi = np.vstack(( (self.Zs[1:,:]-self.Zs[:-1,:]) / (self.y[1:,:]-self.y[:self.ny-1,:]), np.zeros((1,self.nx)) )) # directional ice surface slope at node j+1(?),i [m/m]
    self.dZsdy_jMhi = np.vstack(( np.zeros((1,self.nx)), (self.Zs[1:,:]-self.Zs[:-1,:]) / (self.y[1:,:]-self.y[:-1,:]) )) # directional ice surface slope at node j-1(?),i [m/m]
      
    tauPh = -self.rho*self.g*self.HPh*self.dZsdxPh # driving stress at node j,i+1(?) [Pa] - positive x direction
    tauMh = -self.rho*self.g*self.HMh*self.dZsdxMh # driving stress at node j,i-1(?) [Pa] - negative x direction
    tau_jPhi = -self.rho*self.g*self.H_jPhi*self.dZsdy_jPhi # driving stress at node j+1(?),i [Pa] - positive y direction
    tau_jMhi = -self.rho*self.g*self.H_jMhi*self.dZsdy_jMhi # driving stress at node j-1(?),i [Pa] - negative y direction
    
    qPh = 2*self.A/(self.n+2)*(self.rho*self.g*self.alphaPh)**(self.n-1)*self.HPh**(self.n+1)*tauPh # ice discharge at node j,i+1(?) [m2/s]
    qMh = 2*self.A/(self.n+2)*(self.rho*self.g*self.alphaMh)**(self.n-1)*self.HMh**(self.n+1)*tauMh # ice discharge at node j,i-1(?) [m2/s]
    q_jPhi = 2*self.A/(self.n+2)*(self.rho*self.g*self.alpha_jPhi)**(self.n-1)*self.H_jPhi**(self.n+1)*tau_jPhi # ice discharge at node j+1(?),i [m2/s]
    q_jMhi = 2*self.A/(self.n+2)*(self.rho*self.g*self.alpha_jMhi)**(self.n-1)*self.H_jMhi**(self.n+1)*tau_jMhi # ice discharge at node j-1(?),i [m2/s]
    
    self.uD = ( ((qPh/np.maximum(self.HPh,1E-8) + qMh/np.maximum(self.HMh,1E-8)) / 2.)**2 + ((q_jPhi/np.maximum(self.H_jPhi,1E-8) + q_jMhi/np.maximum(self.H_jMhi,1E-8)) / 2.)**2)**0.5 # absolute depth-averaged deformational velocity at node j,i [m/s]

    if self.sliding:
      uSPh = -self.C0*tauPh # basal sliding velocity at node j,i+1(?) [m/s]
      uSMh = -self.C0*tauMh; 
        # basal sliding velocity at node j,i-1(?) [m/s]
      uS_jPhi = self.C0*tau_jPhi # basal sliding velocity at node j+1(?),i [m/s] # (AW) made these positive, solution looks better now
      uS_jMhi = self.C0*tau_jMhi # basal sliding velocity at node j-1(?),i [m/s] # (AW) but haven't checked why
      self.uS = (((uSPh + uSMh)/2)**2 + ((uS_jPhi + uS_jMhi)/2)**2)**0.5 # absolute basal sliding velocity at node j,i [m/s]       
    else:
      uSPh = uSMh = uS_jPhi = uS_jMhi = self.uS = 0
      
    self.dZsdxPh[self.dZsdxPh == 0] = 1E-6 # directional ice surface slope at node j,i+1(?) [m/m]
    self.dZsdxMh[self.dZsdxMh == 0] = 1E-6 # directional ice surface slope at node j,i-1(?) [m/m]
    self.dZsdy_jPhi[self.dZsdy_jPhi == 0] = 1E-6 # directional ice surface slope at node j+1(?),i [m/m]
    self.dZsdy_jMhi[self.dZsdy_jMhi == 0] = 1E-6 # directional ice surface slope at node j-1(?),i [m/m]
      
    DPh = qPh/self.dZsdxPh; DPh[:,-1] = 0 # diffusion term at node j,i+1(?) [m2/s]
    DMh = qMh/self.dZsdxMh; DMh[:,0] = 0 # diffusion term at node j,i-1(?) [m2/s]     
    D_jPhi = q_jPhi/self.dZsdy_jPhi; D_jPhi[-1,:] = 0 # diffusion term at node j+1(?),i [m2/s]
    D_jMhi = q_jMhi/self.dZsdy_jMhi; D_jMhi[0,:] = 0 # diffusion term at node j-1(?),i [m2/s]  
    
    # create and solve 5-banded matrix:  
    Array_L = (+ DMh*self.dt/self.dx**2 + uSMh*self.dt/2/self.dx) * self.L_term_yes # unknown at node j,i-1 [unitless]
    Array_D = (+ D_jMhi*self.dt/self.dy**2 + uS_jMhi*self.dt/2/self.dy) * self.D_term_yes # unknown at node j-1,i [unitless]
    Array_C = 1 \
      + (- DPh*self.dt/self.dx**2 - uSPh*self.dt/2./self.dx) * self.R_term_yes \
      + (- DMh*self.dt/self.dx**2 + uSMh*self.dt/2./self.dx) * self.L_term_yes \
      + (- D_jPhi*self.dt/self.dy**2 - uS_jPhi*self.dt/2./self.dy) * self.U_term_yes \
      + (- D_jMhi*self.dt/self.dy**2 + uS_jMhi*self.dt/2./self.dy) * self.D_term_yes
      # unknown at node j,i [unitless]
    Array_U = (+ D_jPhi*self.dt/self.dy**2 - uS_jPhi*self.dt/2/self.dy) * self.U_term_yes # unknown at node j+1,i [unitless]
    Array_R = (+ DPh*self.dt/self.dx**2 - uSPh*self.dt/2/self.dx) * self.R_term_yes # unknown at node j,i+1 [unitless]
    Array_rhs = (+ self.b*self.dt + self.H) * self.C_term_yes \
      + DPh*(self.Zb - self.ZbP1)*self.dt/self.dx**2 * self.R_term_yes \
      - DMh*(self.ZbM1 - self.Zb)*self.dt/self.dx**2 * self.L_term_yes \
      + D_jPhi*(self.Zb - self.Zb_jP1i)*self.dt/self.dy**2 * self.U_term_yes \
      - D_jMhi*(self.Zb_jM1i - self.Zb)*self.dt/self.dy**2 * self.D_term_yes # known (right hand side vector) at node j,i [m] 
      
    Vec_D = np.reshape(Array_D, -1, order='F') # reshape unknown at node j-1,i from array to vector [unitless]     
    Vec_L = np.reshape(Array_L, -1, order='F') # reshape unknown at node j,i-1 from array to vector [unitless]   
    Vec_C = np.reshape(Array_C, -1, order='F') # reshape unknown at node j,i from array to vector [unitless]   
    Vec_R = np.reshape(Array_R, -1, order='F') # reshape unknown at node j,i+1 from array to vector [unitless]   
    Vec_U = np.reshape(Array_U, -1, order='F') # reshape unknown at node j+1,1 from array to vector [unitless]   

    # NECESSARY if updating these parameters
    self.ZbP1 = np.hstack((self.Zb[:,1:], np.zeros((self.ny,1)))) # bedrock elevation at node j,i+1 [m]
    self.ZbM1 = np.hstack((np.zeros((self.ny,1)), self.Zb[:,:-1])) # bedrock elevation at node j,i-1 [m]
    self.Zb_jP1i = np.vstack((self.Zb[1:,:], np.zeros((1,self.nx)))) # bedrock elevation at node j+1,i [m]
    self.Zb_jM1i = np.vstack((np.zeros((1,self.nx)), self.Zb[:-1,:])) # bedrock elevation at node j-1,i [m]  

    self.Vec_rhs = np.reshape(Array_rhs, -1, order='F') # reshape known at node j,i from array to vector [unitless]   
    
    # Sparse matrix
    diags = np.vstack(( np.hstack(( Vec_L[self.ny:], np.zeros(self.ny) )), \
                        np.hstack(( Vec_D[1:], 0 )), \
                        Vec_C, \
                        np.hstack(( 0, Vec_U[:-1] )), \
                        np.hstack(( np.zeros(self.ny), Vec_R[:-self.ny] )) ))
    self.Matrix_lhs = scipy.sparse.spdiags(diags, [-self.ny,-1,0,1,self.ny], self.ny*self.nx, self.ny*self.nx, format='csr') # create five-banded sparse matrix [unitless]

  def solve_sparse_equation(self):
    # solve ice thickness, in Matlab using backslash, here using umfpack [m]
    Vec_H = linalg.spsolve(self.Matrix_lhs, self.Vec_rhs, use_umfpack=False)

    self.H = np.reshape(Vec_H, (self.nx, self.ny)).transpose() # reshape ice thickness from vector to array [m]
    self.H[self.H < 1E-8] = 0 # constrain potentially negative ice thicknesses [m]
    
    # potentially implement domain boundary condition:
    H_pre = self.H # note ice thickness before implementing boundary condition [m]
    self.H = self.H*self.BC # implement prescribed boundary condition on ice thickness [m]
    H_post = self.H # note ice thickness again after implementing boundary condition [m]
    
  def update_recorded_mass_balance(self):
    # update mass balance terms between recorded time steps:
    self.outflow_timestep = self.outflow_timestep + sum(sum(H_pre - H_post))* \
      self.dx*self.dy*self.rho # update domain outflow (kg/a)  
    self.b_timestep = self.b_timestep + np.sum(self.b*(self.H>0))*self.dx*self.dy*self.dt*self.rho # update time step surface mass balance (kg/a)
    self.a_timestep = self.a_timestep + np.sum(a*(self.H>0))*self.dx*self.dy*self.dt*self.rho # update time step surface ablation (kg/a)
    self.c_timestep = self.c_timestep + np.sum(c*(self.H>0))*self.dx*self.dy*self.dt*self.rho # update time step surface accumulation (kg/a) 
      
  def record_model_parameters(self, tt, record_frequency):
    self.record_timesteps = np.arange(0, self.run_length+record_frequency/2., record_frequency) # time-steps to be record [a]
    # at selected time steps, record various model parameters:
    if self.t[tt]/self.secyr <= (self.record_timesteps[self.record_index] + self.dt/self.secyr/2.) and self.t[tt]/self.secyr >= (self.record_timesteps[self.record_index] - self.dt/self.secyr/2.):
      print 'model year:', '%10.1f' %(self.t[tt]/self.secyr)
        # display current time step in command window [a]
      self.H_record.append(self.H)
        # record time step ice thickness field [m]
      self.Zb_record.append(self.Zb.copy())
        # record bed elevation [m/a] -- important if isostasy is implemented.
      self.dz_record.append(self.Zb - self.Zb_initial)
        # record bed elevation change from isostasy [m]
      self.uS_record.append(self.uS*self.secyr)
        # record time step basal sliding velocity field [m/a]
      self.uD_record.append(self.uD*self.secyr)
        # record time step deformational velocity field [m/a]
      self.b_record.append(self.b*self.secyr)
        # record mass balance [m/a]
      self.t_record.append(self.t[tt]/self.secyr)
      #self.time_series[self.record_index,:] = np.hstack(( self.record_timesteps[self.record_index], \
      #                                          self.c_timestep/self.record_frequency, \
      #                                          self.a_timestep/self.record_frequency, \
      #                                          self.b_timestep/self.record_frequency, \
      #                                          self.outflow_timestep/self.record_frequency ))
        # update time series of mass balance elements [kg/a]
      self.record_index += 1;
        # update record index [unitless]
      self.b_timestep = 0;
        # reset surface mass balance counter [a]
      self.a_timestep = 0;
        # reset surface ablation counter [a]
      self.c_timestep = 0;
        # reset surface accumulation counter [a]
      self.outflow_timestep = 0;
        # reset domain outflow counter [a]

  def save_output(self, filename='Simulation_Output'):
    out_array = (self.H_record,self.Zb_record,self.uD_record,self.uS_record,self.b_record,self.t_record,
      self.time_series,self.T,self.A,self.C0,self.x,self.y,self.Zb,self.BC,self.T_correction,\
      self.P_factor,self.mu,self.dx,self.dy)
    np.save(filename, out_array) # save simulation output into a single .npy file that can be called on
      # for graphical output at a later time
    

  def plot(self, save=True):
    """
    save == True: save figure images
    save == False: draw the plots on screen
    """
    from matplotlib import pyplot as plt

    plt.figure()
    plt.imshow(self.H, interpolation='nearest')
    plt.colorbar()
    plt.title('Ice Thickness', fontsize=16)
    if save:
      plt.savefig('IceThickness.png')
      plt.close()
    plt.show()

    """
    plt.figure()
    plt.imshow(self.Zb, interpolation='nearest')
    plt.colorbar()
    plt.title('Elevation', fontsize=16)
    if save:
      plt.savefig('Topography.png')
      plt.close()
    else:
      plt.show()
    """
