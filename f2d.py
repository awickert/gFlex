from __future__ import division # No automatic floor division
from base import *
import scipy
from scipy.special import kei

# class F2D inherits Flexure and overrides __init__ therefore setting up the same
# three parameters as class Isostasy; and it then sets up more parameters specific
# to its own type of simulation.
class F2D(Flexure):
  def initialize(self, filename):
    super(F2D, self).initialize(filename)
    if self.Verbose: print 'F2D initialized'

  def run(self):
    if self.method == 'FD':
      # Finite difference
      super(F2D, self).FD()
      self.method_func = self.FD
    elif self.method == 'FFT':
      # Fast Fourier transform
      super(F2D, self).FFT()
      self.method_func = self.FFT
    elif self.method == "SPA":
      # Superposition of analytical solutions
      super(F2D, self).SPA()
      self.method_func = self.SPA   
    elif self.method == "SPA_NG":
      # Superposition of analytical solutions,
      # nonuniform points
      super(F2D, self).SPA_NG()
      self.method_func = self.SPA_NG
    else:
      sys.exit('Error: method must be "FD", "FFT", "SPA", or "SPA_NG"')

    if self.Verbose: print 'F2D run'
    self.method_func()
    #self.imshow(self.w) # debugging

  def finalize(self):
    if self.Verbose: print 'F2D finalized'
    super(F2D, self).finalize()
    
  ########################################
  ## FUNCTIONS FOR EACH SOLUTION METHOD ##
  ########################################

  def FD(self):
    self.elasprep()
    # Only generate coefficient matrix if it is not already provided
    try:
      self.coeff_matrix
    except:
      self.coeff_matrix_creator()
    self.fd_solve()

  def FFT(self):
    sys.exit("The fast Fourier transform solution method is not yet implemented.")

  def SPA(self):
    self.spatialDomainVars()
    self.spatialDomainGridded()

  def SPA_NG(self):
    self.spatialDomainVars()
    self.spatialDomainNoGrid()

  
  ######################################
  ## FUNCTIONS TO SOLVE THE EQUATIONS ##
  ######################################
  
   
  ## SPATIAL DOMAIN SUPERPOSITION OF ANALYTICAL SOLUTIONS
  #########################################################

  # SETUP

  def spatialDomainVars(self):
    self.D = self.E*self.Te**3/(12*(1-self.nu**2)) # Flexural rigidity
    self.alpha = (self.D/(self.drho*self.g))**.25 # 2D flexural parameter
    self.coeff = self.alpha**2/(2*np.pi*self.D)

  # GRIDDED

  def spatialDomainGridded(self):
  
    self.nx = self.q0.shape[1]
    self.x = np.arange(0,self.dx*self.nx,self.dx)
    
    self.ny = self.q0.shape[0]
    self.y = np.arange(0,self.dx*self.nx,self.dx)
    
    # Prepare a large grid of solutions beforehand, so we don't have to
    # keep calculating kei (time-intensive!)
    # This pre-prepared solution will be for a unit load
    bigshape = 2*self.ny+1,2*self.nx+1 # Tuple shape

    dist_ny = np.arange(bigshape[0]) - self.ny
    dist_nx = np.arange(bigshape[1]) - self.nx

    dist_x,dist_y = np.meshgrid(dist_nx*self.dx,dist_ny*self.dy)

    bigdist = np.sqrt(dist_x**2 + dist_y**2) # Distances from center
                                          # Center at [ny,nx]
    
    biggrid = self.coeff * kei(bigdist/self.alpha) # Kelvin fcn solution

    # Now compute the deflections
    self.w = np.zeros((self.ny,self.nx)) # Deflection array
    for i in range(self.nx):
      for j in range(self.ny):
        # Loop over locations that have loads, and sum
        if self.q0[j,i]:
          # Solve by summing portions of "biggrid" while moving origin
          # to location of current cell
          # Load must be multiplied by grid cell size
          self.w += self.q0[j,i] * self.dx * self.dy \
             * biggrid[self.ny-j:2*self.ny-j,self.nx-i:2*self.nx-i]
      # No need to return: w already belongs to "self"

  # NO GRID

  def spatialDomainNoGrid(self):
  
    # Reassign q0 for consistency
    #self.x = self.q0[:,0]
    #self.y = self.q0[:,1]
    #self.q0 = self.q0[:,2]
    
    self.w = np.zeros(self.q0.shape)
    for i in range(len(self.x)):
      # Get the point
      x0 = self.x[i]
      y0 = self.y[i]
      # Create array of distances from point of load
      r = np.sqrt((self.x - x0)**2 + (self.y - y0)**2)
      # Compute and sum deflection
      self.w += self.q0[i] * self.coeff * kei(r/self.alpha)


  ## FINITE DIFFERENCE
  ######################
  
  def elasprep(self):
    """
    dx4, dy4, dx2dy2, D = elasprep(dx,dy,Te,E=1E11,nu=0.25)
    
    Defines the variables that are required to create the 2D finite 
    difference solution coefficient matrix
    """
    self.dx4 = self.dx**4
    self.dy4 = self.dy**4
    self.dx2dy2 = self.dx**2 * self.dy**2
    self.D = self.E*self.Te**3/(12*(1-self.nu**2))
  
  def coeff_matrix_creator(self):
    """
    coeff = coeff_matrix(D,drho,dx4,dy4,dx2dy2,nu=0.25,g=9.8)
    where D is the flexural rigidity, drho is the density difference between
    the mantle and the material filling the depression, nu is Poisson's ratio,
    g is gravitational acceleration at Earth's surface (approx. 9.8 m/s),
    and dx4, dy4, and dx2dy2 are based on the grid dimensions.
    
    All grid parameters except nu and g are generated by the function
    varprep2d, located inside this module

    Calculates the matrix of coefficients that is later used via sparse matrix
    solution techniques (scipy.sparse.linalg.spsolve) to compute the flexural
    response to the load. This step need only be performed once, and the
    coefficient matrix can very rapidly compute flexural solutions to any load.
    This makes this particularly good for probelms with time-variable loads or 
    that require iteration (e.g., water loading, in which additional water 
    causes subsidence, causes additional water detph, etc.).

    This method of coefficient matrix construction utilizes longer-range
    symmetry in the coefficient matrix to build it block-wise, as opposed to
    the much less computationally efficient row-by-row ("serial") method 
    that was previously employed.

    NOTATION FOR COEFFICIENT BIULDING MATRICES (e.g., "cj0i_2"):
    c = "coefficient
    j = columns = x-value
    j0 = no offset: look at center of array
    i = rows = y-value
    i_2 = negative 2 offset (i2 = positive 2 offset)
    """

    self.coeff_start_time = time.time()
    self.BC_selector_and_coeff_matrix_creator()
    self.coeff_creation_time = time.time() - self.coeff_start_time
    print 'Time to construct coefficient (operator) array [s]:', self.coeff_creation_time


  def get_coeff_values(self):
    """
    Build matrices containing all of the values for each of the coefficients
    that must be linearly combined to solve this equation
    13 coefficients: 13 matrices of the same size as the load
    """

    # don't want to keep typing "self." everwhere!
    D = self.D
    drho = self.drho
    dx4 = self.dx4
    dy4 = self.dy4
    dx2dy2 = self.dx2dy2
    nu = self.nu
    g = self.g

    if np.isscalar(self.Te):
      # So much simpler with constant D! And symmetrical stencil
      self.cj2i0 = D/dy4
      self.cj1i_1 = 2*D/dx2dy2
      self.cj1i0 = -4*D/dy4 - 4*D/dx2dy2
      self.cj1i1 = 2*D/dx2dy2
      self.cj0i_2 = D/dx4
      self.cj0i_1 = -4*D/dx4 - 4*D/dx2dy2
      self.cj0i0 = 6*D/dx4 + 6*D/dy4 + 8*D/dx2dy2 + drho*g
      self.cj0i1 = -4*D/dx4 - 4*D/dx2dy2 # Symmetry
      self.cj0i2 = D/dx4 # Symmetry
      self.cj_1i_1 = 2*D/dx2dy2 # Symmetry
      self.cj_1i0 = -4*D/dy4 - 4*D/dx2dy2 # Symmetry
      self.cj_1i1 = 2*D/dx2dy2 # Symmetry
      self.cj_2i0 = D/dy4 # Symmetry
      # Bring up to size
      self.cj2i0 *= np.ones(self.q0.shape)
      self.cj1i_1 *= np.ones(self.q0.shape)
      self.cj1i0 *= np.ones(self.q0.shape)
      self.cj1i1 *= np.ones(self.q0.shape)
      self.cj0i_2 *= np.ones(self.q0.shape)
      self.cj0i_1 *= np.ones(self.q0.shape)
      self.cj0i0 *= np.ones(self.q0.shape)
      self.cj0i1 *= np.ones(self.q0.shape)
      self.cj0i2 *= np.ones(self.q0.shape)
      self.cj_1i_1 *= np.ones(self.q0.shape)
      self.cj_1i0 *= np.ones(self.q0.shape)
      self.cj_1i1 *= np.ones(self.q0.shape)
      self.cj_2i0 *= np.ones(self.q0.shape)
      # Create coefficient arrays to manage boundary conditions
      self.cj2i0_coeff_ij = self.cj2i0.copy()
      self.cj1i_1_coeff_ij = self.cj1i_1.copy()
      self.cj1i0_coeff_ij = self.cj1i0.copy()
      self.cj1i1_coeff_ij = self.cj1i1.copy()
      self.cj0i_2_coeff_ij = self.cj0i_2.copy()
      self.cj0i_1_coeff_ij = self.cj0i_1.copy()
      self.cj0i0_coeff_ij = self.cj0i0.copy()
      self.cj0i1_coeff_ij = self.cj0i1.copy()
      self.cj0i2_coeff_ij = self.cj0i2.copy()
      self.cj_1i_1_coeff_ij = self.cj_1i_1.copy()
      self.cj_1i0_coeff_ij = self.cj_1i0.copy()
      self.cj_1i1_coeff_ij = self.cj_1i1.copy()
      self.cj_2i0_coeff_ij = self.cj_2i0.copy()
      
    elif type(self.Te) == np.ndarray:
    
    #######################################################
    # GENERATE COEFFICIENT VALUES FOR EACH SOLUTION TYPE. #
    #     "THICK" IS THE BEST: LOOSEST ASSUMPTIONS.       #
    #        OTHERS HERE LARGELY FOR COMPARISON           #
    #######################################################
    
      # All derivatives here, to make reading the equations below easier
      D00 = D[1:-1,1:-1]
      D10 = D[1:-1,2:]
      D_10 = D[1:-1,:-2]
      D01 = D[2:,1:-1]
      D0_1 = D[:-2,1:-1]
      D11 = D[2:,2:]
      D_11 = D[2:,:-2]
      D1_1 = D[:-2,2:]
      D_1_1 = D[:-2,:-2]
      print "VAR!"
      # Derivatives of D -- not including /(dx^a dy^b)
      D0  = D00
      Dx  = (-D_10 + D10)/2.
      Dy  = (-D0_1 + D01)/2.
      Dxx = (D_10 - 2.*D00 + D10)
      Dyy = (D0_1 - 2.*D00 + D01)
      Dxy = (D_1_1 - D_11 - D1_1 + D11)/4.
      
      if self.PlateSolutionType == 'Thick':
        # Check that it is thick
        # van Wees and Cloetingh (1994) solution, re-discretized by me
        # using a central difference approx. to 2nd order precision
        # NEW STENCIL
        # x = -2, y = 0
        self.cj_2i0_coeff_ij = (D0 - Dx) / dx4
        # x = 0, y = -2
        self.cj0i_2_coeff_ij = (D0 - Dy) / dy4
        # x = 0, y = 2
        self.cj0i2_coeff_ij = (D0 + Dy) / dy4
        # x = 2, y = 0
        self.cj2i0_coeff_ij = (D0 + Dx) / dx4
        # x = -1, y = -1
        self.cj_1i_1_coeff_ij = (2.*D0 - Dx - Dy + Dxy*(1-nu)/2.) / dx2dy2
        # x = -1, y = 1
        self.cj_1i1_coeff_ij = (2.*D0 - Dx + Dy - Dxy*(1-nu)/2.) / dx2dy2
        # x = 1, y = -1
        self.cj1i_1_coeff_ij = (2.*D0 + Dx - Dy - Dxy*(1-nu)/2.) / dx2dy2
        # x = 1, y = 1
        self.cj1i1_coeff_ij = (2.*D0 + Dx + Dy + Dxy*(1-nu)/2.) / dx2dy2
        # x = -1, y = 0
        self.cj_1i0_coeff_ij = (-4.*D0 + 2.*Dx + Dxx)/dx4 + (-4.*D0 + 2.*Dx + nu*Dyy)/dx2dy2
        # x = 0, y = -1
        self.cj0i_1_coeff_ij = (-4.*D0 + 2.*Dy + Dyy)/dy4 + (-4.*D0 + 2.*Dy + nu*Dxx)/dx2dy2
        # x = 0, y = 1
        self.cj0i1_coeff_ij = (-4.*D0 - 2.*Dy + Dyy)/dy4 + (-4.*D0 - 2.*Dy + nu*Dxx)/dx2dy2
        # x = 1, y = 0
        self.cj1i0_coeff_ij = (-4.*D0 - 2.*Dx + Dxx)/dx4 + (-4.*D0 - 2.*Dx + nu*Dyy)/dx2dy2
        # x = 0, y = 0
        self.cj0i0_coeff_ij = (6.*D0 - 2.*Dxx)/dx4 \
                     + (6.*D0 - 2.*Dyy)/dy4 \
                     + (8.*D0 - 2.*nu*Dxx - 2.*nu*Dyy)/dx2dy2 \
                     + drho*g
                     
      elif self.PlateSolutionType == 'LinearTeVariationsOnly':
        sys.exit("CHECK LATER PARTS OF SOLUTION: NOT SURE IF THEY ARE RIGHT")
        # So check starting with self.c_j1i0
        # These were flipped around in x and y
        # And need a good look over
        # before I will feel OK using them
        # More info here!!!!!!!!!
        # SIMPLER STENCIL: just del**2(D del**2(w)): only linear variations
        # in Te are allowed.
        # x = -2, y = 0
        self.cj_2i0_coeff_ij = D0 / dx4
        # x = 0, y = -2
        self.cj0i_2_coeff_ij = D0 / dy4
        # x = 0, y = 2
        self.cj0i2_coeff_ij = D0 / dy4
        # x = 2, y = 0
        self.cj2i0_coeff_ij = D0 / dx4
        # x = -1, y = -1
        self.cj_1i_1_coeff_ij = 2.*D0 / dx2dy2
        # x = -1, y = 1
        self.cj_1i1_coeff_ij = 2.*D0 / dx2dy2
        # x = 1, y = -1
        self.cj1i_1_coeff_ij = 2.*D0 / dx2dy2
        # x = 1, y = 1
        self.cj1i1_coeff_ij = 2.*D0 / dx2dy2
        # x = -1, y = 0
        self.cj_1i0_coeff_ij = (-4.*D0 + Dxx)/dx4 + (-4.*D0 + Dyy)/dx2dy2
        # x = 0, y = -1
        self.cj0i_1_coeff_ij = (-4.*D0 + Dyy)/dx4 + (-4.*D0 + Dxx)/dx2dy2
        # x = 0, y = 1
        self.cj0i1_coeff_ij = (-4.*D0 + Dyy)/dx4 + (-4.*D0 + Dxx)/dx2dy2
        # x = 1, y = 0
        self.cj1i0_coeff_ij = (-4.*D0 + Dxx)/dx4 + (-4.*D0 + Dyy)/dx2dy2
        # x = 0, y = 0
        self.cj0i0_coeff_ij = (6.*D0 - 2.*Dxx)/dx4 \
                     + (6.*D0 - 2.*Dyy)/dy4 \
                     + (8.*D0 - 2.*Dxx - 2.*Dyy)/dx2dy2 \
                     + drho*g

      elif self.PlateSolutionType == 'Thin':
        # STENCIL FROM GOVERS ET AL. 2009 -- first-order differences
        # x is j and y is i b/c matrix row/column notation
        # Note that this breaks down with b.c.'s that place too much control 
        # on the solution -- harmonic wavetrains
        # x = -2, y = 0
        self.cj_2i0_coeff_ij = D_10/dx4
        # x = -1, y = -1
        self.cj_1i_1_coeff_ij = (D_10 + D0_1)/dx2dy2
        # x = -1, y = 0
        self.cj_1i0_coeff_ij = -2. * ( (D0_1 + D00)/dx2dy2 + (D00 + D_10)/dx4 )
        # x = -1, y = 1
        self.cj_1i1_coeff_ij = (D_10 + D01)/dx2dy2
        # x = 0, y = -2
        self.cj0i_2_coeff_ij = D0_1/dy4
        # x = 0, y = -1
        self.cj0i_1_coeff_ij = -2. * ( (D0_1 + D00)/dx2dy2 + (D00 + D0_1)/dy4)
        # x = 0, y = 0
        self.cj0i0_coeff_ij = (D10 + 4.*D00 + D_10)/dx4 + (D01 + 4.*D00 + D0_1)/dy4 + (8.*D00/dx2dy2) + drho*g
        # x = 0, y = 1
        self.cj0i1_coeff_ij = -2. * ( (D01 + D00)/dy4 + (D00 + D01)/dx2dy2 )
        # x = 0, y = 2
        self.cj0i2_coeff_ij = D0_1/dy4
        # x = 1, y = -1
        self.cj1i_1_coeff_ij = (D10+D0_1)/dx2dy2
        # x = 1, y = 0
        self.cj1i0_coeff_ij = -2. * ( (D10 + D00)/dx4 + (D10 + D00)/dx2dy2 )
        # x = 1, y = 1
        self.cj1i1_coeff_ij = (D10 + D01)/dx2dy2
        # x = 2, y = 0
        self.cj2i0_coeff_ij = D10/dx4

      ################################################################
      # CREATE COEFFICIENT ARRAYS: PLAIN, WITH NO B.C.'S YET APPLIED #
      ################################################################
      # x = -2, y = 0
      self.cj_2i0 = self.cj_2i0_coeff_ij.copy()
      # x = -1, y = -1
      self.cj_1i_1 = self.cj_1i_1_coeff_ij.copy()
      # x = -1, y = 0
      self.cj_1i0 = self.cj_1i0_coeff_ij.copy()
      # x = -1, y = 1
      self.cj_1i1 = self.cj_1i1_coeff_ij.copy()
      # x = 0, y = -2
      self.cj0i_2 = self.cj0i_2_coeff_ij.copy()
      # x = 0, y = -1
      self.cj0i_1 = self.cj0i_1_coeff_ij.copy()
      # x = 0, y = 0
      self.cj0i0 = self.cj0i0_coeff_ij.copy()
      # x = 0, y = 1
      self.cj0i1 = self.cj0i1_coeff_ij.copy()
      # x = 0, y = 2
      self.cj0i2 = self.cj0i2_coeff_ij.copy()
      # x = 1, y = -1
      self.cj1i_1 = self.cj1i_1_coeff_ij.copy()
      # x = 1, y = 0
      self.cj1i0 = self.cj1i0_coeff_ij.copy()
      # x = 1, y = 1
      self.cj1i1 = self.cj1i1_coeff_ij.copy()
      # x = 2, y = 0
      self.cj2i0 = self.cj2i0_coeff_ij.copy()

    # Provide rows and columns in the 2D input to later functions
    self.ncolsx = self.cj0i0.shape[1]
    self.nrowsy = self.cj0i0.shape[0]

  def BC_selector_and_coeff_matrix_creator(self):
    """
    Selects the boundary conditions
    E-W is for inside each panel
    N-S is for the block diagonal matrix ("with fringes")
    Then calls the function to build the diagonal matrix
    """

    # Zeroth, print the boundary conditions to the screen
    print "Boundary condition, West:", self.BC_W, type(self.BC_W)
    print "Boundary condition, East:", self.BC_E, type(self.BC_E)
    print "Boundary condition, North:", self.BC_N, type(self.BC_N)
    print "Boundary condition, South:", self.BC_S, type(self.BC_S)
    
    ################################
    # PERIODIC B.C. VALIDITY CHECK #
    ################################

    # First check to make sure that periodic boundary conditions are applied 
    # properly: and abort if they are not.
    
    # E-W CHECK
    if ((self.BC_W == 'Periodic' and self.BC_E != 'Periodic') \
      or (self.BC_W != 'Periodic' and self.BC_E == 'Periodic')) \
      and (self.BC_W != 'Mirror' and self.BC_E != 'Mirror'):
      # If only one boundary is periodic and the other doesn't implicitly 
      # involve a periodic boundary, this is illegal!
      sys.exit("Having the boundary opposite a periodic boundary condition\n"+
               "be fixed and not include an implicit periodic boundary\n"+
               "condition makes no physical sense.\n"+
               "Please fix the input boundary conditions. Aborting.")

    # N-S CHECK
    if ((self.BC_N == 'Periodic' and self.BC_S != 'Periodic') \
      or (self.BC_N != 'Periodic' and self.BC_S == 'Periodic')) \
      and (self.BC_N != 'Mirror' and self.BC_S != 'Mirror'):
      # If only one boundary is periodic and the other doesn't implicitly 
      # involve a periodic boundary, this is illegal!
      sys.exit("Having the boundary opposite a periodic boundary condition\n"+
               "be fixed and not include an implicit periodic boundary\n"+
               "condition makes no physical sense.\n"+
               "Please fix the input boundary conditions. Aborting.")


    # Need to hold on to q0 array to keep its size; self.q0 will be redefined 
    # for the computation
    self.q0_orig = self.q0.copy()
    
    # Now I can build the coefficient arrays
    self.get_coeff_values()
      
    # Rigidity b.c.?
    self.BC_Rigidity()

    self.build_diags()

  def BC_Rigidity(self):
    """
    Utility function to help implement boundary conditions by specifying 
    them for and applying them to the elastic thickness grid
    """

    ##############################################################
    # AUTOMATICALLY SELECT FLEXURAL RIGIDITY BOUNDARY CONDITIONS #
    ##############################################################
    # West
    if self.BC_W == 'Periodic':
      self.BC_Rigidity_W = 'periodic'
    elif (self.BC_W == np.array(['Dirichlet0', '0Moment0Shear', '0Slope0Shear'])).any():
      self.BC_Rigidity_W = '0 curvature'
    elif self.BC_W == 'Mirror':
      self.BC_Rigidity_W = 'mirror symmetry'
    else:
      sys.exit("Invalid Te B.C. case")
    # East
    if self.BC_E == 'Periodic':
      self.BC_Rigidity_E = 'periodic'
    elif (self.BC_E == np.array(['Dirichlet0', '0Moment0Shear', '0Slope0Shear'])).any():
      self.BC_Rigidity_E = '0 curvature'
    elif self.BC_E == 'Mirror':
      self.BC_Rigidity_E = 'mirror symmetry'
    else:
      sys.exit("Invalid Te B.C. case")
    # North
    if self.BC_N == 'Periodic':
      self.BC_Rigidity_N = 'periodic'
    elif (self.BC_N == np.array(['Dirichlet0', '0Moment0Shear', '0Slope0Shear'])).any():
      self.BC_Rigidity_N = '0 curvature'
    elif self.BC_N == 'Mirror':
      self.BC_Rigidity_N = 'mirror symmetry'
    else:
      sys.exit("Invalid Te B.C. case")
    # South
    if self.BC_S == 'Periodic':
      self.BC_Rigidity_S = 'periodic'
    elif (self.BC_S == np.array(['Dirichlet0', '0Moment0Shear', '0Slope0Shear'])).any():
      self.BC_Rigidity_S = '0 curvature'
    elif self.BC_S == 'Mirror':
      self.BC_Rigidity_S = 'mirror symmetry'
    else:
      sys.exit("Invalid Te B.C. case")
    
    #############
    # PAD ARRAY #
    #############
    #self.D = np.hstack(( np.nan*np.zeros((self.D.shape[0], 1)), self.D, np.nan*np.zeros((self.D.shape[0], 1)) ))
    #self.D = np.vstack(( np.nan*np.zeros(self.D.shape[1]), self.D, np.nan*np.zeros(self.D.shape[1]) ))
    if np.isscalar(self.Te):
      # No longer scalar! But already went through get_coeff_values() step, 
      # which is the other important check
      # Though of course constant Te is just a special case and can be solved
      # by the general solution methods
      self.D *= np.ones(self.q0.shape)
      # THIS WILL GO OUTSIDE IF-STATEMENT ONCE I CHANGE THE SHAPE OF THE
      # Te ARRAYS!!!
      self.D = np.hstack(( np.nan*np.zeros((self.D.shape[0], 1)), self.D, np.nan*np.zeros((self.D.shape[0], 1)) ))
      self.D = np.vstack(( np.nan*np.zeros(self.D.shape[1]), self.D, np.nan*np.zeros(self.D.shape[1]) ))
    else:
      # Temporarily:
      self.D[:,0] = np.nan
      self.D[:,-1] = np.nan
      self.D[0,:] = np.nan
      self.D[-1,:] = np.nan

    ###############################################################
    # APPLY FLEXURAL RIGIDITY BOUNDARY CONDITIONS TO PADDED ARRAY #
    ###############################################################
    if self.BC_Rigidity_W == "0 curvature":
      self.D[:,0] = 2*self.D[:,1] - self.D[:,2]
    if self.BC_Rigidity_E == "0 curvature":
      self.D[:,-1] = 2*self.D[:,-2] - self.D[:,-3]
    if self.BC_Rigidity_N == "0 curvature":
      self.D[0,:] = 2*self.D[1,:] - self.D[2,:]
    if self.BC_Rigidity_S == "0 curvature":
      self.D[-1,:] = 2*self.D[-2,:] - self.D[-3,:]

    if self.BC_Rigidity_W == "mirror symmetry":
      self.D[:,0] = self.D[:,2]
    if self.BC_Rigidity_E == "mirror symmetry":
      self.D[:,-1] = self.D[:,-3]
    if self.BC_Rigidity_N == "mirror symmetry":
      self.D[0,:] = self.D[2,:] # Yes, will work on corners -- double-reflection
    if self.BC_Rigidity_S == "mirror symmetry":
      self.D[-1,:] = self.D[-3,:]
      
    if self.BC_Rigidity_W == "periodic":
      self.D[:,0] = self.D[:,-2]
    if self.BC_Rigidity_E == "periodic":
      self.D[:,-1] = self.D[:,-3]
    if self.BC_Rigidity_N == "periodic":
      self.D[0,:] = self.D[-2,:]
    if self.BC_Rigidity_S == "periodic":
      self.D[-1,:] = self.D[-3,:]
      
### NEW FCN HERE
### AND TRY REMOVING THE REPEAT CODE IN 1D CASE!

    # The next section of code is split over several functions for the 1D 
    # case, but will be all in one function here, at least for now.
    
    # Inf for E-W to separate from nan for N-S. N-S will spill off ends
    # of array (C order, in rows), while E-W will be internal, so I will
    # later change np.inf to 0 to represent where internal boundaries 
    # occur.

    #######################################################################
    # DEFINE COEFFICIENTS TO W_j-2 -- W_j+2 WITH B.C.'S APPLIED (x: W, E) #
    #######################################################################
    
    if self.BC_W == 'Periodic':
      pass
    elif self.BC_W == 'Dirichlet0':
      j = 0
      self.cj_2i0[:,j] += np.inf
      self.cj_1i_1[:,j] += np.inf
      self.cj_1i0[:,j] += np.inf
      self.cj_1i1[:,j] += np.inf
      self.cj0i_2[:,j] += 0
      self.cj0i_1[:,j] += 0
      self.cj0i0[:,j] += 0
      self.cj0i1[:,j] += 0
      self.cj0i2[:,j] += 0
      self.cj1i_1[:,j] += 0
      self.cj1i0[:,j] += 0
      self.cj1i1[:,j] += 0
      self.cj2i0[:,j] += 0
      j = 1
      self.cj_2i0[:,j] += np.inf
      self.cj_1i_1[:,j] += 0
      self.cj_1i0[:,j] += 0
      self.cj_1i1[:,j] += 0
      self.cj0i_2[:,j] += 0
      self.cj0i_1[:,j] += 0
      self.cj0i0[:,j] += 0
      self.cj0i1[:,j] += 0
      self.cj0i2[:,j] += 0
      self.cj1i_1[:,j] += 0
      self.cj1i0[:,j] += 0
      self.cj1i1[:,j] += 0
      self.cj2i0[:,j] += 0
    elif self.BC_W == '0Moment0Shear':
      j = 0
      self.cj_2i0[:,j] += np.inf
      self.cj_1i_1[:,j] += np.inf
      self.cj_1i0[:,j] += np.inf
      self.cj_1i1[:,j] += np.inf
      self.cj0i_2[:,j] += 0
      self.cj0i_1[:,j] += 2*self.cj_1i_1[:,j]
      self.cj0i0[:,j] += 4*self.cj_2i0_coeff_ij[:,j] + 2*self.cj_1i0_coeff_ij[:,j]
      self.cj0i1[:,j] += 2*self.cj_1i1[:,j]
      self.cj0i2[:,j] += 0
      self.cj1i_1[:,j] += -self.cj_1i_1[:,j]
      self.cj1i0[:,j] += -4*self.cj_2i0_coeff_ij[:,j] - self.cj_1i0_coeff_ij[:,j]
      self.cj1i1[:,j] += -self.cj_1i1[:,j]
      self.cj2i0[:,j] += self.cj_2i0_coeff_ij[:,j]
      j = 1
      self.cj_2i0[:,j] += np.inf
      self.cj_1i_1[:,j] += 0
      self.cj_1i0[:,j] += 2*self.cj_2i0_coeff_ij[:,j]
      self.cj_1i1[:,j] += 0
      self.cj0i_2[:,j] += 0
      self.cj0i_1[:,j] += 0
      self.cj0i0[:,j] += 0
      self.cj0i1[:,j] += 0
      self.cj0i2[:,j] += 0
      self.cj1i_1[:,j] += 0
      self.cj1i0[:,j] += -2*self.cj_2i0_coeff_ij[:,j]
      self.cj1i1[:,j] += 0
      self.cj2i0[:,j] += self.cj_2i0_coeff_ij[:,j]
    elif self.BC_W == '0Slope0Shear':
      j = 0
      self.cj_2i0[:,j] += np.inf
      self.cj_1i_1[:,j] += np.inf
      self.cj_1i0[:,j] += np.inf
      self.cj_1i1[:,j] += np.inf
      self.cj0i_2[:,j] += 0
      self.cj0i_1[:,j] += 0
      self.cj0i0[:,j] += 0
      self.cj0i1[:,j] += 0
      self.cj0i2[:,j] += 0
      self.cj1i_1[:,j] += self.cj_1i_1_coeff_ij[:,j] 
      self.cj1i0[:,j] += self.cj_1i0_coeff_ij[:,j]
      self.cj1i1[:,j] += self.cj_1i1_coeff_ij[:,j] #Interference
      self.cj2i0[:,j] += self.cj_2i0_coeff_ij[:,j]
      j = 1
      self.cj_2i0[:,j] += np.inf
      self.cj_1i_1[:,j] += 0
      self.cj_1i0[:,j] += 0
      self.cj_1i1[:,j] += 0
      self.cj0i_2[:,j] += 0
      self.cj0i_1[:,j] += 0
      self.cj0i0[:,j] += 0
      self.cj0i1[:,j] += 0
      self.cj0i2[:,j] += 0
      self.cj1i_1[:,j] += 0
      self.cj1i0[:,j] += 0
      self.cj1i1[:,j] += 0
      self.cj2i0[:,j] += self.cj_2i0_coeff_ij[:,j]
    elif self.BC_W == 'Mirror':
      j = 0
      self.cj_2i0[:,j] += np.inf
      self.cj_1i_1[:,j] += np.inf
      self.cj_1i0[:,j] += np.inf
      self.cj_1i1[:,j] += np.inf
      self.cj0i_2[:,j] += 0
      self.cj0i_1[:,j] += 0
      self.cj0i0[:,j] += 0
      self.cj0i1[:,j] += 0
      self.cj0i2[:,j] += 0
      self.cj1i_1[:,j] += self.cj_1i_1_coeff_ij[:,j] 
      self.cj1i0[:,j] += self.cj_1i0_coeff_ij[:,j]
      self.cj1i1[:,j] += self.cj_1i1_coeff_ij[:,j]
      self.cj2i0[:,j] += self.cj_2i0_coeff_ij[:,j]
      j = 1
      self.cj_2i0[:,j] += np.inf
      self.cj_1i_1[:,j] += 0
      self.cj_1i0[:,j] += 0
      self.cj_1i1[:,j] += 0
      self.cj0i_2[:,j] += 0
      self.cj0i_1[:,j] += 0
      self.cj0i0[:,j] += self.cj_2i0_coeff_ij[:,j]
      self.cj0i1[:,j] += 0
      self.cj0i2[:,j] += 0
      self.cj1i_1[:,j] += 0
      self.cj1i0[:,j] += 0
      self.cj1i1[:,j] += 0
      self.cj2i0[:,j] += 0
    else:
      # Possibly redundant safeguard
      sys.exit("Invalid boundary condition")

    if self.BC_E == 'Periodic':
      pass
    elif self.BC_E == 'Dirichlet0':
      j = -1
      self.cj_2i0[:,j] += 0
      self.cj_1i_1[:,j] += 0
      self.cj_1i0[:,j] += 0
      self.cj_1i1[:,j] += 0
      self.cj0i_2[:,j] += 0
      self.cj0i_1[:,j] += 0
      self.cj0i0[:,j] += 0
      self.cj0i1[:,j] += 0
      self.cj0i2[:,j] += 0
      self.cj1i_1[:,j] += np.inf
      self.cj1i0[:,j] += np.inf
      self.cj1i1[:,j] += np.inf
      self.cj2i0[:,j] += np.inf
      j = -2
      self.cj_2i0[:,j] += 0
      self.cj_1i_1[:,j] += 0
      self.cj_1i0[:,j] += 0
      self.cj_1i1[:,j] += 0
      self.cj0i_2[:,j] += 0
      self.cj0i_1[:,j] += 0
      self.cj0i0[:,j] += 0
      self.cj0i1[:,j] += 0
      self.cj0i2[:,j] += 0
      self.cj1i_1[:,j] += 0
      self.cj1i0[:,j] += 0
      self.cj1i1[:,j] += 0
      self.cj2i0[:,j] += np.inf
    elif self.BC_E == '0Moment0Shear':
      j = -1
      self.cj_2i0[:,j] += self.cj2i0_coeff_ij[:,j]
      self.cj_1i_1[:,j] += -self.cj1i_1_coeff_ij[:,j]
      self.cj_1i0[:,j] += -4*self.cj2i0_coeff_ij[:,j] - self.cj1i0_coeff_ij[:,j]
      self.cj_1i1[:,j] += -self.cj1i1_coeff_ij[:,j]
      self.cj0i_2[:,j] += 0
      self.cj0i_1[:,j] += 2*self.cj1i_1_coeff_ij[:,j]
      self.cj0i0[:,j] += 4*self.cj2i0_coeff_ij[:,j] + 2*self.cj1i0_coeff_ij[:,j]
      self.cj0i1[:,j] += 2*self.cj1i1_coeff_ij[:,j]
      self.cj0i2[:,j] += 0
      self.cj1i_1[:,j] += np.inf
      self.cj1i0[:,j] += np.inf
      self.cj1i1[:,j] += np.inf
      self.cj2i0[:,j] += np.inf
      j = -2
      self.cj_2i0[:,j] += self.cj2i0_coeff_ij[:,j]
      self.cj_1i_1[:,j] += 0
      self.cj_1i0[:,j] += -2*self.cj2i0_coeff_ij[:,j]
      self.cj_1i1[:,j] += 0
      self.cj0i_2[:,j] += 0
      self.cj0i_1[:,j] += 0
      self.cj0i0[:,j] += 0
      self.cj0i1[:,j] += 0
      self.cj0i2[:,j] += 0
      self.cj1i_1[:,j] += 0
      self.cj1i0[:,j] += 2*self.cj2i0_coeff_ij[:,j]
      self.cj1i1[:,j] += 0
      self.cj2i0[:,j] += np.inf
    elif self.BC_E == '0Slope0Shear':
      j = -1
      self.cj_2i0[:,j] += self.cj2i0_coeff_ij[:,j]
      self.cj_1i_1[:,j] += self.cj1i_1_coeff_ij[:,j]
      self.cj_1i0[:,j] += self.cj1i0_coeff_ij[:,j]
      self.cj_1i1[:,j] += self.cj1i1_coeff_ij[:,j]
      self.cj0i_2[:,j] += 0
      self.cj0i_1[:,j] += 0
      self.cj0i0[:,j] += 0
      self.cj0i1[:,j] += 0
      self.cj0i2[:,j] += 0
      self.cj1i_1[:,j] += np.inf
      self.cj1i0[:,j] += np.inf
      self.cj1i1[:,j] += np.inf
      self.cj2i0[:,j] += np.inf
      j = -2
      self.cj_2i0[:,j] += self.cj2i0_coeff_ij[:,j]
      self.cj_1i_1[:,j] += 0
      self.cj_1i0[:,j] += 0
      self.cj_1i1[:,j] += 0
      self.cj0i_2[:,j] += 0
      self.cj0i_1[:,j] += 0
      self.cj0i0[:,j] += 0
      self.cj0i1[:,j] += 0
      self.cj0i2[:,j] += 0
      self.cj1i_1[:,j] += 0
      self.cj1i0[:,j] += 0
      self.cj1i1[:,j] += 0
      self.cj2i0[:,j] += np.inf
    elif self.BC_E == 'Mirror':
      j = -1
      self.cj_2i0[:,j] += self.cj2i0_coeff_ij[:,j]
      self.cj_1i_1[:,j] += self.cj1i_1_coeff_ij[:,j]
      self.cj_1i0[:,j] += self.cj1i0_coeff_ij[:,j]
      self.cj_1i1[:,j] += self.cj1i1_coeff_ij[:,j]
      self.cj0i_2[:,j] += 0
      self.cj0i_1[:,j] += 0
      self.cj0i0[:,j] += 0
      self.cj0i1[:,j] += 0
      self.cj0i2[:,j] += 0
      self.cj1i_1[:,j] += np.inf
      self.cj1i0[:,j] += np.inf
      self.cj1i1[:,j] += np.inf
      self.cj2i0[:,j] += np.inf
      j = -2
      self.cj_2i0[:,j] += 0
      self.cj_1i_1[:,j] += 0
      self.cj_1i0[:,j] += 0
      self.cj_1i1[:,j] += 0
      self.cj0i_2[:,j] += 0
      self.cj0i_1[:,j] += 0
      self.cj0i0[:,j] += self.cj2i0_coeff_ij[:,j]
      self.cj0i1[:,j] += 0
      self.cj0i2[:,j] += 0
      self.cj1i_1[:,j] += 0
      self.cj1i0[:,j] += 0
      self.cj1i1[:,j] += 0
      self.cj2i0[:,j] += np.inf
    else:
      # Possibly redundant safeguard
      sys.exit("Invalid boundary condition")

    #######################################################################
    # DEFINE COEFFICIENTS TO W_i-2 -- W_i+2 WITH B.C.'S APPLIED (y: N, S) #
    #######################################################################
    
    # Infinitiy dominates over nan where they cross, because this is used 
    # to flag places where coeff values should be 0, and this would otherwise
    # cause boundary condition nan's to appear in the cross-derivatives

    if self.BC_N == 'Periodic':
      pass
    elif self.BC_N == 'Dirichlet0':
      i = 0
      self.cj_2i0[i,:] += 0
      self.cj_1i_1[i,:][self.cj_1i_1[i,:] != np.inf] = np.nan
      self.cj_1i0[i,:] += 0
      self.cj_1i1[i,:] += 0
      self.cj0i_2[i,:] += np.nan
      self.cj0i_1[i,:] += np.nan
      self.cj0i0[i,:] += 0
      self.cj0i1[i,:] += 0
      self.cj0i2[i,:] += 0
      self.cj1i_1[i,:][self.cj1i_1[i,:] != np.inf] += np.nan
      self.cj1i0[i,:] += 0
      self.cj1i1[i,:] += 0
      self.cj2i0[i,:] += 0
      i = 1
      self.cj_2i0[i,:] += 0
      self.cj_1i_1[i,:] += 0
      self.cj_1i0[i,:] += 0
      self.cj_1i1[i,:] += 0
      self.cj0i_2[i,:] += np.nan
      self.cj0i_1[i,:] += 0
      self.cj0i0[i,:] += 0
      self.cj0i1[i,:] += 0
      self.cj0i2[i,:] += 0
      self.cj1i_1[i,:] += 0
      self.cj1i0[i,:] += 0
      self.cj1i1[i,:] += 0
      self.cj2i0[i,:] += 0
    elif self.BC_N == '0Moment0Shear':
      i = 0
      self.cj_2i0[i,:] += 0
      self.cj_1i_1[i,:][self.cj_1i_1[i,:] != np.inf] = np.nan
      self.cj_1i0[i,:] += 2*self.cj_1i_1_coeff_ij[i,:]
      self.cj_1i1[i,:] += -self.cj_1i_1_coeff_ij[i,:]
      self.cj0i_2[i,:] += np.nan
      self.cj0i_1[i,:] += np.nan
      self.cj0i0[i,:] += 4*self.cj0i_2_coeff_ij[i,:] + 2*self.cj0i_1_coeff_ij[i,:]
      self.cj0i1[i,:] += -4*self.cj0i_2_coeff_ij[i,:] - self.cj0i_1_coeff_ij[i,:]
      self.cj0i2[i,:] += self.cj0i_2_coeff_ij[i,:]
      self.cj1i_1[i,:][self.cj1i_1[i,:] != np.inf] += np.nan
      self.cj1i0[i,:] += 2*self.cj1i_1_coeff_ij[i,:]
      self.cj1i1[i,:] += -self.cj1i_1_coeff_ij[i,:]
      self.cj2i0[i,:] += 0
      i = 1
      self.cj_2i0[i,:] += 0
      self.cj_1i_1[i,:] += 0
      self.cj_1i0[i,:] += 0
      self.cj_1i1[i,:] += 0
      self.cj0i_2[i,:] += np.nan
      self.cj0i_1[i,:] += 2*self.cj0i_2_coeff_ij[i,:]
      self.cj0i0[i,:] += 0
      self.cj0i1[i,:] += -2*self.cj0i_2_coeff_ij[i,:]
      self.cj0i2[i,:] += self.cj0i_2_coeff_ij[i,:]
      self.cj1i_1[i,:] += 0
      self.cj1i0[i,:] += 0
      self.cj1i1[i,:] += 0
      self.cj2i0[i,:] += 0
    elif self.BC_N == '0Slope0Shear':
      i = 0
      self.cj_2i0[i,:] += 0
      self.cj_1i_1[i,:][self.cj_1i_1[i,:] != np.inf] = np.nan
      self.cj_1i0[i,:] += 0
      self.cj_1i1[i,:] += self.cj_1i_1_coeff_ij[i,:]
      self.cj0i_2[i,:] += np.nan
      self.cj0i_1[i,:] += np.nan
      self.cj0i0[i,:] += 0
      self.cj0i1[i,:] += self.cj0i_1_coeff_ij[i,:]
      self.cj0i2[i,:] += self.cj0i_2_coeff_ij[i,:]
      self.cj1i_1[i,:][self.cj1i_1[i,:] != np.inf] += np.nan
      self.cj1i0[i,:] += 0
      self.cj1i1[i,:] += self.cj1i_1_coeff_ij[i,:]
      self.cj2i0[i,:] += 0
      i = 1
      self.cj_2i0[i,:] += 0
      self.cj_1i_1[i,:] += 0
      self.cj_1i0[i,:] += 0
      self.cj_1i1[i,:] += 0
      self.cj0i_2[i,:] += np.nan
      self.cj0i_1[i,:] += 0
      self.cj0i0[i,:] += 0
      self.cj0i1[i,:] += 0
      self.cj0i2[i,:] += self.cj0i_2_coeff_ij[i,:]
      self.cj1i_1[i,:] += 0
      self.cj1i0[i,:] += 0
      self.cj1i1[i,:] += 0
      self.cj2i0[i,:] += 0
    elif self.BC_N == 'Mirror':
      i = 0
      self.cj_2i0[i,:] += 0
      self.cj_1i_1[i,:][self.cj_1i_1[i,:] != np.inf] = np.nan
      self.cj_1i0[i,:] += 0
      self.cj_1i1[i,:] += self.cj_1i_1_coeff_ij[i,:]
      self.cj0i_2[i,:] += np.nan
      self.cj0i_1[i,:] += np.nan
      self.cj0i0[i,:] += 0
      self.cj0i1[i,:] += self.cj0i_1_coeff_ij[i,:]
      self.cj0i2[i,:] += self.cj0i_2_coeff_ij[i,:]
      self.cj1i_1[i,:][self.cj1i_1[i,:] != np.inf] += np.nan
      self.cj1i0[i,:] += 0
      self.cj1i1[i,:] += self.cj1i_1_coeff_ij[i,:]
      self.cj2i0[i,:] += 0
      i = 1
      self.cj_2i0[i,:] += 0
      self.cj_1i_1[i,:] += 0
      self.cj_1i0[i,:] += 0
      self.cj_1i1[i,:] += 0
      self.cj0i_2[i,:] += np.nan
      self.cj0i_1[i,:] += 0
      self.cj0i0[i,:] += self.cj0i_2_coeff_ij[i,:]
      self.cj0i1[i,:] += 0
      self.cj0i2[i,:] += 0
      self.cj1i_1[i,:] += 0
      self.cj1i0[i,:] += 0
      self.cj1i1[i,:] += 0
      self.cj2i0[i,:] += 0
    else:
      # Possibly redundant safeguard
      sys.exit("Invalid boundary condition")

    if self.BC_S == 'Periodic':
      pass
    elif self.BC_S == 'Dirichlet0':
      i = -2
      self.cj_2i0[i,:] += 0
      self.cj_1i_1[i,:] += 0
      self.cj_1i0[i,:] += 0
      self.cj_1i1[i,:] += 0
      self.cj0i_2[i,:] += 0
      self.cj0i_1[i,:] += 0
      self.cj0i0[i,:] += 0
      self.cj0i1[i,:] += 0
      self.cj0i2[i,:] += np.nan
      self.cj1i1[i,:] += 0
      self.cj1i0[i,:] += 0
      self.cj1i_1[i,:] += 0
      self.cj2i0[i,:] += 0
      i = -1
      self.cj_2i0[i,:] += 0
      self.cj_1i_1[i,:] += 0
      self.cj_1i0[i,:] += 0
      self.cj_1i1[i,:][self.cj_1i1[i,:] != np.inf] += np.nan
      self.cj0i_2[i,:] += 0
      self.cj0i_1[i,:] += 0
      self.cj0i0[i,:] += 0
      self.cj0i1[i,:] += np.nan
      self.cj0i2[i,:] += np.nan
      self.cj1i_1[i,:] += 0
      self.cj1i0[i,:] += 0
      self.cj1i1[i,:][self.cj1i1[i,:] != np.inf] += np.nan
      self.cj2i0[i,:] += 0
    elif self.BC_S == '0Moment0Shear':
      i = -2
      self.cj_2i0[i,:] += 0
      self.cj_1i_1[i,:] += 0
      self.cj_1i0[i,:] += 0
      self.cj_1i1[i,:] += 0
      self.cj0i_2[i,:] += self.cj0i2_coeff_ij[i,:]
      self.cj0i_1[i,:] += -2*self.cj0i2_coeff_ij[i,:]
      self.cj0i0[i,:] += 0
      self.cj0i1[i,:] += 2*self.cj0i2_coeff_ij[i,:]
      self.cj0i2[i,:] += np.nan
      self.cj1i1[i,:] += 0
      self.cj1i0[i,:] += 0
      self.cj1i_1[i,:] += 0
      self.cj2i0[i,:] += 0
      i = -1
      self.cj_2i0[i,:] += 0
      self.cj_1i_1[i,:] += -self.cj1i1_coeff_ij[i,:]
      self.cj_1i0[i,:] += 2*self.cj1i1_coeff_ij[i,:]
      self.cj_1i1[i,:][self.cj_1i1[i,:] != np.inf] += np.nan
      self.cj0i_2[i,:] += self.cj0i2_coeff_ij[i,:]
      self.cj0i_1[i,:] += -4*self.cj0i2_coeff_ij[i,:] - self.cj0i1_coeff_ij[i,:]
      self.cj0i0[i,:] += 4*self.cj0i2_coeff_ij[i,:] + 2*self.cj0i1_coeff_ij[i,:]
      self.cj0i1[i,:] += np.nan
      self.cj0i2[i,:] += np.nan
      self.cj1i_1[i,:] += -self.cj_1i1_coeff_ij[i,:]
      self.cj1i0[i,:] += 2*self.cj_1i1_coeff_ij[i,:]
      self.cj1i1[i,:][self.cj1i1[i,:] != np.inf] += np.nan
      self.cj2i0[i,:] += 0
    elif self.BC_S == '0Slope0Shear':
      i = -2
      self.cj_2i0[i,:] += 0
      self.cj_1i_1[i,:] += 0
      self.cj_1i0[i,:] += 0
      self.cj_1i1[i,:] += 0
      self.cj0i_2[i,:] += self.cj0i2_coeff_ij[i,:]
      self.cj0i_1[i,:] += 0
      self.cj0i0[i,:] += 0
      self.cj0i1[i,:] += 0
      self.cj0i2[i,:] += np.nan
      self.cj1i_1[i,:] += 0
      self.cj1i0[i,:] += 0
      self.cj1i1[i,:] += 0
      self.cj2i0[i,:] += 0
      i = -1
      self.cj_2i0[i,:] += 0
      self.cj_1i_1[i,:] += self.cj_1i1_coeff_ij[i,:]
      self.cj_1i0[i,:] += 0
      self.cj_1i1[i,:][self.cj_1i1[i,:] != np.inf] += np.nan
      self.cj0i_2[i,:] += self.cj0i2_coeff_ij[i,:]
      self.cj0i_1[i,:] += self.cj0i1_coeff_ij[i,:]
      self.cj0i0[i,:] += 0
      self.cj0i1[i,:] += np.nan
      self.cj0i2[i,:] += np.nan
      self.cj1i_1[i,:] += self.cj1i1_coeff_ij[i,:]
      self.cj1i0[i,:] += 0
      self.cj1i1[i,:][self.cj1i1[i,:] != np.inf] += np.nan
      self.cj2i0[i,:] += 0
    elif self.BC_S == 'Mirror':
      i = -2
      self.cj_2i0[i,:] += 0
      self.cj_1i_1[i,:] += 0
      self.cj_1i0[i,:] += 0
      self.cj_1i1[i,:] += 0
      self.cj0i_2[i,:] += 0
      self.cj0i_1[i,:] += 0
      self.cj0i0[i,:] += self.cj0i2_coeff_ij[i,:]
      self.cj0i1[i,:] += 0
      self.cj0i2[i,:] += np.nan
      self.cj1i_1[i,:] += 0
      self.cj1i0[i,:] += 0
      self.cj1i1[i,:] += 0
      self.cj2i0[i,:] += 0
      i = -1
      self.cj_2i0[i,:] += 0
      self.cj_1i_1[i,:] += self.cj_1i1_coeff_ij[i,:]
      self.cj_1i0[i,:] += 0
      self.cj_1i1[i,:][self.cj_1i1[i,:] != np.inf] += np.nan
      self.cj0i_2[i,:] += self.cj0i2_coeff_ij[i,:]
      self.cj0i_1[i,:] += self.cj0i1_coeff_ij[i,:]
      self.cj0i0[i,:] += 0
      self.cj0i1[i,:] += np.nan
      self.cj0i2[i,:] += np.nan
      self.cj1i_1[i,:] += self.cj1i1_coeff_ij[i,:]
      self.cj1i0[i,:] += 0
      self.cj1i1[i,:][self.cj1i1[i,:] != np.inf] += np.nan
      self.cj2i0[i,:] += 0
    else:
      # Possibly redundant safeguard
      sys.exit("Invalid boundary condition")

    #####################################################
    # CORNERS: INTERFERENCE BETWEEN BOUNDARY CONDITIONS #
    #####################################################
    
    # In 2D, have to consider diagonals and interference (additive) among 
    # boundary conditions
    # DIRICHLET -- DO NOTHING.
    # 0Slope0Shear -- is this generic (i.e. should I just change to != Dirichlet0?)
    # How do multiple types of b.c.'s interfere?
    if self.BC_N == '0Slope0Shear' and self.BC_W == '0Slope0Shear':
      self.cj1i1[0,0] += self.cj_1i_1_coeff_ij[0,0]
    if self.BC_N == '0Slope0Shear' and self.BC_E == '0Slope0Shear':
      self.cj_1i1[0,-1] += self.cj1i_1_coeff_ij[0,-1]
    if self.BC_S == '0Slope0Shear' and self.BC_W == '0Slope0Shear':
      self.cj1i_1[-1,0] += self.cj_1i1_coeff_ij[-1,0]
    if self.BC_S == '0Slope0Shear' and self.BC_E == '0Slope0Shear':
      self.cj_1i_1[-1,-1] += self.cj1i1_coeff_ij[-1,-1]
    # 0Moment0Shear
    # Mirror
    # Currently the same as 0Slope0Shear -- a placeholder.
    if self.BC_N == 'Mirror' and self.BC_W == 'Mirror':
      self.cj1i1[0,0] += self.cj_1i_1_coeff_ij[0,0]
    if self.BC_N == 'Mirror' and self.BC_E == 'Mirror':
      self.cj_1i1[0,-1] += self.cj1i_1_coeff_ij[0,-1]
    if self.BC_S == 'Mirror' and self.BC_W == 'Mirror':
      self.cj1i_1[-1,0] += self.cj_1i1_coeff_ij[-1,0]
    if self.BC_S == 'Mirror' and self.BC_E == 'Mirror':
      self.cj_1i_1[-1,-1] += self.cj1i1_coeff_ij[-1,-1]
    # Periodic

    """
    # Template: 1 set
    self.cj_2i0[:,j] += 
    self.cj_1i_1[:,j] += 
    self.cj_1i0[:,j] += 
    self.cj_1i1[:,j] += 
    self.cj0i_2[:,j] += 
    self.cj0i_1[:,j] += 
    self.cj0i0[:,j] += 
    self.cj0i1[:,j] += 
    self.cj0i2[:,j] += 
    self.cj1i_1[:,j] += 
    self.cj1i0[:,j] += 
    self.cj1i1[:,j] += 
    self.cj2i0[:,j] += 

    self.cj_2i0[i,:] += 
    self.cj_1i_1[i,:] += 
    self.cj_1i0[i,:] += 
    self.cj_1i1[i,:] += 
    self.cj0i_2[i,:] += 
    self.cj0i_1[i,:] += 
    self.cj0i0[i,:] += 
    self.cj0i1[i,:] += 
    self.cj0i2[i,:] += 
    self.cj1i_1[i,:] += 
    self.cj1i0[i,:] += 
    self.cj1i1[i,:] += 
    self.cj2i0[i,:] += 

    # Template: All sets
      if self.BC_W == 
        j = -1
        self.cj_2i0[:,j] += np.nan
        self.cj_1i_1[:,j] += np.nan
        self.cj_1i0[:,j] += np.nan
        self.cj_1i1[:,j] += np.nan
        self.cj0i_2[:,j] += 
        self.cj0i_1[:,j] += 
        self.cj0i0[:,j] += 
        self.cj0i1[:,j] += 
        self.cj0i2[:,j] += 
        self.cj1i_1[:,j] += 
        self.cj1i0[:,j] += 
        self.cj1i1[:,j] += 
        self.cj2i0[:,j] += 
        j = -2
        self.cj_2i0[:,j] += np.nan
        self.cj_1i_1[:,j] += 
        self.cj_1i0[:,j] += 
        self.cj_1i1[:,j] += 
        self.cj0i_2[:,j] += 
        self.cj0i_1[:,j] += 
        self.cj0i0[:,j] += 
        self.cj0i1[:,j] += 
        self.cj0i2[:,j] += 
        self.cj1i_1[:,j] += 
        self.cj1i0[:,j] += 
        self.cj1i1[:,j] += 
        self.cj2i0[:,j] += 
      
      if self.BC_E == 
        j = -1
        self.cj_2i0[:,j] += 
        self.cj_1i_1[:,j] += 
        self.cj_1i0[:,j] += 
        self.cj_1i1[:,j] += 
        self.cj0i_2[:,j] += 
        self.cj0i_1[:,j] += 
        self.cj0i0[:,j] += 
        self.cj0i1[:,j] += 
        self.cj0i2[:,j] += 
        self.cj1i_1[:,j] += np.nan
        self.cj1i0[:,j] += np.nan
        self.cj1i1[:,j] += np.nan
        self.cj2i0[:,j] += np.nan
        j = -2
        self.cj_2i0[:,j] += 
        self.cj_1i_1[:,j] += 
        self.cj_1i0[:,j] += 
        self.cj_1i1[:,j] += 
        self.cj0i_2[:,j] += 
        self.cj0i_1[:,j] += 
        self.cj0i0[:,j] += 
        self.cj0i1[:,j] += 
        self.cj0i2[:,j] += 
        self.cj1i_1[:,j] += 
        self.cj1i0[:,j] += 
        self.cj1i1[:,j] += 
        self.cj2i0[:,j] += np.nan

      if self.BC_N == 
        j = -2
        self.cj_2i0[i,:] += np.nan
        self.cj_1i_1[i,:] += 
        self.cj_1i0[i,:] += 
        self.cj_1i1[i,:] += 
        self.cj0i_2[i,:] += 
        self.cj0i_1[i,:] += 
        self.cj0i0[i,:] += 
        self.cj0i1[i,:] += 
        self.cj0i2[i,:] += 
        self.cj1i_1[i,:] += 
        self.cj1i0[i,:] += 
        self.cj1i1[i,:] += 
        self.cj2i0[i,:] += 
        j = -1
        self.cj_2i0[i,:] += np.nan
        self.cj_1i_1[i,:] += np.nan
        self.cj_1i0[i,:] += np.nan
        self.cj_1i1[i,:] += np.nan
        self.cj0i_2[i,:] += 
        self.cj0i_1[i,:] += 
        self.cj0i0[i,:] += 
        self.cj0i1[i,:] += 
        self.cj0i2[i,:] += 
        self.cj1i_1[i,:] += 
        self.cj1i0[i,:] += 
        self.cj1i1[i,:] += 
        self.cj2i0[i,:] += 

      if self.BC_S == 
        i = -2
        self.cj_2i0[i,:] += 
        self.cj_1i_1[i,:] += 
        self.cj_1i0[i,:] += 
        self.cj_1i1[i,:] += 
        self.cj0i_2[i,:] += 
        self.cj0i_1[i,:] += 
        self.cj0i0[i,:] += 
        self.cj0i1[i,:] += 
        self.cj0i2[i,:] += 
        self.cj1i_1[i,:] += 
        self.cj1i0[i,:] += 
        self.cj1i1[i,:] += 
        self.cj2i0[i,:] += np.nan
        i = -1
        self.cj_2i0[i,:] += 
        self.cj_1i_1[i,:] += 
        self.cj_1i0[i,:] += 
        self.cj_1i1[i,:] += 
        self.cj0i_2[i,:] += 
        self.cj0i_1[i,:] += 
        self.cj0i0[i,:] += 
        self.cj0i1[i,:] += 
        self.cj0i2[i,:] += 
        self.cj1i_1[i,:] += np.nan
        self.cj1i0[i,:] += np.nan
        self.cj1i1[i,:] += np.nan
        self.cj2i0[i,:] += np.nan
    """

  def build_diags(self):

    ##########################################################
    # INCORPORATE BOUNDARY CONDITIONS INTO COEFFICIENT ARRAY #
    ##########################################################

    # Roll to keep the proper coefficients at the proper places in the
    # arrays: Python will naturally just do vertical shifts instead of 
    # diagonal shifts, so this takes into account the horizontal compoent 
    # to ensure that boundary values are at the right place.
    
    # Roll x
    self.cj_2i0 = np.roll(self.cj_2i0, -2, 1)
    self.cj_1i0 = np.roll(self.cj_1i0, -1, 1)
    self.cj1i0 = np.roll(self.cj1i0, 1, 1)
    self.cj2i0 = np.roll(self.cj2i0, 2, 1)
    # Roll y
    self.cj0i_2 = np.roll(self.cj0i_2, -2, 0)
    self.cj0i_1 = np.roll(self.cj0i_1, -1, 0)
    self.cj0i1 = np.roll(self.cj0i1, 1, 0)
    self.cj0i2 = np.roll(self.cj0i2, 2, 0)
    # Roll x and y
    self.cj_1i_1 = np.roll(self.cj_1i_1, -1, 1)
    self.cj_1i_1 = np.roll(self.cj_1i_1, -1, 0)
    self.cj_1i1 = np.roll(self.cj_1i1, -1, 1)
    self.cj_1i1 = np.roll(self.cj_1i1, 1, 0)
    self.cj1i_1 = np.roll(self.cj1i_1, 1, 1)
    self.cj1i_1 = np.roll(self.cj1i_1, -1, 0)
    self.cj1i1 = np.roll(self.cj1i1, 1, 1)
    self.cj1i1 = np.roll(self.cj1i1, 1, 0)

    coeff_array_list = [self.cj_2i0, self.cj_1i0, self.cj1i0, self.cj2i0, self.cj0i_2, self.cj0i_1, self.cj0i1, self.cj0i2, self.cj_1i_1, self.cj_1i1, self.cj1i_1, self.cj1i1, self.cj0i0]
    for array in coeff_array_list:
      array[np.isinf(array)] = 0
      #array[np.isnan(array)] = 0 # had been used for testing
    
    # Reshape to put in solver
    vec_cj_2i0 = np.reshape(self.cj_2i0, -1, order='C')
    vec_cj_1i_1 = np.reshape(self.cj_1i_1, -1, order='C')
    vec_cj_1i0 = np.reshape(self.cj_1i0, -1, order='C')
    vec_cj_1i1 = np.reshape(self.cj_1i1, -1, order='C')
    vec_cj0i_2 = np.reshape(self.cj0i_2, -1, order='C')
    vec_cj0i_1 = np.reshape(self.cj0i_1, -1, order='C')
    vec_cj0i0 = np.reshape(self.cj0i0, -1, order='C')
    vec_cj0i1 = np.reshape(self.cj0i1, -1, order='C')
    vec_cj0i2 = np.reshape(self.cj0i2, -1, order='C')
    vec_cj1i_1 = np.reshape(self.cj1i_1, -1, order='C')
    vec_cj1i0 = np.reshape(self.cj1i0, -1, order='C')
    vec_cj1i1 = np.reshape(self.cj1i1, -1, order='C')
    vec_cj2i0 = np.reshape(self.cj2i0, -1, order='C')
    
    # Changed this 6 Nov. 2014 in betahaus Berlin to be x-based
    Up2 = vec_cj0i2
    Up1 = np.vstack(( vec_cj_1i1, vec_cj0i1, vec_cj1i1 ))
    Mid = np.vstack(( vec_cj_2i0, vec_cj_1i0, vec_cj0i0, vec_cj1i0, vec_cj2i0 ))
    Dn1 = np.vstack(( vec_cj_1i_1, vec_cj0i_1, vec_cj1i_1 ))
    Dn2 = vec_cj0i_2
    
    # Arrange in solver
    diags = np.vstack(( Dn2, \
                        Dn1, \
                        Mid, \
                        Up1, \
                        Up2 ))
                        
    # Number of rows and columns for array size and offsets
    self.ny = self.nrowsy
    self.nx = self.ncolsx
                        
    self.coeff_matrix = scipy.sparse.spdiags(diags, [-2*self.nx, -self.nx-1, -self.nx, -self.nx+1, -2, -1, 0, 1, 2, self.nx-1, self.nx, self.nx+1, 2*self.nx], self.ny*self.nx, self.ny*self.nx, format='csr') # create banded sparse matrix

    #self.coeff_matrix = scipy.sparse.spdiags(np.vstack((Up1, Mid)), [-self.nx-1, -self.nx,  -self.nx+1, -2, -1, 0, 1, 2], self.ny*self.nx, self.ny*self.nx, format='csr') # create banded sparse matrix

  def calc_max_flexural_wavelength(self):
    """
    Returns the approximate maximum flexural wavelength
    This is important when padding of the grid is required: in Flexure (this 
    code), grids are padded out to one maximum flexural wavelength, but in any 
    case, the flexural wavelength is a good characteristic distance for any 
    truncation limit
    """
    if np.isscalar(self.D):
      Dmax = self.D
    else:
      Dmax = self.D.max()
    # This is an approximation if there is fill that evolves with iterations 
    # (e.g., water), but should be good enough that this won't do much to it
    alpha = (4*Dmax/(self.drho*self.g))**.25 # 2D flexural parameter
    self.maxFlexuralWavelength = 2*np.pi*alpha
    self.maxFlexuralWavelength_ncells_x = int(np.ceil(self.maxFlexuralWavelength / self.dx))
    self.maxFlexuralWavelength_ncells_y = int(np.ceil(self.maxFlexuralWavelength / self.dy))
    
  def fd_solve(self):
    """
    w = fd_solve()
    Sparse flexural response calculation.
    Can be performed by direct factorization with UMFpack (defuault) 
    or by an iterative minimum residual technique
    These are both the fastest of the standard Scipy builtin techniques in 
    their respective classes 
    Requires the coefficient matrix from "2D.coeff_matrix"
    """
    
    self.solver_start_time = time.time()

    if self.Debug:
      try:
        # Will fail if scalar
        print 'self.Te', self.Te.shape
      except:
        pass
      print 'self.q0', self.q0.shape
      print 'maxFlexuralWavelength_ncells: (x, y):', self.maxFlexuralWavelength_ncells_x, self.maxFlexuralWavelength_ncells_y
    
    q0vector = self.q0.reshape(-1, order='C')
    if self.solver == "iterative" or self.solver == "Iterative":
      if self.Debug:
        print "Using generalized minimal residual method for iterative solution"
      if self.Verbose:
        print "Converging to a tolerance of", self.iterative_ConvergenceTolerance, "m between iterations"
      wvector = scipy.sparse.linalg.isolve.lgmres(self.coeff_matrix, q0vector)#, tol=1E-10)#,x0=woldvector)#,x0=wvector,tol=1E-15)    
      wvector = wvector[0] # Reach into tuple to get my array back
    else:
      if self.solver == "direct" or self.solver == "Direct":
        if self.Debug:
          print "Using direct solution with UMFpack"
      else:
        print "Solution type not understood:"
        print "Defaulting to direct solution with UMFpack"
      wvector = scipy.sparse.linalg.spsolve(self.coeff_matrix, q0vector, use_umfpack=True)

    # Reshape into grid
    self.w = -wvector.reshape(self.q0.shape)
    self.w_padded = self.w.copy() # for troubleshooting

    self.time_to_solve = time.time() - self.solver_start_time
    print 'Time to solve [s]:', self.time_to_solve

