from __future__ import division # No automatic floor division
from base import *
from scipy.sparse.linalg import spsolve
from scipy import sparse
from scipy.special import kei

# class F2D inherits Flexure and overrides __init__ therefore setting up the same
# three parameters as class Isostasy; and it then sets up more parameters specific
# to its own type of simulation.
class F2D(Flexure):
  def initialize(self, filename):
    super(F2D, self).initialize(filename)
    if debug: print 'F2D initialized'

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
      print 'Error: method must be "FD", "FFT", or "SPA"'
      self.abort()

    if debug: print 'F2D run'
    self.method_func()
    #self.imshow(self.w) # debugging

  def finalize(self):
    if debug: print 'F2D finalized'
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
    print "The fast Fourier transform solution method is not yet implemented."
    sys.exit()

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

    # don't want to keep saying "self." everwhere!
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
      
    elif type(self.Te) == np.ndarray:
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
      Dy  = (-D0_1 - D01)/2.
      Dxx = (D_10 - 2.*D00 + D10)
      Dyy = (D0_1 - 2.*D00 + D01)
      Dxy = (D_1_1 - D_11 - D1_1 + D11)/4.
      """
      # NEW STENCIL
      # x = -2, y = 0
      self.cj_2i0 = (D0 - Dx) / dx4
      # x = 0, y = -2
      self.cj0i_2 = (D0 - Dy) / dy4
      # x = 0, y = 2
      self.cj0i2 = (D0 + Dy) / dy4
      # x = 2, y = 0
      self.cj2i0 = (D0 + Dx) / dx4
      # x = -1, y = -1
      self.cj_1i_1 = (2.*D0 - Dx - Dy + Dxy*(1-nu)/2.) / dx2dy2
      # x = -1, y = 1
      self.cj1i_1 = (2.*D0 - Dx + Dy - Dxy*(1-nu)/2.) / dx2dy2
      # x = 1, y = -1
      self.cj_1i1 = (2.*D0 + Dx - Dy - Dxy*(1-nu)/2.) / dx2dy2
      # x = 1, y = 1
      self.cj1i1 = (2.*D0 + Dx + Dy + Dxy*(1-nu)/2.) / dx2dy2
      # x = -1, y = 0
      self.cj0i_1 = (-4.*D0 + 2.*Dx + Dxx)/dx4 + (-4.*D0 + 2.*Dx + nu*Dyy)/dx2dy2
      # x = 0, y = -1
      self.cj_1i0 = (-4.*D0 + 2.*Dy + Dyy)/dy4 + (-4.*D0 + 2.*Dy + nu*Dxx)/dx2dy2
      # x = 0, y = 1
      self.cj1i0 = (-4.*D0 - 2.*Dy + Dyy)/dy4 + (-4.*D0 - 2.*Dy + nu*Dxx)/dx2dy2
      # x = 1, y = 0
      self.cj0i1 = (-4.*D0 - 2.*Dx + Dxx)/dx4 + (-4.*D0 - 2.*Dx + nu*Dyy)/dx2dy2
      # x = 0, y = 0
      self.cj0i0 = (6.*D0 - 2.*Dxx)/dx4 \
                   + (6.*D0 - 2.*Dyy)/dy4 \
                   + (8.*D0 - 2.*nu*Dxx - 2.*nu*Dyy)/dx2dy2 \
                   + drho*g
      """
      """
      # SIMPLER STENCIL: just del**2(D del**2(w)): only linear variations
      # in Te are allowed.
      # x = -2, y = 0
      self.cj_2i0 = D0 / dx4
      # x = 0, y = -2
      self.cj0i_2 = D0 / dy4
      # x = 0, y = 2
      self.cj0i2 = D0 / dy4
      # x = 2, y = 0
      self.cj2i0 = D0 / dx4
      # x = -1, y = -1
      self.cj_1i_1 = 2.*D0 / dx2dy2
      # x = -1, y = 1
      self.cj1i_1 = 2.*D0 / dx2dy2
      # x = 1, y = -1
      self.cj_1i1 = 2.*D0 / dx2dy2
      # x = 1, y = 1
      self.cj1i1 = 2.*D0 / dx2dy2
      # x = -1, y = 0
      self.cj0i_1 = (-4.*D0 + Dxx)/dx4 + (-4.*D0 + Dyy)/dx2dy2
      # x = 0, y = -1
      self.cj_1i0 = (-4.*D0 + Dyy)/dx4 + (-4.*D0 + Dxx)/dx2dy2
      # x = 0, y = 1
      self.cj1i0 = (-4.*D0 + Dyy)/dx4 + (-4.*D0 + Dxx)/dx2dy2
      # x = 1, y = 0
      self.cj0i1 = (-4.*D0 + Dxx)/dx4 + (-4.*D0 + Dyy)/dx2dy2
      # x = 0, y = 0
      self.cj0i0 = (6.*D0 - 2.*Dxx)/dx4 \
                   + (6.*D0 - 2.*Dyy)/dy4 \
                   + (8.*D0 - 2.*Dxx - 2.*Dyy)/dx2dy2 \
                   + drho*g
      """
      #"""
      # STENCIL FROM GOVERS ET AL. 2009 -- first-order differences
      # x is j and y is i b/c matrix row/column notation
      # x = -2, y = 0
      self.cj_2i0 = D_10/dx4
      # x = -1, y = -1
      self.cj_1i_1 = (D_10 + D0_1)/dx2dy2
      # x = -1, y = 0
      self.cj_1i0 = -2. * ( (D0_1 + D00)/dx2dy2 + (D00 + D_10)/dx4 )
      # x = -1, y = 1
      self.cj_1i1 = (D_10 + D01)/dx2dy2
      # x = 0, y = -2
      self.cj0i_2 = D0_1/dy4
      # x = 0, y = -1
      self.cj0i_1 = -2. * ( (D0_1 + D00)/dx2dy2 + (D00 + D0_1)/dy4)
      # x = 0, y = 0
      self.cj0i0 = (D10 + 4.*D00 + D_10)/dx4 + (D01 + 4.*D00 + D0_1)/dy4 + (8.*D00/dx2dy2) + drho*g
      # x = 0, y = 1
      self.cj0i1 = -2. * ( (D01 + D00)/dy4 + (D00 + D01)/dx2dy2 )
      # x = 0, y = 2
      self.cj0i2 = D0_1/dy4
      # x = 1, y = -1
      self.cj1i_1 = (D10+D0_1)/dx2dy2
      # x = 1, y = 0
      self.cj1i0 = -2. * ( (D10 + D00)/dx4 + (D10 + D00)/dx2dy2 )
      # x = 1, y = 1
      self.cj1i1 = (D10 + D01)/dx2dy2
      # x = 2, y = 0
      self.cj2i0 = D10/dx4
      #"""
      """
      # OLD STENCIL -- think I made some mistakes in the discretization
      # (and maybe even in solving the equation!!!)
      self.cj2i0 = (D00 + 0.5*(D01 - D0_1))/dy4
      self.cj1i_1 = (2*D00 + 0.5*(D10 - D_10 + D01 - D0_1) \
        + ((1-nu)/8) * (D11 - D1_1 - D_11 + D_1_1)) / dx2dy2
      self.cj1i0 = (-6*D00 + 2*D0_1)/dy4 \
        + nu*(D10 - 2*D00 + D_10)/dx2dy2 \
        + (-D01 - 4*D00 + D0_1 - D10 + D_10)/dx2dy2
      self.cj1i1 = (2*D00 + 0.5*(D10 - D_10 + D01 - D0_1) \
        + ((1-nu)/8)*(D11 - D1_1 - D_11 + D_1_1)) \
        /dx2dy2
      self.cj0i_2 = (D00 - 0.5*(D10 - D_10)) / dx4
      self.cj0i_1 = (2*D10 - 6*D00) / dx4 \
        + nu*(D01 - 2*D00 + D0_1) / dx2dy2 \
        - 4*D00 / dx2dy2
      self.cj0i0 = (-2*D10 + 10*D00 - 2*D_10) / dx4 \
        + (-2*D01 + 10*D00 - 2*D0_1) / dy4 \
        + (8*D00) / dx2dy2 \
        + nu*(-2*D10 - 2*D_10 + 8*D00 - 2*D01 -2*D0_1) / dx2dy2 \
        + drho*g
      self.cj0i1 = (-6*D00 + 2*D_10) / dx4 \
        + nu*(D01 - 2*D00 + D0_1) / dx2dy2 \
        - 4*D00 / dx2dy2
      self.cj0i2 = (D00 + 0.5*(D10 - D_10)) / dx4
      self.cj_1i_1 = (2*D00 \
        + 0.5*(-D10 + D_10 - D01 + D0_1) \
        + ((1-nu)/8)*(D11 - D1_1 - D_11 + D_1_1)) \
        / dx2dy2
      self.cj_1i0 = (2*D01 - 6*D00) / dy4 \
        + nu*(D10 - 2*D00 + D_10) / dx2dy2 \
        + (D10 + D01 - 4*D00 - D0_1 - D_10) / dx2dy2
      self.cj_1i1 = (2*D00 + 0.5*(-D10 + D_10 - D01 + D0_1) \
        - ((1-nu)/8) * (D11 - D1_1 - D_11 + D_1_1)) / dx2dy2
      self.cj_2i0 = (D00 - 0.5*(D01 - D0_1)) / dy4
      """
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

    ###################################
    # PADDING CAN BE REQUIRED: MIRROR #
    ###################################

    # Go to this step if any of the boundary conditions require padding
    if self.BC_W == 'Mirror' or self.BC_E == 'Mirror' or \
      self.BC_N == 'Mirror' or self.BC_S == 'Mirror':
      # Define an approximate maximum flexural wavelength to obtain
      # required distances to pad the array    
      self.calc_max_flexural_wavelength()
      # Boundaries require padding! Do it!
      self.BCs_that_need_padding()

      # Need to hold on to q0 array to keep its size; self.q0 will be redefined 
      # for the computation
      self.q0_orig = self.q0.copy()
      
      # Then copy over q0
      self.q0 = self.q0pad.copy()

      # Now I can build the coefficient arrays
      self.get_coeff_values()
      
      # Once the array is its definitive size, can start applying real, solid
      # boundary conditions
      # However, because these are related in a non-simple way to the padding 
      # boundary conditions, need to pick these separately
      self.padded_edges_BCs() # Boundaries outside of the padding chosen here

    ####################################################
    # PERIODIC OR DIRICHLET: DECISION MADE IN FUNCTION #
    ####################################################
    else:
      self.get_coeff_values()
      #for i in range(self.nrowsy):
        #self.build_coeff_matrix_nonzero_blocks_1row(i)
        #self.assemble_blocks_sparse_1row(i)
      self.build_diags()

  def BCs_that_need_padding(self):
    """
    This function acts as a main interface for BC_Mirror
    
    It is needed because these functions pad the array, and if one pads it on 
    one side and the other pads it on the other, the final array isn't known 
    until both boundary conditions are evaluated. Because these padding 
    boundary conditions are evaluated outside of the static boundary b.c.'s, 
    it is necessary to have them combined here (instead of above)
    """

    # self.q0 not touched until later; the two functions called from this 
    # one modify self.q0pad
    # This is especially important in 2D because the N-S functions are 
    # dependent on the modified q0pad after the E-W functions are run
    # (if the E-W functions are run)
    self.q0pad = self.q0.copy() # Prep for concatenation

    # Build the proper boundary conditions, E-W first, to flip or pad LR first 
    # (if desired) then flip (or pad) the whole thing N-S
    
    # E-W
    if self.BC_W == 'Mirror' or self.BC_E == 'Mirror':
      print "MIRROR, E-W!"
      self.BC_Mirror_EW()

    # N-S
    if self.BC_N == 'Mirror' or self.BC_S == 'Mirror':
      print "MIRROR, N-S!"
      self.BC_Mirror_NS()

    # Pad Te, if it is an array, to match q0
    self.pad_Te()

  def BC_Mirror_EW(self):
    """
    Mirrors q0 across the boundary on either the west (left) or east (right) 
    side, depending on the selections.
    
    This can, for example, produce a scenario in which you are observing 
    a mountain range up to the range crest (or, more correctly, the halfway 
    point acros the mountain range).
    
    The mirror is run out to one flexural wavelength away from the main 
    part of the grid, after which it is clipped (if longer) or padded with 
    additional zeros (if not).
    
    This has similar rules to the no outside loads condition: if both sides 
    are mirrored, one will be mirrored and a periodic boundary condition will 
    be applied.
    
    BC_Mirror_EW and BC_Mirror_NS are almost identical, but are different 
    enough that I've found it easier to separate them for the time being
    """
        
    # Before starting, make sure that other side isn't "periodic": in that 
    # case, change it to Mirror because the solution is OFTEN the same and 
    # changing it should improve flow control (right-side padding desired 
    # in this case even if the right boundary was periodic)
    # Often but not always b/c if the domain is too short, it won't see its 
    # loads more than once with mirror, which goes against periodic (i.e. 
    # my numerical mirror doesn't work recursively like a real one would)
    if self.BC_E == 'Periodic' or self.BC_W == 'Periodic':
      print "Setting one periodic boundary in conjunction with one Mirror"
      print "boundary is OFTEN but not ALWAYS the same as setting both sides"
      print "to have no outside loads; this is because a periodic boundary"
      print "condition is often used to halve the distance that needs to be"
      print "padded."
      print "DEFAULTING TO BOTH MIRROR BOUNDARIES IN THIS SITUATION!"
      # One is already Mirror, so change both to make sure that we 
      # get both of them
      self.BC_E = 'Mirror'
      self.BC_W = 'Mirror'
    
    # First, see how far to pad: 1 flexural wavelength
    # This is done outside this function
    
    # Second, create the E-W mirrored load grid
    self.q0_mirror_EW = np.fliplr(self.q0pad)
    
    # Third, make padding array (if needed)
    # If doing both sides, just repeat whole grid (if > flexural wavelength 
    # and more efficient than just repeating part of it on both sides)
    if self.q0_mirror_EW.shape[1] < self.maxFlexuralWavelength_ncells_x:
      zeropad = np.zeros( (self.q0_mirror_EW.shape[0], \
        self.maxFlexuralWavelength_ncells_x -  self.q0_mirror_EW.shape[1]) )

    # Fourth, find what may need to be added to each side
    if self.BC_E == 'Mirror':
      if self.q0_mirror_EW.shape[1] < self.maxFlexuralWavelength_ncells_x:
        self.q0_mirror_EW_E = np.concatenate((self.q0_mirror_EW,zeropad), axis=1)
      elif self.q0_mirror_EW.shape[1] > self.maxFlexuralWavelength_ncells_x:
        self.q0_mirror_EW_E = self.q0_mirror_EW[:,:self.maxFlexuralWavelength_ncells_x]
    if self.BC_W == 'Mirror':
      if self.q0_mirror_EW.shape[1] < self.maxFlexuralWavelength_ncells_x:
        self.q0_mirror_EW_W = np.concatenate((zeropad,self.q0_mirror_EW), axis=1)
      elif self.q0_mirror_EW.shape[1] > self.maxFlexuralWavelength_ncells_x:
        self.q0_mirror_EW_W = self.q0_mirror_EW[:,-self.maxFlexuralWavelength_ncells_x:]

    # Fifth, add things properly to each side
    # Starting with both sides being mirror bc's
    if self.BC_E == 'Mirror' and self.BC_W == 'Mirror':
      # Case 1: glom onto both sides because it is too short or too long
      if len(self.q0_mirror_EW) < self.maxFlexuralWavelength_ncells_x \
        or len(self.q0_mirror_EW) > 2*self.maxFlexuralWavelength_ncells_x:
        self.q0pad = np.concatenate((self.q0_mirror_EW_W,self.q0pad,self.q0_mirror_EW_E), axis=1)
      # Case 2: Add to one side and later use a periodic boundary condition  
      # because it is just right and these are more efficient
      else:
        self.q0pad = np.concatenate((self.q0pad,self.q0_mirror_EW), axis=1)        
    # And then if just one side or the other is mirror:
    elif self.BC_E == 'Mirror':
      self.q0pad = np.concatenate((self.q0pad,self.q0_mirror_EW_E), axis=1)
    elif self.BC_W == 'Mirror':
      self.q0pad = np.concatenate((self.q0_mirror_EW_W,self.q0pad), axis=1)
    
  def BC_Mirror_NS(self):
    """
    Mirrors q0 across the boundary on either the north (top) or south (bottom) 
    side, depending on the selections.
    
    This can, for example, produce a scenario in which you are observing 
    a mountain range up to the range crest (or, more correctly, the halfway 
    point acros the mountain range).
    
    The mirror is run out to one flexural wavelength away from the main 
    part of the grid, after which it is clipped (if longer) or padded with 
    additional zeros (if not).
    
    This has similar rules to the no outside loads condition: if both sides 
    are mirrored, one will be mirrored and a periodic boundary condition will 
    be applied.
    
    BC_Mirror_EW and BC_Mirror_NS are almost identical, but are different 
    enough that I've found it easier to separate them for the time being
    """
    
    # Before starting, make sure that other side isn't "periodic": in that 
    # case, change it to Mirror because the solution is OFTEN the same and 
    # changing it should improve flow control (right-side padding desired 
    # in this case even if the right boundary was periodic)
    # Often but not always b/c if the domain is too short, it won't see its 
    # loads more than once with mirror, which goes against periodic (i.e. 
    # my numerical mirror doesn't work recursively like a real one would)
    if self.BC_N == 'Periodic' or self.BC_S == 'Periodic':
      print "Setting one periodic boundary in conjunction with one Mirror"
      print "boundary is OFTEN but not ALWAYS the same as setting both sides"
      print "to have no outside loads; this is because a periodic boundary"
      print "condition is often used to halve the distance that needs to be"
      print "padded."
      print "DEFAULTING TO BOTH MIRROR BOUNDARIES IN THIS SITUATION!"
      # One is already Mirror, so change both to make sure that we 
      # get both of them
      self.BC_N = 'Mirror'
      self.BC_S = 'Mirror'
    
    # First, see how far to pad: 1 flexural wavelength
    # This is done outside this function
    
    # Second, create the N-S mirrored load grid
    # Important to use q0pad here in case the EW padding has been done
    self.q0_mirror_NS = np.flipud(self.q0pad)
    
    # Third, make padding array (if needed)
    # If doing both sides, just repeat whole grid (if > flexural wavelength 
    # and more efficient than just repeating part of it on both sides)
    if self.q0_mirror_NS.shape[0] < self.maxFlexuralWavelength_ncells_y:
      zeropad = np.zeros( (self.maxFlexuralWavelength_ncells_y - \
      self.q0_mirror_NS.shape[0], self.q0_mirror_NS.shape[1],) )

    # Fourth, find what may need to be added to each side
    if self.BC_S == 'Mirror':
      if self.q0_mirror_NS.shape[0] < self.maxFlexuralWavelength_ncells_y:
        self.q0_mirror_NS_S = np.concatenate((self.q0_mirror_NS,zeropad), axis=0)
      elif self.q0_mirror_NS.shape[0] >= self.maxFlexuralWavelength_ncells_y:
        self.q0_mirror_NS_S = self.q0_mirror_NS[:self.maxFlexuralWavelength_ncells_y,:]
    if self.BC_N == 'Mirror':
      if self.q0_mirror_NS.shape[0] < self.maxFlexuralWavelength_ncells_y:
        self.q0_mirror_NS_N = np.concatenate((zeropad,self.q0_mirror_NS), axis=0)
      elif self.q0_mirror_NS.shape[0] >= self.maxFlexuralWavelength_ncells_y:
        self.q0_mirror_NS_N = self.q0_mirror_NS[-self.maxFlexuralWavelength_ncells_y:,:]

    # Fifth, add things properly to each side
    # Starting with both sides being mirror bc's
    if self.BC_S == 'Mirror' and self.BC_N == 'Mirror':
      # Case 1: glom onto both top and bottom because it is too short or too long
      if len(self.q0_mirror_NS) < self.maxFlexuralWavelength_ncells_y \
        or len(self.q0_mirror_NS) > 2*self.maxFlexuralWavelength_ncells_y:
        self.q0pad = np.concatenate((self.q0_mirror_NS_N,self.q0pad,self.q0_mirror_NS_S), axis=0)
      # Case 2: Add to the bottom and later use a periodic boundary condition  
      # because it is just right and these are more efficient
      else:
        self.q0pad = np.concatenate((self.q0pad,self.q0_mirror_NS), axis=0)
    # And then if just one side or the other is mirror:
    elif self.BC_S == 'Mirror':
      self.q0pad = np.concatenate((self.q0pad,self.q0_mirror_NS_S), axis=0)
    elif self.BC_N == 'Mirror':
      self.q0pad = np.concatenate((self.q0_mirror_NS_N,self.q0pad), axis=0)
  
  def pad_Te(self):
    """
    Pad elastic thickness to match padded q0 array.
    This is needed for the Mirror boundary condition
    
    Mirror boundary conditions mirror the elastic thickness array out to the 
    desired distance.

    This function will do nothing if elastic thickness is a scalar... because 
    then it is not a grid that needs to be padded. It will also do nothing if 
    the boundary conditions do not include a "Mirror".

    Use linspace to keep value constant on both sides of the padding seam
    And also update D based on the extended Te array
    """

    if type(self.Te) == np.ndarray:
      self.Te_orig = self.Te.copy() # Save original Te
      
      ###################
      # EAST-WEST FIRST #
      ###################
  
      # Get elastic thicknesses on margins for Mirror b.c.
      if self.BC_W == 'Mirror' or self.BC_E == 'Mirror':
        TeW = self.Te[:,0]
        TeE = self.Te[:,-1]
      # Mirror on both sides
      if self.BC_W == 'Mirror' and self.BC_E == 'Mirror':
        # Case 1: padded on both sides 
        if self.Te.shape[1] < self.maxFlexuralWavelength_ncells_x:
          # Thought I would need extra padding because there are 2 extra cells 
          # in the Te grid when compared to the q0 grid
          # But I guess not - must be taken care of elsewhere - also, it is 
          # very late and I am not thinking so well... :)
          extrapad_length = self.maxFlexuralWavelength_ncells_x - self.Te.shape[1]
          extrapad_W = np.lib.stride_tricks.as_strided(TeW, \
            (extrapad_length, TeW.size), (0, TeW.itemsize)).transpose()
          extrapad_E = np.lib.stride_tricks.as_strided(TeE, \
            (extrapad_length, TeE.size), (0, TeE.itemsize)).transpose()
          padTeW = np.concatenate(( extrapad_W, np.fliplr(self.Te) ), axis=1)
          padTeE = np.concatenate(( np.fliplr(self.Te), extrapad_E ), axis=1)
          self.Te = np.concatenate((padTeW, self.Te, padTeE), axis=1)
        elif self.Te.shape[1] > 2*self.maxFlexuralWavelength_ncells_x:
          padTeW = np.fliplr(self.Te)[:,-self.maxFlexuralWavelength_ncells_x:]
          padTeE = np.fliplr(self.Te)[:,:self.maxFlexuralWavelength_ncells_x]
          self.Te = np.concatenate((padTeW, self.Te, padTeE), axis=1)
        # Case 2: Padded on right
        else:
          # Can't go the whole length because Te_orig is padded
          # SO THIS IS DIFFERENT FROM THE Q0, TE MIRRORING OF THE ENDPOINT 
          # BEFORE!
          padTeE = np.fliplr(self.Te)[:,1:-1]
          self.Te = np.concatenate((self.Te,padTeE), axis=1)
      elif self.BC_W == 'Mirror':
        if self.Te.shape[1] < self.maxFlexuralWavelength_ncells_x:
          extrapad_length = self.maxFlexuralWavelength_ncells_x - self.Te.shape[1]
          extrapad_W = np.lib.stride_tricks.as_strided(TeW, \
            (extrapad_length, TeW.size), (0, TeW.itemsize)).transpose()
          padTeW = np.concatenate(( extrapad_W, np.fliplr(self.Te) ), axis=1)
        else:
          padTeW = np.fliplr(self.Te)[:,-self.maxFlexuralWavelength_ncells_x:]        
        self.Te = np.concatenate((padTeW,self.Te), axis=1)
      elif self.BC_E == 'Mirror':
        if self.Te.shape[1] < self.maxFlexuralWavelength_ncells_x:
          extrapad_length = self.maxFlexuralWavelength_ncells_x - self.Te.shape[1]
          extrapad_W = np.lib.stride_tricks.as_strided(TeW, \
            (extrapad_length, TeW.size), (0, TeW.itemsize)).transpose()
          padTeE = np.concatenate(( np.fliplr(self.Te), extrapad_E ), axis=1)
        else:
          padTeE = np.fliplr(self.Te)[:,:self.maxFlexuralWavelength_ncells_x]
        self.Te = np.concatenate((self.Te,padTeE), axis=1)

      ####################
      # THEN NORTH-SOUTH #
      ####################
      
      # Get elastic thicknesses on margins for Mirror b.c.
      if self.BC_N == 'Mirror' or self.BC_S == 'Mirror':
        TeN = self.Te[0,:]
        TeS = self.Te[-1,:]
      if self.BC_N == 'Mirror' and self.BC_S == 'Mirror':
        # Case 1: padded on both sides 
        if self.Te.shape[0] < self.maxFlexuralWavelength_ncells_y:
          extrapad_length = self.maxFlexuralWavelength_ncells_y - self.Te.shape[0]
          extrapad_N = np.lib.stride_tricks.as_strided(TeN, \
            (extrapad_length, TeN.size), (0, TeN.itemsize))
          extrapad_S = np.lib.stride_tricks.as_strided(TeS, \
            (extrapad_length, TeS.size), (0, TeS.itemsize))
          padTeN = np.concatenate(( extrapad_N, np.fliplr(self.Te) ), axis=0)
          padTeS = np.concatenate(( np.fliplr(self.Te), extrapad_S ), axis=0)
          self.Te = np.concatenate((padTeN, self.Te, padTeS), axis=0)
        elif self.Te.shape[0] > 2*self.maxFlexuralWavelength_ncells_y:
          padTeN = np.flipud(self.Te)[-self.maxFlexuralWavelength_ncells_y:,:]
          padTeS = np.flipud(self.Te)[:self.maxFlexuralWavelength_ncells_y,:]
          self.Te = np.concatenate((padTeN, self.Te, padTeS), axis=0)
        # Case 2: Padded on right
        else:
          # Can't go the whole length because Te_orig is padded
          # SO THIS IS DIFFERENT FROM THE Q0, TE MIRRORING OF THE ENDPOINT 
          # BEFORE!
          padTeS = np.flipud(self.Te)[1:-1,:]
          self.Te = np.concatenate((self.Te,padTeS), axis=0)
      # Combo and mirror-only cases already accounted for, so can just do a 
      # simple if
      elif self.BC_N == 'Mirror':
        if self.Te.shape[0] < self.maxFlexuralWavelength_ncells_y:
          extrapad_length = self.maxFlexuralWavelength_ncells_y - self.Te.shape[0]
          extrapad_N = np.lib.stride_tricks.as_strided(TeN, \
            (extrapad_length, TeN.size), (0, TeN.itemsize))
          padTeN = np.concatenate(( extrapad_N, np.fliplr(self.Te) ), axis=0)
        else:
          padTeN = np.flipud(self.Te)[-self.maxFlexuralWavelength_ncells_y:,:]        
        self.Te = np.concatenate((padTeN,self.Te), axis=0)
      elif self.BC_S == 'Mirror':
        if self.Te.shape[0] < self.maxFlexuralWavelength_ncells_y:
          extrapad_length = self.maxFlexuralWavelength_ncells_y - self.Te.shape[0]
          extrapad_S = np.lib.stride_tricks.as_strided(TeS, \
            (extrapad_length, TeS.size), (0, TeS.itemsize))
          padTeS = np.concatenate(( np.fliplr(self.Te), extrapad_S ), axis=0)
        else:
          padTeS = np.flipud(self.Te)[:self.maxFlexuralWavelength_ncells_y,:]
        self.Te = np.concatenate((self.Te,padTeS), axis=0)
      
      #####################
      # Finally, update D #
      #####################
      
      self.D = self.E*self.Te**3/(12*(1-self.nu**2)) # Flexural rigidity
      print "MOO!"

  def padded_edges_BCs(self):
    """
    Sets the boundary conditions outside of padded edges; this is important 
    for Mirror boundary conditions

    Also makes an archival copy of q0 for while q0 is temporarily replaced 
    by the padded version of it
    """

    ###################
    # East-West first #
    ###################

    # First check if there is an East-West component to consider
    if self.BC_E == 'Mirror' or self.BC_W == 'Mirror':
      if self.BC_E == 'Mirror' and self.BC_W == 'Mirror':
        if len(self.q0_mirror_EW) > self.maxFlexuralWavelength_ncells_x \
          or len(self.q0_mirror_EW) < 2*self.maxFlexuralWavelength_ncells_x:
          BC_EW = 'Periodic'
        else:
          BC_EW = 'Dirichlet'
      else:
        BC_EW = 'Dirichlet'

    ####################
    # Then North-South #
    ####################
    
    # First check if there is a North-South component to consider
    if self.BC_S == 'Mirror' or self.BC_N == 'Mirror':
      if self.BC_S == 'Mirror' and self.BC_N == 'Mirror':
        if len(self.q0_mirror_NS) > self.maxFlexuralWavelength_ncells_y \
          or len(self.q0_mirror_NS) < 2*self.maxFlexuralWavelength_ncells_y:
          BC_NS = 'Periodic'
        else:
          BC_NS = 'Dirichlet'
      else:
        BC_NS = 'Dirichlet'

    ##########################################################################
    # Then take the boundary conditions that the program has decided upon to #
    # build coefficient arrays                                               #
    ##########################################################################
    
    self.build_diags()
    
  
  def build_diags(self):
    # Assign boundary conditions
    # Dirichlet force right now
    self.cj_2i0[:2,:] = 0
    self.cj_2i0[-2:,:] = 0

    self.cj_1i_1[:1,:] = 0
    self.cj_1i_1[-1:,:] = 0
    self.cj_1i_1[:,:1] = 0
    self.cj_1i_1[:,-1:] = 0
    self.cj_1i0[:1,:] = 0
    self.cj_1i0[-1:,:] = 0
    self.cj_1i1[:1,:] = 0
    self.cj_1i1[-1:,:] = 0
    self.cj_1i1[:,:1] = 0
    self.cj_1i1[:,-1:] = 0

    self.cj0i_2[:,:2] = 0
    self.cj0i_2[:,-2:] = 0
    self.cj0i_1[:,:1] = 0
    self.cj0i_1[:,-1:] = 0
    self.cj0i1[:,:1] = 0
    self.cj0i1[:,-1:] = 0
    self.cj0i2[:,:2] = 0
    self.cj0i2[:,-2:] = 0
    
    self.cj1i_1[:1,:] = 0
    self.cj1i_1[-1:,:] = 0
    self.cj1i_1[:,:1] = 0
    self.cj1i_1[:,-1:] = 0
    self.cj_1i0[:1,:] = 0
    self.cj_1i0[-1:,:] = 0
    self.cj1i1[:1,:] = 0
    self.cj1i1[-1:,:] = 0
    self.cj1i1[:,:1] = 0
    self.cj1i1[:,-1:] = 0

    self.cj2i0[:2,:] = 0
    self.cj2i0[-2:,:] = 0

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
    
    Up2 = vec_cj_2i0
    Up1 = np.vstack(( vec_cj_1i_1, vec_cj_1i0, vec_cj_1i1 ))
    Mid = np.vstack(( vec_cj0i_2, vec_cj0i_1, vec_cj0i0, vec_cj0i1, vec_cj0i2 ))
    Dn1 = np.vstack(( vec_cj1i_1, vec_cj1i0, vec_cj1i1 ))
    Dn2 = vec_cj2i0
    
    # Arrange in solver
                        
    diags = np.vstack(( Up2, \
                        Up1, \
                        Mid, \
                        Dn1, \
                        Dn2 ))
                        
    self.ny = self.nrowsy
    self.nx = self.ncolsx
                        
    import scipy
                        
    self.coeff_matrix = scipy.sparse.spdiags(diags, [-2*self.nx, -self.nx-1, -self.nx,  -self.nx+1, -2, -1, 0, 1, 2, self.nx-1, self.nx, self.nx+1, 2*self.nx], self.ny*self.nx, self.ny*self.nx, format='csr') # create banded sparse matrix

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

    try:
      # Will fail if scalar
      print 'self.Te', self.Te.shape
    except:
      pass
    print 'self.q0', self.q0.shape

    #print 'maxFlexuralWavelength_ncells: (x, y):', self.maxFlexuralWavelength_ncells_x, self.maxFlexuralWavelength_ncells_y
    
    if self.solver == "iterative" or self.solver == "Iterative":
      from scipy.sparse.linalg import isolve
      # Set x0 from direct solver
      #coeff_matrix = sparse.csr_matrix(self.coeff_matrix)
      #q0vector = self.q0.reshape(1,np.prod(self.q0.shape))
      #q0vector = sparse.csr_matrix(q0vector)
      # UMFpack is the default direct solver now in python, 
      # but being explicit here
      #wvector = spsolve(coeff_matrix,q0vector,use_umfpack=True)
      # Now iterative
      q0vector = self.q0.reshape(np.prod(self.q0.shape),1, order='F')
      #woldvector = self.wold.reshape(np.prod(self.wold.shape),1)
      #wvector = isolve.minres(self.coeff_matrix,q0vector)    # fast
      wvector = isolve.lgmres(self.coeff_matrix,q0vector)#,x0=woldvector)#,x0=wvector,tol=1E-15)    
      wvector = wvector[0] # Reach into tuple to get my array back

    else:
      if self.solver == "direct" or self.solver == "Direct":
        pass
      else:
        print "Solution type not understood:"
        print "Defaulting to direct solution with UMFpack"
      # Convert coefficient array format to csr for sparse solver
      #coeff_matrix = sparse.csr_matrix(self.coeff_matrix)
      q0vector = self.q0.reshape(-1)
      #q0vector = sparse.csr_matrix(q0vector)
      # UMFpack is the default direct solver now in python, 
      # but being explicit here
      wvector = spsolve(self.coeff_matrix, q0vector, use_umfpack=True)

    # Reshape into grid
    self.w = -wvector.reshape(self.q0.shape)
    self.w_padded = self.w.copy() # for troubleshooting

    self.time_to_solve = time.time() - self.solver_start_time
    print 'Time to solve [s]:', self.time_to_solve
    
    # If needed, revert q0 and Te to original dimensions
    self.back_to_original_q0_Te_w()

    
  def back_to_original_q0_Te_w(self):
    """
    Pull out the parts of q0, Te that we want for special boundary condition
    cases
    """

    # Special case: clipping needed for the Mirror boundary condition.
    
    # IS THIS SOMETHING I NEED TO WORRY ABOUT? LET'S SEE:
    if self.BC_W == 'Mirror' or self.BC_E == 'Mirror' or \
      self.BC_N == 'Mirror' or self.BC_S == 'Mirror':
      # Then return q0 and (if needed) Te to their original states
      self.q0 = self.q0_orig
      if type(self.Te) == np.ndarray:
        self.Te_padded = self.Te.copy() # For diagnostics
        self.Te = self.Te_orig

    # STARTING WITH EAST-WEST: THIS IS THE SAME AS FOR THE 1D CASE:
    # Mirror
    if self.BC_E == 'Mirror' and self.BC_W == 'Mirror':
      # Case 1: padded on both sides 
      if self.q0_mirror_EW.shape[1] < self.maxFlexuralWavelength_ncells_x \
        or self.q0_mirror_EW.shape[1] > 2*self.maxFlexuralWavelength_ncells_x:
        self.w = self.w[:,self.maxFlexuralWavelength_ncells_x:-self.maxFlexuralWavelength_ncells_x]
      # Case 2: Padded on right
      else:
        self.w = self.w[:,:self.q0_orig.shape[1]]
    # Combo and mirror-only cases already accounted for, so can just do a 
    # simple if
    elif self.BC_W == 'Mirror':
      self.w = self.w[:,self.maxFlexuralWavelength_ncells_x:]
    elif self.BC_E == 'Mirror':
      self.w = self.w[:,:-self.maxFlexuralWavelength_ncells_x]
  
    # NOW FOR NORTH-SOUTH: NOT MUCH DIFFERENT FROM EAST-WEST IN PRINCIPLE:
    # Mirror
    if self.BC_S == 'Mirror' and self.BC_N == 'Mirror':
      # Case 1: padded on both sides 
      if self.q0_mirror_NS.shape[0] < self.maxFlexuralWavelength_ncells_y \
        or self.q0_mirror_NS.shape[0] > 2*self.maxFlexuralWavelength_ncells_y:
        self.w = self.w[self.maxFlexuralWavelength_ncells_y:-self.maxFlexuralWavelength_ncells_y,:]
      # Case 2: Padded on right
      else:
        self.w = self.w[:self.q0_orig.shape[0],:]
    # Combo and mirror-only cases already accounted for, so can just do a 
    # simple if
    elif self.BC_N == 'Mirror':
      self.w = self.w[self.maxFlexuralWavelength_ncells_y:,:]
    elif self.BC_S == 'Mirror':
      self.w = self.w[:-self.maxFlexuralWavelength_ncells_y,:]

