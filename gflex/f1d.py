"""
This file is part of gFlex.
gFlex computes lithospheric flexural isostasy with heterogeneous rigidity
Copyright (C) 2010-2018 Andrew D. Wickert

gFlex is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

gFlex is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with gFlex.  If not, see <http://www.gnu.org/licenses/>.
"""

from __future__ import division, print_function # No automatic floor division
from base import *
from scipy.sparse import spdiags
from scipy.sparse.linalg import spsolve, isolve

class F1D(Flexure):
  def initialize(self, filename=None):
    self.dimension = 1 # Set it here in case it wasn't set for selection before
    super(F1D, self).initialize()
    if self.Verbose: print("F1D initialized")

  def run(self):
    self.bc_check()
    self.solver_start_time = time.time()
    if self.Method == 'FD':
      # Finite difference
      super(F1D, self).FD()
      self.method_func = self.FD
    elif self.Method == 'FFT':
      # Fast Fourier transform
      super(F1D, self).FFT()
      self.method_func = self.FFT
    elif self.Method == "SAS":
      # Superposition of analytical solutions
      super(F1D, self).SAS()
      self.method_func = self.SAS
    elif self.Method == "SAS_NG":
      # Superposition of analytical solutions,
      # nonuniform points
      super(F1D, self).SAS_NG()
      self.method_func = self.SAS_NG
    else:
      sys.exit('Error: method must be "FD", "FFT", "SAS", or "SAS_NG"')

    if self.Verbose: print("F1D run")
    self.method_func()

    self.time_to_solve = time.time() - self.solver_start_time
    if self.Quiet == False:
      print("Time to solve [s]:", self.time_to_solve)

  def finalize(self):
    # If elastic thickness has been padded, return it to its original
    # value, so this is not messed up for repeat operations in a 
    # model-coupling exercise
    try:
      self.Te = self.Te_unpadded
    except:
      pass
    if self.Verbose: print("F1D finalized")
    super(F1D, self).finalize()   
    
  ########################################
  ## FUNCTIONS FOR EACH SOLUTION METHOD ##
  ########################################
  
  def FD(self):
    self.gridded_x()
    # Only generate coefficient matrix if it is not already provided
    if self.coeff_matrix is not None:
      pass
    else:
      self.elasprepFD() # define dx4 and D within self
      self.BC_selector_and_coeff_matrix_creator()
    self.fd_solve() # Get the deflection, "w"

  def FFT(self):
    if self.plotChoice:
      self.gridded_x()
    sys.exit("The fast Fourier transform solution method is not yet implemented.")
    
  def SAS(self):
    self.gridded_x()
    self.spatialDomainVarsSAS()
    self.spatialDomainGridded()

  def SAS_NG(self):
    self.spatialDomainVarsSAS()
    self.spatialDomainNoGrid()

  ######################################
  ## FUNCTIONS TO SOLVE THE EQUATIONS ##
  ######################################


  ## UTILITY
  ############

  def gridded_x(self):
    self.nx = self.qs.shape[0]
    self._x_local = np.arange(0,self.dx*self.nx,self.dx)
    
  
  ## SPATIAL DOMAIN SUPERPOSITION OF ANALYTICAL SOLUTIONS
  #########################################################

  # SETUP

  def spatialDomainVarsSAS(self):
    self.D = self.E*self.Te**3/(12*(1-self.nu**2)) # Flexural rigidity
    self.alpha = (4*self.D/(self.drho*self.g))**.25 # 1D flexural parameter
    self.coeff = self.alpha**3/(8*self.D)

  # UNIFORM DX ("GRIDDED"): LOADS PROVIDED AS AN ARRAY WITH KNOWN DX TO
  # CONVERT LOAD MAGNITUDE AT A POINT INTO MASS INTEGRATED ACROSS DX

  def spatialDomainGridded(self):
  
    self.w = np.zeros(self.nx) # Deflection array
    
    for i in range(self.nx):
      # Loop over locations that have loads, and sum
      if self.qs[i]:
        dist = abs(self._x_local[i]-self._x_local)
        # -= b/c pos load leads to neg (downward) deflection
        self.w -= self.qs[i] * self.coeff * self.dx * np.exp(-dist/self.alpha) * \
          (np.cos(dist/self.alpha) + np.sin(dist/self.alpha))
    # No need to return: w already belongs to "self"
    

  # NONUNIFORM DX (NO GRID): ARBITRARILY-SPACED POINT LOADS
  # So essentially a sum of Green's functions for flexural response

  def spatialDomainNoGrid(self):
    """
    Superposition of analytical solutions without a gridded domain
    """
    self.w = np.zeros(self.xw.shape)

    if self.Debug:
      print("w = ")
      print(self.w.shape)
    
    for i in range(len(self.q)):
      # More efficient if we have created some 0-load points
      # (e.g., for where we want output)
      if self.q[i] != 0:
        dist = np.abs(self.xw - self.x[i])
        self.w -= self.q[i] * self.coeff * np.exp(-dist/self.alpha) * \
          ( np.cos(dist/self.alpha) + np.sin(dist/self.alpha) )

  ## FINITE DIFFERENCE
  ######################
  
  def elasprepFD(self):
    """
    dx4, D = elasprepFD(dx,Te,E=1E11,nu=0.25)
    
    Defines the variables (except for the subset flexural rigidity) that are
    needed to run "coeff_matrix_1d"
    """
    self.dx4 = self.dx**4
    self.D = self.E*self.Te**3/(12*(1-self.nu**2))

  def BC_selector_and_coeff_matrix_creator(self):
    """
    Selects the boundary conditions
    Then calls the function to build the pentadiagonal matrix to solve 
    1D flexure with variable (or constant) elsatic thickness
    """
    
    # Zeroth, start the timer and print the boundary conditions to the screen
    self.coeff_start_time = time.time()
    if self.Verbose:
      print("Boundary condition, West:", self.BC_W, type(self.BC_W))
      print("Boundary condition, East:", self.BC_E, type(self.BC_E))

    # First, set flexural rigidity boundary conditions to flesh out this padded
    # array
    self.BC_Rigidity()
    
    # Second, build the coefficient arrays -- with the rigidity b.c.'s
    self.get_coeff_values()

    # Third, apply boundary conditions to the coeff_arrays to create the 
    # flexural solution
    self.BC_Flexure()
    
    # Fourth, construct the sparse diagonal array
    self.build_diagonals()
    
    # Finally, compute the total time this process took    
    self.coeff_creation_time = time.time() - self.coeff_start_time
    if self.Quiet == False:
      print("Time to construct coefficient (operator) array [s]:", self.coeff_creation_time)

  def BC_Rigidity(self):
    """
    Utility function to help implement boundary conditions by specifying 
    them for and applying them to the elastic thickness grid
    """

    #########################################
    # FLEXURAL RIGIDITY BOUNDARY CONDITIONS #
    #########################################
    # West
    if self.BC_W == 'Periodic':
      self.BC_Rigidity_W = 'periodic'
    elif (self.BC_W == np.array(['0Displacement0Slope', '0Moment0Shear', '0Slope0Shear'])).any():
      self.BC_Rigidity_W = '0 curvature'
    elif self.BC_W == 'Mirror':
      self.BC_Rigidity_W = 'mirror symmetry'
    else:
      sys.exit("Invalid Te B.C. case")
    # East
    if self.BC_E == 'Periodic':
      self.BC_Rigidity_E = 'periodic'
    elif (self.BC_E == np.array(['0Displacement0Slope', '0Moment0Shear', '0Slope0Shear'])).any():
      self.BC_Rigidity_E = '0 curvature'
    elif self.BC_E == 'Mirror':
      self.BC_Rigidity_E = 'mirror symmetry'
    else:
      sys.exit("Invalid Te B.C. case")
    
    #############
    # PAD ARRAY #
    #############
    if np.isscalar(self.Te):
      self.D *= np.ones(self.qs.shape) # And leave Te as a scalar for checks
    else:
      self.Te_unpadded = self.Te.copy()
    # F2D keeps this inside the "else" and handles this differently, 
    # largely because it has different ways of computing the flexural
    # response with variable Te. We'll keep everything simpler here and 
    # just pad this array so it can be sent through the same process
    # to create the coefficient arrays.
    self.D = np.hstack([np.nan, self.D, np.nan])

    ###############################################################
    # APPLY FLEXURAL RIGIDITY BOUNDARY CONDITIONS TO PADDED ARRAY #
    ###############################################################
    if self.BC_Rigidity_W == "0 curvature":
      self.D[0] = 2*self.D[1] - self.D[2]
    if self.BC_Rigidity_E == "0 curvature":
      self.D[-1] = 2*self.D[-2] - self.D[-3]
    if self.BC_Rigidity_W == "mirror symmetry":
      self.D[0] = self.D[2]
    if self.BC_Rigidity_E == "mirror symmetry":
      self.D[-1] = self.D[-3]
    if self.BC_Rigidity_W == "periodic":
      self.D[0] = self.D[-2]
    if self.BC_Rigidity_E == "periodic":
      self.D[-1] = self.D[-3]
    
  def get_coeff_values(self):
  
    ##############################
    # BUILD GENERAL COEFFICIENTS #
    ##############################

    # l2 corresponds to top value in solution vector, so to the left (-) side
    # Good reference for how to determine central difference (and other) coefficients is:
    # Fornberg, 1998: Generation of Finite Difference Formulas on Arbitrarily Spaced Grids

    ###################################################
    # DEFINE SUB-ARRAYS FOR DERIVATIVE DISCRETIZATION #
    ###################################################
    Dm1 = self.D[:-2]
    D0  = self.D[1:-1]
    Dp1 = self.D[2:]

    ###########################################################
    # DEFINE COEFFICIENTS TO W_-2 -- W_+2 WITH B.C.'S APPLIED #
    ###########################################################
    self.l2_coeff_i = ( Dm1/2. + D0 - Dp1/2. ) / self.dx4
    self.l1_coeff_i = ( -6.*D0 + 2.*Dp1 ) / self.dx4
    self.c0_coeff_i = ( -2.*Dm1 + 10.*D0 - 2.*Dp1 ) / self.dx4 + self.drho*self.g
    self.r1_coeff_i = ( 2.*Dm1 - 6.*D0 ) / self.dx4
    self.r2_coeff_i = ( -Dm1/2. + D0 + Dp1/2. ) / self.dx4
    # These will be just the 1, -4, 6, -4, 1 for constant Te

    ###################################################################
    # START DIAGONALS AS SIMPLY THE BASE COEFFICIENTS, WITH NO B.C.'S #
    ###################################################################
    self.l2 = self.l2_coeff_i.copy()
    self.l1 = self.l1_coeff_i.copy()
    self.c0 = self.c0_coeff_i.copy()
    self.r1 = self.r1_coeff_i.copy()
    self.r2 = self.r2_coeff_i.copy()

    # Number of columns; equals number of rows too - square coeff matrix
    self.ncolsx = self.c0.shape[0]
    
    # Either way, the way that Scipy stacks is not the same way that I calculate
    # the rows. It runs offsets down the column instead of across the row. So
    # to simulate this, I need to re-zero everything. To do so, I use 
    # numpy.roll. (See self.build_diagonals.)
    
  def BC_Flexure(self):

    # Some links that helped me teach myself how to set up the boundary conditions
    # in the matrix for the flexure problem:
    # 
    # Good explanation of and examples of boundary conditions
    # https://en.wikipedia.org/wiki/Euler%E2%80%93Bernoulli_beam_theory#Boundary_considerations
    # 
    # Copy of Fornberg table:
    # https://en.wikipedia.org/wiki/Finite_difference_coefficient
    # 
    # Implementing b.c.'s:
    # http://scicomp.stackexchange.com/questions/5355/writing-the-poisson-equation-finite-difference-matrix-with-neumann-boundary-cond
    # http://scicomp.stackexchange.com/questions/7175/trouble-implementing-neumann-boundary-conditions-because-the-ghost-points-cannot
    
    if self.Verbose:
      print("Boundary condition, West:", self.BC_W, type(self.BC_W))
      print("Boundary condition, East:", self.BC_E, type(self.BC_E))

    # In 2D, these are handled inside the function; in 1D, there are separate
    # defined functions. Keeping these due to inertia and fear of cut/paste
    # mistakes
    if self.BC_E == '0Displacement0Slope' or self.BC_W == '0Displacement0Slope':
      self.BC_0Displacement0Slope()
    if self.BC_E == '0Slope0Shear' or self.BC_W == '0Slope0Shear':
      self.BC_0Slope0Shear()
    if self.BC_E == '0Moment0Shear' or self.BC_W == '0Moment0Shear':
      self.BC_0Moment0Shear()
    if self.BC_E == 'Mirror' or self.BC_W == 'Mirror':
      self.BC_Mirror()
    if self.BC_E == 'Periodic' and self.BC_W == 'Periodic':
      self.BC_Periodic()
    if self.BC_E == 'Sandbox' or self.BC_W == 'Sandbox':
      # Sandbox is the developer's testing ground
      sys.exit("Sandbox Closed")

  def build_diagonals(self):
    """
    Builds the diagonals for the coefficient array
    """

    ##########################################################
    # INCORPORATE BOUNDARY CONDITIONS INTO COEFFICIENT ARRAY #
    ##########################################################

    # Roll to keep the proper coefficients at the proper places in the
    # arrays: Python will naturally just do vertical shifts instead of 
    # diagonal shifts, so this takes into account the horizontal compoent 
    # to ensure that boundary values are at the right place.
    self.l2 = np.roll(self.l2, -2)
    self.l1 = np.roll(self.l1, -1)
    self.r1 = np.roll(self.r1, 1)
    self.r2 = np.roll(self.r2, 2)

    # Then assemble these rows: this is where the periodic boundary condition 
    # can matter.
    if self.coeff_matrix is not None:
      pass
    elif self.BC_E == 'Periodic' and self.BC_W == 'Periodic':
      # In this case, the boundary-condition-related stacking has already 
      # happened inside b.c.-handling function. This is because periodic
      # boundary conditions require extra diagonals to exist on the edges of 
      # the solution array
      pass
    else:
      self.diags = np.vstack((self.l2,self.l1,self.c0,self.r1,self.r2))
      self.offsets = np.array([-2,-1,0,1,2])

    # Everybody now (including periodic b.c. cases)
    self.coeff_matrix = spdiags(self.diags, self.offsets, self.nx, self.nx, format='csr')
  
  def BC_Periodic(self):
    """
    Periodic boundary conditions: wraparound to the other side.
    """
    if self.BC_E == 'Periodic' and self.BC_W == 'Periodic':
      # If both boundaries are periodic, we are good to go (and self-consistent)
      pass # It is just a shift in the coeff. matrix creation.
    else:
      # If only one boundary is periodic and the other doesn't implicitly 
      # involve a periodic boundary, this is illegal!
      # I could allow it, but would have to rewrite the Periodic b.c. case,
      # which I don't want to do to allow something that doesn't make 
      # physical sense... so if anyone wants to do this for some unforeseen 
      # reason, they can just split my function into two pieces themselves.i
      sys.exit("Having the boundary opposite a periodic boundary condition\n"+
               "be fixed and not include an implicit periodic boundary\n"+
               "condition makes no physical sense.\n"+
               "Please fix the input boundary conditions. Aborting.")
    self.diags = np.vstack((self.r1,self.r2,self.l2,self.l1,self.c0,self.r1,self.r2,self.l2,self.l1))
    self.offsets = np.array([1-self.ncolsx,2-self.ncolsx,-2,-1,0,1,2,self.ncolsx-2,self.ncolsx-1])

  def BC_0Displacement0Slope(self):
    """
    0Displacement0Slope boundary condition for 0 deflection.
    This requires that nothing be done to the edges of the solution array, 
    because the lack of the off-grid terms implies that they go to 0
    Here we just turn the cells outside the array into nan, to ensure that 
    we are not accidentally including the wrong cells here (and for consistency 
    with the other solution types -- this takes negligible time)
    """
    if self.BC_W == '0Displacement0Slope':
      i=0
      self.l2[i] = np.nan
      self.l1[i] = np.nan
      self.c0[i] += 0
      self.r1[i] += 0
      self.r2[i] += 0
      i=1
      self.l2[i] = np.nan
      self.l1[i] += 0
      self.c0[i] += 0
      self.r1[i] += 0
      self.r2[i] += 0
    if self.BC_E == '0Displacement0Slope':
      i=-2
      self.l2[i] += 0
      self.l1[i] += 0
      self.c0[i] += 0
      self.r1[i] += 0
      self.r2[i] = np.nan
      i=-1
      self.l2[i] += 0
      self.l1[i] += 0
      self.c0[i] += 0
      self.r1[i] = np.nan
      self.r2[i] = np.nan

  def BC_0Slope0Shear(self):
    i=0
    """
    This boundary condition is esentially a Neumann 0-gradient boundary 
    condition with that 0-gradient state extended over a longer part of 
    the grid such that the third derivative also equals 0.
    
    This boundary condition has more of a geometric meaning than a physical 
    meaning. It produces a state in which the boundaries have to have all 
    gradients in deflection go to 0 (i.e. approach constant values) while 
    not specifying what those values must be.
    
    This uses a 0-curvature boundary condition for elastic thickness 
    that extends outside of the computational domain.
    """
    
    if self.BC_W == '0Slope0Shear':
      i=0
      self.l2[i] = np.nan
      self.l1[i] = np.nan
      self.c0[i] += 0
      self.r1[i] += self.l1_coeff_i[i]
      self.r2[i] += self.l2_coeff_i[i]
      i=1
      self.l2[i] = np.nan
      self.l1[i] += 0
      self.c0[i] += 0
      self.r1[i] += 0
      self.r2[i] += self.l2_coeff_i[i]
    if self.BC_E == '0Slope0Shear':
      i=-2
      self.l2[i] += self.r2_coeff_i[i]
      self.l1[i] += 0
      self.c0[i] += 0
      self.r1[i] += 0
      self.r2[i] = np.nan
      i=-1
      self.l2[i] += self.r2_coeff_i[i]
      self.l1[i] += self.r1_coeff_i[i]
      self.c0[i] += 0
      self.r1[i] = np.nan
      self.r2[i] = np.nan

  def BC_0Moment0Shear(self):
    """
    d2w/dx2 = d3w/dx3 = 0
    (no moment or shear)
    This simulates a free end (broken plate, end of a cantilevered beam: 
    think diving board tip)
    It is *not* yet set up to have loads placed on the ends themselves: 
    (look up how to do this, thought Wikipdia has some info, but can't find
    it... what I read said something about generalizing)
    """

    # First, just define coefficients for each of the positions in the array
    # These will be added in code instead of being directly combined by 
    # the programmer (as I did above (now deleted) for constant Te), which might add 
    # rather negligibly to the compute time but save a bunch of possibility 
    # for unfortunate typos!

    # Also using 0-curvature boundary condition for D (i.e. Te)
    if self.BC_W == '0Moment0Shear':
      i=0
      self.l2[i] += np.nan
      self.l1[i] += np.nan
      self.c0[i] += 4*self.l2_coeff_i[i] + 2*self.l1_coeff_i[i]
      self.r1[i] += -4*self.l2_coeff_i[i] - self.l1_coeff_i[i]
      self.r2[i] += self.l2_coeff_i[i]
      i=1
      self.l2[i] += np.nan
      self.l1[i] += 2*self.l2_coeff_i[i]
      self.c0[i] += 0
      self.r1[i] += -2*self.l2_coeff_i[i]
      self.r2[i] += self.l2_coeff_i[i]
    
    if self.BC_E == '0Moment0Shear':
      i=-2
      self.l2[i] += self.r2_coeff_i[i]
      self.l1[i] += -2*self.r2_coeff_i[i]
      self.c0[i] += 0
      self.r1[i] += 2*self.r2_coeff_i[i]
      self.r2[i] += np.nan
      i=-1
      self.l2[i] += self.r2_coeff_i[i]
      self.l1[i] += -4*self.r2_coeff_i[i] - self.r1_coeff_i[i]
      self.c0[i] += 4*self.r2_coeff_i[i] + + 2*self.r1_coeff_i[i]
      self.r1[i] += np.nan
      self.r2[i] += np.nan

  def BC_Mirror(self):
    """
    Mirrors qs across the boundary on either the west (left) or east (right) 
    side, depending on the selections.
    
    This can, for example, produce a scenario in which you are observing 
    a mountain range up to the range crest (or, more correctly, the halfway 
    point across the mountain range).
    """
    if self.BC_W == 'Mirror':
      i=0
      #self.l2[i] += np.nan
      #self.l1[i] += np.nan
      self.c0[i] += 0
      self.r1[i] += self.l1_coeff_i[i]
      self.r2[i] += self.l2_coeff_i[i]
      i=1
      #self.l2[i] += np.nan
      self.l1[i] += 0
      self.c0[i] += self.l2_coeff_i[i]
      self.r1[i] += 0
      self.r2[i] += 0
    
    if self.BC_E == 'Mirror':
      i=-2
      self.l2[i] += 0
      self.l1[i] += 0
      self.c0[i] += self.r2_coeff_i[i]
      self.r1[i] += 0
      #self.r2[i] += np.nan
      i=-1
      self.l2[i] += self.r2_coeff_i[i]
      self.l1[i] += self.r1_coeff_i[i]
      self.c0[i] += 0
      #self.r1[i] += np.nan
      #self.r2[i] += np.nan
    
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
    self.maxFlexuralWavelength_ncells = int(np.ceil(self.maxFlexuralWavelength / self.dx))
    
  def fd_solve(self):
    """
    w = fd_solve()
      where coeff is the sparse coefficient matrix output from function
      coeff_matrix and qs is the array of loads (stresses)

    Sparse solver for one-dimensional flexure of an elastic plate
    """
    
    if self.Debug:
      print("qs", self.qs.shape)
      print("Te", self.Te.shape)
      self.calc_max_flexural_wavelength()
      print("maxFlexuralWavelength_ncells', self.maxFlexuralWavelength_ncells")
    
    if self.Solver == "iterative" or self.Solver == "Iterative":
      if self.Debug:
        print("Using generalized minimal residual method for iterative solution")
      if self.Verbose:
        print("Converging to a tolerance of", self.iterative_ConvergenceTolerance, "m between iterations")
      # qs negative so bends down with positive load, bends up with neative load 
      # (i.e. material removed)
      w = isolve.lgmres(self.coeff_matrix, -self.qs, tol=self.iterative_ConvergenceTolerance)  
      self.w = w[0] # Reach into tuple to get my array back
    else:
      if self.Solver == 'direct' or self.Solver == 'Direct':
        if self.Debug:
          print("Using direct solution with UMFpack")
      else:
        print("Solution type not understood:")
        print("Defaulting to direct solution with UMFpack")
      # UMFpack is now the default, but setting true just to be sure in case
      # anything changes
      # qs negative so bends down with positive load, bends up with neative load 
      # (i.e. material removed)
      self.w = spsolve(self.coeff_matrix, -self.qs, use_umfpack=True)
    
    if self.Debug:
      print("w.shape:")
      print(self.w.shape)
      print("w:")
      print(self.w)
    
