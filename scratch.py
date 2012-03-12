    # And point to all of the class variables so I don't have to say "self"
    # constantly... and mostly so I don't have to change a lot of local code
    ncolsx = self.ncolsx
    nrowsy = self.nrowsy
    cj2i0 = self.cj2i0
    cj1i_1 = self.cj1i_1
    cj1i0 = self.cj1i0
    cj1i1 = self.cj1i1
    cj0i_2 = self.cj0i_2
    cj0i_1 = self.cj0i_1
    cj0i0 = self.cj0i0
    cj0i1 = self.cj0i1
    cj0i2 = self.cj0i2
    cj_1i_1 = self.cj_1i_1
    cj_1i0 = self.cj_1i0
    cj_1i1 = self.cj_1i1
    cj_2i0 = self.cj_2i0
    

  def pad_for_NoOutsideLoads(self):
    """
    Creates an array, q0pad, that is padded to the East, Southeast, and
    South. This makes it easy to re-extract the domain of interest. 
    It then is given periodic boundary conditions on all sides to 
    more efficiently share the region of padding
    """

    # First, see how far to pad: 1 flexural wavelength
    if np.isscalar(self.D):
      Dmax = self.D
    else:
      Dmax = self.D.max()

    # This is an approximation if there is fill that evolves with iterations 
    # (e.g., water), but should be good enough that this won't do much to it
    alpha = (Dmax/(self.drho*self.g))**.25 # 2D flexural parameter
    flexuralWavelength = 2*np.pi*alpha

    # Then calculate the components for padding
    # Distance
    q0pad_nx = int(np.ceil(flexuralWavelength / self.dx))
    q0pad_ny = int(np.ceil(flexuralWavelength / self.dy))
    # Component
    q0padE = np.zeros((self.q0.shape[0],q0pad_nx))
    q0padS = np.zeros((q0pad_ny,self.q0.shape[1]))
    q0padSE = np.zeros((q0pad_ny,q0pad_nx))
    # Further component building
    q0west = np.concatenate((self.q0,q0padS))
    q0east = np.concatenate((q0padE,q0padSE))

  def assemble_blocks_sparse_1row(self,i):
    """
    Build the coefficient matrix by combining array entries and blank arrays
    into one gigantic sparse matrix

    Loop over rows, with each block being a square with both side
    lengths being the number of columns (i.e. the number of entries 
    in each row)

    Adding blocks of zeros of size ncolsx to each side of the array
    in order to pad it for later concatenation (i.e. making all same size
    and putting nonzero blocks in correct places)
    There has got to be a more elegant way to do this, but I don't know it
    and it will take longer to figure out than to just write it all out now.
    """

    # Boundary conditions never change interior cells
    if i>=3 and i<=self.nrowsy-4: # If no truncation at edges (normal case)
      leftzeros = sparse.dia_matrix((self.ncolsx,self.ncolsx*(i-2)))
      rightzeros = sparse.dia_matrix((self.ncolsx,self.ncolsx*(self.nrowsy-3-i)))
      coeff_row = sparse.hstack( [leftzeros,self.l2,self.l1,self.c0,self.r1,self.r2,rightzeros] )
    # Shouldn't have to be a special case - but have to b/c no leftzeros
    elif i==2:
      rightzeros = sparse.dia_matrix((self.ncolsx,self.ncolsx*(self.nrowsy-3-i)))
      coeff_row = sparse.hstack( [self.l2,self.l1,self.c0,self.r1,self.r2,rightzeros] )
    # Shouldn't have to be a special case - but have to b/c no rightzeros
    elif i==self.nrowsy-3:
      leftzeros = sparse.dia_matrix((self.ncolsx,self.ncolsx*(i-2)))
      coeff_row = sparse.hstack( [leftzeros,self.l2,self.l1,self.c0,self.r1,self.r2] )
    else:
      # NEED TO FIX THIS: DIRICHLET BC SUPERCEDES ANYTHING PERIODIC
      if self.BC_N == 'Dirichlet' or self.BC_S == 'Dirichlet':
        if self.BC_N == 'Dirichlet':
          print 'DN'
          if i==0:
            print 'DN'
            rightzeros = sparse.dia_matrix((self.ncolsx,self.ncolsx*(self.nrowsy-3-i)))
            coeff_row = sparse.hstack( [self.c0,self.r1,self.r2,rightzeros] )
          elif i==1:
            print 'DN'
            rightzeros = sparse.dia_matrix((self.ncolsx,self.ncolsx*(self.nrowsy-3-i)))
            coeff_row = sparse.hstack( [self.l1,self.c0,self.r1,self.r2,rightzeros] )
        if self.BC_S == 'Dirichlet':
          if i==self.nrowsy-1:
            print 'DS'
            leftzeros = sparse.dia_matrix((self.ncolsx,self.ncolsx*(i-2)))
            coeff_row = sparse.hstack( [leftzeros,self.l2,self.l1,self.c0] )
          elif i==self.nrowsy-2:
            print 'DS'
            leftzeros = sparse.dia_matrix((self.ncolsx,self.ncolsx*(i-2)))
            coeff_row = sparse.hstack( [leftzeros,self.l2,self.l1,self.c0,self.r1] )
      if self.BC_N == 'Periodic' and self.BC_S == 'Periodic':




    # E-W
    if self.BC_W == 'Periodic' and self.BC_E == 'Periodic':
      # If both boundaries are periodic, we are good to go (and self-consistent)
      self.BC_Periodic_EW()

    # N-S
    if self.BC_N == 'Periodic' and self.BC_S == 'Periodic':
      # If both boundaries are periodic, we are good to go (and self-consistent)
      self.BC_Periodic_NS()


      


    #############
    # DIRICHLET #
    #############

    # DEFAULTS TO 0-DISPLACEMENT DIRICHLET BOUNDARY CONDITIONS IF NOTHING
    # ELSE IS SET - SO DIRICHLET HERE IS BLANK!




  ##############
  ## PLOTTING ##
  ##############

  def imshow(self,image):
    # Plot if you want to - for troubleshooting
    from matplotlib.pyplot import imshow, show, colorbar, figure
    figure()
    imshow(image,interpolation='nearest') #,interpolation='nearest'
    colorbar()
    show()
    
