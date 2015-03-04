#! /usr/bin/env python

import gflex
import numpy as np
from matplotlib import pyplot as plt

flex = gflex.F1D()

flex.Quiet = True

flex.Method = 'FD'
flex.Solver = 'direct'

flex.g = 9.8 # acceleration due to gravity
flex.E = 65E10 # Young's Modulus
flex.nu = 0.25 # Poisson's Ratio
flex.rho_m = 3300. # MantleDensity
flex.rho_fill = 0. # InfiillMaterialDensity

flex.Te = 35000. # Elastic thickness -- scalar but may be an array
flex.qs = np.zeros(50); flex.qs[10:40] += 1E6 # surface load stresses
flex.dx = 5000.
flex.BC_W = 'Dirichlet0' # west boundary condition
flex. BC_E = '0Moment0Shear' # east boundary condition

flex.initialize()
flex.run()
flex.finalize()

# If you want to plot the output
flex.plotChoice='combo'
# An output file could also be defined here
# flex.wOutFile = 
flex.output() # Plots and/or saves output, or does nothing, depending on
              # whether flex.plotChoice and/or flex.wOutFile have been set
