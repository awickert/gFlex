#! /usr/bin/env python

import gflex
import numpy as np
from matplotlib import pyplot as plt

flex = gflex.F2D()

flex.Quiet = False

flex.Method = 'FD'
flex.PlateSolutionType = 'vWC1994'
flex.Solver = 'direct'

flex.g = 9.8 # acceleration due to gravity
flex.E = 65E9 # Young's Modulus
flex.nu = 0.25 # Poisson's Ratio
flex.rho_m = 3300. # MantleDensity
flex.rho_fill = 0. # InfiillMaterialDensity

flex.Te = 80000. # Elastic thickness -- scalar but may be an array
flex.qs = np.zeros((720, 360)) # Template array for surface load stresses
flex.qs[100:150, 100:150] += 1E6 # Populating this template
flex.dx = 80000.
flex.dy = 111000.
flex.BC_W = 'Periodic' # west boundary condition
flex.BC_E = 'Periodic' # east boundary condition
flex.BC_S = 'Periodic' # south boundary condition
flex.BC_N = 'Periodic' # north boundary condition

flex.initialize()
flex.run()
flex.finalize()

# If you want to plot the output
flex.plotChoice='both'
# An output file could also be defined here
# flex.wOutFile =
flex.output() # Plots and/or saves output, or does nothing, depending on
              # whether flex.plotChoice and/or flex.wOutFile have been set
