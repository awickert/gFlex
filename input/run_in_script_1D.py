#! /usr/bin/env python

import numpy as np

import gflex

flex = gflex.F1D()

flex.quiet = True

flex.method = "FD"  # Solution method: * FD (finite difference)
#                  * FFT (spectral)
#                  * SAS (superposition of analytical solutions)
#                  * SAS_NG (ungridded SAS)


flex.g = 9.8  # acceleration due to gravity
flex.E = 65e9  # Young's Modulus
flex.nu = 0.25  # Poisson's Ratio
flex.rho_m = 3300.0  # MantleDensity
flex.rho_fill = 1000.0  # InfillMaterialDensity

flex.te = 30000.0  # *np.ones(500) # Elastic thickness -- scalar but may be an array
# flex.te[-3:] = 0
flex.qs = np.zeros(300)
flex.qs[100:200] += 1e6  # surface load stresses
flex.dx = 4000.0  # grid cell size [m]
flex.bc_west = "zero_displacement_zero_slope"  # west boundary condition
flex.bc_east = "zero_moment_zero_shear"  # east boundary condition

flex.sigma_xx = 100.0  # Normal stress on the edge of the plate

flex.initialize()
flex.run()
flex.finalize()

# If you want to plot the output
flex.plotChoice = "combo"
# An output file for deflections could also be defined here
# flex.wOutFile =
flex.output()  # Plots and/or saves output, or does nothing, depending on
# whether flex.plotChoice and/or flex.wOutFile have been set
# TO OBTAIN OUTPUT DIRECTLY IN PYTHON, you can assign the internal variable,
# flex.w, to another variable -- or as an element in a list if you are looping
# over many runs of gFlex:
deflection = flex.w
