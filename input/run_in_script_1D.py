#! /usr/bin/env python

import numpy as np

import gflex

flex = gflex.F1D()

flex.quiet = True

flex.method = "fd"  # Solution method: * fd (finite difference)
#                  * fft (spectral)
#                  * sas (superposition of analytical solutions)
#                  * sas_ng (ungridded SAS)


flex.g = 9.8  # acceleration due to gravity
flex.E = 65e9  # Young's Modulus
flex.nu = 0.25  # Poisson's Ratio
flex.rho_m = 3300.0  # mantle_density
flex.rho_fill = 1000.0  # infill_material_density

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

# Read deflection BEFORE calling finalize (finalize clears w and qs)
deflection = flex.w

# If you want to plot the output
flex.plot_choice = "combo"
# An output file for deflections could also be defined here
# flex.w_out_file =
flex.output()  # Plots and/or saves output, or does nothing, depending on
# whether flex.plot_choice and/or flex.w_out_file have been set

flex.finalize()
