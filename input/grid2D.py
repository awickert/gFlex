#! /usr/bin/env python

import numpy as np

import gflex

flex = gflex.F2D()

flex.quiet = False

flex.method = "fd"

flex.g = 9.8  # acceleration due to gravity
flex.E = 65e9  # Young's Modulus
flex.nu = 0.25  # Poisson's Ratio
flex.rho_m = 3300.0  # mantle_density
flex.rho_fill = 0.0  # infill_material_density

flex.te = 80000.0  # Elastic thickness -- scalar but may be an array
flex.qs = np.zeros((720, 360))  # Template array for surface load stresses
flex.qs[100:150, 100:150] += 1e6  # Populating this template
flex.dx = 80000.0
flex.dy = 111000.0
flex.bc_west = "periodic"  # west boundary condition
flex.bc_east = "periodic"  # east boundary condition
flex.bc_south = "periodic"  # south boundary condition
flex.bc_north = "periodic"  # north boundary condition

flex.initialize()
flex.run()
flex.finalize()

# If you want to plot the output
flex.plot_choice = "both"
# An output file could also be defined here
# flex.w_out_file =
flex.output()  # Plots and/or saves output, or does nothing, depending on
# whether flex.plot_choice and/or flex.w_out_file have been set
