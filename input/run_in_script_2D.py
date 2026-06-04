#! /usr/bin/env python

import numpy as np

import gflex

flex = gflex.F2D()

flex.quiet = False

flex.method = "fd"  # Solution method: * fd (finite difference)
#                  * fft (spectral)
#                  * sas (superposition of analytical solutions)
#                  * sas_ng (ungridded SAS)

flex.g = 9.8  # acceleration due to gravity
flex.E = 65e9  # Young's Modulus
flex.nu = 0.25  # Poisson's Ratio
flex.rho_m = 3300.0  # mantle_density
flex.rho_fill = 0.0  # infill_material_density

flex.T_e = 35000.0 * np.ones(
    (50, 50)
)  # Elastic thickness [m] -- scalar but may be an array
flex.T_e[:, -3:] = 0.0
flex.qs = np.zeros((50, 50))  # Template array for surface load stresses
flex.qs[10:40, 10:40] += 1e6  # Populating this template
flex.dx = 5000.0  # grid cell size, x-oriented [m]
flex.dy = 5000.0  # grid cell size, y-oriented [m]
# Boundary conditions can be:
# (FD): zero_displacement_zero_slope, zero_displacement_zero_moment, zero_moment_zero_shear, zero_slope_zero_shear, mirror, or periodic
# For SAS or SAS_NG, no_outside_loads is valid, and no entry defaults to this
flex.bc_west = "zero_displacement_zero_slope"  # west boundary condition
flex.bc_east = "zero_moment_zero_shear"  # east boundary condition
flex.bc_south = "zero_displacement_zero_slope"  # south boundary condition
flex.bc_north = "zero_displacement_zero_slope"  # north boundary condition

# latitude/longitude solutions are exact for SAS, approximate otherwise
# latlon = # true/false: flag to enable lat/lon input. Defaults False.
# PlanetaryRadius = # radius of planet [m], for lat/lon solutions

flex.initialize()
flex.run()

# Read deflection BEFORE calling finalize (finalize clears w and qs)
deflection = flex.w

# If you want to plot the output
flex.plot_choice = "both"
# An output file for deflections could also be defined here
# flex.w_out_file =
flex.output()  # Plots and/or saves output, or does nothing, depending on
# whether flex.plot_choice and/or flex.w_out_file have been set

flex.finalize()
