#! /usr/bin/env python

import numpy as np

import gflex

flex = gflex.F2D()

flex.quiet = False

flex.method = "FD"  # Solution method: * FD (finite difference)
#                  * FFT (spectral)
#                  * SAS (superposition of analytical solutions)
#                  * SAS_NG (ungridded SAS)

flex.g = 9.8  # acceleration due to gravity
flex.E = 65e9  # Young's Modulus
flex.nu = 0.25  # Poisson's Ratio
flex.rho_m = 3300.0  # MantleDensity
flex.rho_fill = 0.0  # InfillMaterialDensity

flex.te = 35000.0 * np.ones(
    (50, 50)
)  # Elastic thickness [m] -- scalar but may be an array
flex.te[:, -3:] = 0.0
flex.qs = np.zeros((50, 50))  # Template array for surface load stresses
flex.qs[10:40, 10:40] += 1e6  # Populating this template
flex.dx = 5000.0  # grid cell size, x-oriented [m]
flex.dy = 5000.0  # grid cell size, y-oriented [m]
# Boundary conditions can be:
# (FD): 0Displacement0Slope, 0Displacement0Moment, 0Moment0Shear, 0Slope0Shear, Mirror, or Periodic
# For SAS or SAS_NG, NoOutsideLoads is valid, and no entry defaults to this
flex.bc_west = "0Displacement0Slope"  # west boundary condition
flex.bc_east = "0Moment0Shear"  # east boundary condition
flex.bc_south = "0Displacement0Slope"  # south boundary condition
flex.bc_north = "0Displacement0Slope"  # north boundary condition

# latitude/longitude solutions are exact for SAS, approximate otherwise
# latlon = # true/false: flag to enable lat/lon input. Defaults False.
# PlanetaryRadius = # radius of planet [m], for lat/lon solutions

flex.initialize()
flex.run()
flex.finalize()

# If you want to plot the output
flex.plotChoice = "both"
# An output file for deflections could also be defined here
# flex.wOutFile =
flex.output()  # Plots and/or saves output, or does nothing, depending on
# whether flex.plotChoice and/or flex.wOutFile have been set
# TO OBTAIN OUTPUT DIRECTLY IN PYTHON, you can assign the internal variable,
# flex.w, to another variable -- or as an element in a list if you are looping
# over many runs of gFlex:
deflection = flex.w
