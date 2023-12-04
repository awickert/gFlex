#! /usr/bin/env python

import numpy as np

from gflex.f2d import F2D


def test_main():
    flex = F2D()

    flex.Quiet = False

    flex.Method = "FD"  # Solution method: * FD (finite difference)
    #                  * SAS (superposition of analytical solutions)
    #                  * SAS_NG (ungridded SAS)
    flex.PlateSolutionType = "vWC1994"  # van Wees and Cloetingh (1994)
    # The other option is 'G2009': Govers et al. (2009)
    flex.Solver = "direct"  # direct or iterative
    # convergence = 1E-3 # convergence between iterations, if an iterative solution
    # method is chosen

    flex.g = 9.8  # acceleration due to gravity
    flex.E = 65e9  # Young's Modulus
    flex.nu = 0.25  # Poisson's Ratio
    flex.rho_m = 3300.0  # MantleDensity
    flex.rho_fill = 0.0  # InfiillMaterialDensity

    flex.Te = 35000.0 * np.ones(
        (50, 50)
    )  # Elastic thickness [m] -- scalar but may be an array
    flex.Te[:, -3:] = 0.0
    flex.qs = np.zeros((50, 50))  # Template array for surface load stresses
    flex.qs[10:40, 10:40] += 1e6  # Populating this template
    flex.dx = 5000.0  # grid cell size, x-oriented [m]
    flex.dy = 5000.0  # grid cell size, y-oriented [m]
    # Boundary conditions can be:
    # (FD): 0Slope0Shear, 0Moment0Shear, 0Displacement0Slope, Mirror, or Periodic
    # For SAS or SAS_NG, NoOutsideLoads is valid, and no entry defaults to this
    flex.BC_W = "0Displacement0Slope"  # west boundary condition
    flex.BC_E = "0Moment0Shear"  # east boundary condition
    flex.BC_S = "0Displacement0Slope"  # south boundary condition
    flex.BC_N = "0Displacement0Slope"  # north boundary condition

    # latitude/longitude solutions are exact for SAS, approximate otherwise
    # latlon = # true/false: flag to enable lat/lon input. Defaults False.
    # PlanetaryRadius = # radius of planet [m], for lat/lon solutions

    flex.initialize()
    flex.run()
    flex.finalize()

    # If you want to plot the output
    # flex.plotChoice='both'
    # An output file for deflections could also be defined here
    # flex.wOutFile =
    flex.output()  # Plots and/or saves output, or does nothing, depending on
    # whether flex.plotChoice and/or flex.wOutFile have been set


if __name__ == "__main__":
    test_main()
