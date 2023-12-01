#! /usr/bin/env python

import numpy as np
from matplotlib import pyplot as plt

import gflex


def test_main():
    flex = gflex.f2d.F1D()

    flex.Quiet = True

    flex.Method = 'SAS' # Solution method: * FD (finite difference)
                       #                  * SAS (superposition of analytical solutions)
                       #                  * SAS_NG (ungridded SAS)

    flex.g = 9.8 # acceleration due to gravity
    flex.E = 65E9 # Young's Modulus
    flex.nu = 0.25 # Poisson's Ratio
    flex.rho_m = 3300. # MantleDensity
    flex.rho_fill = 1000. # InfiillMaterialDensity

    flex.Te = 30000.#*np.ones(500) # Elastic thickness -- scalar but may be an array
    flex.qs = np.zeros(10); flex.qs[5] += 1E6 # surface load stresses
    flex.dx = 4000. # grid cell size [m]
    flex.BC_W = '0Displacement0Slope' # west boundary condition
    flex. BC_E = '0Moment0Shear' # east boundary condition

    flex.initialize()
    flex.run()
    flex.finalize()

    # If you want to plot the output
    #flex.plotChoice='combo'
    # An output file for deflections could also be defined here
    # flex.wOutFile =
    flex.output() # Plots and/or saves output, or does nothing, depending on
                  # whether flex.plotChoice and/or flex.wOutFile have been set
    # TO OBTAIN OUTPUT DIRECTLY IN PYTHON, you can assign the internal variable,
    # flex.w, to another variable -- or as an element in a list if you are looping
    # over many runs of gFlex:
    deflection = flex.w

if __name__ == '__main__':
    test_main()
