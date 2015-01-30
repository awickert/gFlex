#! /usr/bin/python

import gflex
import numpy as np
from matplotlib import pyplot as plt

obj = gflex.F1D()

obj.Quiet = True

obj.Method = 'FD'
obj.PlateSolutionType = 'vWC1994'
obj.Solver = 'direct'

obj.g = 9.8 # acceleration due to gravity
obj.E = 65E10 # Young's Modulus
obj.nu = 0.25 # Poisson's Ratio
obj.rho_m = 3300. # MantleDensity
obj.rho_fill = 0. # InfiillMaterialDensity

obj.Te = 35000. # Elastic thickness
obj.qs = 1E6*np.ones(50) # surface load stresses
obj.dx = 5000.
obj.BC_W = '0Slope0Shear' # west boundary condition
obj. BC_E = 'Dirichlet0' # east boundary condition
obj.AAA = 'GridSpacing_x', 5000.

obj.initialize()
obj.run()
obj.finalize()
