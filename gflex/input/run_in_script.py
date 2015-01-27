#! /usr/bin/python

import gflex
import numpy as np
from matplotlib import pyplot as plt

obj = gflex.F1D()
obj.set_value('Quiet', True)

obj.set_value('Method', 'FD')
obj.set_value('PlateSolutionType', 'vWC1994')
obj.set_value('Solver', 'direct')

obj.set_value('GravAccel', 9.8)
obj.set_value('YoungsModulus', 65E10)
obj.set_value('PoissonsRatio', .25)
obj.set_value('MantleDensity', 3300.)
obj.set_value('InfiillMaterialDensity', 0.)

obj.set_value('ElasticThickness', 35000.)
obj.set_value('Loads_grid_stress', 1E6*np.ones(50))
obj.set_value('GridSpacing_x', 5000.)
obj.set_value('BoundaryCondition_East', 'Dirichlet0')
obj.set_value('BoundaryCondition_West', '0Slope0Shear')

obj.initialize()
obj.run()
obj.finalize()
