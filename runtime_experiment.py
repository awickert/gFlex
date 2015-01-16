#! /usr/bin/python

# GFLEX
from base import *
from f1d import *
from f2d import *
from prattairy import *

# PYTHON
import numpy as np
import time



# This code is for 2D flexural isostasy (and instantiate)
obj = F2D()
obj.set_value('model', 'flexure')
# Set verbosity
obj.set_value('Quiet', True)
# Always use the van Wees and Cloetingh (1994) solution type.
# It is the best.
obj.set_value('PlateSolutionType', 'vWC1994')
 
# Make a bunch of standard selections
obj.set_value('GravAccel', 9.8)
obj.set_value('YoungsModulus', 65E10)
obj.set_value('PoissonsRatio', .25)
obj.set_value('MantleDensity', 3300.)
obj.set_value('InfillMaterialDensity', 0.)

# And solver / iterations (if needed)
obj.set_value('Solver', 'direct')
#obj.set_value('ConvergenceTolerance', ConvergenceTolerance)

# Elastic thickness
obj.set_value('ElasticThickness', FlexureTe)

for 

  # Grid size and spacing
  # SIZE -- MAINTAIN BLOCKS INSIDE, TOO -- START WORK HERE
  obj.set_value('GridSpacing_x', dx)
  obj.set_value('GridSpacing_y', dx)

  if method == 'SAS':
    obj.set_value('method', 'SAS')
  elif method == 'FD':
    obj.set_value('method', 'FD')
    obj.set_value('Solver', 'direct')
    # Always use the van Wees and Cloetingh (1994) solution type.
    # It is the best.
    obj.set_value('PlateSolutionType', 'vWC1994')
  # No need for "else" here:
  # Will automatically fail via parser if value is out of range

  

  # Set all boundary conditions
  obj.set_value('BoundaryCondition_East', bc)
  obj.set_value('BoundaryCondition_West', bc)
  obj.set_value('BoundaryCondition_North', bc)
  obj.set_value('BoundaryCondition_South', bc)



# CALCULATE!
obj.initialize()
obj.run()
obj.finalize()










