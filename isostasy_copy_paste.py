try:
  already_run
except:
  already_run = False

import sys
if already_run:
  reload(sys.modules['base'])
  reload(sys.modules['f1d'])
  reload(sys.modules['f2d'])
  reload(sys.modules['prattairy'])

already_run = True

from base import *
from f1d import *
from f2d import *
from prattairy import *

from matplotlib import pyplot as plt

obj = Isostasy()

# Looks like it wants to be an input file!
filename = 'input/input_f1d' # it works for usage (1) and (2)
obj.whichModel(filename)

obj.model = 'flexure'
obj.dimension = 1

## SET MODEL TYPE AND DIMENSIONS HERE ##
########################################
if obj.model == 'flexure':
  if obj.dimension == 1:
    obj = F1D()
  elif obj.dimension == 2:
    obj = F2D()
elif obj.model == 'PrattAiry':
  obj = PrattAiry()

obj.initialize(filename)
## SET ALL OTHER MODEL PARAMETERS HERE ##
# obj.set_value('method','FD')
#########################################
obj.run()
obj.output() # Not part of IRF or BMI: Does standalone plotting and file output
## GET VALUES HERE ##
#wout = obj.get_value('Deflection')
#print wout
#####################
obj.finalize()


#plt.imshow(obj.Te_padded); plt.show(block=False)

#plt.figure()
#plt.imshow(obj.Te_orig); plt.show(block=False)

#plt.imshow(obj.w_padded); plt.show()
