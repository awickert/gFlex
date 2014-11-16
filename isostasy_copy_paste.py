import sys

try:
  reload(sys.modules['base'])
  reload(sys.modules['f1d'])
  reload(sys.modules['f2d'])
  reload(sys.modules['prattairy'])
  print "RELOADING"
except:
  pass # 1st time running

from base import *
from f1d import *
from f2d import *
from prattairy import *

from matplotlib import pyplot as plt

# Looks like it wants to be an input file!
filename = 'input/input_f2d' # it works for usage (1) and (2)
obj = WhichModel(filename)

## SET MODEL TYPE AND DIMENSIONS HERE ##
########################################
if obj.model == 'flexure':
  if obj.dimension == 1:
    obj = F1D(filename)
  elif obj.dimension == 2:
    obj = F2D(filename)
elif obj.model == 'PrattAiry':
  obj = PrattAiry(filename)

self = obj # easier interaction
obj.initialize(filename)

## SET ALL OTHER MODEL PARAMETERS HERE ##
# obj.set_value('method','FD')
#########################################
obj.run()
obj.finalize()
#####################
obj.output() # Not part of IRF or BMI: Does standalone plotting and file output
## GET VALUES HERE ##
#wout = obj.get_value('Deflection')
#print wout

#plt.imshow(obj.Te_padded); plt.show(block=False)

#plt.figure()
#plt.imshow(obj.Te_orig); plt.show(block=False)

#plt.imshow(obj.w_padded); plt.show()
