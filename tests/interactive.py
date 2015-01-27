try:
  reload(gflex)
  print "RELOADING"
except:
  import gflex # first time running

from matplotlib import pyplot as plt
import numpy as np

# Looks like it wants to be an input file!
filename = '../input/input_f2d' # it works for usage (1) and (2)
obj = gflex.WhichModel(filename)

## SET MODEL TYPE AND DIMENSIONS HERE ##
########################################
if obj.dimension == 1:
  obj = gflex.F1D(filename)
elif obj.dimension == 2:
  obj = gflex.F2D(filename)

self = obj # easier interaction

#cm = obj.coeff_matrix.copy()
obj.set_value('CoeffArray', cm)

obj.initialize(filename)

#obj.set_value('GridSpacing_x', 50000)

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
