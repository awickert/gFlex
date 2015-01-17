#! /usr/bin/python

"""
Code for testing gFlex

It is run in a semi-automated way to obtain results from different scenarios
that are then manually re-combined, plotted, and interpreted.
"""

# GFLEX
from base import *
from f1d import *
from f2d import *
from prattairy import *

# PYTHON
import numpy as np
import time

# plotting
from matplotlib import pyplot as plt

# This code is for 2D flexural isostasy (and instantiate)
obj = F2D()
obj.set_value('dimension', 2) # MEMO WAS NOT RECEIVED BY THE REST OF THE PROGRAM!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

# Always use the van Wees and Cloetingh (1994) solution type.
# It is the best.
obj.set_value('PlateSolutionType', 'vWC1994')

# Domain length
L = 1000E3 # domain length
q_load = 9702000

nelements = []
nloadelements = []
solvetime = []
i = 0

# load length
for l in [200E3, 100E3, 200E3, 500E3, 1000E3]:
  #for Te in [25000, GRID]:
  for Te in [25000]:
    #for bc in ["0Slope0Shear", 'Dirichlet0', 'Mirror', '0Slope0Shear']:
    #for bc in ['Periodic']:
    for bc in ['NoOutsideLoads']:
      #ne = []
      #nle = []
      #st = []
      #nelements = .append([])
      #nloadelements = []
      #solvetime = []
      #for method in ['SAS', 'SAS_NG', 'FD']:
      for method in ['SAS']:
        #for dx in [100, 500, 1000, 2000, 2500, 5000, 10000, 20000, 25000, 50000]:
        #for dx in [2000, 2500, 4000, 5000, 10000, 20000, 25000, 40000, 50000]: #1000, 2000, 2500, 4000, 5000, 10000, 
        for dx in [10000, 20000, 25000, 40000, 50000]: #1000, 2000, 2500, 4000, 5000, 10000, 

            obj.set_value('method', method)
            # Grid size and spacing
            # SIZE -- MAINTAIN BLOCKS INSIDE, TOO -- START WORK HERE
            obj.set_value('GridSpacing_x', dx)
            obj.set_value('GridSpacing_y', dx)
            
            
            # Make the array
            if L % dx != 0 or l % dx !=0:
              pass
            else:
              print dx 
              q = np.zeros((L/dx, L/dx))
              q[(L/2-l/2)/dx : (L/2+l/2)/dx+.0001, (L/2-l/2)/dx : (L/2+l/2)/dx+.0001] = q_load
              obj.set_value('Loads', q)
              
              #plt.imshow(obj.q0); plt.show()

              # Set all boundary conditions
              obj.set_value('BoundaryCondition_East', bc)
              obj.set_value('BoundaryCondition_West', bc)
              obj.set_value('BoundaryCondition_North', bc)
              obj.set_value('BoundaryCondition_South', bc)

              # Elastic thickness
              obj.set_value('ElasticThickness', Te)

              # CALCULATE!
              obj.initialize()
              obj.run()
              obj.finalize()
              del obj.D
              try:
                del obj.coeff_matrix # This was not being properly overwritten!
              except:
                pass

              nelements.append(np.prod(q.shape))
              nloadelements.append(np.sum(q > 0))
              solvetime.append(obj.time_to_solve)
              
              print 'Time to solve [s]:', obj.get_value('SolverTime')
              
              #ne.append(np.prod(q.shape))
              #nle.append(np.sum(q > 0))
              #st.append(obj.time_to_solve)
              
      #nelements.append(ne)
      #nloadelements.append(nle)
      #solvetime.append(st)



"""
# For analytical, should be number of cells TIMES number of times to calculate
# At least in some way, even with the stencil that is moved around and placed over
# -- well, but it will always be a power law relationship, so slope of relation
# to ncells or ncells on one side will just change
from scipy.optimize import curve_fit
def func(x, a, b):
  return a * x**b

popt, pcov = curve_fit(func, np.array(nelements)*np.array(nelements)/100., solvetime)

fig = plt.figure()
ax = fig.add_subplot(111)
plt.loglog(nelements, solvetime, 'ko')
x_cf = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 100)
plt.loglog(x_cf, func(x_cf, popt[0], popt[1]), 'k-')
plt.show()
"""

nelements = np.array(nelements)
nloadelements = np.array(nloadelements)
solvetime = np.array(solvetime)

# For Finite Difference or both, I think

# Curve fit
from scipy.optimize import curve_fit
def func(x, a, b):
  return a * x**b

popt, pcov = curve_fit(func, nelements, solvetime)

fig = plt.figure()
ax = fig.add_subplot(111)
plt.loglog(nelements, solvetime, 'ko')
x_cf = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 100)
plt.loglog(x_cf, func(x_cf, popt[0], popt[1]), 'k-')
plt.show()

outarray = np.vstack((nelements, nloadelements, nelements*nloadelements, solvetime)).transpose()

#np.savetxt('benchmark/FD_Periodic_variable_load_constant_Te.csv.csv', outarray, delimiter=',',)
#np.savetxt('benchmark/FD_not_Periodic_200km_load_constant_Te.csv', outarray, delimiter=',',)

#np.savetxt('benchmark/SAS_NG.csv', outarray, delimiter=',',)
#np.savetxt('benchmark/SAS_NG.csv', outarray, delimiter=',',)


