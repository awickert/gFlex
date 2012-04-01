#! /usr/bin/python

"""
Started 13 March 2011 to vary dx in 1D flexure and see where the solution starts to break down when a gradient in elastic thickness is involved
In this way, I can define a phase space of stable solutions
"""

# FLEXURE
from base import *
from f1d import *
from f2d import *
from prattairy import *

# Start up grids to store compensation and set dx
#compensation = []
dx = np.linspace(10,20000,100)
#dx = np.array([1000,3000])

# Set up other grids
amplitudes = np.linspace(0,20000,100)
lambda_min = [] # Approx. cutoff for good solutions

# Compute grids for solutions
minTe = 20000
maxmin_ratio = 2*amplitudes / minTe

# Option to use a pre-made input file as a shortcut to keep from needing all getters and setters
#Filename = 'input/input_f1d'
# Instantiate
obj = F1D()
obj.set_value('model', 'flexure')
obj.set_value('dimension', 1)
obj.set_value('GravAccel', 9.8)
obj.set_value('method', 'FD')
obj.set_value('Solver', 'direct')

# Make a bunch of standard selections
obj.set_value('YoungsModulus', 1.0E11)
obj.set_value('PoissonsRatio', 0.25)
obj.set_value('GravAccel', 9.8)
obj.set_value('MantleDensity', 3300)

# Set all boundary conditions to Mirror
obj.set_value('BoundaryCondition_East', 'Periodic')
obj.set_value('BoundaryCondition_West', 'Periodic')


###############################
# Loads and elastic thickness #
###############################

ncells = 51 # for q0, Te is 2 greater than this

q0 = 3300*9.8*100*np.ones(ncells)
#q0[:ncells/3] = 0
#q0[-ncells/3:] = 0

obj.set_value('Loads', q0)

print q0.shape

# Always keep drho = 0
obj.set_value('InfillMaterialDensity', 0)
obj.drho = obj.rho_m - obj.rho_fill


# Plotting
obj.plotChoice = None


# Initialize
# No filename: getter/setter interface isn't totally worked out
obj.filename = None
Filename = None
obj.initialize(Filename)

for amplitude in amplitudes:

  compensation = []

  # Sinusoidal flexural wavelengths that allow periodic b.c.'s:
  #amplitude = 5000 # 1/2 of total variation in Te
  #minTe = 20000
  xi_periodic = np.arange(0,ncells+2) * 2 * np.pi / (ncells+2)
  Te = amplitude * (np.cos(xi_periodic)) + amplitude + minTe

  obj.set_value('ElasticThickness', Te)# np.ones(Te.shape)*20000)#

  for dxi in dx:

    # Vary the dx
    obj.set_value('GridSpacing_x', dxi)

  #  obj.drho = obj.rho_m - obj.rho_fill
    
    # Run
    obj.run()

    # Get deflection array and calculate the total area deflected compared to what 
    # the load should be doing at full compensation
    w = obj.get_value('Deflection')
    q0 = obj.get_value('Loads')
    h_q0_mantle_equiv = q0 / (obj.g * obj.rho_m)
    compensation.append(-1 * np.mean(w) / np.mean(h_q0_mantle_equiv))

  # Finalize
  obj.finalize()

  cdiff = np.diff(compensation)

  # cdiff[0] happens at x[.5]
  # calculating the maximum negative lambda, so ... actually make it good later
  # just test for now
  try:
    minvalid_index = (cdiff<0).nonzero()[0].max()
    lambda_min.append(minvalid_index / float(ncells))
  except:
    # If there is no clear minumum
    lambda_min.append(np.nan)

# Asymptote of hyperbola


#dxmean = (dx[:-1] + dx[1:]) / 2

from matplotlib.pyplot import *

lambda_min = lambda_min[3:]
maxmin_ratio = maxmin_ratio[3:]

figure(1)
plot(dx, compensation,'ko-')
xlabel('dx')
ylabel('compensation')

figure(2)
plot(maxmin_ratio, lambda_min,'ko-')
xlabel('Maximum Te / Minimum Te')
ylabel('Min cell size as as fraction\nof Te variability wavelength')

show()

"""
plot(dxmean, cdiff, 'k-')
#ylim((-1,1))
xlabel('dx')
ylabel('compensation')
show()
"""

print maxmin_ratio
print lambda_min

from scipy.optimize import curve_fit

logm = np.log10(maxmin_ratio)
logl = np.log10(lambda_min)

def linfit(x,a,b):
  return a*x+b

popt, pcov = curve_fit(linfit, logm, logl)


"""
def sqrtfunc(x,a,b,c):
  return a * (x-b)**.5 + c

popt, pcov = curve_fit(sqrtfunc, np.array(lambda_min), maxmin_ratio)


def powerfunc(x,a,b,c):
  return a * x**2 + c

popt, pcov = curve_fit(sqrtfunc, maxmin_ratio, np.array(lambda_min),p0=(1,.5,.2))

"""
#print obj.Te

#plot(Te); show()
