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
amplitudes = np.linspace(0,10000,600)
lambda_min = [] # Approx. cutoff for good solutions
teslope_max = []

# Compute Te ratio grids for solutions
minTe = 20000.
max_Tes = (2.*amplitudes + minTe)
maxmin_ratio_Te = max_Tes / minTe

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

ncells = 201 # for q0, Te is 2 greater than this

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


#------------------------
# Compute flexural wavelengths (no infill!)
# MINIMUM
minD = minTe**3 * obj.E / (12.*(1.-obj.nu))
min_alpha = (4 * minD / (obj.drho * obj.g))**.25 # 4 b/c 1D
min_lambdaf = 2. * np.pi * min_alpha
# MAXIMUMS
maxDs = max_Tes**3 * obj.E / (12.*(1-obj.nu))
max_alphas = (4 * maxDs / (obj.drho * obj.g))**.25 # 4 b/c 1D
max_lambdafs = 2. * np.pi * max_alphas
# RATIO
maxmin_ratio_lambdaf = max_lambdafs / min_lambdaf
#------------------------


# problem setup
lambda_Tevar = 500000. # wavelength of Te variations
ncells_range = np.arange(501,8,-10)
dx_range = lambda_Tevar / ncells_range

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
    # Run
    obj.run()
    #
    # Get deflection array and calculate the total area deflected compared to what 
    # the load should be doing at full compensation
    w = obj.get_value('Deflection')
    q0 = obj.get_value('Loads')
    h_q0_mantle_equiv = q0 / (obj.g * obj.rho_m)
    compensation.append(-1 * np.mean(w) / np.mean(h_q0_mantle_equiv))
  
  """
  compensation = []
  # Vary the number of cells per flexural wavelength
  counter = 0
  for ncells in ncells_range:
    #
    q0 = 3300*9.8*100*np.ones(ncells)
    obj.set_value('Loads', q0)
    #
    xi_periodic = np.arange(0,ncells+2) * 2 * np.pi / (ncells+2)
    Te = amplitude * (np.cos(xi_periodic)) + amplitude + minTe
    obj.set_value('ElasticThickness', Te)# np.ones(Te.shape)*20000)#
    #
    obj.set_value('GridSpacing_x',dx_range[counter])
    counter += 1
    # Run
    obj.run()
    #
    # Get deflection array and calculate the total area deflected compared to what 
    # the load should be doing at full compensation
    w = obj.get_value('Deflection')
    q0 = obj.get_value('Loads')
    h_q0_mantle_equiv = q0 / (obj.g * obj.rho_m)
    compensation.append(-1 * np.mean(w) / np.mean(h_q0_mantle_equiv))
  """
  
  # Finalize
  obj.finalize()

  cdiff = np.diff(compensation)

  #plot(dx_range,compensation)
  #show()

  # cdiff[0] happens at x[.5]
  # calculating the maximum negative lambda, so ... actually make it good later
  # just test for now
  try:
    minvalid_index = (cdiff<0).nonzero()[0].max()
    # WRONG - mixing dx variability with ncells for a given dx
    # ORIGINAL INDEX RIGHT BELOW THIS COMMENT:
    #lambda_min.append(minvalid_index / float(ncells)) # b/c ncells = 1 wavelength of Te variability # ORIGINAL
    #lambda_min.append(minvalid_index * dxi)
    # RIGHT: minimum valid wavelength of dx changes
    lambda_min.append(dx[minvalid_index] * ncells / 1000.) # wavelength [km]
    # RIGHT: minimum valid dx normalized to whole flexural wavelength = ncells
    #lambda_min.append(dx[minvalid_index]
    teslope_max.append(np.max(np.diff(Te)) / dx[minvalid_index])
  except:
    # If there is no clear minumum
    lambda_min.append(np.nan)
    teslope_max.append(np.nan)

# Asymptote of hyperbola

#dxmean = (dx[:-1] + dx[1:]) / 2

from matplotlib.pyplot import *

# Get rid of first element - doesn't follow trend and doesn't make sense - 
# possibly a fluke of my analysis
maxmin_ratio_lambdaf = np.array(maxmin_ratio_lambdaf[1:])
maxmin_ratio_Te = np.array(maxmin_ratio_Te[1:])
lambda_min = np.array(lambda_min[1:])
teslope_max = np.array(teslope_max[1:])

# Get rid of NaNs
maxmin_ratio_lambdaf = maxmin_ratio_lambdaf[np.isnan(lambda_min) == 0]
maxmin_ratio_Te = maxmin_ratio_Te[np.isnan(lambda_min) == 0]
teslope_max = teslope_max[np.isnan(lambda_min) == 0]
lambda_min = lambda_min[np.isnan(lambda_min) == 0]


"""
figure(1)
plot(dx, compensation,'ko-')
xlabel('dx')
ylabel('compensation')

figure(2)
plot(maxmin_ratio_Te, lambda_min,'ko-')
xlabel('Maximum Te / Minimum Te')
#ylabel('Min cell size as as fraction\nof Te variability wavelength')
ylabel('Minimum wavelength of Te variability to be solvable [km]')
"""
figure(3)
plot(maxmin_ratio_lambdaf, lambda_min,'ko-')
xlabel('Maximum Flexural Wavelength / Minimum Flexural Wavelength')
#ylabel('Min cell size as as fraction\nof Te variability wavelength')
ylabel('Minimum wavelength of Te variability to be solvable [km]')

"""
figure(4)
plot(maxmin_ratio_lambdaf, teslope_max,'ko-')
xlabel('Maximum Flexural Wavelength / Minimum Flexural Wavelength')
ylabel('Maximum solvable Te slope')
"""

figure(5)
plot(teslope_max, lambda_min, 'ko-')
xlabel('Maximum solvable Te slope')
ylabel('Minimum wavelength of Te variability to be solvable [km]')

show()

"""
plot(dxmean, cdiff, 'k-')
#ylim((-1,1))
xlabel('dx')
ylabel('compensation')
show()
"""

print maxmin_ratio_Te
print lambda_min

from scipy.optimize import curve_fit

def powerfunc(x,a,b,c):
  return a * (x-b)**c + d

#popt, pcov = curve_fit(sqrtfunc, np.array(lambda_min), maxmin_ratio, p0=[1, min(maxmin_ratio)+1)

"""
def powerfunc(x,a,b,c):
  return a * x**2 + c

popt, pcov = curve_fit(sqrtfunc, maxmin_ratio, np.array(lambda_min),p0=(1,.5,.2))

"""
#print obj.Te

#plot(Te); show()
