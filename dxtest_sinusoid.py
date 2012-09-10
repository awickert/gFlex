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

# Plotting
from matplotlib.pyplot import *

# Start up grids to store compensation and set dx
compensation = []
dx = np.linspace(10,20000,100)
#dx = np.array([1000,3000])

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

# Sinusoidal flexural wavelengths that allow periodic b.c.'s:
amplitude = 30000 # 1/2 of total variation in Te
minTe = 5000
xi_periodic = np.arange(0,ncells+2) * 2 * np.pi / (ncells+2)
Te = amplitude * (np.cos(xi_periodic)) + amplitude + minTe

print Te.shape

q0 = 3300*9.8*100*np.ones(ncells)
#q0[:ncells/3] = 0
#q0[-ncells/3:] = 0

print q0.shape

obj.set_value('Loads', q0)
obj.set_value('ElasticThickness', Te)# np.ones(Te.shape)*20000)#


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

counter=0
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

  # PLOT
  #if counter % (len(dx)/20) == 0:
  if counter == 17:
    figure(1)
    plot(h_q0_mantle_equiv,'k--')
    plot(w, label=r'\lambda_{Te} = ' + str(ncells*dxi/(2*np.pi)/1000) + ' km')

  counter+=1

legend()
# Finalize
obj.finalize()

cdiff = np.diff(compensation)
dxmean = (dx[:-1] + dx[1:]) / 2

figure(2)
plot(dx, compensation,'ko-')
xlabel('dx')
ylabel('compensation')
show()

"""
plot(dxmean, cdiff, 'k-')
#ylim((-1,1))
xlabel('dx')
ylabel('compensation')
show()
"""

#print obj.Te

#plot(Te); show()
