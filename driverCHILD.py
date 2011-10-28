#! /usr/local/python/bin/python

# Import 2D stuff and dependencies from Flexure, so we can talk to it.
# Unnecessary imports commented out
#import sys
#import sys
from base import *
#from f1d import *
from f2d import *
#from prattairy import *

from numpy import vstack

# First, we initialize

# CHILD
child_port.initialize( "/home/beach/faculty/gtucker/Runs/Biggertest/biggertest.in --silent-mode" )

# Flexure
flexobj = Isostasy()
# Which solution type: this will not change so long as we are using constant Te with CHILD.
# But I am setting these values and running the flow control just in case we do want to
# eventually.
flexobj.set_value('model','flexure')
flexobj.set_value('dimension',2)
flexobj.set_value('method','SPA_NG')
if flexobj.model == 'flexure':
  if flexobj.dimension == 1:
    flexobj = F1D()
  elif flexobj.dimension == 2:
    flexobj = F2D()
elif flexobj.model == 'PrattAiry':
  flexobj = PrattAiry()
# Variables
flexobj.initialize() # INITIALIZE, allowing filename to default to None
flexobj.set_value('YoungsModulus',1.0E11)
flexobj.set_value('PoissonsRatio',0.25)
flexobj.set_value('GravAccel',9.8)
flexobj.set_value('MantleDensity',3300)
flexobj.set_value('InfillMaterialDensity',0)
flexobj.set_value('ElasticThickness',10000)
# Now initialized except for loading

# set up, get results of initial loads
x = child_port.get_value_set('XCOORDS')
y = child_port.get_value_set('XCOORDS')
loads = child_port.get_value_set('LOADS')
q0 = vstack(x,y,loads)
flexobj.set_value( 'Loads', q0 )
flexobj.run()
w = flexobj.get_value( 'Deflection' )
w0 = w.copy()
flexobj.finalize()

# set run parameters
number_of_iterations = 1000
dt = 1000.0

# main loop
for i in range( 0, number_of_iterations ):

  print '***** ITERATION',i+1,'*****'

  # run child and get its new loads
  child_port.run( dt )
  print 'Done with run'
  x = child_port. get_value_set('XCOORDS')
  y = child_port. get_value_set('XCOORDS')
  loads = child_port.get_value_set('LOADS')
  q0 = vstack(x,y,loads)

  # set loads in flex and re-calculate
  flexobj.set_value( 'Loads', q0 )
  flexobj.run()
  w = flexobj.get_value( 'Deflection' )

  # calculate the net deflection and send to child
  dw = w - w0
  w0 = w.copy()
  child_port.adjust_elevations( dw )

# clean up
child_port.finalize()
flexobj.finalize()
print 'Done with go!'
