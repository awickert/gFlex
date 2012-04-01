#! /usr/bin/python

"""
Started 11 March 2012 as a GRASS interface for Flexure
"""

# FLEXURE
from base import *
from f1d import *
from f2d import *
from prattairy import *

# GRASS
import numpy as np
from grass.script import core as grass
from grass.script import mapcalc
from grass.script import db as db
import grass.script.array as garray
import time

# PARSER
import argparse

parser = argparse.ArgumentParser(description='Communicate between the user, GRASS, and Flexure')
parser.add_argument('q0', type=str, help='input load raster (q0 = rho * g * h)')
parser.add_argument('Te', type=str, help='elastic thickness raster [km]')
parser.add_argument('output', type=str, help='name of raster of deflections')
parser.add_argument('rho_fill', type=float, help='density of material that fills flexural depressions', default=0)
parser.add_argument('resolution', type=float, help='density of material that fills flexural depressions', default=grass.region()['nsres'])

args = parser.parse_args()

# Change grid resolution to solver resolution
grass.run_command('g.region', rast=args.Te) # Start out with some default resolution in case something gets messed up
grass.run_command('g.region', nsres=args.resolution, ewres=args.resolution)

grass.run_command('r.resamp.interp', input=args.q0, output='q0resamp', method='bicubic', overwrite=True)
grass.run_command('r.resamp.interp', input=args.Te, output='Teresamp', method='bicubic', overwrite=True)

# Do a better job of interpolating
#grass.run_command('r.resamp.rst', input=args.q0, elev='q0resamp', ns_res=args.resolution, ew_res=args.resolution, overwrite=True)
#grass.run_command('r.resamp.rst', input=args.Te, output='Teresamp', overwrite=True)

# Automatically decide that we are doing 2D finite difference flexural isostasy
# with a direct(?) solution method
obj = F2D()
obj.set_value('model', 'flexure')
obj.set_value('dimension', 2)
obj.set_value('GravAccel', 9.8)
obj.set_value('method', 'FD')
obj.set_value('Solver', 'direct')

# No filename: getter/setter interface isn't totally worked out
obj.filename = None

# Make a bunch of standard selections
obj.set_value('YoungsModulus', 65E9)#70E6/(600/3300.))#
obj.set_value('PoissonsRatio', 0.25)
obj.set_value('GravAccel', 9.8)
obj.set_value('MantleDensity', 3300)

# Set all boundary conditions to Mirror
obj.set_value('BoundaryCondition_East', 'Mirror')
obj.set_value('BoundaryCondition_West', 'Mirror')
obj.set_value('BoundaryCondition_North', 'NoOutsideLoads')
obj.set_value('BoundaryCondition_South', 'NoOutsideLoads')

# Get grid spacing from GRASS
obj.set_value('GridSpacing_x', grass.region()['ewres'])
obj.set_value('GridSpacing_y', grass.region()['nsres'])
#obj.set_value('GridSpacing_x', 50000)
#obj.set_value('GridSpacing_y', 50000)


# Get raster grids from GRASS
q0 = garray.array()
#q0.read(args.q0)
q0.read('q0resamp')
Te = garray.array()
#Te.read(args.Te)
Te.read('Teresamp')

# Change these grids into basic numpy arrays for easier use with Flexure
q0 = np.array(q0)
Te = np.array(Te * 1000) # *1000 for km-->m
#Te[Te>35000] = 35000
#Te[Te<30000] = 30000

#q0 = np.zeros(q0.shape)
#q0[q0.shape[0]/4:3*q0.shape[0]/4,:] = 1000

#sh = Te.shape
#gr = np.linspace(1000,60000,Te.shape[0])
#Te = np.lib.stride_tricks.as_strided(gr, (Te.shape[1], gr.size), (0, gr.itemsize)).transpose()

# Try to "guess" a known good low-res solution to imprpove higher-res solutions 

#wold = garray.array()
#wold.read('w')
#wold = np.array(wold)
#wold = wold[1:-1,1:-1]
#obj.wold = wold
#print obj.wold.shape
#print Te.shape
#print q0.shape

#Te /= 2 # Fudge for stability to see how cell size vs. Te works

from matplotlib.pyplot import *
#figure(1); imshow(q0, interpolation='nearest'); colorbar()
#figure(2); imshow(Te, interpolation='nearest'); colorbar()
#show()

# Values set by user
obj.set_value('Loads', q0[1:-1,1:-1]) # Te needs to be 1 cell bigger on each edge
obj.set_value('ElasticThickness', Te)# np.ones(Te.shape)*20000)#
obj.set_value('InfillMaterialDensity', args.rho_fill) # defaults to 0

# Calculated values
obj.drho = obj.rho_m - obj.rho_fill

# CALCULATE!
#obj.initialize(None)
obj.run()
obj.finalize()

# Write to GRASS
# First, shrink the region by 1 cell so it accepts the flexural solution
n = grass.region()['n'] - grass.region()['nsres']
s = grass.region()['s'] + grass.region()['nsres']
e = grass.region()['e'] - grass.region()['ewres']
w = grass.region()['w'] + grass.region()['ewres']
nrows = grass.region()['rows']-2
ncols = grass.region()['cols']-2
grass.run_command('g.region', n=n, s=s, w=w, e=e, rows=nrows, cols=ncols) 
# Then create a new garray buffer and write to it
outbuffer = garray.array() # Instantiate output buffer
outbuffer[...] = obj.w
outbuffer.write(args.output, overwrite=True) # Write it with the desired name
# And create a nice colormap!
grass.run_command('r.colors', map=args.output, color='rainbow')#, flags='e')
# Then revert to the old region
grass.run_command('g.region', n=n, s=s, w=w, e=e) 
n = grass.region()['n'] + grass.region()['ewres']
s = grass.region()['s'] - grass.region()['ewres']
e = grass.region()['e'] + grass.region()['ewres']
w = grass.region()['w'] - grass.region()['ewres']
grass.run_command('g.region', n=n, s=s, w=w, e=e)

# Finally, return to original resolution (overwrites previous region selection)
grass.run_command('g.region', rast=args.Te, flags='p')
grass.run_command('r.resamp.interp', input=args.output, output=args.output + '_interp', method='lanczos', overwrite=True)
grass.run_command('r.colors', map=args.output + '_interp', color='rainbow')#, flags='e')

imshow(obj.w, interpolation='nearest'), show()
