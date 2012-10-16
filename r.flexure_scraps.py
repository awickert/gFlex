
"""
Can't implement this quite right (i.e. to pull default resolution)
#%option
#%  key: resolution
#%  type: string
#%  description: Resolution of flexural calculations [m]
#%  answer: resolution of Te grid
#%  required: no
#%end
"""

# Too long text:
# Te: Raster map of elastic thickness [km] - constant or subtle variability, or else the code will produce misleading results
# l: Allows the code to be run with geographic (lat/lon) coordinates with the assumption that there are 111 km betwen lines of latitude and longitude (appropriate only for near the equator)
# Interface to Andy Wickert's lithospheric flexure model (see http://csdms.colorado.edu/wiki/Model:Flexure)

  """
  # Change grid resolution to solver resolution
  grass.run_command('g.region', rast=Te) # Start out with some default resolution in case something gets messed up
  grass.run_command('g.region', nsres=resolution, ewres=resolution)
  """

  # Now the real code

  grass.run_command('r.resamp.interp', input=q, output='q0resamp', method='bicubic', overwrite=True, quiet=True)
  grass.run_command('r.resamp.interp', input=Te, output='Teresamp', method='bicubic', overwrite=True, quiet=True)

  # Do this just in case we lose some area in the Te grid
  grass.run_command('g.region', rast='Teresamp') # Start out with some default resolution in case something gets messed up


