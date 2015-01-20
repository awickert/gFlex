#! /usr/local/python/bin/python

# isostasy.py
# Written by Andrew Wickert with help from Greg Tucker and Eric Hutton
# This is the standalone driver for the flexural istostasy program
# 2010-2011

"""
Solves flexural and non-flexural (Airy & Pratt) isostasy.

Non-flexural isostasy solved analytically on 1D and 2D grids.

Flexural isostasy is solved both analytically (for constant flexural rigidity)
and numerically (for either variable or constant flexural rigidity); the former
uses Green's functions and the latter uses a direct sparse matrix solution.
"""

from socket import gethostname
import sys
from base import *
from f1d import *
from f2d import *

def displayUsage():
  print ""
  print "Command line usage help:"
  print 'Usage: (1) python isostasy.py path_to_input_file'
  print '       (2) ./isostasy.py path_to_input_file'
  print 'Other arguments: (3) --help or -h: this menu'
  print '                 (4) --getset: force to run with only getters and setters'
  print 'If no arguments are provided:'
  print '     If on Beach, will run with getters and setters'
  print '     Otherwise, this menu will appear'

def main():
  # Choose how to instantiate
  if len(sys.argv) == 2:
    if sys.argv[1] == '--help' or sys.argv[1] == '-h':
      displayUsage()
      return
    elif sys.argv[1] == '--getset':
      print ""
      print "No input file: running entirely with getters and setters."
      filename = None
    else:
      # Looks like it wants to be an input file!
      filename = sys.argv[1] # it works for usage (1) and (2)
      obj = WhichModel(filename)
  elif len(sys.argv) == 1:
    print ""
    print "No input file: running entirely with getters and setters."
    filename = None
    print ""
    if gethostname()=='beach':
      print ""
      print "Welcome to Beach!"
      print ""
    else:
      print "You are not running on Beach; are you sure you want to do this?"
      print 'Add "--getset" as an argument (python isostasy.py --getset)'
      print 'when you run isostasy to confirm that you did not just forget'
      print 'to set an input file'
      displayUsage()
      print ""
      sys.exit()
  else:
    print "Too many input parameters provided; exiting."
    displayUsage()
    sys.exit()
  
  ## SET MODEL TYPE AND DIMENSIONS HERE ##
  ########################################
  if obj.dimension == 1:
    obj = F1D(filename)
  elif obj.dimension == 2:
    obj = F2D(filename)

  obj.initialize(filename) # Does nothing
  
  if obj.Debug: print 'Command line:', sys.argv

  ####################################
  ##   SET MODEL PARAMETERS HERE    ##
  ## (if not defined in input file) ##
  #################################### 
  # obj.set_value('method','FD') # for example

  obj.run()
  obj.finalize()

  obj.output() # Not part of IRF or BMI: Does standalone plotting and file output

  #####################
  ## GET VALUES HERE ##
  ##   (if desired)  ##
  ##################### 
  #wout = obj.get_value('Deflection') # for example


def supercomputer_or_standalone():
# Not used right now; holding onto this if it is useful for the BMI/CMI
# implementation
  """
  Find if the program is running on Beach (and therefore we need to load the 
  CSDMS utility modules)
  """
  if gethostname()=='beach':
    import os.system
    if obj.Verbose: print "Running on Beach; loading CSDMS utilities"
    os.system("module load internal")
    # CSDMS_base used to manage getters/setters to avoid typing issues w/ bocca
    from csdms_utils.CSDMS_base import CSDMS_component as cc
  else:
    if obj.Verbose: print "Running in standalone mode"
    from CSDMS_base import CSDMS_component as cc



if __name__ == '__main__':
  main()
