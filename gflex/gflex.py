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
  print "***********************************************"
  print "*** WELCOME to the gFlex development branch ***"
  print "***********************************************"
  print ""
  print "Open-source licensed under GNU GPL v3"
  print ""
  print 'Usage: gflex.py path_to_input_file'
  print 'Other arguments: --help or -h: this menu'
  print '     Other options for running gflex include:'
  print '      * Typing "import gflex" in a Python script to run with getters and setters'
  print '      * Using the GRASS GIS interfaces, r.flexure and v.flexure'

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

if __name__ == '__main__':
  main()
