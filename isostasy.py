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
from prattairy import *

# Not used right now; holding onto this if it is useful for the BMI/CMI
# implementation
def supercomputer_or_standalone():
  """
  Find if the program is running on Beach (and therefore we need to load the 
  CSDMS utility modules)
  """
  if gethostname()=='beach':
    import os.system
    if debug: print "Running on Beach; loading CSDMS utilities"
    os.system("module load internal")
    # CSDMS_base used to manage getters/setters to avoid typing issues w/ bocca
    from csdms_utils.CSDMS_base import CSDMS_component as cc
  else:
    if debug: print "Running in standalone mode"
    from CSDMS_base import CSDMS_component as cc

def displayUsage():
  print 'Usage: (1) python isostasy.py path_to_input_file'
  print '       (2) ./isostasy.py path_to_input_file'
  print 'Other arguments: (3) --help or -h: this menu'
  print '                 (4) --getset: force to run with only getters and setters'
  print 'If no arguments are provided:'
  print '     If on Beach, will run with getters and setters'
  print '     Otherwise, this menu will appear'

def main():
  infile = 1 # start by assuming that there is an input file

  print ""

  if len(sys.argv) > 1:
    if sys.argv[1] == '--help' or sys.argv[1] == '-h':
      displayUsage()
      return
    elif sys.argv[1] == '--getset':
      print "No input file: running entirely with getters and setters."
    else:
      displayUsage()
  elif len(sys.argv) == 1:
    infile = None
    print "No input file: running entirely with getters and setters."
    print ""
    if gethostname()=='beach':
      print "Welcome to Beach!"
      print ""
    else:
      print "You are not running on Beach; are you sure you want to do this?"
      print 'Add "--getset" as an argument (python isostasy.py --getset)'
      print 'when you run isostasy to confirm that you did not just forget'
      print 'to set an input file'
      sys.exit()      
  
  if debug: print 'Command line: ',sys.argv
  
  if infile:
    filename = sys.argv[1] # it works for usage (1) and (2)
  else:
    filename = None
  obj = Isostasy()
  obj.whichModel(filename)
  ## SET MODEL TYPE AND DIMENSIONS HERE ##
  ########################################
  if obj.model == 'flexure':
    if obj.dimension == 1:
      obj = F1D()
    elif obj.dimension == 2:
      obj = F2D()
  elif obj.model == 'PrattAiry':
    obj = PrattAiry()

  obj.initialize(filename)
  ## SET ALL OTHER MODEL PARAMETERS HERE ##
  # obj.set_value('method','FD')
  #########################################
  obj.run()
  obj.output() # Not part of IRF: Does standalone plotting and file output
  ## GET VALUES HERE ##
  #wout = obj.get_value('Deflection')
  #print wout
  #####################
  obj.finalize()

if __name__ == '__main__':
  main()
