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


import sys
from base import *
from f1d import *
from f2d import *
from prattairy import *

def supercomputer_or_standalone():
  """
  Find if the program is running on Beach (and therefore we need to load the 
  CSDMS utility modules)
  """
  from socket import gethostname
  if gethostname()==beach:
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

def main():
  infile = 1 # start by assuming that there is an input file

  if len(sys.argv) > 1:
    if sys.argv[1] == '--help' or sys.argv[1] == '-h':
      displayUsage()
      return
  elif len(sys.argv) == 1:
    infile = None
    if debug:
      print "No input file: running entirely with getters and setters."
      print 'Comment out the "return # TEMPORARY" immediately after this'
      print 'in "isostasy.py" to start working on the get/set routine.'
    #return # TEMPORARY
  
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
