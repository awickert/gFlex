#! /usr/bin/env python

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

import os.path
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
  print 'Usage: gflex.py <<path_to_input_file>> [-h OR --help for more information]'
  print ""
  
def displayScriptInclusionInstructions():
  print ""
  print "USAGE NOTE FOR NO CONFIGURATION FILE"
  print "--------------------------------------------------------------------"
  print "No configuration file: to run entirely with getters and setters, it is"
  print "not posslible to simply run 'gflex.py'. Instead one must write a script"
  print "in Python or a compatible language that includes 'import gflex' and"
  print "then defines a flexure object (F1D or F2D), like:"
  print ""
  print "import gflex"
  print "model_object = gflex.F1D()"
  print "model_object.set_value(VALUE_KEY, VALUE)"
  print "#...more..."
  print ""

def main():
  # Choose how to instantiate
  if len(sys.argv) == 2:
    if sys.argv[1] == '--help' or sys.argv[1] == '-h':
      displayUsage()
      displayScriptInclusionInstructions()
      return
    else:
      # Looks like it wants to be an configuration file!
      filename = sys.argv[1] # it works for usage (1) and (2)
      # Let's see if there is a file there
      try:
        obj = WhichModel(filename)
      except:
        displayUsage()
        if os.path.isfile(filename):
          print ">>>> Error: configuration file contains an error <<<<"
          sys.exit("")
        else:
          print ">>>> Error: can't locate specified configuration file. <<<<"
          sys.exit("")

  elif len(sys.argv) == 1:
    displayUsage()
    displayScriptInclusionInstructions()
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

  obj.initialize(filename)
  
  if obj.Debug: print 'Command line:', sys.argv

  ############################################
  ##       SET MODEL PARAMETERS HERE        ##
  ## (if not defined in configuration file) ##
  ############################################
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
