#! /usr/bin/env python

"""
This file is part of gFlex.
gFlex computes lithospheric flexural isostasy with heterogeneous rigidity
Copyright (C) 2010-2018 Andrew D. Wickert

gFlex is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

gFlex is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with gFlex.  If not, see <http://www.gnu.org/licenses/>.
"""

from __future__ import print_function
import os.path
import sys
from base import *
from f1d import *
from f2d import *

try:
  from _version import __version__
except:
  __version__ = '? version file missing ?'

def welcome():
  print("")
  print("**************************"+"*"*len(__version__))
  print("*** WELCOME to gFlex v"+__version__+" ***")
  print("**************************"+"*"*len(__version__))
  print("")

def displayUsage():
  print("Open-source licensed under GNU GPL v3")
  print("")
  print("Usage:")
  print("gflex <<path_to_configuration_file>>  # TO RUN STANDALONE")
  print("gflex -h  *OR*  gflex --help          # DISPLAY ADDITIONAL HELP")
  print("gflex -v  *OR*  gflex --version       # DISPLAY VERSION NUMBER")
  print("import gflex                          # WITHIN PYTHON SHELL OR SCRIPT")
  print("")
  
def furtherHelp():
  print("")
  print("ADDITIONAL HELP:")
  print("--------------- ")
  print("")
  print('To generate an input file, please see the examples in the "input"')
  print("directory of this install.")
  print("")
  print("To run in a Python script or shell, follow this general pattern:")
  print("import gflex")
  print("flex = gflex.F1D()")
  print("flex.method = ...")
  print("# ...more variable setting...")
  print("# see the 'input' directory for examples")
  print("")

def main():
  # Choose how to instantiate
  if len(sys.argv) == 2:
    if sys.argv[1] == '--help' or sys.argv[1] == '-h':
      welcome()
      displayUsage()
      furtherHelp()
      return
    if sys.argv[1] == '--version' or sys.argv[1] == '-v':
      print("gFlex v"+__version__)
      return
    else:
      # Looks like it wants to be an configuration file!
      filename = sys.argv[1] # it works for usage (1) and (2)
      obj = WhichModel(filename)

  elif len(sys.argv) == 1:
    welcome()
    displayUsage()
    print("")
    sys.exit()
  else:
    welcome()
    print(">>>> ERROR: Too many input parameters provided; exiting. <<<<")
    print("")
    displayUsage()
    print("")
    sys.exit()
  
  
  ########################################
  ## SET MODEL TYPE AND DIMENSIONS HERE ##
  ########################################
  
  if obj.dimension == 1:
    obj = F1D(filename)
  elif obj.dimension == 2:
    obj = F2D(filename)

  obj.initialize(filename)
  
  if obj.Debug: print("Command line:", sys.argv)


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

