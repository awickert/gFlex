import sys
import os
# If run from up one level
sys.path.append( os.path.dirname(os.path.realpath(__file__)) )
os.chdir( os.path.dirname(os.path.realpath(__file__)) )
print(os.getcwd())
from gflex import *
from f1d import *
from f2d import *
from base import *
