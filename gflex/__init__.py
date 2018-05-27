import sys
import os
# If run from up one level
sys.path.append( os.path.dirname(os.path.realpath(__file__)) )
os.chdir( os.path.dirname(os.path.realpath(__file__)) )
from gflex import *
