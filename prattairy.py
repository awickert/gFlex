from __future__ import division # No automatic floor division
from base import *

"""
This code is here as a vestige of a former time.
It is really just a whole lot of architecture to solve a trivial equation.
It was somewhat useful to figure out how to build classes way back when...
This will be its final update before it is deleted sometime before the 1.0
release, or so I plan.
-- ADW, 19 January 2015
"""

class PrattAiry(Isostasy):
  def initialize(self, filename):
    super(PrattAiry, self).whichModel(filename)
    super(PrattAiry, self).initialize(filename)
    if debug: print 'PrattAiry initialized'

  def run(self):
    self.noflexure ()
    if debug: print 'PrattAiry run'

  def finalize(self):
    if debug: print 'PrattAiry finalized'
    
    
  ######################################
  ## FUNCTIONS TO SOLVE THE EQUATIONS ##
  ######################################

  def noflexure(self):
    """
    q0 must be rho_load * g * h_load
    rho --> Pratt, h --> Airy
    Nothing to return: all calls to self
    """
    self.w = - self.q0 / (self.rho_m * self.g)

