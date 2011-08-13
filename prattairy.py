from __future__ import division # No automatic floor division
from base import *


class PrattAiry(Isostasy):
  def initialize(self, filename):
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

