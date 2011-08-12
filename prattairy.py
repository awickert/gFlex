from __future__ import division # No automatic floor division
from base import *


class PrattAiry(Isostasy):
  def initialize(self, filename):
    super(PrattAiry, self).initialize(filename)
    if debug: print 'PrattAiry initialized'

  def run(self):
    self.noflexure ()
    if debug: print 'PrattAiry run'
    #self.imshow(self.w,1) # in here temporarily
    #self.imshow(self.q0/(self.rho_m * self.g),2)

  def finalize(self):
    if debug: print 'PrattAiry finalized'
    
    
  ######################################
  ## FUNCTIONS TO SOLVE THE EQUATIONS ##
  ######################################

  def noflexure(self):
    # q0 must be rho_load * g * h_load
    # rho --> Pratt, h --> Airy
    # Nothing to return: all calls to self
    self.w = - self.q0 / (self.rho_m * self.g)


  ##############
  ## PLOTTING ##
  ##############
  
  def plot(self):
    """
    Plot if you want to - for troubleshooting
    """
    from matplotlib.pyplot import plot, show, figure
    from numpy import arange
    figure()
    plot(arange(0,self.dx*self.w.shape[0],self.dx),self.w)
    show()

  def imshow(self,image,fignum):
    # Plot if you want to - for troubleshooting
    from matplotlib.pyplot import imshow, show, colorbar, figure
    figure(fignum)
    imshow(image,interpolation='nearest') #,interpolation='nearest'
    colorbar()
    show()
