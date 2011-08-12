# Using FFT for flexure
# See p.178 (start), Watts Flexure and Isostasy of the Lithosphere
# Figure out how density should figure in

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

def next_pwr_of_2(n):
  """
  Finds next higher power of 2 by finding the highest nonzero bit, moving it 
  up by "1", and setting all lower bits to 0.
  """
  Ln = len(bin(n))-2 # -2 b/c it counts the "0b"
  result = 1 << Ln
  return result
  

h0 = np.zeros(100)
h0[20:80] = 100

dx = 2000 # meters
L = h0.shape[0]*dx # Length of domain
k = np.arange(1/L,1/dx+1/L,1/L) # wavenumber range: domain to cell
NFFT = next_pwr_of_2(L)
x = np.arange(0,L)*dx # x values
ks = 1/dx # sample wavenumber
#k = ks/2*np.linspace(0,1,L.shape[0]/2+1)
#k = ks/2*np.linspace(0,1,NFFT/2+1);

D = 1E22 # Flexural rigidity
rho_m = 3300
rho_fill = 0
g = 9.8

rho_c = 2700

Phi = ( (D*k**4) / ( (rho_m - rho_fill) * g) + 1 )**-1 # FlexResponseFcn
# Fourier transformed result
W = rho_m / (rho_m - rho_fill) * np.fft.fft(h0) * Phi
w = np.fft.ifft(W)

plt.plot(w)
plt.show()

