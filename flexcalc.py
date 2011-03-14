"""
Functions for finite-difference solutions to the thin-plate flexure equation. \\
\\
Written by Andrew D. Wickert. \\
Started 1 November 2010, Preliminary version finished 6 November 2010. \\
More edits on 7, 8 Nov. 2010; changed sign of "w" (flex down w/ positive load) \\
Based on MATLAB code started in late Spring / early Summer 2010 \\
\\
Originally called fvet = "Flexure: Variable Elastic Thickness"
Planned expansion to analytical solutions required new name: "flexcalc"
"""

from __future__ import division # No floor division here

def importascii(f,d=' '):
  """
  imported_array_varname = importarray(filename,delimiter=' ') \\
  \\
  Both filename and delimiter must be strings (grep commands work, \\
  e.g., use \t as delim for tab)
  """
  from numpy import loadtxt
  imported_array = loadtxt(fname=f,delimiter=d)
  return imported_array
  
def importsparse(filename):
  """
  varname = importsparse(filename): \\
  \\
  "filename" must be an object or string
  """
  from scipy.io import mmread
#  from scipy.sparse import lil_matrix
  a = mmread(filename)
#  a = lil_matrix(a)
  return a
  
def varprep2d(dx,dy,Te,E=1E11,nu=0.25,rho_m=3300,rho_fill=0):
  """
  dx4, dy4, dx2dy2, D, drho = varprep2d(dx,dy,Te,E=,nu=0.25,rho_m=3300,rho_fill=0) \\
  \\
  Defines the variables (except for the subset flexural rigidity) that are needed \\
  to run "coeff_matrix_2d"
  """
  dx4 = dx**4
  dy4 = dy**4
  dx2dy2 = dx**2 * dy**2
  D = E*Te**3/(12*(1-nu**2))
  drho = rho_m - rho_fill
  output = dx4, dy4, dx2dy2, D, drho
  return output

def subset_2d(a):
  """
  b = subset_2d(a) \\
  \\
  Takes only the interior elements of the array (skipping first and last two). \\
  This should be used for the elastic thickness array, which should be imported \\
  as two elements wider on each edge as the load array (because the 4th-order \\
  finite difference requires this for a solution within the area of the load).
  """
  b = a[2:-2,2:-2]
  return b

def coeff_matrix_2d(D,D_subset,drho,dx4,dy4,dx2dy2,nu=0.25,g=9.8):
  """
  coeff = coeff_matrix_2d(D=D,D_subset,nu=0.25,drho=drho,g=9.8,dx4=dx4,dy4=dy4,dx2dy2=dx2dy2) \\
  \\
  Calculates the matrix of coefficients that is later used via a Thomas algorithm to compute \\
  the flexural response to the load. This step need only be performed once, and the coefficient \\
  matrix can very rapidly compute flexural solutions to any load. This makes this particularly \\
  good for probelms with time-variable loads or that require iteration (e.g., water loading, in \\
  which additional water causes subsidence, causes additional water detph, etc.).
  """

  from numpy import array, prod
  from scipy.sparse import lil_matrix
  
  coeff = lil_matrix((prod(D_subset.shape),prod(D_subset.shape)))
  
  row = 0
  
  c = lil_matrix(D.shape) # Define it here as well for matrix size calc.
  endlength = prod(c[2:-2,2:-2].shape) # Reshaped matrix size
  
  for i in range(2,D.shape[1]-2):
    for j in range(2,D.shape[0]-2):

      # Should add some rapid matrix creator in here for the case in which all or most
      # of the D's are the same, esp. to speed up benchmarking against analytical
      # solutions, which necessarily have constant D.

      c = lil_matrix(D.shape) # Coefficient matrix before reshaping into single row

      c[j+2,i] = (D[j,i] + 0.5*(D[j+1,i]-D[j-1,i])) / dy4
      c[j+1,i-1] = (2*D[j,i] \
          + 0.5*(-D[j,i+1]+D[j,i-1]+D[j+1,i]-D[j-1,i]) \
          + ((1-nu)/8)*(D[j+1,i+1]-D[j-1,i+1]-D[j+1,i-1]+D[j-1,i-1])) \
          / dx2dy2
      c[j+1,i] = (-6*D[j,i]+2*D[j-1,i])/dy4 \
          + nu*(D[j,i+1]-2*D[j,i]+D[j,i-1])/dx2dy2 \
          + (-D[j+1,i]-4*D[j,i]+D[j-1,i])/dx2dy2
      c[j+1,i+1] = (2*D[j,i] \
          + 0.5*(D[j,i+1]-D[j,i-1]+D[j+1,i]-D[j-1,i]) \
          + ((1-nu)/8)*(D[j+1,i+1]-D[j-1,i+1]-D[j+1,i-1]+D[j-1,i-1])) \
          / dx2dy2
      c[j,i-2] = (D[j,i] - 0.5*(D[j,i+1]-D[j,i-1])) / dx4
      c[j,i-1] = (2*D[j,i+1]-6*D[j,i])/dx4 \
          + nu*(D[j+1,i]-2*D[j,i]+D[j-1,i])/dx2dy2 \
          + (D[j,i+1]-4*D[j,i]-D[j,i-1])/dx2dy2
      c[j,i] = (-2*D[j,i+1]+10*D[j,i]-2*D[j,i-1])/dx4 \
          + (-2*D[j+1,i]+10*D[j,i]-2*D[j-1,i])/dy4 \
          + (8*D[j,i])/dx2dy2 \
          + nu*(-2*D[j,i+1]-2*D[j,i-1]+8*D[j,i]-2*D[j+1,i]-2*D[j-1,i])/dx2dy2 \
          + drho*g
      c[j,i+1] = (-6*D[j,i]+2*D[j,i-1])/dx4 \
          + nu*(D[j+1,i]-2*D[j,i]+D[j-1,i])/dx2dy2 \
          + (-D[j,i+1]-4*D[j,i]+D[j,i-1])/dx2dy2
      c[j,i+2] = (D[j,i] + 0.5*(D[j,i+1]-D[j,i-1])) / dx4
      c[j-1,i-1] = (2*D[j,i] \
          + 0.5*(-D[j,i+1]+D[j,i-1]-D[j+1,i]+D[j-1,i]) \
          + ((1-nu)/8)*(D[j+1,i+1]-D[j-1,i+1]-D[j+1,i-1]+D[j-1,i-1])) \
          / dx2dy2
      c[j-1,i] = (2*D[j+1,i]-6*D[j,i])/dy4 \
          + nu*(D[j,i+1]-2*D[j,i]+D[j,i-1])/dx2dy2 \
          + (D[j+1,i]-4*D[j,i]-D[j-1,i])/dx2dy2
      c[j-1,i+1] = (2*D[j,i] \
          + 0.5*(D[j,i+1]-D[j,i-1]-D[j+1,i]+D[j-1,i]) \
          - ((1-nu)/8)*(D[j+1,i+1]-D[j-1,i+1]-D[j+1,i-1]+D[j-1,i-1])) \
          / dx2dy2
      c[j-2,i] = (D[j,i] - 0.5*(D[j+1,i]-D[j-1,i])) / dy4
      
      #ccol = (c[2:-2,2:-2],endlength)
      # This conversion is REALLY inefficient; how to reshape sparse matrices?
      ccol = c.todense()
      # ccol = reshape
      ccol = ccol[2:-2,2:-2].reshape(1,endlength)
      coeff[row,:] = lil_matrix(ccol)
      
      row += 1 # Increment
      
  return coeff
  
def coeff_matrix_1D(D,nu,drho,g,dx4,dx2):
  """
  1D pentadiagonal matrix to solve 1D flexure with variable elastic thickness \\
  via a Thomas algorithm. \\
  """
  
  # Based on and translated in part from "flexvar.m", written in 2001 by Bob
  # Anderson, which in turn is based on a code written by Laura Wallace
  
  from numpy import array, prod
  from scipy.sparse import lil_matrix
  
  coeff = lil_matrix((prod(D_subset.shape),prod(D_subset.shape)))
  
  row = 0
  
  c = lil_matrix(D.shape) # Define it here as well for matrix size calc.
  endlength = prod(c[2:-2,2:-2].shape) # Reshaped matrix size
  
  
  
  
def direct_fd_solve_2d(coeff,q0):
  """
  w = flex_calc(coeff,q0) \\
  \\
  Sparse Thomas algorithm-based flexural response calculation.
  Requires the coefficient matrix from "coeff_matrix_2d"
  """
  from scipy.sparse.linalg import spsolve
  from scipy.sparse import lil_matrix
  from numpy import prod
  q0t = q0.transpose() # Need to do this to make the "reshape" commands
                       # work correctly
  q0vector = q0t.reshape(prod(q0.shape),1)
  q0vector = lil_matrix(q0vector)
  wvector = spsolve(coeff,q0vector)
  w = -wvector.reshape(q0.shape)
  return w

def direct_fd_solve_1d(coeff,q0):
  """
  In the works.
  """

def adi():
  """
  Alternating direction implicit: not written yet
  """
  
def fft_1D():
  """
  1D fast Fourier transform: not written yet
  """

def fft_2D():
  """
  2D fast Fourier transform: not written yet
  """
  
def analytical_1D():
  """
  T&S solution, necessarily constant D
  Sum line loads
  """
  
def analytical_2D():
  """
  Kelvin-Bessel fcn. solution
  Sum point loads
  """

def savesparse(filename,varname):
  """
  savesparse(filename,varname) \\
  \\
  "filename" must be a string with a ".mtx" extension. \\
  Will save file as a Matrix Market file. \\
  Good ASCII I/O option for the sparse matrix.
  """
  from scipy.io import mmwrite
  mmwrite(filename,varname)
  
def saveascii(filename,varname):
  """
  saveascii(filename,varname) \\
  \\
  "filename" must be a string and requires an extension (e.g., ".txt") \\
  Will save file as a space-delimited ASCII text file. \\
  Good for easy universal readability of arrays that should be small (e.g., the \\
  output flexural response).
  """
  from numpy import savetxt
  savetxt(fname=filename,X=varname,delimiter=' ')

def implot(a,x_label='',y_label='',fig_title=''):
  """
  plot(a,x_label,y_label,fig_title) \\
  \\
  Plots an image-stype figure with automatically generates axis labels and titles, \\
  colorbars, and properly-scaled dimensions.
  """
  from matplotlib.pylab import imshow, colorbar, figure, axis, xlabel, ylabel, title
  figure(1)
  imshow(Te)
  colorbar()
  axis('image')
  xlabel(x_label)
  ylabel(y_label)
  title(fig_title)
  
def implotall(Te,w,q0,x_label='',y_label='',q0_title='Load [kg/m$^2$ = Pa]',w_title='Flexure [m]',Te_title='Elastic thickness [km]'):
  """
  plotall(Te,w,q0,x_label,y_label,q0_title,w_title,Te_title) \\
  \\
  Plots figures for Te, w, and q0. Automatically generates axis labels, titles, \\
  and colorbars, and causes the plots to have properly-scaled dimensions. \\
  Of course, you don't have to make these figures be of Te, w, and q0... but this \\
  is what the variable names will say that they are.
  """
  from matplotlib.pylab import imshow, colorbar, figure, axis, xlabel, ylabel, title, show
  figure(1)
  imshow(Te/1000)
  colorbar()
  axis('image')
  xlabel(x_label)
  ylabel(y_label)
  title(Te_title)
  figure(2)
  imshow(w)
  colorbar()
  axis('image')
  xlabel(x_label)
  ylabel(y_label)
  title(w_title)
  figure(3)
  imshow(q0)
  colorbar()
  axis('image')
  xlabel(x_label)
  ylabel(y_label)
  title(q0_title)
  show()
  
  
#def varprep1d(dx,Te=Te,E=1E11,nu=0.25,rho_m=3300,rho_fill=0):
#  dx4 = dx**4
  
  
#def pad2d_const_value(a):
#  """
#  apad = pad2d_cost_value(a)
#  """
#  from numpy import nan
#  from numpy import array
#  apad = array([[nan,nan,a[0,:],nan,nan],[nan,a[0,0],a[0,:],a[0,-1],nan], \
#  [a[:,0],a[:,0],a,a[:,-1],a[:,-1]],[nan,a[-1,0],a[-1,:],a[-1,-1],nan], \
#  [nan,nan,a[-1,:],nan,nan]])
#  return apad
  
#def pad2d_const_gradient(a):
  
#def pad1d_const_value(a):

#def pad2d_const_gradient(a):

