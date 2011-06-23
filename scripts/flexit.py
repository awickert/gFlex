# Frontend to flexure calculations with the flexcalc library (created by ADW, 1-9 Nov. 2010,
# based on MATLAB script written sometime around June
# 
# Started by ADW on 7 NOV 2010
# Added lots more options (rho, Te, nu) on 8-9 Nov.
# Changed "subset" to "subset

from flexure import flexcalc
import optparse


###########################
## PARSE TO COMMAND LINE ##
###########################

parser = optparse.OptionParser ()

# Flags
parser.add_option("-v", action="store_true", dest="verbose", help="If selected, turns on verbosity")
parser.add_option("-c", action="store_true", dest="calc_coeff", help="Runs the coefficient matrix creator")
parser.add_option("-r", action="store_true", dest="flex_response", help="Runs the flexural response calculator")
parser.add_option("-p", action="store_true", dest="plot_input_and_response", help="Plots Te, q0, and w (and therefore requires both Te and q0 to be defined in the input parameters)")

# Make coefficient matrix
parser.add_option ("--dx", action="store", type="float", dest="dx", help="Sets dx [m] (dy = dx).")
parser.add_option ("--Te", action="store", type="string", dest="Te", help="Filename of array of elastic thicknesses to import (must be 2 cells wider on each side than the load matrix). Default =  %default", default="Te.txt")
parser.add_option ("--E", action="store", type="float", dest="E", help="Young's modulus [Pa]. Default = %default",  default = 1E11)
parser.add_option ("--rho_m", action="store", type="float", dest="rho_m", help="Mantle density [kg/m^3]. Default = %default",  default = 3300)
parser.add_option ("--rho_fill", action="store", type="float", dest="rho_fill", help="Density of material that fills depressions [kg/m^3]. Default = %default",  default = 0)
parser.add_option ("--nu", action="store", type="float", dest="nu", help="Poisson's ratio [-]. Default = %default",  default = 1E11)
parser.add_option ("--delimiter", action="store", type="string", dest="delimiter", help="Delimiter of arrays that are imported (default is space, and exported arrays are automatically space-delimited)", default=" ")
parser.add_option ("--coeff-out", action="store", type="string", dest="coeff_out", help="Filename of sparse coefficient matrix to output. Default =  %default", default="coefficient_matrix.mtx")

# Flexural response calculation
parser.add_option ("--coeff-file", action="store", type="string", dest="coeff", help="Filename of coefficient matrix for Thomas algorithm to import. Default =  %default", default="coefficient_matrix.mtx")
parser.add_option ("--q0", action="store", type="string", dest="q0", help="Array of loads to import. Default =  %default", default="q0.txt")
parser.add_option ("--w-out", action="store", type="string", dest="w_out", help="Filename of flexural response matrix to output (as ASCII). Default =  %default", default="w.txt")

# Other loading commands (for plotting)
parser.add_option ("--w", action="store", type="string", dest="w", help="Filename of flexural response matrix to import (i.e., for plotting). Default =  %default", default="w.txt")

# Plotting
#parser.add_option

(options, args) = parser.parse_args ()


#############
## EXECUTE ##
#############

###############################
## CREATE COEFFICIENT MATRIX ##
###############################

if options.calc_coeff==True:
#  if options.Te!=None:
#    if options.dx!=None:
    if options.verbose==True:
      print('Creating coefficient matrix')
    
    #################################
    ## DEFINE AND IMPORT VARIABLES ##
    #################################

    Te = flexcalc.importascii(options.Te)

    dx = options.dx
    dy = dx # My solution method requires that dy and dx are equal
            # Though I bet you could get away with small differences
            # (In fact, when I tried it, it looked convincing)
            # So probably it is just rigorously correct for dx=dy
            
    E = options.E
    rho_m = options.rho_m
    rho_fill = options.rho_fill
    nu = options.nu

    # Optional variables to define: defaults given here:
    # E = 1E11,rho_m = 3300,rho_fill=0,nu=0.25
    dx4, dy4, dx2dy2, D, drho = flexcalc.varprep2d(dx=dx,dy=dy,Te=Te,E=E,rho_m=rho_m,rho_fill=rho_fill,nu=0.25)

    D_subset = flexcalc.subset_2d(D)

    ##################################
    ## CALCULATE COEFFICIENT MATRIX ##
    ##################################

    coeff = flexcalc.coeff_matrix_2d(D=D,D_subset=D_subset,drho=drho,dx4=dx4,dy4=dy4,dx2dy2=dx2dy2)

    #############################
    ## SAVE COEFFICIENT MATRIX ##
    #############################

    flexcalc.savesparse(options.coeff_out,coeff)

if options.flex_response==True:
  if options.verbose==True:
    print('Calculating flexural response')

  #################################
  ## DEFINE AND IMPORT VARIABLES ##
  #################################

  q0 = flexcalc.importascii(options.q0,options.delimiter) # Load matrix
  coeff = flexcalc.importsparse(options.coeff) # Coefficient matrix to solve
                                                      # Thomas algorithm

  #####################################
  ## SOLVE FOR FLEXURAL RESPONSE (w) ##
  #####################################

  w = flexcalc.direct_fd_solve_2d(coeff,q0)

  ################################
  ## SAVE FLEXURAL RESPONSE (w) ##
  ################################

  flexcalc.saveascii(options.w_out,w) 
  
##########
## PLOT ##
##########

if options.plot_input_and_response==True:

  Te = flexcalc.importascii(options.Te,options.delimiter)
  w = flexcalc.importascii(options.w,options.delimiter)
  q0 = flexcalc.importascii(options.q0,options.delimiter)
  flexcalc.implotall(Te,w,q0)
  
