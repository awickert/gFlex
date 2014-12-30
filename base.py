import sys, ConfigParser, os
import numpy as np
import time # For efficiency counting
import types # For flow control
#import CSDMS_base

# BMI Interface
# edited by Andy
class BMI (object):
  # Replaced "file" with "filename", since "file" should be reserved for
  # the Python object
  def initialize (self, filename):
      pass
  def update (self):
    self.run()
  def finalize (self):
    pass
  def run_model (self):
    self.run()
  def run (self):
    # Added by Andy: self.update() isn't so meaningful becasue this is a
    # 1-time-step sort of model... unless I add some sort of mantle 
    # response code
    pass

  # NONE OF THESE ARE IMPLEMENTED
  def get_var_type (self, var):
      pass
  def get_var_units (self, var):
      pass
  def get_var_rank (self, var):
      pass
  def get_var_name (self, var):
      pass

  """
  # JUST CAUSING PROBLEMS BY PREVENTING CODE FROM RUNNING NORMALLY
  # EVEN IF TEMPLATE SEEMED LIKE A GOOD IDEA WHEN BEICHUAN INTRODUCED IT TO ME
  def get_value (self, long_var_name):
      pass
  # replaced "src" with "value" to match my code as it stands
  def set_value (self, long_var_name, value):
      pass
  """

  # NONE OF THESE ARE IMPLEMENTED
  def get_component_name (self):
      pass
  def get_input_var_names (self):
      pass
  def get_output_var_names (self):
      pass

  # NONE OF THESE ARE IMPLEMENTED
  def get_grid_dimen (self, long_var_name):
      pass
  def get_grid_res (self, long_var_name):
      pass
  def get_grid_corner (self, long_var_name):
      pass

  # NOT IMPLEMENTED
  def is_raster_grid (self, long_var_name):
      pass

class Utility(object):
  """
  Generic utility functions
  """
  def configGet(self,vartype,category,name,optional=False,specialReturnMessage=None):
    """
    Wraps a try / except and a check for self.filename around ConfigParser
    as it talks to the input file.
    Also, checks for existence of input file so this won't execute (and fail) 
    when no input file is provided (e.g., running in coupled mode with CSDMS 
    entirely with getters and setters)

    vartype can be 'float', 'str' or 'string' (str and string are the same),
    or 'int' or 'integer' (also the same).
    
    "Optional" determines whether or not the program will exit if the variable 
    fails to load. Set it to "True" if you don't want it to exit. In this case,
    the variable will be set to "None". Otherwise, it defaults to "False".
    
    "specialReturnMessage" is something that you would like to add at the end 
    of a failure to execute message. By default it does not print.
    """
    
    try:
      if vartype == 'float':
        var = self.config.getfloat(category,name)
      elif vartype == 'string' or vartype == 'str':
        var = self.config.get(category,name)
        if var == ""  and optional == False:
          print "Input strings cannot be empty unless they are optional."
          print name, "is not optional."
          print "Program crash likely to occur."
      elif vartype == 'integer' or vartype == 'int':
        var = self.config.getint(category,name)
      elif vartype == 'boolean' or vartype == 'bool':
        var = self.config.getboolean(category,name)
      else:
        print "Please enter 'float', 'string' (or 'str'), 'integer' (or 'int'), or 'boolean (or 'bool') for vartype"
        sys.exit() # Won't exit, but will lead to exception
      return var
    except:
      if optional:
        # Carry on if the variable is optional
        var = None
        if self.Verbose or self.Debug:
          print ""
          print 'No value entered for optional parameter "' + name + '"'
          print 'in category "' + category + '" in input file.'
          print 'No action related to this optional parameter will be taken.'
          print ""
      else:
        print 'Problem loading ' + vartype + ' "' + name + '" in category "' + category + '" from input file.'
        if specialReturnMessage:
          print specialReturnMessage
        sys.exit("Exiting.")

  def set_value(self, value_key, value):
    """
    Universal setter
    """
    # FIRST, VALUES THAT EVERYONE NEEDS::   
    
    # [mode]
    # Model type
    if value_key == 'model':
      self.model = value
    elif value_key =='dimension':
      self.dimension = value

    # [parameter]
    elif value_key == 'GravAccel':
      self.g = value
    elif value_key == 'MantleDensity':
      self.rho_m = value
    elif value_key == 'InfillMaterialDensity':
      self.rho_fill = value

    # [input]
    # Loading grid
    elif value_key == 'Loads':
      self.q0 = value

    # [numerical]
    # Grid spacing
    elif value_key == 'GridSpacing_x':
      self.dx = value
    elif value_key == 'GridSpacing_y':
      if self.dimension == 1:
        print "No dy in 1D problems; doing nothing"
      else:
        self.dy = value
    # Boundary conditions
    # "Dirichlet0" - 0-displacement at edges) # TO DO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # "0Slope0Shear" - First and third derivatives are 0: not so physical, but basically means that dw/dx_i at boundaries is flat and at a value that is not externally imposed
    # "0Moment0Shear" - second and third derivatives are 0: free cantilever edge (same as Stewart and Watts (1997) used, and standard in CivE)
    # "Periodic" - wraparound on edges
    # "Mirror" - reflects on edges
    elif value_key == 'BoundaryCondition_East':
      self.BC_E = value
    elif value_key == 'BoundaryCondition_West':
      self.BC_W = value
    elif value_key == 'BoundaryCondition_North':
      self.BC_N = value
    elif value_key == 'BoundaryCondition_South':
      self.BC_S = value

    # [verbosity]
    elif value_key == 'Verbosity' or value_key == 'Verbose':
      if self.Quiet:
        self.Verbose = False
      else:
        self.Verbose = value
    elif value_key == 'Debug':
      if self.Quiet:
        self.Debug = False
      else:
        self.Debug = value
    elif value_key == 'Quiet':
      self.Quiet = True
      self.Debug = False
      self.Verbose = False

    # Output
    elif value_key == 'DeflectionOut':
      # Output file name (string)
      self.wOutFile = value
    elif value_key == 'Plot':
      # 'q0', 'w', 'both', or (1D) 'combo'
      self.plotChoice = value
    
    # THEN, LITHOSPHERIC ELASTIC PROPERTIES AND SOLUTION METHOD
    # FOR FLEXURAL ISOSTASY

    # [Mode]
    # The lowercase version is here from earlier work; should phase it out
    elif value_key == 'method' or value_key == 'Method':
      self.method = value
    elif value_key == 'PlateSolutionType':
      self.PlateSolutionType = value

    # [Parameters]
    elif value_key == 'YoungsModulus':
      self.E  = value
    elif value_key == 'PoissonsRatio':
      self.nu = value
    elif value_key == 'ElasticThickness':
      self.Te = value
      
    # [Input]
    elif value_key == 'CoeffArray':
      # This coefficient array is what is used with the UMFPACK direct solver
      # or the iterative solver
      self.coeff_matrix = CoeffArray
      self.readyCoeff() # See if this is a sparse or something that needs to be loaded
      coeffArraySizeCheck() # Make sure that array size is all right
      # if so, let everyone know
      if self.Verbose:
        print "LOADING COEFFICIENT ARRAY"
        print "Elastic thickness maps will not be used for the solution."

    # [Numerical]
    elif value_key == 'Solver':
      self.solver = value # Direct or iterative
    elif value_key == 'ConvergenceTolerance':
      self.iterative_ConvergenceTolerance = value
    # None of the above?
    else:
      sys.exit('Error setting value: "'+value_key+'"')
    # Inherit from higher-level setter, if not one of these
    # Currently not inheriting from the whichModel() setter, as I
    # figure that has to be done right away or not at all
    
  def readyCoeff(self):
    from scipy import sparse
    if sparse.issparse(self.coeff_matrix):
      pass # Good type
    # Otherwise, try to load from file
    elif type(self.coeff_matrix) is types.StringType:
      pass
      print "Loading sparse coefficient arrays is not yet implemented."
      print "This must be done soon."
      print "Exiting."
      sys.exit()
    else:
      try:
        self.coeff_matrix = sparse.dia_matrix(self.coeff_matrix)
      except:
        "Failed to make a sparse array or load a sparse matrix from the input."
  
  # UNIVERSAL GETTER
  def get_value(self, val_string):
    if val_string=='Deflection':
      # This is the primary model output
      return self.w
    elif val_string=='CoeffMatrix':
      # This is to hold onto the coefficient matrix in memory so it doesn't 
      # need to be reloaded or recreated
      return self.coeff_matrix
    # Times to calculate efficiency of solution method
    elif val_string=='CoeffMatrixCreationTime':
      # Time to generate the sparse coefficient / operator matrix
      return self.coeff_creation_time
    elif val_string=='SolverTime':
      # Amount of time taken by the solver (direct or iterative)
      return self.time_to_solve
    if val_string=='Loads':
      # This is the model input
      return self.q0
    if val_string=='ElasticThickness':
      # This is the model input
      return self.Te
    if val_string=='Verbosity' or val_string=='Verbose':
      return self.Verbose


class Plotting(object):
  # Plot, if desired
  def plotting(self):
    # try:
    #   self.plotChoice
    # except:
    #   self.plotChoice = None
    if self.plotChoice:
      if self.Verbose: print "Starting to plot " + self.plotChoice
      if self.dimension==1:
        if self.plotChoice == 'q0':
          self.lineplot(self.q0/(self.rho_m*self.g),
            'Load thickness, mantle equivalent [m]')
        elif self.plotChoice == 'w':
          self.lineplot(self.w,'Deflection [m]')
        elif self.plotChoice == 'both':
          self.linesubplots()
        elif self.plotChoice == 'combo':
          self.plotTogether()
        else:
          if self.Quiet == False:
            print 'Incorrect plotChoice input, "' + self.plotChoice + '" provided.'
            print "Possible input strings are: q0, w, both, and (for 1D) combo"
            print "Unable to produce plot."
      elif self.dimension==2:
        if self.plotChoice == 'q0':
          self.surfplot(self.q0/(self.rho_m*self.g),
            'Load thickness, mantle equivalent [m]')
        elif self.plotChoice == 'w':
          self.surfplot(self.w,'Deflection [m]')
        elif self.plotChoice == 'both':
          self.surfsubplots()
        else:
          if self.Quiet == False:
            print 'Incorrect plotChoice input, "' + self.plotChoice + '" provided.'
            print "Possible input strings are: q0, w, both, and (for 1D) combo"
            print "Unable to produce plot."

  def linesubplots(self,figNum=1):
    from matplotlib.pyplot import plot, show, figure, subplot, xlabel, \
                                  ylabel, title
    
    figure(figNum)

    subplot(211)
    title('Loads and Lithospheric Deflections',fontsize=20)
    if self.method == "SAS_NG":
      plot(self.x,self.q0/(self.rho_m*self.g),'o')
    else:
      plot(self.x,self.q0/(self.rho_m*self.g))
    ylabel('Load thickness, mantle equivalent [m]')

    subplot(212)
    if self.method == "SAS_NG":
      plot(self.x,self.w,'o')
    else:
      plot(self.x,self.w)
    xlabel('Distance along profile [m]',fontsize=16)
    ylabel('Deflection [m]')

    show()

  def lineplot(self,data,ytext,xtext='Distance along profile [m]',titletext='',
      fontsize=16,figNum=1):
    """
    Plot if you want to - for troubleshooting
    """
    from matplotlib.pyplot import plot, show, figure, xlabel, ylabel, title
    
    figure(figNum)

    plot(self.x,data)

    xlabel(xtext,fontsize=16)
    ylabel(ytext,fontsize=16)
    title(titletext,fontsize=16)
    
    show()

  def plotTogether(self,figNum=1,titletext='Loads and Lithospheric Deflections'):
    from matplotlib.pyplot import plot, show, figure, xlabel, \
                                  ylabel, title, axis, ylim, legend
    
    fig = figure(figNum,figsize=(10,6))
    ax = fig.add_subplot(1,1,1)
    
    xkm = self.x/1000

    # Plot undeflected load
    if self.method == "SAS_NG":
      ax.plot(xkm,self.q0/(self.rho_m*self.g),'go',linewidth=2,label="Load thickness [m mantle equivalent]")
    else:
      ax.plot(xkm,self.q0/(self.rho_m*self.g),'g--',linewidth=2,label="Load thickness [m mantle equivalent]")
    
    # Plot deflected load
    if self.method == "SAS_NG":
      ax.plot(xkm,self.q0/(self.rho_m*self.g) + self.w,'go-',linewidth=2,label="Load thickness [m mantle equivalent]")
    else:
      ax.plot(xkm,self.q0/(self.rho_m*self.g) + self.w,'g',linewidth=2,label="Deflection [m] + load thickness [m mantle equivalent]")

    # Plot deflection
    if self.method == "SAS_NG":
      ax.plot(xkm,self.w,'ko-',linewidth=2,label="Deflection [m mantle equivalent]")
    else:
      ax.plot(xkm,self.w,'k',linewidth=2,label="Deflection [m mantle equivalent]")
    
    # Set y min to equal to the absolute value maximum of y max and y min
    # (and therefore show isostasy better)
    yabsmax = max(abs(np.array(ylim())))
    ylim((-yabsmax,yabsmax))

    if self.method == "FD":
      if type(self.Te) is np.ndarray:
        if (self.Te != (self.Te).mean()).any():
          title(titletext,fontsize=20)       
        else:
          title(titletext + ', $T_e$ = ' + str((self.Te / 1000).mean()) + " km",fontsize=20)
      else:
        title(titletext + ', $T_e$ = ' + str(self.Te / 1000) + " km",fontsize=20)
    else:
      title(titletext + ', $T_e$ = ' + str(self.Te / 1000) + " km",fontsize=20)
      
    ylabel('Loads and flexural response [m]',fontsize=16)
    xlabel('Distance along profile [km]',fontsize=16)
    
    legend(loc=0,numpoints=1,fancybox=True)
    
    show()


  def surfplot(self,data,titletext,figNum=1):
    """
    Plot if you want to - for troubleshooting
    """
    from matplotlib.pyplot import imshow, show, figure, colorbar, title
    
    figure(figNum)

    imshow(data, extent=(0, self.dx/1000.*data.shape[0], self.dy/1000.*data.shape[1], 0)) #,interpolation='nearest'
    xlabel('x [km]')
    ylabel('y [km]')
    colorbar()

    title(titletext,fontsize=16)
    
    show()

  def surfsubplots(self,figNum=1):
    from matplotlib.pyplot import imshow, show, figure, subplot, xlabel, \
                                  ylabel, title, colorbar
    
    figure(figNum,figsize=(6,9))

    subplot(211)
    title('Load thickness, mantle equivalent [m]',fontsize=16)
    imshow(self.q0/(self.rho_m*self.g), extent=(0, self.dx/1000.*self.q0.shape[0], self.dy/1000.*self.q0.shape[1], 0))
    xlabel('x [km]')
    ylabel('y [km]')
    colorbar()

    subplot(212)
    title('Deflection [m]')
    imshow(self.w, extent=(0, self.dx/1000.*self.w.shape[0], self.dy/1000.*self.w.shape[1], 0))
    xlabel('x [km]')
    ylabel('y [km]')
    colorbar()

    show()


class WhichModel(Utility):
  def __init__(self, filename=None):
    self.filename = filename
    if self.filename:
      try:
        # only let this function imoprt things once
        self.whichModel_AlreadyRun
      except:
        # Open parser and get what kind of model
        self.config = ConfigParser.ConfigParser()
        try:
          self.config.read(filename)
          self.inpath = os.path.dirname(os.path.realpath(filename)) + '/'
          # Need to have these guys inside "try" to make sure it is set up OK
          # (at least for them)
          self.model     = self.configGet("string", "mode", "model")
          self.dimension = self.configGet("integer", "mode", "dimension")
          self.whichModel_AlreadyRun = True
        except:
          sys.exit("No input file at specified path, or input file configured incorrectly")

# class Isostasy inherits IRF interface, and it determines the simulation type
# by reading three parameters from input file, but it does not set up other
# parameters, which is the responsibility of derived concrete classes.
class Isostasy(BMI, Utility, Plotting):

  def __init__(self, filename=None):
    # 17 Nov 2014: Splitting out initialize from __init__ to allow space
    # to use getters and setters to define values
    
    # Use standard routine to pull out values
    # If no filename provided, will not initialize input file.
    self.filename = filename
    # Set default "quiet" to False, unless set by setter
    # This will not work with input file; is designed for GRASS GIS
    self.Quiet = False
    # And also set default verbosity
    self.Verbose = True
    self.Debug = False

  def initialize(self, filename=None):
    # Values from input file
    
    # Introduce model
    if self.Quiet == False:
      print "" # Blank line at start of run
      print "*********************************************"
      print "*** Initializing gFlex development branch ***"
      print "*********************************************"
      print ""
      print "Open-source licensed under GNU GPL v3"
      print ""

    # If a filename is provided here, overwrite any prior value
    if filename:
      if self.filename:
        pass # Don't overwrite if filename is None-type
        # "Debug" not yet defined.
        #if self.Debug:
        #  print "Overwriting filename from '__init__' step with that from\n"+\
        #        "initialize step."
      else:
        # Update to new filename
        self.filename = filename

    if self.filename:
      self.config = ConfigParser.ConfigParser()
      try:
        self.config.read(self.filename)
        self.inpath = os.path.dirname(os.path.realpath(self.filename)) + '/'
        # Need to have these guys inside "try" to make sure it is set up OK
        # (at least for them)
        self.model     = self.configGet("string", "mode", "model")
        self.dimension = self.configGet("integer", "mode", "dimension")
        self.whichModel_AlreadyRun = True
      except:
        sys.exit("No input file at specified path, or input file configured incorrectly")

    if self.filename:
      # Set clocks to None so if they are called by the getter before the 
      # calculation is performed, there won't be an error
      self.coeff_creation_time = None
      self.time_to_solve = None
      
      # Boundary conditions
      # Not optional: flexural solutions can be very sensitive to b.c.'s
      self.BC_E = self.configGet('string', 'numerical', 'BoundaryCondition_East', optional=False)
      self.BC_W = self.configGet('string', 'numerical', 'BoundaryCondition_West', optional=False)
      self.bclist = [self.BC_E, self.BC_W]
      if self.dimension == 2:
        self.BC_N = self.configGet('string', 'numerical2D', 'BoundaryCondition_North', optional=False)
        self.BC_S = self.configGet('string', 'numerical2D', 'BoundaryCondition_South', optional=False)
        self.bclist += [self.BC_N, self.BC_S]
      # Check that boundary conditions are acceptable with code implementation
      # Acceptable b.c.'s
      self.bc1D = np.array(['Dirichlet0', 'Periodic', 'Mirror', '0Moment0Shear', '0Slope0Shear'])
      self.bc2D = np.array(['Dirichlet0', 'Periodic', 'Mirror', '0Moment0Shear', '0Slope0Shear'])
      for bc in self.bclist:
        if self.dimension == 1:
          if (bc == self.bc1D).any():
            pass
          else:
            sys.exit("'"+bc+"'"+ " is not an acceptable 1D boundary condition and/or\n"\
                     +"is not yet implement in the code. Acceptable boundary conditions\n"\
                     +"are:\n"\
                     +str(self.bc1D)+"\n"\
                     +"Exiting.")
        elif self.dimension == 2:
          if (bc == self.bc2D).any():
            pass
          else:
            sys.exit("'"+bc+"'"+ " is not an acceptable 2D boundary condition and/or\n"\
                     +"is not yet implement in the code. Exiting.")
        else:
          sys.exit("For a flexural solution, grid must be 1D or 2D. Exiting.")

      # Parameters
      self.g = self.configGet('float', "parameter", "GravAccel")
      self.rho_m = self.configGet('float', "parameter", "MantleDensity")
      self.rho_fill = self.configGet('float', "parameter", "InfillMaterialDensity")

      # Grid spacing
      # Unnecessary for PrattAiry, but good to keep along, I think, for use 
      # in model output and plotting.
      # No meaning for ungridded superimposed analytical solutions
      # From input file
      self.dx = self.configGet("float", "numerical", "GridSpacing_x")
      if self.dimension == 2:
        self.dy = self.configGet("float", "numerical2D", "GridSpacing_y")
      # Loading grid
      self.q0 = self.configGet('string', "input", "Loads")
      
      # Verbosity
      self.Verbose = self.configGet("bool", "verbosity", "Verbose")
      # Deebug means that whole arrays, etc., can be printed
      self.Debug = self.configGet("bool", "verbosity", "Debug")

    # Stop program if there is no q0 defined or if it is None-type
    try:
      self.q0
      # Stop program if q0 is None-type
      if type(self.q0) == None: # if is None type, just be patient
        sys.exit("Must define non-None-type q0 by this stage in the initialization step\n"+\
                 "from either input file (string) or direct array import")
    except:
      sys.exit("Must define q0 by this stage in the initialization step\n"+\
               "from either input file (string) or direct array import")

    # If a q0 is a string (i.e. we need to load something)
    if type(self.q0) == str:
      try:
        # First see if it is a full path or directly links from the current
        # working directory
        self.q0 = np.load(self.q0)
        if self.Verbose: print "Loading q0 from numpy binary"
      except:
        try:
          self.q0 = np.loadtxt(self.q0)
          if self.Verbose: print "Loading q0 ASCII"
        except:
          # Then see if it is relative to the location of the input file
          try:
            self.q0 = load(self.inpath + self.q0)
            if self.Verbose: print "Loading q0 from numpy binary"
          except:
            try:
              self.q0 = np.loadtxt(self.inpath + self.q0)
              if self.Verbose: print "Loading q0 ASCII"
            # If failure
            except:
              print "Cannot find q0 file"
              print "q0path = " + self.q0
              print "Looked relative to model python files."
              print "Also looked relative to input file path, " + self.inpath
              print "Exiting."
              sys.exit()
    
    # Check consistency of dimensions
    if self.q0.ndim != self.dimension:
      print "Number of dimensions in loads file is inconsistent with"
      print "number of dimensions in solution technique."
      print "Loads", self.q0.ndim
      print "Dimensions", self.dimension
      print self.q0
      print "Exiting."
      sys.exit()

    # Plotting selection
    self.plotChoice = self.configGet("string", "output", "Plot", optional=True)

  # Finalize
  def finalize(self):
    # Just print a line to stdout
    if self.Quiet==False:
      print ""

  # SAVING TO FILE AND PLOTTING STEPS

  # Output: One of the functions run by isostasy.py; not part of IRF
  # (for standalone model use)
  def output(self):
    if self.Verbose: print 'Output step'
    self.outputDeflections()
    self.plotting()

  # Save output deflections to file, if desired
  def outputDeflections(self):
    """
    Outputs a grid of deflections if an output directory is defined in the 
    input file
    
    If the filename given in the input file ends in ".npy", then a binary 
    numpy grid will be exported.
    
    Otherwise, an ASCII grid will be exported.
    """
    try:
      # If wOutFile exists, has already been set by a setter
      self.wOutFile
      if self.Verbose:
        print "Output filename provided by setter"
        print "Not saving file with this code; that should be handled by the driver"
        
    # Otherwise, it needs to be set by an input file
    except:
      try:
        self.wOutFile = self.configGet("string", "output", "DeflectionOut", optional=True)
        # If this exists and is a string, write output to a file
        if self.wOutFile[-4:] == '.npy':
          from numpy import save
          save(self.wOutFile,self.w)
        else:
          from numpy import savetxt
          # Shouldn't need more than mm precision, at very most
          savetxt(self.wOutFile,self.w,fmt='%.3f')
          if self.Verbose:
            print 'Saving deflections --> ' + self.wOutFile
      except:
        # if there is no parsable output string, do not generate output;
        # this allows the user to leave the line blank and produce no output
        if self.Debug:
          print 'Not writing any deflection output to file'

# class Flexure inherits Isostay and it overrides the __init__ method. It also
# define three different solution methods, which are implemented by its subclass.
class Flexure(Isostasy):
  
  def coeffArraySizeCheck(self):
    """
    Make sure that q0 and coefficient array are the right size compared to 
    each other; otherwise, exit.
    """
    if prod(self.coeff_matrix.shape) != long(prod(np.array(self.q0.shape,dtype=int64)+2)**2):
      print "Inconsistent size of q0 array and coefficient mattrix"
      print "Exiting."
      sys.exit()
      
  def TeArraySizeCheck(self):
    """
    Checks that Te and q0 array sizes are compatible
    """
    # Only if they are both defined and are arrays
    # Both being arrays is a possible bug in this check routine that I have 
    # intentionally introduced
    if type(self.Te) == np.ndarray and type(self.q0) == np.ndarray:
      if self.Te.any() and self.q0.any():
        # Doesn't touch non-arrays or 1D arrays
        if type(self.Te) is np.ndarray:
          if (np.array(self.Te.shape) - 2 != np.array(self.q0.shape)).any():
            sys.exit("q0 and Te arrays have incompatible shapes. Exiting.")
        else:
          if self.Debug: print "Te and q0 array sizes pass consistency check"

  def initialize(self, filename=None):
    super(Flexure, self).initialize(filename)

    # Mode: solution method and type of plate solution (if applicable)
    if self.filename:
      self.method = self.configGet("string", "mode", "method")
      if self.dimension == 2:
        self.PlateSolutionType = self.configGet("string", "mode", "PlateSolutionType")
    
    # Parameters
    self.drho = self.rho_m - self.rho_fill
    if self.filename:
      self.E  = self.configGet("float", "parameter", "YoungsModulus")
      self.nu = self.configGet("float", "parameter", "PoissonsRatio")
    
  ### need to determine its interface, it is best to have a uniform interface
  ### no matter it is 1D or 2D; but if it can't be that way, we can set up a
  ### variable-length arguments, which is the way how Python overloads functions.
  def FD(self):
    if self.Verbose:
      print "Finite Difference Solution Technique"
    # Find the solver
    try:
      self.solver # See if it exists already
    except:
      self.solver = self.configGet("string", "numerical", "Solver")
    if self.filename:
      # In the case that it is iterative, find the convergence criterion
      self.iterative_ConvergenceTolerance = self.configGet("float", "numerical", "ConvergenceTolerance")    
      # Try to import Te grid or scalar for the finite difference solution
      Tepath = self.configGet("string", "input", "ElasticThickness", optional=True)
      # See if there is a pre-made coefficient matrix to import
      coeffPath = self.configGet("string", "input", "CoeffMatrix", optional=True)
      # If there is, import it.
    # or if getter/setter
    elif type(self.Te) == str: # Not super stable here either
      # Try to import Te grid or scalar for the finite difference solution
      Tepath = self.Te
      # If there is, import it.
    try:
      coeffPath
      if coeffPath:
        try:
          self.coeff_matrix = np.load(coeffPath)
          if self.Verbose: print "Loading coefficient array as numpy array binary"
        except:
          try:
            self.coeff_matrix = np.loadtxt(coeffPath)
            if self.Verbose: print "Loading coefficient array as ASCII grid"
          except:
            print "Could not load coefficient array; check filename provided"
            print "Exiting."
            sys.exit()
        print "Any coefficient matrix provided in input file has been ignored,"
        print "as a pre-provided coefficient matrix array is available"

      # Check consistency of size if coeff array was loaded
      coeffArraySizeCheck()

    except:
      coeffPath = None
      
    try:
      Tepath
    except:
      Tepath = None

    # Only get Te if you aren't loading a pre-made coefficient matrix
    if coeffPath == None:
      # No grid?
      if Tepath == None or Tepath == '':
        # Go through this only if using an input file
        if self.filename:
          if self.Verbose: print "Trying to use the scalar elastic thickness"
          # Is there are scalar file?
          try:
            # No Te grid path defined, so going for scalar Te
            self.Te = self.config.getfloat("parameter", "ElasticThickness")
          except:
            try:
              self.Te = float(self.Te)
            except:
              # No Te input provided - scalar or array path
              print "Requested Te file cannot be found, and no scalar elastic"
              print "thickness is provided in input file"
              print "You should add one or the other to the input file."
              print "Exiting."
              sys.exit()
      else:
        # In this case, Tepath has been defined as something by the input file
        try:
          # First see if it is a full path or directly links from the current
          # working directory
          self.Te = np.loadtxt(Tepath)
          print "Loading elastic thickness array from provided file"
        except:
          try:
            # Then see if it is relative to the location of the input file
            self.Te = np.loadtxt(self.inpath + Tepath)
          except:
            # At this point, a Tepath has been provided, but has failed.
            # No unambiguous way to find what input is desired
            # 2 options: (1) Te scalar exists, (2) nothing exists
            # Either way, need to exit
            try:
              TeScalar = self.config.getfloat("parameter", "ElasticThickness")
              print "Requested Te file cannot be found, but a scalar elastic"
              print "thickness is provided in input file."
              print "Ambiguous as to whether a Te grid or a scalar Te value was"
              print "desired."
              print "(Typo in path to input Te grid?)"
            except:
              # No Te input provided - scalar or array path
              print "Requested Te file is provided but cannot be located."
              print "No scalar elastic thickness is provided in input file"
              print "(Typo in path to input Te grid?)"
            print "Exiting."
            sys.exit()
        # At this point, Te from a grid has been successfully loaded.
        # See if there is an ambiguity with a scalar Te also defined
        try:
          TeScalar = self.config.getfloat("parameter", "ElasticThickness")
        except:
          TeScalar = None
        if TeScalar:
          # If this works, need to exit - ambiguous
          print "Both an elastic thickness array and a scalar elastic thickness"
          print "have been loaded - ambiguity; cannot continue"
          print "Exiting."
          sys.exit()
    # Otherwise, all set!
    
  ### need work
  def FFT(self):
    pass

  # SAS and SAS_NG are the exact same here; leaving separate just for symmetry 
  # with other functions

  def SAS(self):
    if self.filename:
      # Define the (scalar) elastic thickness
      self.Te = self.configGet("float", "parameter", "ElasticThickness")
      if self.dimension == 2:
        from scipy.special import kei

  def SAS_NG(self):
    if self.filename:
      # Define the (scalar) elastic thickness
      self.Te = self.configGet("float", "parameter", "ElasticThickness")
      if self.dimension == 2:
        from scipy.special import kei
    
class PrattAiry(Isostasy):
  pass
