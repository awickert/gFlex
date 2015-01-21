import sys, ConfigParser, os
import numpy as np
import time # For efficiency counting
import types # For flow control
from matplotlib import pyplot as plt

class Utility(object):
  """
  Generic utility functions
  """
  def configGet(self, vartype, category, name, optional=False, specialReturnMessage=None):
    """
    Wraps a try / except and a check for self.filename around ConfigParser
    as it talks to the configuration file.
    Also, checks for existence of configuration file so this won't execute (and fail) 
    when no configuration file is provided (e.g., running in coupled mode with CSDMS 
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
        var = self.config.getfloat(category, name)
      elif vartype == 'string' or vartype == 'str':
        var = self.config.get(category, name)
        if var == ""  and optional == False:
          # but "" is acceptable for boundary conditions
          if name[:17] != 'BoundaryCondition':
            if self.Quiet != True:
              print "An empty input string here is not an acceptable option."
              print name, "is not optional."
              print "Program crash likely to occur."
      elif vartype == 'integer' or vartype == 'int':
        var = self.config.getint(category, name)
      elif vartype == 'boolean' or vartype == 'bool':
        var = self.config.getboolean(category, name)
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
          print 'in category "' + category + '" in configuration file.'
          print 'No action related to this optional parameter will be taken.'
          print ""
      else:
        print 'Problem loading ' + vartype + ' "' + name + '" in category "' + category + '" from configuration file.'
        if specialReturnMessage:
          print specialReturnMessage
        sys.exit("Exiting.")

  def set_value(self, value_key, value):
    """
    Universal setter
    """
    # FIRST, VALUES THAT EVERYONE NEEDS::   
    
    # [mode]
    # Dimensions -- 1D or 2D solution.
    if value_key =='dimension':
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
    # NOT IN CONFIGURATION FILE
    elif value_key == 'Loads_grid_stress':
      self.qs = value
    elif value_key == 'Loads_force':
      self.q = value

    # [numerical]
    # Grid spacing
    elif value_key == 'GridSpacing_x':
      self.dx = value
    elif value_key == 'GridSpacing_y':
      if self.Debug:
        print "Setting y-value; should be done only for 2D problems"
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
    # NOT IN CONFIGURATION FILE -- FOR POINT LOADS
    # vectors of x and y values
    elif value_key == 'x':
      self.x = value
    elif value_key == 'y':
      self.y = value

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
    elif value_key == value_key == 'Method':
      self.method = value
    elif value_key == 'PlateSolutionType':
      self.PlateSolutionType = value

    # [Parameters]
    elif value_key == 'YoungsModulus':
      self.E  = value
    elif value_key == 'PoissonsRatio':
      self.nu = value
      
    # [Input]
    elif value_key == 'CoeffArray':
      # This coefficient array is what is used with the UMFPACK direct solver
      # or the iterative solver
      self.coeff_matrix = value
      self.readyCoeff() # See if this is a sparse or something that needs to be loaded
      # if so, let everyone know
      if self.Verbose:
        print "LOADING COEFFICIENT ARRAY"
        print "Elastic thickness maps will not be used for the solution."
    elif value_key == 'ElasticThickness':
      self.Te = value

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
    elif val_string=='SolverTime':
      # Amount of time taken by the solver (direct or iterative)
      return self.time_to_solve
    if val_string=='Loads':
      # This is the model input for the gridded case
      return self.q0
    if val_string=='x':
      # This is a component of the ungridded model input
      # It is also produced during the gridded model run
      # But will be overwritten in those cases
      return self.x
    if val_string=='y':
      # This is a component of the ungridded model input
      # It is also produced during the gridded model run
      # But will be overwritten in those cases
      return self.y
    if val_string=='LoadVerticalNormalForces':
      # This is a component of the ungridded model input
      return self.q
    if val_string=='ElasticThickness':
      # This is the model input
      return self.Te
    if val_string=='Verbosity' or val_string=='Verbose':
      return self.Verbose


class Plotting(object):
  # Plot, if desired
  # 1D all here, 2D in functions
  # Just because there is often more code in 2D plotting functions
  # Also, yes, this portion of the code is NOT efficient or elegant in how it
  # handles functions. But it's just a simple way to visualize results 
  # easily! And not too hard to improve with a bit of time. Anyway, the main
  # goal here is the visualization, not the beauty of the code : )
  def plotting(self):
    #try:
    #  self.plotChoice
    #except:
    #  self.plotChoice = None
    if self.plotChoice:
      if self.Verbose: print "Starting to plot " + self.plotChoice
      if self.dimension==1:
        if self.plotChoice == 'q':
          plt.figure(1)
          if self.method == 'SAS_NG':
            plt.plot(self.x/1000., self.q/(self.rho_m*self.g), 'ko-')
            plt.ylabel('Load volume, mantle equivalent [m$^3$]', fontsize=12, fontweight='bold')
          else:
            plt.plot(self.x/1000., self.qs/(self.rho_m*self.g), 'k-')
            plt.ylabel('Load thickness, mantle equivalent [km]', fontsize=12, fontweight='bold')
          plt.xlabel('Distance along profile [km]', fontsize=12, fontweight='bold')
          plt.tight_layout()
          plt.show()
        elif self.plotChoice == 'w':
          plt.figure(1)
          if self.method == 'SAS_NG':
            plt.plot(self.x/1000., self.w, 'ko-')
          else:
            plt.plot(self.x/1000., self.w, 'k-')
          plt.ylabel('Deflection [m]', fontsize=12, fontweight='bold')
          plt.xlabel('Distance along profile [km]', fontsize=12, fontweight='bold')
          plt.tight_layout()
          plt.show()
        elif self.plotChoice == 'both':
          plt.figure(1,figsize=(6,9))
          plt.subplot(211)
          plt.title('Loads and Lithospheric Deflections', fontsize=16)
          if self.method == 'SAS_NG':
            plt.plot(self.x/1000., self.q/(self.rho_m*self.g), 'ko-')
            plt.ylabel('Load volume, mantle equivalent [m$^3$]', fontsize=12, fontweight='bold')
          else:
            plt.plot(self.x/1000., self.qs/(self.rho_m*self.g), 'k-')
            plt.ylabel('Load thickness, mantle equivalent [m]', fontsize=12, fontweight='bold')
          plt.xlabel('Distance along profile [km]', fontsize=12, fontweight='bold')
          plt.subplot(212)
          if self.method == "SAS_NG":
            plt.plot(self.x, self.w, 'ko-')
          else:
            plt.plot(self.x, self.w, 'k-')
          plt.ylabel('Deflection [m]', fontsize=12, fontweight='bold')
          plt.xlabel('Distance along profile [m]', fontsize=12, fontweight='bold')
          plt.tight_layout()
          plt.show()
        elif self.plotChoice == 'combo':
          fig = plt.figure(1,figsize=(10,6))
          titletext='Loads and Lithospheric Deflections'
          xkm = self.x/1000
          ax = fig.add_subplot(1,1,1)
          # Plot undeflected load
          if self.method == "SAS_NG":
            if self.Quiet == False:
              print "Combo plot can't work with SAS_NG! Don't have mechanism in place\nto calculate load width."
              print "Big problem -- what is the area represented by the loads at the\nextreme ends of the array?"
            #ax.plot(xkm, self.q/(self.rho_m*self.g),'go',linewidth=2,label="Load volume [m^3 mantle equivalent]") # MUST FIX!!!! Turn into m mantle equivalent
          else:
            ax.plot(xkm,self.qs/(self.rho_m*self.g),'g--',linewidth=2,label="Load thickness [m mantle equivalent]")
          # Plot deflected load
          if self.method == "SAS_NG":
            pass
            #ax.plot(xkm,self.q/(self.rho_m*self.g) + self.w,'go-',linewidth=2,label="Load volume [m^3] mantle equivalent]")
          else:
            ax.plot(xkm,self.qs/(self.rho_m*self.g) + self.w,'g-',linewidth=2,label="Deflection [m] + load thickness [m mantle equivalent]")
          # Plot deflection
          if self.method == "SAS_NG":
            ax.plot(xkm, self.w, 'ko-', linewidth=2, label="Deflection [m mantle equivalent]")
          else:
            ax.plot(xkm,self.w, 'k-', linewidth=2, label="Deflection [m mantle equivalent]")
          # Set y min to equal to the absolute value maximum of y max and y min
          # (and therefore show isostasy better)
          yabsmax = max(abs(np.array(plt.ylim())))
          # Y axis label
          plt.ylim((-yabsmax,yabsmax))
          # Plot title selector -- be infomrative
          try:
            self.Te
            if self.method == "FD":
              if type(self.Te) is np.ndarray:
                if (self.Te != (self.Te).mean()).any():
                  plt.title(titletext,fontsize=16)       
                else:
                  plt.title(titletext + ', $T_e$ = ' + str((self.Te / 1000).mean()) + " km", fontsize=16)
              else:
                plt.title(titletext + ', $T_e$ = ' + str(self.Te / 1000) + " km", fontsize=16)
            else:
              plt.title(titletext + ', $T_e$ = ' + str(self.Te / 1000) + " km", fontsize=16)
          except:
            plt.title(titletext,fontsize=16)       
          # x and y labels
          plt.ylabel('Loads and flexural response [m]',fontsize=16)
          plt.xlabel('Distance along profile [km]',fontsize=16)
          # legend -- based on lables
          plt.legend(loc=0,numpoints=1,fancybox=True)
          plt.tight_layout()
          plt.show()
        else:
          if self.Quiet == False:
            print 'Incorrect plotChoice input, "' + self.plotChoice + '" provided.'
            print "Possible input strings are: q, w, both, and (for 1D) combo"
            print "Unable to produce plot."
      elif self.dimension==2:
        if self.plotChoice == 'q':
          fig = plt.figure(1, figsize=(8,6))
          if self.method != 'SAS_NG':
            self.surfplot(self.qs/(self.rho_m*self.g), 'Load thickness, mantle equivalent [m]')
            plt.show()
          else:
            self.xyzinterp(self.q, 'Load volume, mantle equivalent [m$^3$]')
          plt.tight_layout()
          plt.show()
        elif self.plotChoice == 'w':
          fig = plt.figure(1, figsize=(8,6))
          if self.method != 'SAS_NG':
            self.surfplot(self.w, 'Deflection [m]')
            plt.show()
          else:
            self.xyzinterp(self.w, 'Deflection [m]')
          plt.tight_layout()
          plt.show()
        elif self.plotChoice == 'both':
          plt.figure(1,figsize=(6,9))
          if self.method != 'SAS_NG':
            self.twoSurfplots()
            plt.show()
          else:
            plt.subplot(211)
            self.xyzinterp(self.q, 'Load volume, mantle equivalent [m$^3$]')
            plt.subplot(212)
            self.xyzinterp(self.w, 'Deflection [m]')
            plt.tight_layout()
            plt.show()
        else:
          if self.Quiet == False:
            print 'Incorrect plotChoice input, "' + self.plotChoice + '" provided.'
            print "Possible input strings are: q, w, both, and (for 1D) combo"
            print "Unable to produce plot."

  def surfplot(self, z, titletext):
    """
    Plot if you want to - for troubleshooting - 1 figure
    """
    plt.imshow(z, extent=(0, self.dx/1000.*z.shape[0], self.dy/1000.*z.shape[1], 0)) #,interpolation='nearest'
    plt.xlabel('x [km]', fontsize=12)
    plt.ylabel('y [km]', fontsize=12)
    plt.colorbar()

    plt.title(titletext,fontsize=16)

  def twoSurfplots(self):
    """
    Plot multiple subplot figure for 2D array
    """
    # Could more elegantly just call surfplot twice
    # And also could include xyzinterp as an option inside surfplot.
    # Noted here in case anyone wants to take that on in the future...

    plt.subplot(211)
    plt.title('Load thickness, mantle equivalent [m]',fontsize=16)
    plt.imshow(self.qs/(self.rho_m*self.g), extent=(0, self.dx/1000.*self.qs.shape[0], self.dy/1000.*self.qs.shape[1], 0))
    plt.xlabel('x [km]', fontsize=12, fontweight='bold')
    plt.ylabel('y [km]', fontsize=12, fontweight='bold')
    plt.colorbar()

    plt.subplot(212)
    plt.title('Deflection [m]')
    plt.imshow(self.w, extent=(0, self.dx/1000.*self.w.shape[0], self.dy/1000.*self.w.shape[1], 0))
    plt.xlabel('x [km]', fontsize=12, fontweight='bold')
    plt.ylabel('y [km]', fontsize=12, fontweight='bold')
    plt.colorbar()
  
  def xyzinterp(self, z, titletext):
    """
    Interpolates and plots ungridded model outputs from SAS_NG solution
    """
    # Help from http://wiki.scipy.org/Cookbook/Matplotlib/Gridding_irregularly_spaced_data
    
    if self.Verbose:
      print "Starting to interpolate grid -- can be a slow process!"
    
    from scipy.interpolate import griddata
    import numpy.ma as ma
    
    # define grid.
    xmin = np.min(self.x)
    xmean = np.mean(self.x) # not used right now
    xmax = np.max(self.x)
    ymin = np.min(self.y)
    ymean = np.mean(self.y) # not used right now
    ymax = np.max(self.y)
    x_range = xmax - xmin
    y_range = ymax - ymin
    
    # x and y grids
    # 100 cells on each side -- just for plotting, not so important
    # to optimize with how many points are plotted
    #xi = np.linspace(xmin-.05*x_range, xmax+.05*x_range, 200)
    #yi = np.linspace(ymin-.05*y_range, ymax+.05*y_range, 200)
    xi = np.linspace(xmin, xmax, 200)
    yi = np.linspace(ymin, ymax, 200)
    # grid the z-axis
    zi = griddata((self.x, self.y), z, (xi[None,:], yi[:,None]), method='cubic')
    # contour the gridded outputs, plotting dots at the randomly spaced data points.
    #CS = plt.contour(xi,yi,zi,15,linewidths=0.5,colors='k') -- don't need lines
    CS = plt.contourf(xi/1000.,yi/1000.,zi,100,cmap=plt.cm.jet)
    plt.colorbar() # draw colorbar
    # plot model points.
    plt.scatter(self.x/1000., self.y/1000., marker='.', c='k', s=1)
    #plt.hexbin(self.x, self.y, C=self.w) -- show colors on points -- harder to see
    plt.xlabel('x [km]')
    plt.ylabel('y [km]')
    # Limits -- to not get messed up by points (view wants to be wider so whole point visible)
    plt.xlim( (xi[0]/1000., xi[-1]/1000.) )
    plt.ylim( (yi[0]/1000., yi[-1]/1000.) )
    # Title
    plt.title(titletext, fontsize=16)

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
          self.dimension = self.configGet("integer", "mode", "dimension")
          self.whichModel_AlreadyRun = True
        except:
          sys.exit("No configuration file at specified path, or configuration file configured incorrectly")

# class Isostasy inherits IRF interface, and it determines the simulation type
# by reading three parameters from configuration file, but it does not set up other
# parameters, which is the responsibility of derived concrete classes.
class Isostasy(Utility, Plotting):

  def __init__(self, filename=None):
    # 17 Nov 2014: Splitting out initialize from __init__ to allow space
    # to use getters and setters to define values
    
    # Use standard routine to pull out values
    # If no filename provided, will not initialize configuration file.
    self.filename = filename
    
    # DEFAULT VERBOSITY
    # Set default "quiet" to False, unless set by setter or overwritten by 
    # the configuration file.
    self.Quiet = False
    # And also set default verbosity
    self.Verbose = True
    self.Debug = False

  def initialize(self, filename=None):
    # Values from configuration file

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
      # Set up ConfigParser
      self.config = ConfigParser.ConfigParser()
      try:
        self.config.read(self.filename)
        self.inpath = os.path.dirname(os.path.realpath(self.filename)) + '/'
        # Need to have these guys inside "try" to make sure it is set up OK
        # (at least for them)
        self.dimension = self.configGet("integer", "mode", "dimension")
        self.whichModel_AlreadyRun = True
      except:
        sys.exit("No configuration file at specified path, or configuration file configured incorrectly")

      # Set verbosity for model run
      # Default is "verbose" with no debug or quiet
      # Verbose
      try:
        self.Verbose = self.configGet("bool", "verbosity", "Verbose", optional=False)
      except:
        pass
      # Deebug means that whole arrays, etc., can be printed
      try:
        self.Debug = self.configGet("bool", "verbosity", "Debug", optional=False)
      except:
        pass
      # Deebug means that whole arrays, etc., can be printed
      try:
        self.Quiet = self.configGet("bool", "verbosity", "Quiet", optional=False)
      except:
        pass
      # Quiet overrides all others
      if self.Quiet:
        self.Debug = False
        self.Verbose = False
    
    # Introduce model
    # After configuration file can define "Quiet", and getter/setter should be done
    # by this point if we are going that way.
    if self.Quiet == False:
      print "" # Blank line at start of run
      print "*********************************************"
      print "*** Initializing gFlex development branch ***"
      print "*********************************************"
      print ""
      print "Open-source licensed under GNU GPL v3"
      print ""

    if self.filename:
      # Set clocks to None so if they are called by the getter before the 
      # calculation is performed, there won't be an error
      self.coeff_creation_time = None
      self.time_to_solve = None
      
      self.method = self.configGet("string", "mode", "method")
      # Boundary conditions
      # This used to be nested inside an "if self.method == 'FD'", but it seems 
      # better to define these to ensure there aren't mistaken impressions
      # about what they do for the SAS case
      # Not optional: flexural solutions can be very sensitive to b.c.'s
      self.BC_E = self.configGet('string', 'numerical', 'BoundaryCondition_East', optional=False)
      self.BC_W = self.configGet('string', 'numerical', 'BoundaryCondition_West', optional=False)
      self.bclist = [self.BC_E, self.BC_W]
      if self.dimension == 2:
        self.BC_N = self.configGet('string', 'numerical2D', 'BoundaryCondition_North', optional=False)
        self.BC_S = self.configGet('string', 'numerical2D', 'BoundaryCondition_South', optional=False)
        self.bclist += [self.BC_N, self.BC_S]

      # Parameters
      self.g = self.configGet('float', "parameter", "GravAccel")
      self.rho_m = self.configGet('float', "parameter", "MantleDensity")
      self.rho_fill = self.configGet('float', "parameter", "InfillMaterialDensity")

      # Grid spacing
      if self.method != 'SAS_NG':
        # No meaning for ungridded superimposed analytical solutions
        # From configuration file
        self.dx = self.configGet("float", "numerical", "GridSpacing_x")
        if self.dimension == 2:
          self.dy = self.configGet("float", "numerical2D", "GridSpacing_y")

      # Loading grid
      # q0 is either a load array or an x,y,q array.
      # Therefore q_0, initial q, before figuring out what it really is
      # for grid, q0 could also be written as $q_\sigma$ or q/(dx*(dy))
      # it is a surface normal stress that is h_load * rho_load * g
      # it later is combined with dx and (if 2D) dy for FD cases
      # for point loads, need mass: q0 should be written as [x, (y), force])
      self.q0 = self.configGet('string', "input", "Loads")
      
    # Stop program if there is no q0 defined or if it is None-type
    try:
      self.q0
      # Stop program if q0 is None-type
      if type(self.q0) == None: # if is None type, just be patient
        sys.exit("Must define non-None-type q0 by this stage in the initialization step\n"+\
                 "from either configuration file (string) or direct array import")
    except:
      try:
        self.q
      except:
        try:
          self.qs
        except:
          sys.exit("Must define q0, q, or qs by this stage in the initialization step\n"+\
                   "from either configuration file (string) or direct array import")

    # Ignore this if no q0 set
    try:
      self.q0
    except:
      self.q0 = None
    if self.q0:
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
            # Then see if it is relative to the location of the configuration file
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
                print "Also looked relative to configuration file path, " + self.inpath
                print "Exiting."
                sys.exit()
      
      # Check consistency of dimensions
      if self.method != 'SAS_NG':
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
    configuration file
    
    If the filename given in the configuration file ends in ".npy", then a binary 
    numpy grid will be exported.
    
    Otherwise, an ASCII grid will be exported.
    """
    try:
      # If wOutFile exists, has already been set by a setter
      self.wOutFile
      if self.Verbose:
        print "Output filename provided by setter"
        print "Not saving file with this code; that should be handled by the driver"
        
    # Otherwise, it needs to be set by an configuration file
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

  def bc_check(self):

    # Check that boundary conditions are acceptable with code implementation
    # Acceptable b.c.'s
    if self.method == 'FD':
      self.bc1D = np.array(['Dirichlet0', 'Periodic', 'Mirror', '0Moment0Shear', '0Slope0Shear'])
      self.bc2D = np.array(['Dirichlet0', 'Periodic', 'Mirror', '0Moment0Shear', '0Slope0Shear'])
      for bc in self.bclist:
        if self.dimension == 1:
          if (bc == self.bc1D).any():
            pass
          else:
            sys.exit("'"+bc+"'"+ " is not an acceptable 1D finite difference boundary condition\n"\
                     +"and/or is not yet implement in the code. Acceptable boundary conditions\n"\
                     +"are:\n"\
                     +str(self.bc1D)+"\n"\
                     +"Exiting.")
        elif self.dimension == 2:
          if (bc == self.bc2D).any():
            pass
          else:
            sys.exit("'"+bc+"'"+ " is not an acceptable 2D finite difference boundary condition\n"\
                     +"and/or is not yet implement in the code. Acceptable boundary conditions\n"\
                     +"are:\n"\
                     +str(self.bc2D)+"\n"\
                     +"Exiting.")
        else:
          sys.exit("For a flexural solution, grid must be 1D or 2D. Exiting.")
    else:
      if self.BC_E == 'NoOutsideLoads' or self.BC_E == '' \
         and self.BC_W == 'NoOutsideLoads' or self.BC_W == '':
        if self.BC_E == '' or self.BC_W == '':
          if self.Verbose:
            print "Assuming NoOutsideLoads boundary condition, as this is implicit in the " 
            print "  superposition-based analytical solution"
      else:
        if self.Quiet == False:
          print ""
          print ">>> BOUNDARY CONDITIONS IMPROPERLY DEFINED <<<"
          print ""
          print "For analytical solutions the boundaries must be either:"
          print ""
          print "* NoOutsideLoads (explicitly)"
          print "* <left blank>"
          print ""
          print "The latter is to implictly indicate a desire to use the only"
          print "boundary condition available for the superposition-based"
          print "analytical solutions."
          print "This check is in place to ensure that the user does not apply"
          print "boundary conditions for finite difference solutions to the"
          print "analytical solutions and expect them to work."
          print ""
          sys.exit()
  
# class Flexure inherits Isostay and it overrides the __init__ method. It also
# define three different solution methods, which are implemented by its subclass.
class Flexure(Isostasy):
  
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
    
  def coeffArraySizeCheck(self):
    """
    Make sure that q0 and coefficient array are the right size compared to 
    each other (for finite difference if loading a pre-build coefficient
    array). Otherwise, exit.
    """
    if prod(self.coeff_matrix.shape) != long(prod(np.array(self.qs.shape,dtype=int64)+2)**2):
      print "Inconsistent size of q0 array and coefficient mattrix"
      print "Exiting."
      sys.exit()
      
  def TeArraySizeCheck(self):
    """
    Checks that Te and q0 array sizes are compatible
    For finite difference solution.
    """
    # Only if they are both defined and are arrays
    # Both being arrays is a possible bug in this check routine that I have 
    # intentionally introduced
    if type(self.Te) == np.ndarray and type(self.q0) == np.ndarray:
      if self.Te.any() and self.qs.any():
        # Doesn't touch non-arrays or 1D arrays
        if type(self.Te) is np.ndarray:
          if (np.array(self.Te.shape) - 2 != np.array(self.q0.shape)).any():
            sys.exit("q0 and Te arrays have incompatible shapes. Exiting.")
        else:
          if self.Debug: print "Te and qs array sizes pass consistency check"

  ### need to determine its interface, it is best to have a uniform interface
  ### no matter it is 1D or 2D; but if it can't be that way, we can set up a
  ### variable-length arguments, which is the way how Python overloads functions.
  def FD(self):
    """
    Handles set-up for the finite difference solution method
    """
    if self.Verbose:
      print "Finite Difference Solution Technique"
    # Define a stress-based qs = q0
    self.qs = self.q0.copy()
    # Remove self.q0 to avoid issues with multiply-defined inputs
    # q0 is the parsable input to either a qs grid or contains (x,(y),q)
    del self.q0
    # Is there a slver defined?
    try:
      self.solver # See if it exists already
    except:
      # Well, will fail if it doesn't see this, maybe not the most reasonable
      # error message.
      if self.filename:
        self.solver = self.configGet("string", "numerical", "Solver")
      else:
        sys.exit("No solver defined!")
    # Check if a coefficient array has been defined
    # It would only be by a getter or setter;
    # no way to do I/O with this with present configuration files
    try:
      self.coeff_matrix
    except:
      self.coeff_matrix = None
    # Check consistency of size if coeff array was loaded
    if self.filename:
      # In the case that it is iterative, find the convergence criterion
      self.iterative_ConvergenceTolerance = self.configGet("float", "numerical", "ConvergenceTolerance")    
      # Try to import Te grid or scalar for the finite difference solution
      try:
        self.Te = self.configGet("float", "input", "ElasticThickness", optional=True)
        Tepath = None
      except:
        Tepath = self.configGet("string", "input", "ElasticThickness", optional=True)
      if self.Te is None:
        if self.coeff_matrix is not None:
          pass
        else:
          # Have to bring this out here in case it was discovered in the 
          # try statement that there is no value given
          sys.exit("No input elastic thickness supplied.")
    # or if getter/setter
    if type(self.Te) == str: 
      # Try to import Te grid or scalar for the finite difference solution
      Tepath = self.Te
    else:
      Tepath = None # in case no self.filename present (like for GRASS GIS)
    # If there is a Tepath, import Te
    # Assume that even if a coeff_matrix is defined
    # That the user wants Te if they gave the path
    if Tepath:
      try:
        # First see if it is a full path or directly links from the current
        # working directory
        self.Te = np.loadtxt(Tepath)
        if self.Verbose:
          print "Loading elastic thickness array from provided file path"
      except:
        try:
          # Then see if it is relative to the location of the configuration file
            self.Te = np.loadtxt(self.inpath + Tepath)
            if self.Verbose:
                print "Elastic thickness array loaded from provided filename"
        except:
          if quiet == False:
            print "Requested Te file is provided but cannot be located."
            print "No scalar elastic thickness is provided in configuration file"
            print "(Typo in path to input Te grid?)"
          if self.coeff_matrix is not None:
            if quiet == False:
              print "But a coefficient matrix has been found."
              print "Calculations will be carried forward using it."
          else:
            if quiet == False:
              print "Exiting."
            sys.exit()

      # Check that Te is the proper size if it was loaded
      if self.Te:
        self.TeArraySizeCheck()
    
  ### need work
  def FFT(self):
    pass

  # SAS and SAS_NG are the exact same here; leaving separate just for symmetry 
  # with other functions

  def SAS(self):
    if self.filename:
      # Define the (scalar) elastic thickness
      self.Te = self.configGet("float", "input", "ElasticThickness")
      if self.dimension == 2:
        from scipy.special import kei
    # Define a stress-based qs = q0
    self.qs = self.q0.copy()
    # Remove self.q0 to avoid issues with multiply-defined inputs
    # q0 is the parsable input to either a qs grid or contains (x,(y),q)
    del self.q0

  def SAS_NG(self):
    if self.filename:
      # Define the (scalar) elastic thickness
      self.Te = self.configGet("float", "input", "ElasticThickness")
      if self.dimension == 2:
        from scipy.special import kei
    # Parse out input q0 into variables of imoprtance for solution
    if self.dimension == 1:
      try:
        # If these have already been set, e.g., by getters/setters, great!
        self.x
        self.q
      except:
        # Using [x, y, w] configuration file
        if self.q0.shape[1] == 2:
          self.x = self.q0[:,0]
          self.q = self.q0[:,1]
      else:
        sys.exit("For 1D (ungridded) SAS_NG configuration file, need [x,w] array. Your dimensions are: "+str(self.q0.shape))
    else:
      try:
        # If these have already been set, e.g., by getters/setters, great!
        self.x
        self.y
        self.q
      except:
        # Using [x, y, w] configuration file
        if self.q0.shape[1] == 3:
          self.x = self.q0[:,0]
          self.y = self.q0[:,1]
          self.q = self.q0[:,2]
        else:
          sys.exit("For 2D (ungridded) SAS_NG configuration file, need [x,y,w] array. Your dimensions are: "+str(self.q0.shape))
    # Remove self.q0 to avoid issues with multiply-defined inputs
    # q0 is the parsable input to either a qs grid or contains (x,(y),q)
    del self.q0
    
