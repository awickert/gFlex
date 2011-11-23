import sys, ConfigParser, os
import numpy as np
import types # For flow control
#import CSDMS_base

debug = True

# IRF interface
class IRF(object):
  def initialize(self, file):
    pass

  def run(self):
    pass
  
  def run_step(self, i):
    # Maybe meant to run one step?
    # Right now, does nothing.
    # I should make something that does an iterative solution... hmm...
    pass
  
  def finalize(self):
    pass
  
  def get(self, *value):
    pass

  def set(self, *value):
    pass


# class Isostasy inherits IRF interface, and it determines the simulation type
# by reading three parameters from input file, but it does not set up other
# parameters, which is the responsibility of derived concrete classes.
class Isostasy(IRF):

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
        if var == "":
          print "Input strings cannot be empty unless they are optional."
          sys.exit()
      elif vartype == 'integer' or vartype == 'int':
        var = self.config.getint(category,name)
      else:
        print "Please enter 'float', 'string' (or 'str'), or 'integer' (or 'int') for vartype"
        sys.exit() # Won't exit, but will lead to exception
      return var
    except:
      if optional:
        # Carry on if the variable is optional
        var = None
        print 'No value entered for optional parameter "' + name + '" in category "' + category + '" in input file.'
        print 'No action related to this optional parameter will be taken.'
      else:
        print 'Problem loading ' + vartype + ' "' + name + '" in category "' + category + '" from input file.'
        if specialReturnMessage:
          print specialReturnMessage
        print "Exiting."
        sys.exit()

  def whichModel(self, filename=None):
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
          print "Successfully loaded input file and obtained model mode"
          self.whichModel_AlreadyRun = True
        except:
          print "No input file at specified path, or input file configured incorrectly"
          sys.exit()

  def initialize(self, filename=None):
  
    if debug: print "" # Blank line at start of run
    if debug: print "Starting to initialize..."

    self.whichModel(filename) # Use standard routine to pull out values
                         # If no filename provided, will not initialize input 
                         # file. If input file provided this should have 
                         # already run and raised that flag.
                         # So really, this should never be executed.

    # Parameters
    # From input file
    self.g = self.configGet('float', "parameter", "GravAccel")
    self.rho_m = self.configGet('float', "parameter", "MantleDensity")
    self.rho_fill = self.configGet('float', "parameter", "InfillMaterialDensity")

    # Grid spacing
    # Unnecessary for PrattAiry, but good to keep along, I think, for use 
    # in model output and plotting.
    # No meaning for ungridded superimposed analytical solutions
    # From input file
    self.dx = self.configGet("float", "numerical", "GridSpacing")
    
    # Loading grid
    q0path = self.configGet('string', "input", "Loads")

    try:
      # First see if it is a full path or directly links from the current
      # working directory
      try:
        self.q0 = load(q0path)
        if debug: print "Loading q0 from numpy binary"
      except:
        self.q0 = np.loadtxt(q0path)
        if debug: print "Loading q0 ASCII"
    except:
      try:
        # Then see if it is relative to the location of the input file
        try:
          self.q0 = load(self.inpath + q0path)
          if debug: print "Loading q0 from numpy binary"
        except:
          self.q0 = np.loadtxt(self.inpath + q0path)
          if debug: print "Loading q0 ASCII"
      except:
        print "Cannot find q0 file"
        print "q0path = " + q0path
        print "Looked relative to model python files."
        print "Also looked relative to input file path, " + self.inpath
        print "Exiting."
        sys.exit()
      
    # Check consistency of dimensions
    if self.q0.ndim != self.dimension:
      print "Number of dimensions in loads file is inconsistent with"
      print "number of dimensions in solution technique."
      print "Exiting."
      sys.exit()

    # Plotting selection
    self.plotChoice = self.configGet("string","output","Plot",optional=True)

  # Finalize: just print a line to stdout
  def finalize(self):
    print ""

  # UNIVERSAL SETTER: VALUES THAT EVERYONE NEEDS
  def set_value(self, value_key, value):
    # Model type
    if value_key == 'model':
      self.model = value
    elif value_key =='dimension':
      self.dimension = value
    # Parameters
    if value_key == 'GravAccel':
      self.g = value
    elif value_key == 'MantleDensity':
      self.rho_m = value
    elif value_key == 'InfillMaterialDensity':
      self.rho_fill = value
    # Grid spacing
    elif value_key == 'GridSpacing':
      self.dx = value
    # Loading grid
    elif value_key == 'Loads':
      self.q0 = value
      # Te setting
      try:
        # If self.q0 == None, then we know we have tried to make Te first and  
        # need a grid for some reason
        # "Try" b/c self.q0 might be undefined
        if self.q0 == None:
          readyElasticThickness() # Elastic thicnkess may need to be modified
      except:
        pass # If this doesn't succeed, everything should be fine.
             # And if it does, program will try to make everything work
    # Dimensions
    elif value_key == "x":
      self.x = value
    elif value_key == "y":
      self.y = value
    # Output
    elif value_key == 'DeflectionOut':
      # Output file name (string)
      self.wOutFile = value
    elif value_key == 'Plot':
      # 'q0', 'w', 'both', or (1D) 'combo'
      self.plotChoice = value
  
  # UNIVERSAL GETTER
  def get_value(self, val_string):
    if val_string=='Deflection':
      # This is the primary model output
      return self.w
    elif val_string=='CoeffMatrix':
      # This is to hold onto the coefficient matrix in memory so it doesn't 
      # need to be reloaded or recreated
      return self.coeff


  # SAVING TO FILE AND PLOTTING STEPS

  # Output: One of the functions run by isostasy.py; not part of IRF
  # (for standalone model use)
  def output(self):
    if debug: print 'Output step'
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
      if debug:
        print "Output filename provided by setter"
        print "Not saving file with this code; that should be handled by the driver"
        
    # Otherwise, it needs to be set by an input file
    except:
      try:
        self.wOutFile = self.configGet("string", "output", "DeflectionOut",optional=True)
        # If this exists and is a string, write output to a file
        if self.WoutFile[-4:] == '.npy':
          from numpy import save
          save(self.wOutFile,self.w)
        else:
          from numpy import savetxt
          # Shouldn't need more than mm precision, at very most
          savetxt(self.wOutFile,self.w,fmt='%.3f')
        if debug: print 'Saving deflections --> ' + self.wOutFile
      except:
        # if there is no parsable output string, do not generate output;
        # this allows the user to leave the line blank and produce no output
        print 'Not writing any deflection output to file'

  # Plot, if desired
  def plotting(self):
    if self.plotChoice:
      if debug: print "Starting to plot " + self.plotChoice
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
          print 'Incorrect plotChoice input, "' + self.plotChoice + '" provided.'
          print "Possible input strings are: q0, w, both, and (for 1D) combo"
          print "Unable to produce plot."

  def linesubplots(self,figNum=1):
    from matplotlib.pyplot import plot, show, figure, subplot, xlabel, \
                                  ylabel, title
    
    figure(figNum)

    subplot(211)
    title('Loads and Lithospheric Deflections',fontsize=20)
    if self.method == "SPA_NG":
      plot(self.x,self.q0/(self.rho_m*self.g),'o')
    else:
      plot(self.x,self.q0/(self.rho_m*self.g))
    ylabel('Load thickness, mantle equivalent [m]')

    subplot(212)
    if self.method == "SPA_NG":
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
                                  ylabel, title, axis, ylim
    
    fig = figure(figNum,figsize=(10,6))
    ax = fig.add_subplot(1,1,1)
    
    xkm = self.x/1000

    # Plot undeflected load
    if self.method == "SPA_NG":
      ax.plot(xkm,self.q0/(self.rho_m*self.g),'go',linewidth=2)
    else:
      ax.plot(xkm,self.q0/(self.rho_m*self.g),'g--',linewidth=2)
    
    # Plot deflected load
    if self.method == "SPA_NG":
      ax.plot(xkm,self.q0/(self.rho_m*self.g) + self.w,'go-',linewidth=2)
    else:
      ax.plot(xkm,self.q0/(self.rho_m*self.g) + self.w,'g',linewidth=2)

    # Plot deflection
    if self.method == "SPA_NG":
      ax.plot(xkm,self.w,'ko-',linewidth=2)
    else:
      ax.plot(xkm,self.w,'k',linewidth=2)
    
    # Set y min to equal y max (and therefore show isostasy better)
    ymax = axis()[-1]
    ylim((-ymax,ymax))

    if self.method == "FD":
      if (self.Te != (self.Te).mean()).any():
        title(titletext,fontsize=20)       
      else:
        title(titletext + ', $T_e$ = ' + str((self.Te / 1000).mean()) + " km",fontsize=20)
    else:
      title(titletext + ', $T_e$ = ' + str(self.Te / 1000) + " km",fontsize=20)
      
    ylabel('Load thickness, mantle equivalent [m]',fontsize=16)
    xlabel('Distance along profile [km]',fontsize=16)
    
    show()


  def surfplot(self,data,titletext,figNum=1):
    """
    Plot if you want to - for troubleshooting
    """
    from matplotlib.pyplot import imshow, show, figure, colorbar, title
    
    figure(figNum)

    imshow(data) #,interpolation='nearest'
    colorbar()

    title(titletext,fontsize=16)
    
    show()

  def surfsubplots(self,figNum=1):
    from matplotlib.pyplot import imshow, show, figure, subplot, xlabel, \
                                  ylabel, title, colorbar
    
    figure(figNum,figsize=(6,9))

    subplot(211)
    title('Load thickness, mantle equivalent [m]',fontsize=16)
    imshow(self.q0/(self.rho_m*self.g))
    colorbar()

    subplot(212)
    title('Deflection [m]')
    imshow(self.w)
    colorbar()

    show()

  # CODE TO ABORT RUN AFTER ERROR
  # mainly to avoid all of the gobbedygook after I print my error message
  # that I've tailored to this code
  def abort(self):
    print("Aborting.")
    sys.exit() # Stop program from running after printing only
               # my error message (calls "raise systemExit")

# class Flexure inherits Isostay and it overrides the __init__ method. It also
# define three different solution methods, which are implemented by its subclass.
class Flexure(Isostasy):
  def coeffArraySizeCheck(self):
    """
    Make sure that q0 and coefficient array are the right size compared to 
    each other; otherwise, exit.
    """
    if prod(self.coeff.shape) != long(prod(np.array(self.q0.shape,dtype=int64)+2)**2):
      print "Inconsistent size of q0 array and coefficient mattrix"
      print "Exiting."
      sys.exit()
      
  def TeArraySizeCheck(self):
    """
    Checks that Te and q0 array sizes are compatible
    """
    # Only if they are both defined
    if self.Te and self.q0:
      # Doesn't touch non-arrays or 1D arrays
      if type(self.Te) is np.ndarray:
        if prod(self.Te.shape) > 1:
          # +2 for fringes in Te required for finite difference solution
          if prod(self.Te.shape) != long(prod(np.array(self.q0.shape,dtype=int64)+2)):
            print "Inconsistent Te and q0 array shapes."
            print "Exiting."
            sys.exit()
      else:
        if debug: print "Te and q0 array sizes pass consistency check"


  def initialize(self, filename=None):
    super(Flexure, self).initialize(filename)

    # Solution method
    if self.filename:
      self.method = self.configGet("string", "mode", "method")
    
    # Parameters
    self.drho = self.rho_m - self.rho_fill
    if self.filename:
      self.E  = self.configGet("float", "parameter", "YoungsModulus")
      self.nu = self.configGet("float", "parameter", "PoissonsRatio")
    
  ### need to determine its interface, it is best to have a uniform interface
  ### no matter it is 1D or 2D; but if it can't be that way, we can set up a
  ### variable-length arguments, which is the way how Python overloads functions.
  def FD(self):
    print "Finite Difference Solution Technique"
    if self.filename:
      # Try to import Te grid or scalar for the finite difference solution
      Tepath = self.configGet("string", "input", "ElasticThickness",optional=True)
      # See if there is a pre-made coefficient matrix to import
      coeffPath = self.configGet("string", "input", "CoeffMatrix",optional=True)
      # If there is, import it.
      if coeffPath:
        try:
          self.coeff = np.load(coeffPath)
          if debug: print "Loading coefficient array as numpy array binary"
        except:
          try:
            self.coeff = np.loadtxt(coeffPath)
            if debug: print "Loading coefficient array as ASCII grid"
          except:
            print "Could not load coefficient array; check filename provided"
            print "Exiting."
            sys.exit()
        print "Any coefficient matrix provided in input file has been ignored,"
        print "as a pre-provided coefficient matrix array is available"

        # Check consistency of size
        coeffArraySizeCheck()

      # Only get Te if you aren't loading a pre-made coefficient matrix
      if coeffPath == None:
        # No grid?
        if Tepath == None:
          if debug: print "Trying to use the scalar elastic thickness"
          # Is there are scalar file?
          try:
            # No Te grid path defined, so going for scalar Te
            TeScalar = self.config.getfloat("parameter", "ElasticThickness")
            q0shape = np.array(self.q0.shape)
            for i in range(len(q0shape)):
              q0shape[i] += 2 # padding for numerical solution
            self.Te = TeScalar*np.ones(q0shape)
            print "Using constant elastic thickness at provided value"
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

  def SPA(self):
    if self.filename:
      # Define the (scalar) elastic thickness
      self.Te = self.configGet("float", "parameter", "ElasticThickness")

  def SPA_NG(self):
    if self.filename:
      # Define the (scalar) elastic thickness
      self.Te = self.configGet("float", "parameter", "ElasticThickness")

    
  # UNIVERSAL SETTER: LITHOSPHERIC ELASTIC PROPERTIES AND SOLUTION METHOD
  # FOR FLEXURAL ISOSTASY
  def set_value(self, value_key, value):
    # Parameters
    if value_key == 'YoungsModulus':
      self.E  = value
    elif value_key == 'PoissonsRatio':
      self.nu = value
    elif value_key == 'ElasticThickness':
      self.Te = value # How to dynamically type for scalar or array?
                      # Python is smart enough to handle it.
      self.readyElasticThickness() # But need a program to handle converting 
                                   # scalars and arrays, as well as potentially 
                                   # needing to load a Te file
    elif value_key == 'CoeffArray':
      # This coefficient array is what is used with the UMFPACK direct solver
      self.coeff = CoeffArray
      self.readyCoeff() # See if this is a sparse or something that needs to be loaded
      coeffArraySizeCheck() # Make sure that array size is all right
      # if so, let everyone know
      print "LOADING COEFFICIENT ARRAY"
      print "Elastic thickness maps will not be used for the solution."
    elif value_key == 'method':
      self.method = value
      print "method set"
    # Inherit from higher-level setter, if not one of these
    # Currently not inheriting from the whichModel() setter, as I
    # figure that has to be done right away or not at all
    super(Flexure, self).set_value(value_key, value)
    
  def readyCoeff(self):
    from scipy import sparse
    if sparse.issparse(self.coeff):
      pass # Good type
    # Otherwise, try to load from file
    elif type(self.coeff) is types.StringType:
      pass
      print "Loading sparse coefficient arrays is not yet implemented."
      print "This must be done soon."
      print "Exiting."
      sys.exit()
    else:
      try:
        self.coeff = sparse.dia_matrix(self.coeff)
      except:
        "Failed to make a sparse array or load a sparse matrix from the input."
    
  def readyElasticThickness(self):
    """
    Is run at the right time to either import elastic thickness values from 
    an array, format them as a grid, or keep them as a scalar
    """
    # Need method to know whether you need a scalar or grid
    try:
      self.method
    except:
      self.method = None
      
    # If method is defined as something real
    if self.method:
      if self.method is 'FD':
        if type(self.coeff) is np.ndarray:
          if prod(array.shape) == 1:
            # 0D array only works for a constant Te value
            self.coeff = self.coeff[0]
          else:
            # Make sure array is proper size
            try:
              # requires self.q0 to be defined
              self.q0
            except:
              self.q0 = None
            # Only do this check if self.q0 is defined, so you can get its 
            # error check message
            if self.q0:
              # Checks that coeff array is of a workable size
              TeArraySizeCheck(self)
        elif np.isscalar(self.Te):
          try:
            self.q0 # need a grid size to be defined
          except:
            # Just exit and try again if this is not defined yet
            # But warn just in case
            self.q0 = None
            print "Attempting to define a scalar Te for a finite difference solution"
            print "Before the size of the domain (via q0) is defined."
            print "Unless q0 is defined before Te is first needed, program will crash."
          # If self.q0 is defined, convert the scalar to an array
          if self.q0:
            q0s = np.array(q0.shape)
            TeShapeTup = (())
            for i in len(q0s): # work for the number of dimensions you have
              # Make shape tuple
              TeShapeTup += q0s[i]+2
            # And apply it to self.Te to make it be a grid
            self.Te = self.Te * np.ones(TeShapeTup)
        # Otherwise, check if it is a string to a file path
        elif type(self.Te) is types.StringType:
          if self.Te[-4:] == '.npy':
            # Load binary array
            self.Te = np.load(self.Te)
          else:
            # Otherwise, assume ASCII with nothing special
            self.Te = np.genfromtxt(self.Te)
        else:
          print "Can't recognize the input Te type as being able to be imported"
          print "into code or read as a filename."
          print "Exiting."
          sys.exit()
      else:
        # Any other method requires a scalar Te
        if type(self.coeff) is np.ndarray:
          print "Converting numpy array to scalar elastic thickness for your solution method."
          # Check if array is really a scalar
          if prod(array.shape) == 1:
            self.coeff = self.coeff[0]
          # Or if array is uniform
          elif (a == a.mean()).all():
            self.coeff = a.mean() # Should find out how to take the mean only once
          else:
            print "Cannot figure out how to make your elastic thickness into a scalar."
            print "Exiting."
            sys.exit()
          print "Te format conversion complete."
        elif np.isscalar(self.Te):
          # All good here
          pass
    else:
      pass # Will just try again next time, when self.method becomes defined (hopefully)

    
class PrattAiry(Isostasy):
  pass
