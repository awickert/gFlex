import sys, ConfigParser, os
from numpy import loadtxt, ones, array
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
      if optional:
        # Must be first, so we don't exit program if optional flag is raised
        var = None
        print 'No value entered for optional parameter "' + name + '" in category "' + category + '" in input file.'
        print 'No action related to this optional parameter will be taken.'
      elif vartype == 'float':
        var = self.config.getfloat(category,name)
      elif vartype == 'string' or vartype == 'str':
        var = self.config.get(category,name)
        if var == "":
          print "Input strings cannot be empty."
          sys.exit()
      elif vartype == 'integer' or vartype == 'int':
        var = self.config.getint(category,name)
      else:
        print "Please enter 'float' or 'string' (or 'str') for vartype"
        sys.exit() # Won't exit, but will lead to exception
      return var
    except:
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
  
    print "" # Blank line at start of run

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
      self.q0 = loadtxt(q0path)
    except:
      try:
        # Then see if it is relative to the location of the input file
        self.q0 = loadtxt(self.inpath + q0path)
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
      print self.rho_m
    elif value_key == 'InfillMaterialDensity':
      self.rho_fill = value
      print self.rho_fill
    # Grid spacing
    elif value_key == 'GridSpacing':
      self.dx = value
    # Loading grid
    elif value_key == 'Loads':
      self.q0 = value # had [np.]float64(value) before
    # Dimensions
    elif value_key == "x":
      print "about to set x"
      self.x = value
      print "x set"
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

  # Save to file, if desired
  def outputDeflections(self):
    """
    Outputs a space-delimited ASCII grid of deflections if an output 
    directory is defined in the input file
    """
    try:
      # If wOutFile exists, has already been set by a setter
      self.wOutFile
      print "Output filename provided by setter"
    # Otherwise, it needs to be set by an input file
    except:
      try:
        self.wOutFile = self.configGet("string", "output", "DeflectionOut",optional=True)
        # If this exists and is a string, write output to a file
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
      print "Plotting " + self.plotChoice
      if debug: print 'Plotting...'
      if self.dimension==1:
        if self.plotChoice == 'q0':
          self.lineplot(self.q0/(self.rho_m*self.g),
            'Load thickness, mantle equivalent [m]')
        elif self.plotChoice == 'w':
          self.lineplot(self.w,'Deflection [m]')
        elif self.plotChoice == 'both':
          self.linesubplots()
        elif self.plotChoice == 'combo':
          #titletext=self.config.get(
          self.plotTogether()
      elif self.dimension==2:
        if self.plotChoice == 'q0':
          self.surfplot(self.q0/(self.rho_m*self.g),
            'Load thickness, mantle equivalent [m]')
        elif self.plotChoice == 'w':
          self.surfplot(self.w,'Deflection [m]')
        elif self.plotChoice == 'both':
          self.surfsubplots()
      else:
        print "Incorrect plotChoice input provided; proper input strings are:"
        print "q0, w, both, combo"
    else:
      try:
        # Double check that we didn't forget to include plotting selection
        self.plotChoice = self.config.get("output", "Plot")
        self.plotting()
      except:
        self.plotChoice = None
        print "No plotting today"

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
      Tepath = self.config.get("input", "ElasticThickness")
      
      # No grid?
      if len(Tepath) == 0:
        if debug: print "Trying to use the scalar elastic thickness"
        # Is there are scalar file?
        try:
          # No Te grid path defined, so going for scalar Te
          TeScalar = self.config.getfloat("parameter", "ElasticThickness")
          q0shape = array(self.q0.shape)
          for i in range(len(q0shape)):
            q0shape[i] += 2 # padding for numerical solution
          self.Te = TeScalar*ones(q0shape)
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
          self.Te = loadtxt(Tepath)
          print "Loading elastic thickness array from provided file"
        except:
          try:
            # Then see if it is relative to the location of the input file
            self.Te = loadtxt(self.inpath + Tepath)
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
          print "Exiting"
          print TeScalar
          sys.exit()
        # Otherwise, all set!
                  
  ### need work
  def FFT(self):
    pass

  def SPA(self):
    if self.filename:
      # Define the (scalar) elastic thickness
      self.Te = self.config.getfloat("parameter", "ElasticThickness")

  def SPA_NG(self):
    if self.filename:
      # Define the (scalar) elastic thickness
      self.Te = self.config.getfloat("parameter", "ElasticThickness")

    
  # UNIVERSAL SETTER: LITHOSPHERIC ELASTIC PROPERTIES AND SOLUTION METHOD
  def set_value(self, value_key, value):
    # Parameters
    if value_key == 'YoungsModulus':
      self.E  = value
    elif value_key == 'PoissonsRatio':
      self.nu = value
    # Elastic thickness: array or scalar  
    elif value_key == 'ElasticThickness':
      self.Te = value # How to dynamically type for scalar or array?
                      # Python is smart enough to handle it.
    elif value_key == 'method':
      self.method = value
      print "method set"
    # Inherit from higher-level setter, if not one of these
    # Currently not inheriting from the whichModel() setter, as I
    # figure that has to be done right away or not at all
    super(Flexure, self).set_value(value_key, value)
    
class PrattAiry(Isostasy):
  pass
