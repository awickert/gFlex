import sys, ConfigParser
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
  def whichModel(self, filename=None):
    self.filename = filename
    if self.filename:
      # Open parser and get what kind of model
      self.config = ConfigParser.ConfigParser()
      self.config.read(filename)
      self.model     = self.config.get("mode", "model")
      self.dimension = self.config.getint("mode", "dimension")

  def initialize(self, filename=None):
    self.filename = filename # Redundant with whichModel()?
    self.plotChoice = None # Default for no plotting
    self.wOutFile = None # Default with no output
    if self.filename:
      # Open parser if needed (whichModel() should have run already, so
      # shouldn't be needed.)
#      if self.config:
#        pass
#      else:
      self.config = ConfigParser.ConfigParser()
      self.config.read(filename)

      # Parameters
      # From input file
      self.g = self.config.getfloat("parameter", "GravAccel")
      self.rho_m = self.config.getfloat("parameter", "MantleDensity")
      self.rho_fill = self.config.getfloat("parameter", "InfillMaterialDensity")

      # Grid spacing
      # Unnecessary for PrattAiry, but good to keep along, I think, for use 
      # in model output and plotting.
      # No meaning for ungridded superimposed analytical solutions
      # From input file
      self.dx = self.config.getfloat("numerical", "GridSpacing")
      
      # Loading grid
      q0path = self.config.get("input", "Loads")
      self.q0 = loadtxt(q0path)

      # Plotting selection
      self.plotChoice = self.config.get("output","Plot")

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
      if self.wOutFile:
        print "Output filename provided by setter"
      else:  
        self.wOutFile = self.config.get("output", "DeflectionOut")
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
    print self.plotChoice
    if self.plotChoice:
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
        self.plotChoice = self.config.get("output", "Plot")
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
        # No "/1000" needed for FD implementation Te; already in km
        title(titletext + ', $T_e$ = ' + str((self.Te).mean()) + " km",fontsize=20)
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
    raise SystemExit # Stop program from running after printing only
                     # my error message

# class Flexure inherits Isostay and it overrides the __init__ method. It also
# define three different solution methods, which are implemented by its subclass.
class Flexure(Isostasy):
  def initialize(self, filename=None):
    super(Flexure, self).whichModel(filename)
    super(Flexure, self).initialize(filename)

    # Solution method
    if self.filename:
      self.method = self.config.get("mode", "method")
    
    # Parameters
    self.drho = self.rho_m - self.rho_fill
    if self.filename:
      self.E  = self.config.getfloat("parameter", "YoungsModulus")
      self.nu = self.config.getfloat("parameter", "PoissonsRatio")
    
  ### need to determine its interface, it is best to have a uniform interface
  ### no matter it is 1D or 2D; but if it can't be that way, we can set up a
  ### variable-length arguments, which is the way how Python overloads functions.
  def FD(self):
    print "Finite Difference Solution Technique"
    if self.filename:
    # Import Te grid for the finite difference solution
      Tepath = self.config.get("input", "ElasticThickness")
      if len(Tepath) > 0:
        try:
          self.Te = loadtxt(Tepath)
          print "Loading elastic thickness array from provided file"
        except:
          print "Cannot find an ASCII elastic thickness file at Tepath; aborting"
          raise SystemExit #sys.exit() (sys not imported, I think...)
      else:
        q0shape = array(self.q0.shape)
        for i in range(len(q0shape)):
          q0shape[i] += 2 # padding for numerical solution
        self.Te = ones(q0shape)
        print "Using constant elastic thickness at provided value"
        
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
