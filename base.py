import sys, ConfigParser
from numpy import loadtxt

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
  def initialize(self, filename):
    # Open parser and get what kind of model
    self.config = ConfigParser.ConfigParser()
    self.config.read(filename)
    self.model     = self.config.get("mode", "model")
    self.dimension = self.config.getint("mode", "dimension")
    
    # Parameters
    # From input file
    self.g = self.config.getfloat("parameter", "GravAccel")
    self.rho_m = self.config.getfloat("parameter", "MantleDensity")
    self.rho_fill = self.config.getfloat("parameter", "InfillMaterialDensity")
    self.drho = self.rho_m - self.rho_fill # Move to Flexure class?

    # Grid spacing
    # Unnecessary for PrattAiry, but good to keep along, I think, for use 
    # in model output and plotting.
    # From input file
    self.dx = self.config.getfloat("numerical", "GridSpacing")
    
    
    # Loading grid
    q0path = self.config.get("input", "Loads")
    self.q0 = loadtxt(q0path)


    # UNIVERSAL SETTER: START WITH WHAT EVERYONE NEEDS
    def set_value(self, value_key, value):
      # Model type
      if value_key == 'model':
        self.model = str(value)
      elif value_key =='dimension':
        self.dimension = int(value)
      # Parameters
      elif value_key == 'GravAccel':
        self.g = float(value)
      elif value_key == 'MantleDensity':
        self.rho_m = float(value)
        # Update drho
        if self.rho_fill:
          drho = self.rho_m - self.rho_fill # Move to Flexure class?
      elif value_key == 'InfillMaterialDensity':
        self.rho_fill = float(value)
        # Update drho
        if self.rho_m:
          drho = self.rho_m - self.rho_fill # Move to Flexure class?
      # Grid spacing
      elif value_key == 'GridSpacing':
        self.dx = float(value) # Never tried it with anything but an int, but
                               # should work
      # Loading grid
      elif value_key == 'Loads':
        self.q0 = value # Anything special needed to pass arrays?
    
    
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
    plotChoice = self.config.get("output", "Plot")
    if plotChoice == 'q0' or plotChoice == 'w' or plotChoice == 'both':
      if debug: print 'Plotting...'
      if self.dimension==1:
        if plotChoice == 'q0':
          self.lineplot(self.q0/(self.rho_m*self.g),
            'Load thickness, mantle equivalent [m]')
        elif plotChoice == 'w':
          self.lineplot(self.w,'Deflection [m]')
        elif plotChoice == 'both':
          self.linesubplots()
      elif self.dimension==2:
        if plotChoice == 'q0':
          self.surfplot(self.q0/(self.rho_m*self.g),
            'Load thickness, mantle equivalent [m]')
        elif plotChoice == 'w':
          self.surfplot(self.w,'Deflection [m]')
        elif plotChoice == 'both':
          self.surfsubplots()

  def linesubplots(self,figNum=1):
    from matplotlib.pyplot import plot, show, figure, subplot, xlabel, \
                                  ylabel, title
    from numpy import arange
    
    figure(figNum)

    subplot(211)
    title('Loads and Lithospheric Deflections',fontsize=20)
    plot(arange(0,self.dx*self.q0.shape[0],self.dx),self.q0/(self.rho_m*self.g))
    ylabel('Load thickness, mantle equivalent [m]')

    subplot(212)
    plot(arange(0,self.dx*self.w.shape[0],self.dx),self.w)
    xlabel('Distance along profile [m]',fontsize=16)
    ylabel('Deflection [m]')

    show()

  def lineplot(self,data,ytext,xtext='Distance along profile [m]',titletext='',
      fontsize=16,figNum=1):
    """
    Plot if you want to - for troubleshooting
    """
    from matplotlib.pyplot import plot, show, figure, xlabel, ylabel, title
    from numpy import arange
    
    figure(figNum)

    plot(arange(0,self.dx*data.shape[0],self.dx),data)

    xlabel(xtext,fontsize=16)
    ylabel(ytext,fontsize=16)
    title(titletext,fontsize=16)
    
    show()

  def surfplot(self,data,titletext,figNum=1):
    """
    Plot if you want to - for troubleshooting
    """
    from matplotlib.pyplot import imshow, show, figure, colorbar, title
    from numpy import arange
    
    figure(figNum)

    imshow(data) #,interpolation='nearest'
    colorbar()

    title(titletext,fontsize=16)
    
    show()

  def surfsubplots(self,figNum=1):
    from matplotlib.pyplot import imshow, show, figure, subplot, xlabel, \
                                  ylabel, title, colorbar
    from numpy import arange
    
    figure(figNum)

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
  def initialize(self, filename):
    super(Flexure, self).initialize(filename)
    
    # Parameters
    self.E  = self.config.getfloat("parameter", "YoungsModulus")
    self.nu = self.config.getfloat("parameter", "PoissonsRatio")

  ### need to determine its interface, it is best to have a uniform interface
  ### no matter it is 1D or 2D; but if it can't be that way, we can set up a
  ### variable-length arguments, which is the way how Python overloads functions.
  def FD(self):
    # Import Te grid for the finite difference solution
    Tepath = self.config.get("input", "ElasticThickness")
    self.Te = loadtxt(Tepath)

  ### need work
  def FFT(self):
    pass

  def SPA(self):
    # Define the (scalar) elastic thickness
    self.Te = self.config.getfloat("parameter", "ElasticThickness")

  def SPA_NG(self):
    # Define the (scalar) elastic thickness
    self.Te = self.config.getfloat("parameter", "ElasticThickness")

    
  # UNIVERSAL SETTER: LITHOSPHERIC ELASTIC PROPERTIES ADDED
  def set_value(self, value_key, value):
    # Inherit from higher-level setter
    super(Flexure, self).set_value(value_key, value)
    # Parameters
    if value_key == 'YoungsModulus':
      self.E  = value
    elif value_key == 'PoissonsRatio':
      self.nu = value
    # Elastic thickness: array or scalar  
    elif value_key == 'ElasticThickness':
      self.Te = value # How to dynamically type for scalar or array?
    
class PrattAiry(Isostasy):
  pass
