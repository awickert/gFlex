gFlex
=====

***Multiple methods to solve elastic plate flexure, designed for applications to Earth's lithosphere.***

These instructions are menat to take an user familiar with computers but new to (or a beginner with) Python through the basics of how to get gFlex to work. The Python scripting part towards the end should be pretty straightforward as well, insofar as information is provided on how to get and set the chosen values inside gFlex. *Please leave a message if you have trouble working with gFlex; your comments could assist both you and the more general improvement of this documentation.*

## Installation

#### Python

For users who are new to Python, follow these directions to install the Python interpreters onto your computer. It might not be necessary to install all of the packages used by gFlex right away, as Python's "setuptools" does a good job of finding and installing the package dependencies as well, but these instructions are for a full initial install of everything.

###### Linux

Use your package manager to download and install the required Python packages. For Debian/Ubuntu, it will be something like:

```
# Basic packages
sudo apt-get install python python-numpy python-scipy python-setuptools python-configparser python-matplotlib

# iPython console -- very useful (optional)
sudo apt-get install ipython

# Sypder IDE (I don't personally use it but many others like it: optional)
sudo apt-get install spyder
```

You may not need all of these packages depending on how much is packaged within setuptools... 

###### Windows

Download [**pythonxy**](https://code.google.com/p/pythonxy/wiki/Downloads) or another full-featured distribution. Pythonxy and several others also contain the required packages (including the numerical libraries), the iPython console, and the Spyder IDE; [**Spyder**](https://code.google.com/p/spyderlib/) is a nice IDE that will provide a familiar-looking interface for users accustomed to Matlab.

###### Mac

The current recommendation is to use a package manager like [**homebrew**](http://brew.sh/) and follow the Linux instructions; recent efforts to download Python distributions (both **Anaconda** and **Enthought**) have not met with success with both gFlex and GRASS, though **Anaconda** has been tested successfully with Windows. As a result, it should be more successful to keep the Python packages managed better by something like **homebrew**.

#### gFlex

Install gFlex at the command prompt using [setuptools](https://pypi.python.org/pypi/setuptools). If you have administrator privileges, which *is often also the case when doing this install under Windows*, you may drop the "sudo". For standard Linux or Mac users, the "sudo" will remain necessary, and you will have to enter your administrator password for the program to be added to your local set of applications (e.g., as "/usr/local/bin/gflex").

```
# For standard Linux/Mac users:
sudo python setup.py install

# For Windows users or Unix-type users with SuperUser privileges:
python setup.py install
```

## Running

Once gFlex is installed, it is possible to run it in four ways:
 1. With a configuration file
 2. Within a Python script
 3. Within GRASS GIS
 4. As a coupled component of a set of models via the Community Surface Dynamics Modeling System [Component Model Interface (CMI)](http://csdms.colorado.edu/wiki/CMI_Description)

For options 1 and 2, there are pre-built methods that can be selected along the way to visualize results. These use Python's Matplotlib plotting library. For option 3, GRASS GIS is used for visualization. In Option 4, output from CSDMS sets of models can be visualized using tools such as [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit/) ([CSDMS page about VisIt](http://csdms.colorado.edu/wiki/CMT_visualization)) and [ParaView](http://www.paraview.org/). Paraview also now has [Python bindings](http://www.paraview.org/python/), which can further be used to visualize outputs produced wiht any of these methods.

#### With configuration file

A configuration file can be generated to run gFlex; see examples in the **input/** directory. To run gFlex using this file, one simply opens a terminal window and types:

```
# run like this:
gflex <path-to-configuration-file>
```

For help constructing configuration files, see the blank template files **input/template1D** and **input/template2D**, as well as the other examples found in the **input/** directory. The **input/** directory also contains **input/README.md**, which provides a further local description of the files available. **input/input_help** provides a longer explanation of what the parameters are, and is therefore reproduced immediately below for reference:

```
; input_help
; All units are SI. Not all entries are needed.
; Standard parameter values for Earth are included.

[mode]
dimension=2 ; 1 (line) or 2 (surface) dimensions
method=SPA ; Solution method: FD (Finite Difference), FFT (Fast Fourier 
;          ; Transform, not yet implemented), SAS (Spatial domain analytical 
;          ; solutions), or SAS_NG (SPA, but do not require a uniform grid
;          ; - NG = "no grid")
;          ; For SAS_NG, 1D data must be provided and will be returned in 
;          ; two columns: (x,q0) --> (x,w). 2D data are similar, except
;          ; will be of the form (x,y,[q0/in or w/out])
;          ; I am working on gridded output for these, so this might change
;          ; in the future.
;          ; Both the FFT and SPA techniques rely on superposition 
;          ; of solutions, because they can be combined linearly, whether in 
;          ; the spectral or the spatial domain)
;
PlateSolutionType=vWC1994 ; Plate solutions can be:
;                         ; * vWC1994 (best), or
;                         ; * G2009 (from Govers et al., 2009; not bad, but not 
;                         ;          as robust as vWC1994)

[parameter]
YoungsModulus=6.5E10
PoissonsRatio=0.25
GravAccel=9.8
MantleDensity=3300
InfillMaterialDensity=0 ; This is the density of material (e.g., air, water) 
;                       ; that is filling (or leaving) the hole that was 
;                       ; created by flexure. If you do not have a constant 
;                       ; density of infilling material, for example, at a 
;                       ; subsiding shoreline, you must instead iterate (see
;                       ; [numerical], below).

[input]
; space-delimited array of loads
; stresses (rho*g*h) if gridded (dx (and if applicable, dy) will be applied
;   to convert them into masses
; forces (rho*g*h*Area) if not gridded (SAS_NG)
; If the solution method (above) is selected as "SAS_NG", then this file
; will actually be of the format (x,[y],q0) and the code will sort it out.
; (Once again, working on a gridded output option for ungridded inputs)
Loads=q0_sample/2D/central_square_load.txt
;
; scalar value or space-delimited array of elastic thickness(es)
; array used for finite difference solutions
ElasticThickness=Te_sample/2D/10km_const.txt
;
; xw and yw are vectors of desired output points for the SAS_NG method.
; If they are not specified and a SAS_NG solution is run, the solution will be 
; calculated at the points with the loads.
; they are ignored if a different solution method is chosen.
xw=
yw=

[output]
; DeflectionOut is for writing an output file. 
; If this is blank, no output is printed.
; Otherwise, a space-delimited ASCII file of 
; outputs is with this file name (and path).
DeflectionOut=tmpout.txt
;
; Acceptable inputs to "Plot" are q0 (loads), w (deflection), or both; any 
; other entry here will result in no plotting.
; Automatically plots a 1D line or 2D surface based on the choice 
; of "dimension" variable in [mode]
Plot=both

[numerical]
GridSpacing_x= ; dx [m]
;
; Boundary conditions can be:
; (FD): 0Slope0Shear, 0Moment0Shear, Dirichlet0, Mirror, or Periodic
; For SAS or SAS_NG, NoOutsideLoads is valid, and no entry defaults to this
BoundaryCondition_West=
BoundaryCondition_East=
;
; Solver can be direct or iterative
Solver=
;
; If you have chosen an iterative solution type ("Solver"), it will iterate
; until this is the difference between two subsequent iterations.
; Set as 0 if you don't want to iterate
convergence=1E-3 ; Tolerance between iterations [m]

[numerical2D]
GridSpacing_y= ; dy [m]
;
; Boundary conditions can be:
; (FD): 0Slope0Shear, 0Moment0Shear, Dirichlet0, Mirror, or Periodic
; For SAS or SAS_NG, NoOutsideLoads is valid, and no entry defaults to this
BoundaryCondition_North=
BoundaryCondition_South=
; 
; Flag to enable lat/lon input. By default, this is false
latlon= ; true/false: flag to enable lat/lon input. Defaults False.
PlanetaryRadius= ; radius of planet [m], for lat/lon solutions

[verbosity]
Verbose= ; true/false. Defaults to True.
Debug= ; true/false. Defaults to False.
Quiet= ; true/false -- total silence if True. Defaults to False.
```

#### Within a Python script (with or without a configuration file)

You may run gFlex from other Python programs. When you install it (above), this also produces a Python module that you may import to access it while scripting.

##### With no configuration file
**gflex/input/run_in_script_2D.py**, reproduced below, is a good example of how to set the variables and run the model. This method requires no input file, as all of the values are set inside the Python script that imports gflex. This is essentially how the GRASS GIS interface was written, and is a way to embed the abilities of gFlex into another model. A one-dimensional example, **gflex/input/run_in_script_1D.py**, is also available.

```
import gflex
import numpy as np
from matplotlib import pyplot as plt

flex = gflex.F2D()

flex.Quiet = False

flex.Method = 'FD'
flex.PlateSolutionType = 'vWC1994'
flex.Solver = 'direct'

flex.g = 9.8 # acceleration due to gravity
flex.E = 65E10 # Young's Modulus
flex.nu = 0.25 # Poisson's Ratio
flex.rho_m = 3300. # MantleDensity
flex.rho_fill = 0. # InfiillMaterialDensity

flex.Te = 35000. # Elastic thickness -- scalar but may be an array
flex.qs = np.zeros((50, 50)) # Template array for surface load stresses
flex.qs[10:40, 10:40] += 1E6 # Populating this template
flex.dx = 5000.
flex.dy = 5000.
flex.BC_W = 'Dirichlet0' # west boundary condition
flex.BC_E = '0Moment0Shear' # east boundary condition
flex.BC_S = 'Periodic' # south boundary condition
flex.BC_N = 'Periodic' # north boundary condition

flex.initialize()
flex.run()
flex.finalize()

# If you want to plot the output
flex.plotChoice='both'
# An output file could also be defined here
# flex.wOutFile = 
flex.output() # Plots and/or saves output, or does nothing, depending on
              # whether flex.plotChoice and/or flex.wOutFile have been set
```

##### With a configuration file

If you would like to use a Python script with a configuration file, this is also possible.

```
import gflex

# To use a configuration file:
filename = '../gflex/input/input_f1d_test' # it works for usage (1) and (2)
obj = gflex.WhichModel(filename)

## SET MODEL TYPE AND DIMENSIONS HERE ##
########################################
if obj.dimension == 1:
  obj = gflex.F1D(filename)
elif obj.dimension == 2:
  obj = gflex.F2D(filename)

# Then run the code!
obj.initialize(filename)
obj.run()
obj.finalize()

# Standalone plotting output if you so desire
flex.plotChoice='w'
obj.output()
```


#### Within GRASS GIS

(Outdated)
To run gFlex inside of GRASS GIS, you may copy **r.flexure** to the **scripts** directory of your GRASS GIS installation. This isn't the real way to install GRASS GIS add-ons, but it works for the moment until gFlex is complete enough to be submitted to the GRASS GIS add-ons repository.

When running **r.flexure**, it is important to ensure that the elastic thickness map is at or properly interpolated to the computational region (**g.region**) resolution before solving. A nearest-neighbor interpolated Te map will cause perceived gradients in elastic thickness to be very sharp, and this will strongly affect (and misdirect) the flexural solutions.

#### As part of the CSDMS CMI

To run gFlex within the CSDMS Component Model Interface (CMI) environment, see **gflex_bmi.py**.
