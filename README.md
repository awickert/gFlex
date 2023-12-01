[![CSDMS Component](https://custom-icon-badges.demolab.com/badge/CSDMS-Component-2473c2?logo=csdms&style=for-the-badge)](https://csdms.colorado.edu/wiki/Model:GFlex)

[![Build Status](https://travis-ci.org/awickert/gFlex.svg?branch=master)](https://travis-ci.org/awickert/gFlex)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5034652.svg)](https://doi.org/10.5281/zenodo.5034652)

# gFlex

***Multiple methods to solve elastic plate flexure, designed for applications to Earth's lithosphere.***

These instructions are meant to take an user familiar with computers but new to (or a beginner with) Python through the basics of how to get gFlex to work. The Python scripting part towards the end should be pretty straightforward as well, insofar as information is provided on how to get and set the chosen values inside gFlex. *Please leave a message if you have trouble working with gFlex; your comments could assist both you and the more general improvement of this documentation.*

When you use gFlex, please cite:

**Wickert, A. D. (2016), [Open-source modular solutions for flexural isostasy: gFlex v1.0](https://www.geosci-model-dev.net/9/997/2016/gmd-9-997-2016.html), *Geosci. Model Dev.*, *9*(3), 997–1017, doi:10.5194/gmd-9-997-2016.**

## Download and Installation

#### Python

gFlex has been tested on **Python 3.10+**.

In order to run properly, gFlex requires the following Python dependencies:
* numpy
* scipy
* matplotlib

*For users who are new to Python, follow these directions to install the Python interpreters onto your computer.*

###### Linux

Use your package manager to download and install the required Python packages. For Debian/Ubuntu, it will be something like:

```bash
# Basic packages
sudo apt-get install \
python python-numpy python-scipy \
python-setuptools python-matplotlib

# pip (recommended for automatic installs via setuptools)
sudo apt-get install python-pip

# iPython console -- very useful (optional)
sudo apt-get install ipython

# Sypder IDE (I don't personally use it but many others like it: optional)
sudo apt-get install spyder
```

###### Windows

Download [**python(x,y)**](https://code.google.com/p/pythonxy/wiki/Downloads) or another full-featured distribution such as **Anaconda**; both of these distributions have been tested successfully with gFlex. Python(x,y) and several others also contain the required packages (including the numerical libraries), the iPython console, and the Spyder IDE; [**Spyder**](https://code.google.com/p/spyderlib/) is a nice IDE that will provide a familiar-looking interface for users accustomed to Matlab.

###### Mac

The current recommendation is to use a package manager like [**homebrew**](http://brew.sh/). With this you can install Python, and then move on to using **pip** (or **homebrew**) to install the Python modules. A good introduction to this can be found here: http://www.thisisthegreenroom.com/2011/installing-python-numpy-scipy-matplotlib-and-ipython-on-lion. See the **Linux** instructions for the list of packages that you will need; after installing pip, these commands can be substituted as follows, e.g.,
```bash
# Homebrew
sudo brew install python-numpy
# Pip
pip install -r requirements.txt
```

Recent efforts to download Python distributions (both **Anaconda** and **Enthought**) have not met with success with both gFlex and GRASS, though **Anaconda** has been tested successfully with Windows. As a result, it should be more successful to keep the Python packages managed better by something like **homebrew** with **pip**.

#### gFlex

##### Downloading and Installing in One Step from PyPI using pip

gFlex is downloadable from the Python Package Index ([PyPI](https://pypi.python.org/pypi)); see https://pypi.python.org/pypi/gFlex.

If you have **pip**, you may simply type:
```bash
pip install
pip install gflex
# Or if the destination install folder requires sudo access
# (for UNIX-like systems)
sudo pip install gflex
# pip install gFlex works too -- install is caps-insensitive
```
and you will have a full, running copy of the latest release version of gFlex.

##### Downloading

gFlex may be downloaded here at GitHub, by either:
* Copying the link at right and pasting it into the command prompt as follows:
```bash
git clone <LINK>
```
* Downloading and extracting the compressed ZIP file (link at right)
* Clicking on the link to add gFlex to your local GitHub desktop app (for Windows or Mac)

# Installing

Install gFlex at the command prompt using [setuptools](https://pypi.python.org/pypi/setuptools). If you have administrator privileges, which *is often also the case when doing this install under Windows*, you may drop the "sudo". For standard Linux or Mac users, the "sudo" will remain necessary, and you will have to enter your administrator password for the program to be added to your local set of applications (e.g., as "/usr/local/bin/gflex").

```bash
# For standard Linux/Mac users:
sudo python setup.py install
# OR
sudo python setup.py develop # If you want the install to see instantly
                             # any changes made in the source repository

# For Windows users or Unix-type users with SuperUser privileges:
python setup.py install
# OR
python setup.py develop # If you want the install to see instantly
                        # any changes made in the source repository
```

## Running

Once gFlex is installed, it is possible to run it in four ways:
 1. With a configuration file
 2. Within a Python script
 3. Within GRASS GIS
 4. As part of the Landlab Earth-surface modeling framework, including an interface to the the Community Surface Dynamics Modeling System [Component Model Interface (CMI)](http://csdms.colorado.edu/wiki/CMI_Description)

For options 1 and 2, there are pre-built methods that can be selected along the way to visualize results. These use Python's Matplotlib plotting library. For option 3, GRASS GIS is used for visualization. In Option 4, output from Landlab can be visualized with Matplotlib, and output from CSDMS sets of models can be visualized using tools such as [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit/) ([CSDMS page about VisIt](http://csdms.colorado.edu/wiki/CMT_visualization)) and [ParaView](http://www.paraview.org/). ParaView also now has [Python bindings](http://www.paraview.org/python/), which can further be used to visualize outputs produced with any of these methods.

#### With configuration file

A configuration file can be generated to run gFlex; see examples in the **input/** directory. To run gFlex using this file, one simply opens a terminal window and types:

```bash
# run like this:
gflex <path-to-configuration-file>
```

This can be run from any directory, as the installation of gFlex adds the program "gflex" to the system path.

For help constructing configuration files, see the blank template files **input/template1D** and **input/template2D**, as well as the other examples found in the **input/** directory. The **input/** directory also contains **input/README.md**, which provides a further local description of the files available. **input/input_help** provides a longer explanation of what the parameters are, and is therefore reproduced immediately below for reference:

```Lisp
; input_help
; All units are SI. Not all entries are needed.
; Standard parameter values for Earth are included.

[mode]
; 1 (line) or 2 (surface) dimensions
dimension=2
; Solution method: FD (Finite Difference), FFT (Fast Fourier
; Transform, not yet implemented), SAS (Spatial domain analytical
; solutions), or SAS_NG (SPA, but do not require a uniform grid
; - NG = "no grid")
; For SAS_NG, 1D data must be provided and will be returned in
; two columns: (x,q0) --> (x,w). 2D data are similar, except
; will be of the form (x,y,[q0/in or w/out])
; I am working on gridded output for these, so this might change
; in the future.
; Both the FFT and SPA techniques rely on superposition
; of solutions, because they can be combined linearly, whether in
; the spectral or the spatial domain)
method=SPA
; Plate solutions can be:
;  * vWC1994 (best), or
;  * G2009 (from Govers et al., 2009; not bad, but not
;           as robust as vWC1994)
PlateSolutionType=vWC1994

[parameter]
YoungsModulus=65E9
PoissonsRatio=0.25
GravAccel=9.8
MantleDensity=3300
; This is the density of material (e.g., air, water)
; that is filling (or leaving) the hole that was
; created by flexure. If you do not have a constant
; density of infilling material, for example, at a
; subsiding shoreline, you must instead iterate (see
; [numerical], below).
InfillMaterialDensity=0
[input]
; space-delimited array of loads
; stresses (rho*g*h) if gridded (dx (and if applicable, dy)) will be applied
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
; dx [m]
GridSpacing_x=
;
; Boundary conditions can be:
; (FD): 0Slope0Shear, 0Moment0Shear, 0Displacement0Slope, Mirror, or Periodic
; For SAS or SAS_NG, NoOutsideLoads is valid, and no entry defaults to this
BoundaryCondition_West=
BoundaryCondition_East=
;
; Solver can be direct or iterative
Solver=
; Tolerance between iterations [m]
; If you have chosen an iterative solution type ("Solver"), it will iterate
; until this is the difference between two subsequent iterations.
; Set as 0 if you don't want to iterate
convergence=1E-3

[numerical2D]
; dy [m]
GridSpacing_y=
;
; Boundary conditions can be:
; (FD): 0Slope0Shear, 0Moment0Shear, 0Displacement0Slope, Mirror, or Periodic
; For SAS or SAS_NG, NoOutsideLoads is valid, and no entry defaults to this
BoundaryCondition_North=
BoundaryCondition_South=
;
; Flag to enable lat/lon input (true/false). By default, this is false
latlon=
; radius of planet [m], for lat/lon solutions
PlanetaryRadius=

[verbosity]
; true/false. Defaults to true.
Verbose=
; true/false. Defaults to false.
Debug=
; true/false -- total silence if true. Defaults to false.
Quiet=
```

#### Within a Python script (with or without a configuration file)

You may run gFlex from other Python programs. When you install it (above), this also produces a Python module that you may import to access it while scripting.

##### With no configuration file (recommended)
**input/run_in_script_2D.py**, reproduced below, is a good example of how to set the variables and run the model. This method requires no input file, as all of the values are set inside the Python script that imports gflex. This is essentially how the GRASS GIS interface was written, and is a way to embed the abilities of gFlex into another model. A one-dimensional example, **input/run_in_script_1D.py**, is also available.

```python
#! /usr/bin/env python

import gflex
import numpy as np
from matplotlib import pyplot as plt

flex = gflex.F2D()

flex.Quiet = False

flex.Method = 'FD' # Solution method: * FD (finite difference)
                   #                  * SAS (superposition of analytical solutions)
                   #                  * SAS_NG (ungridded SAS)
flex.PlateSolutionType = 'vWC1994' # van Wees and Cloetingh (1994)
                                   # The other option is 'G2009': Govers et al. (2009)
flex.Solver = 'direct' # direct or iterative
# convergence = 1E-3 # convergence between iterations, if an iterative solution
                     # method is chosen

flex.g = 9.8 # acceleration due to gravity
flex.E = 65E9 # Young's Modulus
flex.nu = 0.25 # Poisson's Ratio
flex.rho_m = 3300. # MantleDensity
flex.rho_fill = 0. # InfiillMaterialDensity

flex.Te = 35000.*np.ones((50, 50)) # Elastic thickness [m] -- scalar but may be an array
flex.Te[:,-3:] = 0.
flex.qs = np.zeros((50, 50)) # Template array for surface load stresses
flex.qs[10:40, 10:40] += 1E6 # Populating this template
flex.dx = 5000. # grid cell size, x-oriented [m]
flex.dy = 5000. # grid cell size, y-oriented [m]
# Boundary conditions can be:
# (FD): 0Slope0Shear, 0Moment0Shear, 0Displacement0Slope, Mirror, or Periodic
# For SAS or SAS_NG, NoOutsideLoads is valid, and no entry defaults to this
flex.BC_W = '0Displacement0Slope' # west boundary condition
flex.BC_E = '0Moment0Shear' # east boundary condition
flex.BC_S = '0Displacement0Slope' # south boundary condition
flex.BC_N = '0Displacement0Slope' # north boundary condition

# latitude/longitude solutions are exact for SAS, approximate otherwise
#latlon = # true/false: flag to enable lat/lon input. Defaults False.
#PlanetaryRadius = # radius of planet [m], for lat/lon solutions

flex.initialize()
flex.run()
flex.finalize()

# If you want to plot the output
flex.plotChoice='both'
# An output file for deflections could also be defined here
# flex.wOutFile =
flex.output() # Plots and/or saves output, or does nothing, depending on
              # whether flex.plotChoice and/or flex.wOutFile have been set
# TO OBTAIN OUTPUT DIRECTLY IN PYTHON, you can assign the internal variable,
# flex.w, to another variable -- or as an element in a list if you are looping
# over many runs of gFlex:
deflection = flex.w
```

##### With a configuration file

If you would like to use a Python script with a configuration file, this is also possible.

```python
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

To run gFlex inside of GRASS GIS 7, run the following commands from within a GRASS GIS session:

```bash
g.extension r.flexure
g.extension v.flexure
```

This will reach into the GRASS GIS subversion repository, download the source code, and install the packages. **r.flexure** is used for raster grids by either finite difference or analytical methods. **v.flexure** takes advantage of the ungridded analytical method to solve for flexure at an aribtrary set of load points, albeit limited to cases with constant elastic thickness. These are stored at and have help files located at, respectively:

* **r.flexure**
 * Source: http://trac.osgeo.org/grass/browser/grass-addons/grass7/raster/r.flexure
 * Manual page (HTML): http://grass.osgeo.org/grass70/manuals/addons/r.flexure.html
* **v.flexure**
 * Source: http://trac.osgeo.org/grass/browser/grass-addons/grass7/vector/v.flexure
 * Manual page (HTML): http://grass.osgeo.org/grass70/manuals/addons/v.flexure.html

When running **r.flexure**, it is important to ensure that the elastic thickness map is at or properly interpolated to the computational region (**g.region**) resolution before solving. A nearest-neighbor interpolated Te map will cause perceived gradients in elastic thickness to be very sharp, and this will strongly affect (and misdirect) the flexural solutions.

#### As part of Landlab and the CSDMS CMI

Landlab is an in-development (but nearing release) Earth-surface modeling framework built to facilitate easy integration of geomorphic, ecological, hydrological, geological, etc. Earth-surface related models to simulate and investigate the links between multiple processes. gFlex can be linked with Landlab, and the code to do this is available within the Landlab repository at https://github.com/landlab/landlab/tree/master/landlab/components/gFlex.

The Landlab interface to gFlex also provides gFlex with the Community Surface Dynamics Modeling System (CSDMS) [Component Model Interface (CMI)](http://csdms.colorado.edu/wiki/CMI_Description) interface. This allows it to be run as a coupled component across multiple programming languages and paradigms as part of the CSDMS community of models. For more information on model coupling with CSDMS, see the example presentation at http://csdms.colorado.edu/w/images/CSDMS_lecture7.pdf and the paper on the model coupling published by [Peckham et al., "A component-based approach to integrated modeling in the geosciences: The design of CSDMS"](http://www.sciencedirect.com/science/article/pii/S0098300412001252).

### Plotting

There are four plot choices, defined via `self.plotChoice`:
* `'q'`: plots the load in mantle-density-equivalent units of length
* `'w'`: plots the deflection in units of length
* `'both'`: plots both deflection and loads in separate panels of a 2-subplot figure
* `'combo'`: (1D only): plots lithospheric deflections and the deflected mantle-density-equivalent load atop it.
  * Note that the load does not affect the area above/below the datum filled when `rho_fill != 0`. This affects the buoyant balance associated with the motion of the plate, with no additional considerations for topogrpahy. If you would like to include topogrpahy, an iterative approach (e.g., finding areas below sea level, filling them, flexing, finding new areas below sea level, and so on) is recommended.

## Utilities

The "utilities" folder currently contains only one program, `flexural_wavelength_calculator.py`. Operating it is simple and fairly rudimentary: just edit the input variables directly in the calculator Python file, and then run it to see what the flexural parameter, first zero-crossing point (on the load-side of the forebulge), and the flexural wavelength.
