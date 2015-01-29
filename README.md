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

Use a package manager like [**homebrew**](http://brew.sh/) and follow the Linux instructions, or install a distribution following the Windows instructions. Either option should work. gFlex has not yet been tested on Mac, so please leave a message if you do so.

#### gFlex

Install gFlex at the command prompt using [setuptools](https://pypi.python.org/pypi/setuptools). If you have administrator privileges, which *is often also the case when doing this install under Windows*, you may drop the "sudo". For standard Linux or Mac users, the "sudo" will remain necessary, and you will have to enter your administrator password for the program to be added to your local set of applications (e.g., as "/usr/local/bin/gflex").

```
# For standard Linux/Mac users:
sudo python setup.py install

# For Windows users or Unix-type users with SuperUser privileges:
python setup.py install
```

## Running

#### With configuration file

Then you can generate a configuration file (see examples in the **input** directory) and input data and run gFlex. It can produce plots and ASCII files with model outputs. Good files to look up include **input/template1D**, **input/template2D**, and **input/input_help**. The last of these, **input_help**, provides a longer explanation of what the parameters are.

```
# run like this:
gflex <path-to-configuration-file>
```

To run gFlex inside of GRASS GIS, you may copy **r.flexure** to the **scripts** directory of your GRASS GIS installation. This isn't the real way to install GRASS GIS add-ons, but it works for the moment until gFlex is complete enough to be submitted to the GRASS GIS add-ons repository.

To run gFlex within the CSDMS environment, see **gflex_bmi.py**.

#### Within a Python script (with or without a configuration file)

To run gFlex from other Python programs, simply add code like you will find in **gflex_copy_paste.py**, and see **gflex/input/run_in_script.py** as an example, and my scratch file (if it is tracked on git) is **gflex/tests/interactive.py**:

```
import gflex

# If you want to use an input file:
filename = '../gflex/input/input_f1d_test' # it works for usage (1) and (2)
obj = gflex.WhichModel(filename)

## SET MODEL TYPE AND DIMENSIONS HERE ##
########################################
if obj.dimension == 1:
  obj = gflex.F1D(filename)
elif obj.dimension == 2:
  obj = gflex.F2D(filename)

# Othrwise just initialize with getters and setters
# e.g.,
# obj.set_value('GridSpacing_y', 50000)

# Then run the code!
obj.initialize(filename)
# You can use getters/setters here too
obj.run()
obj.finalize()

# Standalone plotting output if you so desire
obj.output()
```
