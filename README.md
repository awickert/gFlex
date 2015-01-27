gFlex
=====

Multiple methods to solve elastic beam and plate flexure, designed for applications to Earth's lithosphere.

To run gflex, first install the program at the command prompt using setuptools.

```
sudo python setup.py install
```

Then you can generate a configuration file (see examples in the **gflex/input** directory) and input data and run gFlex. It can produce plots and ASCII files with model outputs. Good files to look up include **gflex/input/template1D**, **gflex/input/template2D**, and **gflex/input/input_help**. The last of these, **input_help**, provides a longer explanation of what the parameters are.

```
# run like this:
gflex <path-to-configuration-file>
```

To run gFlex inside of GRASS GIS, you may copy **r.flexure** to the **scripts** directory of your GRASS GIS installation. This isn't the real way to install GRASS GIS add-ons, but it works for the moment until gFlex is complete enough to be submitted to the GRASS GIS add-ons repository.

To run gFlex within the CSDMS environment, see **gflex_bmi.py**.

To run gFlex from other Python programs, simply add code like you will find in **gflex_copy_paste.py**, and see **gflex/input/run_in_script** as an example:

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

