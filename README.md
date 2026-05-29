[![CSDMS Component][csdms_badge]][csdms_gflex] [![DOI][doi_badge]][doi_link] [![Test][test_badge]][test_workflow]

# gFlex

***Multiple methods to solve elastic plate flexure, designed for applications to Earth's lithosphere.***

These instructions are meant to take a user familiar with computers but new to (or a beginner with) Python through the basics of how to get gFlex to work. The Python scripting part towards the end should be pretty straightforward as well, insofar as information is provided on how to get and set the chosen values inside gFlex. *Please leave a message if you have trouble working with gFlex; your comments could assist both you and the more general improvement of this documentation.*

When you use gFlex, please cite:

**Wickert, A. D. (2016), [Open-source modular solutions for flexural isostasy: gFlex v1.0][paper_doi], *Geosci. Model Dev.*, *9*(3), 997–1017, doi:10.5194/gmd-9-997-2016.**

If you additionally want an up-to-date citation for the latest source-code release, please see that given in [CITATION.cff](CITATION.cff).

## Documentation

Full documentation, including a configuration file parameter reference, accuracy benchmarks, and the full API, is available at **[gflex.readthedocs.io](https://gflex.readthedocs.io)**.

## Installation

```bash
pip install gflex
```

gFlex requires **Python ≥ 3.10**. Dependencies (numpy, scipy, matplotlib, pyyaml) are installed automatically.

For a development install from source:

```bash
git clone https://github.com/awickert/gFlex.git
cd gFlex
pip install -e .
```

## Running

Once gFlex is installed, it is possible to run it in four ways:
 1. With a configuration file
 2. Within a Python script
 3. Within GRASS GIS
 4. As part of the Landlab Earth-surface modeling framework, with a CSDMS Basic Model Interface (BMI)

For options 1 and 2, there are pre-built methods that can be selected along the way to visualize results. These use Python's Matplotlib plotting library. For option 3, GRASS GIS is used for visualization. In option 4, output from Landlab can be visualized with Matplotlib.

#### With configuration file

A configuration file can be generated to run gFlex; see examples in the **input/** directory. To run gFlex using this file, one simply opens a terminal window and types:

```bash
gflex <path-to-configuration-file>
```

This can be run from any directory, as the installation of gFlex adds the program "gflex" to the system path.

Two configuration file formats are supported:

* **YAML** (`.yaml` / `.yml` extension) — recommended for new workflows. See **input/input_f1d.yaml** and **input/input_f2d.yaml** for complete examples.
* **INI** (any extension) — the original legacy format.

A minimal 2-D YAML configuration file looks like:

```yaml
mode:
  dimension: 2
  method: FD
  PlateSolutionType: vWC1994
parameter:
  YoungsModulus: 6.5e10
  PoissonsRatio: 0.25
  GravAccel: 9.8
  MantleDensity: 3300
  InfillMaterialDensity: 0
input:
  Loads: path/to/loads.txt
  ElasticThickness: path/to/Te.txt
output:
  Plot: both
numerical:
  GridSpacing_x: 4000
  BoundaryCondition_West: 0Moment0Shear
  BoundaryCondition_East: 0Displacement0Slope
  Solver: direct
  ConvergenceTolerance: 1.0e-3
numerical2D:
  GridSpacing_y: 4000
  BoundaryCondition_North: Mirror
  BoundaryCondition_South: 0Slope0Shear
```

For a full parameter reference, see the [Configuration Files](https://gflex.readthedocs.io/en/latest/configuration.html) page on ReadTheDocs. The annotated INI template **input/input_help** is also reproduced below for quick reference:

```ini
; input_help
; All units are SI. Not all entries are needed.
; Standard parameter values for Earth are included.

[mode]
; 1 (line) or 2 (surface) dimensions
dimension=2
; Solution method: FD (Finite Difference), FFT (spectral, 2-D only),
; SAS (Superposition of Analytical Solutions), or SAS_NG (SAS, but
; on an unstructured grid — NG = "no grid").
; For SAS_NG, 1D data must be provided and will be returned in
; two columns: (x,q0) --> (x,w). 2D data are similar, except
; will be of the form (x,y,[q0/in or w/out]).
method=SAS
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
; subsiding shoreline, you must instead iterate.
InfillMaterialDensity=0

[input]
; space-delimited array of loads
; stresses (rho*g*h) if gridded (dx (and if applicable, dy)) will be applied
;   to convert them into forces
; forces (rho*g*h*Area) if not gridded (SAS_NG)
; If the solution method (above) is selected as "SAS_NG", then this file
; will actually be of the format (x,[y],q0) and the code will sort it out.
Loads=q0_sample/2D/central_square_load.txt
;
; scalar value or space-delimited array of elastic thickness(es) [m]
; array required for finite difference solutions with variable Te
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
; Set as 0 if you don't want to iterate.
ConvergenceTolerance=1E-3

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

#### Finite-difference boundary conditions

Five boundary conditions are available for FD solutions (see also Table 1 in Wickert, 2016):

| Name | Condition | Physical interpretation |
|------|-----------|------------------------|
| `0Displacement0Slope` | w = 0, dw/dx = 0 | Plate is pinned to zero deflection and zero slope at the boundary |
| `0Moment0Shear` | d²w/dx² = d³w/dx³ = 0 | Broken plate: free cantilever end with no moment or shear |
| `0Slope0Shear` | dw/dx = d³w/dx³ = 0 | Plate is level at the boundary but free to deflect there; no shear transmitted |
| `Mirror` | w(b − x) = w(b + x) | Mirror-symmetry plane — model only half of a symmetric system |
| `Periodic` | w(0) = w(L) | Wrap-around: the domain tiles infinitely in both directions |

For SAS and SAS_NG, `NoOutsideLoads` (or a blank entry) is used instead; the plate is assumed undeflected at infinity.

**FD boundary-condition warnings:** when running F1D or F2D with the finite-difference solver, gFlex issues `UserWarning` messages for `'0Moment0Shear'` (free broken plate end — verify a rifted margin is intended), `'0Slope0Shear'` (no clear geological analog), and when the nearest loaded cell is within one flexural wavelength of a `'0Displacement0Slope'` boundary (the forebulge would be suppressed). See the [API reference](https://gflex.readthedocs.io/en/latest/api.html#fd-boundary-condition-warnings) for how to suppress or re-enable these warnings.

**A note on `0Slope0Shear`:** the label in Wickert (2016) is "free displacement of a horizontally clamped boundary." The plate is forced to be exactly level at the boundary (dw/dx = 0) — as if it were clamped against rotation — while its vertical position is unconstrained and no shear force is transmitted (d³w/dx³ = 0). This is superficially similar to `Mirror` (both enforce zero slope at the boundary), but `0Slope0Shear` uses a different finite-difference stencil and the two produce noticeably different solutions even far from the boundary. For symmetry problems — e.g., modelling half of a symmetric mountain range or ice sheet — `Mirror` is the more accurate choice.

#### Within a Python script (with or without a configuration file)

You may run gFlex from other Python programs. When you install it (above), this also produces a Python module that you may import to access it while scripting.

##### With no configuration file (recommended)
**input/run_in_script_2D.py**, reproduced below, is a good example of how to set the variables and run the model. This method requires no input file, as all of the values are set inside the Python script that imports gflex. This is essentially how the GRASS GIS interface was written, and is a way to embed the abilities of gFlex into another model. A one-dimensional example, **input/run_in_script_1D.py**, is also available.

```python
#! /usr/bin/env python

import gflex
import numpy as np

flex = gflex.F2D()

flex.Quiet = False

flex.Method = 'FD' # Solution method: * FD (finite difference)
                   #                  * FFT (spectral, 2-D only)
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
flex.rho_fill = 0. # InfillMaterialDensity

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

# Optional: in-plane stresses [Pa] (supported by FD and FFT solvers)
#flex.sigma_xx = 0.  # east–west compression/tension
#flex.sigma_yy = 0.  # north–south compression/tension
#flex.sigma_xy = 0.  # shear — couples x and y deflection

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

# To use a configuration file (INI or YAML):
filename = 'input/input_f1d.yaml'
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
obj.plotChoice='w'
obj.output()
```


#### Within GRASS GIS

To run gFlex inside of GRASS GIS 8, install the addons from within a GRASS GIS session:

```bash
g.extension r.flexure
g.extension v.flexure
```

**r.flexure** is used for raster grids by either finite difference or analytical methods. **v.flexure** takes advantage of the ungridded analytical method to solve for flexure at an arbitrary set of load points, albeit limited to cases with constant elastic thickness. The source code and manual pages are available in the [GRASS GIS addons repository](https://github.com/OSGeo/grass-addons).

When running **r.flexure**, it is important to ensure that the elastic thickness map is at or properly interpolated to the computational region (**g.region**) resolution before solving. A nearest-neighbor interpolated Te map will cause perceived gradients in elastic thickness to be very sharp, and this will strongly affect (and misdirect) the flexural solutions.

#### As part of Landlab and CSDMS

[Landlab](https://landlab.github.io) is an Earth-surface modeling framework built to facilitate easy integration of geomorphic, ecological, hydrological, geological, and other Earth-surface models. gFlex can be used as a Landlab component; see the [Landlab repository](https://github.com/landlab/landlab) for details.

gFlex also implements the [CSDMS Basic Model Interface (BMI)](https://bmi.readthedocs.io), enabling it to be coupled with other models in the CSDMS framework. The BMI wrapper is available as `gflex.BmiGflex` and requires the optional `bmipy` dependency:

```bash
pip install gflex[bmi]
```

### Plotting

There are four plot choices, defined via `self.plotChoice`:
* `'q'`: plots the load in mantle-density-equivalent units of length
* `'w'`: plots the deflection in units of length
* `'both'`: plots both deflection and loads in separate panels of a 2-subplot figure
* `'combo'`: (1D only): plots lithospheric deflections and the deflected mantle-density-equivalent load atop it.
  * Note that the load does not affect the area above/below the datum filled when `rho_fill != 0`. This affects the buoyant balance associated with the motion of the plate, with no additional considerations for topography. If you would like to include topography, an iterative approach (e.g., finding areas below sea level, filling them, flexing, finding new areas below sea level, and so on) is recommended.

## Utilities

The **input/** directory contains several example scripts. The public API also provides standalone utility functions importable from `gflex`:

* `flexural_wavelengths(Te, ...)` — computes the flexural parameter α, first zero-crossing, and flexural wavelength for a given elastic thickness; useful for choosing grid spacing and domain size.

Domain-padding utilities (reduce spurious boundary effects when using variable *Te*):

*2-D (F2D):*
* `recommended_pad_width(Te, dx, ...)` — returns the recommended padding width (in cells).
* `smooth_pad_Te(Te, pad_width, ...)` — extends a 2-D variable-*Te* array with a smooth linear taper.
* `pad_domain(Te, qs, dx, ...)` — pads both the 2-D elastic thickness and load arrays and returns the padding width.

*1-D (F1D):*
* `recommended_pad_width_1d(Te, dx, ...)` — returns the recommended 1-D padding width (in cells).
* `smooth_pad_Te_1d(Te, pad_width, ...)` — extends a 1-D variable-*Te* array with a smooth linear taper.
* `pad_domain_1d(Te, qs, dx, ...)` — pads both the 1-D elastic thickness and load arrays and returns the padding width.

See the [API reference](https://gflex.readthedocs.io/en/latest/api.html) for full documentation of these functions.


[csdms_badge]: https://custom-icon-badges.demolab.com/badge/CSDMS-Component-2473c2?logo=csdms&style=for-the-badge
[csdms_gflex]: https://csdms.colorado.edu/wiki/Model:GFlex
[doi_badge]: https://zenodo.org/badge/DOI/10.5281/zenodo.5034651.svg
[doi_link]: https://doi.org/10.5281/zenodo.5034651
[test_badge]: https://github.com/awickert/gflex/actions/workflows/test.yml/badge.svg
[test_workflow]: https://github.com/awickert/gflex/actions/workflows/test.yml
[paper_doi]: https://www.geosci-model-dev.net/9/997/2016/gmd-9-997-2016.html
