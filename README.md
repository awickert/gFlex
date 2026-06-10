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

gFlex requires **Python ≥ 3.11**. Dependencies (numpy, scipy, matplotlib, pyyaml) are installed automatically.

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

Configuration files must be YAML (`.yaml` / `.yml` extension). See **input/input_f1d.yaml** and **input/input_f2d.yaml** for complete examples.

A minimal 2-D YAML configuration file looks like:

```yaml
mode:
  dimension: 2
  method: fd
parameter:
  youngs_modulus: 6.5e10
  poissons_ratio: 0.25
  gravitational_acceleration: 9.8
  mantle_density: 3300
  infill_material_density: 0
input:
  loads: path/to/loads.txt
  elastic_thickness: path/to/Te.txt
output:
  plot: both
numerical:
  grid_spacing_x: 4000
  boundary_condition_west: zero_moment_zero_shear
  boundary_condition_east: zero_displacement_zero_slope
numerical2D:
  grid_spacing_y: 4000
  boundary_condition_north: mirror
  boundary_condition_south: zero_slope_zero_shear
```

For a full parameter reference, see the [Configuration Files](https://gflex.readthedocs.io/en/latest/configuration.html) page on ReadTheDocs.

#### Finite-difference boundary conditions

Six boundary conditions are available for FD solutions (see also Table 1 in Wickert, 2016):

| Name | Aliases | Condition | Physical interpretation |
|------|---------|-----------|------------------------|
| `zero_displacement_zero_slope` | `clamped` | w = 0, dw/dx = 0 | Clamped end: zero deflection and zero slope (no rotation) |
| `zero_displacement_zero_moment` | `pinned` | w = 0, d²w/dx² = 0 | Simply supported (pinned) end: zero deflection, free to rotate |
| `zero_moment_zero_shear` | `free` | d²w/dx² = d³w/dx³ = 0 | Free (broken) plate end: no moment or shear force transmitted |
| `zero_slope_zero_shear` | `mirror` | dw/dx = d³w/dx³ = 0 | Plate is level at the boundary but free to deflect there; no shear transmitted — use for symmetry planes |
| `no_outside_loads` | `infinite` | w → 0 at ∞ | Infinite plate: the solver auto-pads the domain by one flexural wavelength |
| `periodic` | — | w(0) = w(L) | Wrap-around: the domain tiles infinitely in both directions |

For SAS and SAS_NG, `no_outside_loads` (or a blank entry) is always assumed; the analytical solution is inherently infinite-plate.

**FD boundary-condition warnings:** when running F1D or F2D with the finite-difference solver, gFlex issues `UserWarning` messages for `'zero_moment_zero_shear'` (free broken plate end — verify a rifted margin is intended), `'zero_slope_zero_shear'` (no clear geological analog), and when the nearest loaded cell is within one flexural wavelength of a `'zero_displacement_zero_slope'` boundary (the forebulge would be suppressed). See the [API reference](https://gflex.readthedocs.io/en/latest/api.html#fd-boundary-condition-warnings) for how to suppress or re-enable these warnings.

**A note on `zero_slope_zero_shear` / `mirror`:** the label in Wickert (2016) is "free displacement of a horizontally clamped boundary." The plate is forced to be exactly level at the boundary (dw/dx = 0) — as if it were clamped against rotation — while its vertical position is unconstrained and no shear force is transmitted (d³w/dx³ = 0). `mirror` is an alias for `zero_slope_zero_shear`; both produce identical solutions and are interchangeable. This BC is the natural choice for symmetry planes — e.g., modelling half of a symmetric mountain range or ice sheet.

#### Within a Python script (with or without a configuration file)

You may run gFlex from other Python programs. When you install it (above), this also produces a Python module that you may import to access it while scripting.

##### With no configuration file (recommended)
**input/run_in_script_2D.py**, reproduced below, is a good example of how to set the variables and run the model. This method requires no input file, as all of the values are set inside the Python script that imports gflex. This is essentially how the GRASS GIS interface was written, and is a way to embed the abilities of gFlex into another model. A one-dimensional example, **input/run_in_script_1D.py**, is also available.

```python
#! /usr/bin/env python

import gflex
import numpy as np

flex = gflex.F2D()

flex.quiet = False

flex.method = 'fd' # Solution method: * fd (finite difference)
                   #                  * fft (spectral)
                   #                  * sas (superposition of analytical solutions)
                   #                  * sas_ng (ungridded SAS)

flex.g = 9.8 # acceleration due to gravity
flex.E = 65E9 # Young's Modulus
flex.nu = 0.25 # Poisson's Ratio
flex.rho_m = 3300. # mantle_density
flex.rho_fill = 0. # infill_material_density

flex.T_e = 35000.*np.ones((50, 50)) # Elastic thickness [m] -- scalar but may be an array
flex.T_e[:,-3:] = 0.
flex.qs = np.zeros((50, 50)) # Template array for surface load stresses
flex.qs[10:40, 10:40] += 1E6 # Populating this template
flex.dx = 5000. # grid cell size, x-oriented [m]
flex.dy = 5000. # grid cell size, y-oriented [m]
# Boundary conditions (FD): zero_displacement_zero_slope, zero_displacement_zero_moment,
#   zero_moment_zero_shear, zero_slope_zero_shear, no_outside_loads, periodic
#   (aliases: clamped, pinned, free, mirror, infinite)
# For SAS or SAS_NG, no_outside_loads is always assumed; no entry defaults to this
flex.bc_west = 'zero_displacement_zero_slope' # west boundary condition
flex.bc_east = 'zero_moment_zero_shear' # east boundary condition
flex.bc_south = 'zero_displacement_zero_slope' # south boundary condition
flex.bc_north = 'zero_displacement_zero_slope' # north boundary condition

# latitude/longitude solutions are exact for SAS, approximate otherwise
#latlon = # true/false: flag to enable lat/lon input. Defaults False.
#planetary_radius = # radius of planet [m], for lat/lon solutions

# Optional: in-plane stresses [Pa] (supported by FD and FFT solvers)
#flex.sigma_xx = 0.  # east–west compression/tension
#flex.sigma_yy = 0.  # north–south compression/tension
#flex.sigma_xy = 0.  # shear — couples x and y deflection

flex.initialize()
flex.run()

# Read deflection BEFORE calling finalize (finalize clears w and qs)
deflection = flex.w  # assign to a variable for later use

# If you want to plot or save output, do it before finalize:
flex.plot_choice='both'
# An output file for deflections could also be defined here
# flex.w_out_file =
flex.output() # Plots and/or saves output, or does nothing, depending on
              # whether flex.plot_choice and/or flex.w_out_file have been set

flex.finalize()
```

##### With a configuration file

If you would like to use a Python script with a configuration file, this is also possible.

```python
import gflex

# Pass the YAML file to the constructor; F1D for 1-D, F2D for 2-D problems.
filename = 'input/input_f1d.yaml'
flex = gflex.F1D(filename)

flex.initialize()
flex.run()

# Read deflection before finalize (finalize clears w and qs)
deflection = flex.w

# Standalone plotting output if you so desire
flex.plot_choice = 'w'
flex.output()

flex.finalize()
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

There are four plot choices, defined via `self.plot_choice`:
* `'q'`: plots the load in mantle-density-equivalent units of length
* `'w'`: plots the deflection in units of length
* `'both'`: plots both deflection and loads in separate panels of a 2-subplot figure
* `'combo'`: (1D only): plots lithospheric deflections and the deflected mantle-density-equivalent load atop it.
  * Note that the load does not affect the area above/below the datum filled when `rho_fill != 0`. This affects the buoyant balance associated with the motion of the plate, with no additional considerations for topography. If you would like to include topography, an iterative approach (e.g., finding areas below sea level, filling them, flexing, finding new areas below sea level, and so on) is recommended.

## Utilities

The **input/** directory contains several example scripts. The public API also provides standalone utility functions importable from `gflex`:

* `flexural_wavelengths(Te, ...)` — computes the flexural parameter α, first zero-crossing, and flexural wavelength for a given elastic thickness; useful for choosing grid spacing and domain size.

Domain-padding utilities (reduce spurious boundary effects when using variable *Te*):

* `pad_domain(Te, qs, dx, ...)` — unified helper for both 1-D and 2-D; pads the elastic thickness and load arrays and returns the padding width.

Lower-level building blocks:

*2-D (F2D):*
* `recommended_pad_width(Te, dx, ...)` — returns the recommended padding width (in cells).
* `smooth_pad_Te(Te, pad_width, ...)` — extends a 2-D variable-*Te* array with a smooth linear taper.

*1-D (F1D):*
* `recommended_pad_width_1d(Te, dx, ...)` — returns the recommended 1-D padding width (in cells).
* `smooth_pad_Te_1d(Te, pad_width, ...)` — extends a 1-D variable-*Te* array with a smooth linear taper.

See the [API reference](https://gflex.readthedocs.io/en/latest/api.html) for full documentation of these functions.


[csdms_badge]: https://custom-icon-badges.demolab.com/badge/CSDMS-Component-2473c2?logo=csdms&style=for-the-badge
[csdms_gflex]: https://csdms.colorado.edu/wiki/Model:GFlex
[doi_badge]: https://zenodo.org/badge/DOI/10.5281/zenodo.10471939.svg
[doi_link]: https://doi.org/10.5281/zenodo.10471939
[test_badge]: https://github.com/awickert/gflex/actions/workflows/test.yml/badge.svg
[test_workflow]: https://github.com/awickert/gflex/actions/workflows/test.yml
[paper_doi]: https://www.geosci-model-dev.net/9/997/2016/gmd-9-997-2016.html
