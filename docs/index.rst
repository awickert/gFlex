gFlex
=====

**Multiple methods to solve elastic plate flexure, designed for applications
to Earth's lithosphere.**

gFlex computes lithospheric flexural isostasy — the bending of Earth's
elastic outer shell under surface loads such as ice sheets, sediment, lava
flows, or water.  Both one-dimensional (profile) and two-dimensional (map
view) solutions are supported, using either a finite-difference method (which
handles spatially variable elastic thickness) or superposition of analytical
solutions (fast, for constant elastic thickness).

The deflection :math:`w` satisfies

.. math::

   \nabla^2(D\,\nabla^2 w) - T_e\,\boldsymbol{\sigma} : \nabla\nabla w + \Delta\rho\, g\, w = q,

where :math:`\boldsymbol{\sigma} : \nabla\nabla w` is the double contraction of the
in-plane stress tensor with the Hessian of :math:`w`, :math:`q` [Pa] is the applied
surface normal stress, :math:`\Delta\rho = \rho_m - \rho_\text{fill}` [kg m⁻³] is
the mantle minus infill density, :math:`g` [m s⁻²] is gravitational acceleration,
and :math:`D = E T_e^3 / \bigl[12(1 - \nu^2)\bigr]` is the flexural rigidity
(:math:`E` = Young's modulus, :math:`T_e` = elastic thickness,
:math:`\nu` = Poisson's ratio).  See :doc:`theory` for the full expanded equations
and physical interpretation of each term.

.. note::

   When you use gFlex, please cite:

   Wickert, A. D. (2016), `Open-source modular solutions for flexural
   isostasy: gFlex v1.0 <https://doi.org/10.5194/gmd-9-997-2016>`_,
   *Geosci. Model Dev.*, *9*\(3), 997–1017.

Workflow
--------

.. figure:: _static/flowchart.*
   :alt: gFlex computational workflow
   :align: center
   :width: 90%

   gFlex computational workflow.  Bold-bordered nodes are CSDMS Basic Model
   Interface (BMI) methods.  Boundary conditions apply to finite-difference
   solutions; padding-based approximations are also shown.

Installation
------------

.. code-block:: bash

   pip install gflex

gFlex requires Python ≥ 3.10 and depends on NumPy, SciPy, and Matplotlib.

Quick start
-----------

2-D finite-difference deflection under a rectangular load::

   import numpy as np
   from gflex import F2D

   flex = F2D()
   flex.quiet = True
   flex.method = 'FD'
   flex.g = 9.8
   flex.E = 65e9
   flex.nu = 0.25
   flex.rho_m = 3300.
   flex.rho_fill = 0.
   flex.te = 35e3 * np.ones((50, 50))   # uniform 35 km elastic thickness
   flex.qs = np.zeros((50, 50))
   flex.qs[10:40, 10:40] = 1e6          # 150 × 150 km load at 1 MPa
   flex.dx = flex.dy = 5000.            # 5 km grid
   flex.bc_west = flex.bc_south = flex.bc_north = 'zero_displacement_zero_slope'
   flex.bc_east = 'zero_moment_zero_shear'
   flex.initialize()
   flex.run()
   flex.finalize()

   deflection = flex.w   # (50, 50) array; negative values = downward

Configuration files
-------------------

As an alternative to the programmatic API, gFlex can be driven by a
configuration file — passed to the ``gflex`` CLI or to the
:class:`~gflex.F1D` / :class:`~gflex.F2D` constructor.  The file must
be YAML (``.yaml`` / ``.yml``); see ``input/input_f1d.yaml`` and
``input/input_f2d.yaml`` for complete 1-D and 2-D examples.

A minimal 2-D YAML configuration:

.. code-block:: yaml

   mode:
     dimension: 2
     method: FD
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

Run from the command line::

   gflex path/to/config.yaml

See :doc:`configuration` for a full parameter reference and annotated
examples.

The package also provides domain-padding utilities for variable-*Te* grids
(see :func:`~gflex.pad_domain` for 2-D and :func:`~gflex.pad_domain_1d`
for 1-D), a flexural wavelength calculator
(see :func:`~gflex.flexural_wavelengths`), a Landlab component, a CSDMS
Basic Model Interface (:class:`~gflex.BmiGflex`), and a command-line entry
point (``gflex <config_file>``).

----

.. toctree::
   :maxdepth: 2
   :caption: Contents

   theory
   numerical_methods
   boundary_conditions
   example
   api
   configuration
   accuracy
   references
   changelog

----

Documentation prepared with AI assistance (Claude, Anthropic) drawing on
the gFlex source code and Wickert (2016), reviewed by A. D. Wickert.
