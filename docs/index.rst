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
:math:`\nu` = Poisson's ratio).  See :doc:`theory_and_numerics` for the full expanded equations
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

gFlex requires **Python ≥ 3.11** and depends on NumPy, SciPy, and Matplotlib.

.. note::

   This documentation describes gFlex **2.0.0**.  If ``pip install gflex``
   installs an older release (e.g. 2.0.0b1), some features documented here
   may not be available.  Install the latest release explicitly::

      pip install "gflex>=2.0.0"

   Downstream tools that embed gFlex (QGIS Processing provider, GRASS GIS
   addons) should pin ``gflex>=2.0.0`` in their requirements.

Quick start
-----------

2-D finite-difference deflection under a rectangular load::

   import numpy as np
   from gflex import F2D

   flex = F2D()
   flex.quiet = True
   flex.method = 'fd'
   flex.g = 9.8
   flex.E = 65e9
   flex.nu = 0.25
   flex.rho_m = 3300.
   flex.rho_fill = 0.
   flex.T_e = 35e3 * np.ones((50, 50))   # uniform 35 km elastic thickness
   flex.qs = np.zeros((50, 50))
   flex.qs[10:40, 10:40] = 1e6          # 150 × 150 km load at 1 MPa
   flex.dx = flex.dy = 5000.            # 5 km grid
   flex.bc_west = flex.bc_south = flex.bc_north = 'zero_displacement_zero_slope'
   flex.bc_east = 'zero_moment_zero_shear'
   flex.initialize()
   flex.run()

   deflection = flex.w   # (50, 50) array; negative values = downward
   flex.finalize()        # releases w, qs, and the coefficient matrix

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
     boundary_condition_north: zero_slope_zero_shear
     boundary_condition_south: zero_slope_zero_shear

Run from the command line::

   gflex path/to/config.yaml

See :doc:`configuration` for a full parameter reference and annotated
examples.

Interfaces
----------

gFlex can be accessed through several front-ends depending on your workflow:

**Python and command-line**

.. list-table::
   :header-rows: 1
   :widths: 24 12 64

   * - Interface
     - :math:`T_e`
     - Description
   * - **Python API** (:class:`~gflex.F1D` / :class:`~gflex.F2D`)
     - scalar or array
     - Full programmatic control.  The primary interface; all other
       front-ends call into it.  See the :doc:`tutorial` and :doc:`api`.
   * - **CLI** (``gflex <config.yaml>``)
     - scalar or array
     - Drive gFlex from a YAML configuration file with no Python code.
       See :doc:`configuration`.

**Modelling frameworks**

.. list-table::
   :header-rows: 1
   :widths: 24 12 64

   * - Interface
     - :math:`T_e`
     - Description
   * - **CSDMS BMI** (:class:`~gflex.BmiGflex`)
     - scalar or array
     - CSDMS Basic Model Interface for coupling in the CSDMS framework.
       Requires ``pip install gflex[bmi]``.  See :doc:`api`.
   * - **Landlab component** (``landlab.components.gFlex``)
     - scalar or array
     - Landlab Earth-surface modelling framework component.  Uses
       ``grid.at_node`` fields; compatible with the CSDMS Standard Names
       used by the BMI.  Install with ``pip install landlab``.

**GIS**

.. list-table::
   :header-rows: 1
   :widths: 24 12 64

   * - Interface
     - :math:`T_e`
     - Description
   * - **GRASS GIS** (``r.flexure``, ``v.flexure``)
     - scalar or array
     - Raster and vector interfaces for use inside a GRASS GIS session.
       ``r.flexure`` uses FD, FFT, or SAS with optional variable
       :math:`T_e`; ``v.flexure`` uses SAS_NG for scattered point loads.
       Install with ``g.extension``.
       `r.flexure <https://github.com/awickert/r.flexure>`_ —
       `v.flexure <https://github.com/awickert/v.flexure>`_.
   * - **QGIS Processing provider** (``processing_gflex``)
     - scalar or array
     - No-code access from the QGIS Processing Toolbox and Graphical
       Modeler.  Supports all 2-D methods, variable :math:`T_e`, all
       boundary conditions, and in-plane stresses.
       `processing_gflex <https://github.com/awickert/processing_gflex>`_.

----

.. toctree::
   :maxdepth: 2
   :caption: Contents

   theory_and_numerics
   boundary_conditions
   tutorial
   greenland_example
   api
   configuration
   accuracy
   benchmarks
   references
   changelog

----

Documentation prepared with AI assistance (Claude, Anthropic) drawing on
the gFlex source code and Wickert (2016), reviewed by A. D. Wickert.
