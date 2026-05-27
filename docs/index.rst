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

.. note::

   When you use gFlex, please cite:

   Wickert, A. D. (2016), `Open-source modular solutions for flexural
   isostasy: gFlex v1.0 <https://doi.org/10.5194/gmd-9-997-2016>`_,
   *Geosci. Model Dev.*, *9*\(3), 997–1017.

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
   flex.Quiet = True
   flex.Method = 'FD'
   flex.PlateSolutionType = 'vWC1994'
   flex.Solver = 'direct'
   flex.g = 9.8
   flex.E = 65e9
   flex.nu = 0.25
   flex.rho_m = 3300.
   flex.rho_fill = 0.
   flex.Te = 35e3 * np.ones((50, 50))   # uniform 35 km elastic thickness
   flex.qs = np.zeros((50, 50))
   flex.qs[10:40, 10:40] = 1e6          # 150 × 150 km load at 1 MPa
   flex.dx = flex.dy = 5000.            # 5 km grid
   flex.BC_W = flex.BC_S = flex.BC_N = '0Displacement0Slope'
   flex.BC_E = '0Moment0Shear'
   flex.initialize()
   flex.run()
   flex.finalize()

   deflection = flex.w   # (50, 50) array; negative values = downward

The package also provides domain-padding utilities for variable-*Te* grids
(see :func:`~gflex.pad_domain`), a flexural wavelength calculator
(see :func:`~gflex.flexural_wavelengths`), a Landlab component, a CSDMS
Basic Model Interface (:class:`~gflex.BmiGflex`), and a command-line entry
point (``gflex <config_file>``).

----

.. toctree::
   :maxdepth: 2
   :caption: Contents

   api
   changelog
