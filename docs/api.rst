API Reference
=============

Primary solvers
---------------

Both solvers follow the same lifecycle::

   flex.initialize()
   flex.run()
   w = flex.w          # read deflection before finalize clears it
   flex.output()       # optional: save to file or display plots
   flex.finalize()     # releases w, qs, and the coefficient matrix

.. warning::

   :meth:`finalize` deletes ``w``, ``qs``, and the cached coefficient
   matrix.  Read ``w`` (and call :meth:`~gflex.base.Flexure.output` if
   needed) **before** calling ``finalize``.  Accessing ``w`` afterwards
   raises :exc:`AttributeError`.

.. autoclass:: gflex.F2D
   :members: initialize, run, finalize

.. autoclass:: gflex.F1D
   :members: initialize, run, finalize

Boundary-condition properties
-----------------------------

.. versionadded:: 2.0.0

The boundary conditions are set as properties on the solver instance before
calling :meth:`initialize`.  Each accepts a canonical BC string, a short
alias, or a ``dict`` for inhomogeneous (prescribed-value) conditions:

.. code-block:: python

   flex.bc_west  = "zero_displacement_zero_slope"   # or "clamped"
   flex.bc_east  = "zero_moment_zero_shear"          # or "free"
   flex.bc_south = "zero_slope_zero_shear"           # or "mirror"
   flex.bc_north = {"displacement": w_arr, "slope": dw_arr}

All valid strings are listed in :data:`gflex.VALID_BC_STRINGS_1D` and
:data:`gflex.VALID_BC_STRINGS_2D`.

.. data:: gflex.VALID_BC_STRINGS_1D

   .. versionadded:: 2.0.0

   :class:`frozenset` of every accepted BC string for :class:`~gflex.F1D`
   (canonical names and aliases).  Use this to validate user input without
   maintaining a parallel copy that may drift with new releases.

.. data:: gflex.VALID_BC_STRINGS_2D

   .. versionadded:: 2.0.0

   :class:`frozenset` of every accepted BC string for :class:`~gflex.F2D`.

Output
------

.. automethod:: gflex.base.Flexure.output

In-plane stresses
-----------------

.. versionadded:: 1.4.0

In-plane stresses are set as attributes directly on the solver instance
before calling :meth:`initialize`.  They are not available as configuration
file keys.

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Attribute
     - Solvers
     - Description
   * - ``sigma_xx``
     - FD, FFT (1-D and 2-D)
     - Normal stress in the x-direction :math:`\sigma_{xx}` [Pa].
       Default ``0``.
   * - ``sigma_yy``
     - FD, FFT (2-D only)
     - Normal stress in the y-direction :math:`\sigma_{yy}` [Pa].
       Default ``0``.
   * - ``sigma_xy``
     - FD, FFT (2-D only)
     - Shear stress :math:`\sigma_{xy}` [Pa].  Default ``0``.

All three default to zero if not assigned; setting any of them with
``SAS`` or ``SAS_NG`` raises a :exc:`RuntimeWarning` and has no effect.
See :doc:`theory_and_numerics` for the governing equations that include these terms.

Domain-padding utilities
------------------------

These functions help when running :class:`~gflex.F1D` or :class:`~gflex.F2D`
with a spatially variable elastic thickness grid.  A smooth padding zone
reduces spurious deflections at the domain boundary caused by sharp rigidity
gradients, and ensures that the flexural forebulge can develop freely before
reaching the boundary.

**All-in-one helper**

:func:`~gflex.pad_domain` handles both 1-D and 2-D grids and both scalar
and array elastic thickness.  It calls the lower-level helpers below and
is the recommended starting point.

.. versionadded:: 1.4.0
   1-D support (dispatching on ``qs.ndim``) added in 2.0.0.

.. autofunction:: gflex.pad_domain

**Lower-level building blocks**

These two pairs of functions are the building blocks used by
:func:`~gflex.pad_domain`.  Use them directly when you need finer control —
for example, to compute the pad width once and apply it to multiple arrays,
or to inspect the tapered Te grid before running the solver.

*Recommended pad width* (number of cells):

.. versionadded:: 1.4.0

.. autofunction:: gflex.recommended_pad_width

.. autofunction:: gflex.recommended_pad_width_1d

*Smooth Te taper* (extends a Te array into the padding zone):

.. versionadded:: 1.4.0

.. autofunction:: gflex.smooth_pad_Te

.. autofunction:: gflex.smooth_pad_Te_1d

FD boundary-condition warnings
------------------------------

When running :class:`~gflex.F1D` or :class:`~gflex.F2D` with the finite-difference
solver, gFlex issues :exc:`UserWarning` messages for two categories of potentially
problematic boundary conditions.

**BC-type warnings** fire whenever a side carries a BC whose physical
interpretation deserves verification:

* ``'zero_moment_zero_shear'`` (alias ``'free'``) — assumes a free broken plate end (zero moment and shear
  force).  Physically appropriate for rifted or passive continental margins,
  subduction trenches with an applied edge load (slab pull), and broken-plate
  flexure (Turcotte & Schubert).  Often applied uncritically elsewhere in the
  literature — verify that one of these settings applies.

**Proximity warnings** fire for ``'zero_displacement_zero_slope'`` (alias ``'clamped'``) boundaries when the
nearest loaded cell is within one flexural wavelength
(:math:`\lambda = 2\pi\alpha`, where
:math:`\alpha = (4D / \Delta\rho g)^{1/4}`)
of the boundary.  Within this distance the flexural forebulge — which peaks at
:math:`\approx \pi\alpha` from the load — will be suppressed by the zero-displacement
condition, contaminating the solution.  The warning message reports the distance as
a fraction of the local flexural wavelength and directs you to the domain-padding
utilities.

**Warning deduplication in model-coupling loops**

Python's default warning filter shows each unique warning *once per call site*
per interpreter session.  In a time-stepping or iterative-coupling loop such as::

   for load in loads:
       flex.qs = load
       flex.run()          # warning fires on first iteration, silenced thereafter

the proximity warning fires on the first iteration and is not repeated, even if
the load subsequently moves closer to the boundary.  To re-enable the warning on
every call::

   import warnings
   warnings.filterwarnings("always", category=UserWarning, module="gflex")

**Suppressing warnings you have verified**

Once you have confirmed that a boundary condition is appropriate for your setup,
suppress the corresponding warning by message text::

   import warnings
   warnings.filterwarnings("ignore", message=".*zero_moment_zero_shear.*")

To suppress all gFlex warnings::

   warnings.filterwarnings("ignore", module="gflex")

Or use a context manager for a single run::

   with warnings.catch_warnings():
       warnings.simplefilter("ignore")
       flex.run()

LU factorization cache
----------------------

.. versionadded:: 2.0.0

For coupling workflows that call :meth:`~gflex.F1D.run` (or
:meth:`~gflex.F2D.run`) repeatedly with the same grid, elastic thickness,
and boundary conditions, the sparse-LU factorization of the coefficient
matrix can be cached to avoid re-factorizing on every call.  Set the
attribute **before** calling :meth:`~gflex.base.Flexure.initialize`:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Value
     - Behaviour
   * - ``False`` *(default)*
     - No caching.  The matrix is factorized on every :meth:`run` call.
   * - ``True``
     - Cache with hash check.  The factorization is reused when a
       hash of the coefficient matrix matches the stored hash; it is
       recomputed when any of Te, dx/dy, BCs, or physical parameters change.
   * - ``"no_check"``
     - Cache without hash check.  The stored factorization is reused on
       every call without computing a hash.  Gives the maximum performance
       benefit.  The coefficient matrix is freed from memory immediately
       after factorization; only the LU factorization is retained.  Smart
       invalidation (see below) still applies: reassigning a
       matrix-determining input clears the cache and triggers a rebuild on
       the next call.

**Smart invalidation**

Reassigning any matrix-determining attribute — ``T_e``, ``E``, ``nu``,
``g``, ``rho_m``, ``rho_fill``, ``dx``, ``dy``, boundary conditions, or
in-plane stresses — automatically clears the cached coefficient matrix and
LU factorization.  No explicit cache management is needed between solves
when only ``qs`` changes.

.. note::

   Smart invalidation is triggered by *assignment* (``flex.T_e = new_array``),
   not by in-place mutation of a NumPy array (``flex.T_e[5] = 40e3``).  If
   you mutate an array in place, reassign it afterwards to ensure the cache
   is correctly invalidated::

      flex.T_e[5] = 40e3
      flex.T_e = flex.T_e   # trigger invalidation

Example (coupling loop)::

   flex.cache_factorization = True
   flex.initialize()

   for load in load_sequence:
       flex.qs = load
       flex.run()
       w = flex.w
       # ... process w ...

   flex.finalize()

The cache is cleared by :meth:`~gflex.base.Flexure.finalize`.

Timing attributes
-----------------

.. versionadded:: 2.0.0

After each :meth:`~gflex.base.Flexure.run` call, the following read-only
attributes report wall-clock times measured with :func:`time.perf_counter`:

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Attribute
     - Solvers
     - Description
   * - ``time_to_solve``
     - All
     - Total solve wall time [s] from the start of :meth:`run` to the end.
   * - ``coeff_creation_time``
     - FD only
     - Time [s] to construct the sparse coefficient matrix.  Not set for
       SAS or FFT solvers.  In a coupling loop with
       ``cache_factorization`` enabled, this is zero on cache-hit calls.
   * - ``linear_solve_time``
     - FD only
     - Time [s] for the LU backsolve (triangular solve).  Not set for
       SAS or FFT solvers.

``coeff_creation_time`` and ``linear_solve_time`` together account for most
of ``time_to_solve``; the remainder is boundary-condition setup and array
housekeeping.  Both are most useful in coupling loops:
``coeff_creation_time`` confirms the coefficient matrix was reused (value ≈ 0)
and ``linear_solve_time`` shows the marginal cost per :meth:`run` call.

Flexural wavelengths
--------------------

.. versionadded:: 1.4.0

.. autofunction:: gflex.flexural_wavelengths

Coupling guide
--------------

gFlex is material-agnostic: it receives a surface-normal stress [Pa] and
returns a deflection [m], regardless of whether the load comes from ice,
water, sediment, volcanic edifices, or any combination.  The caller is
responsible for converting source-specific quantities into Pa before passing
them to gFlex.

**Load conversion**

.. code-block:: python

   # Glacial isostasy
   qs = rho_ice * g * ice_thickness        # kg m⁻³ × m s⁻² × m → Pa

   # Sediment or volcanic load
   qs = rho_sediment * g * sediment_thickness

   # Multiple sources: sum them
   qs = rho_ice * g * h_ice + rho_sed * g * h_sed + rho_water * g * h_water

The ``rho_fill`` parameter should be set to the density of the material that
*replaces* the load inside the flexural depression — ``rho_fill=0`` for
subaerially exposed basins, ``rho_fill=1030`` for submarine basins,
``rho_fill=2000``–``2700`` for sediment-filled basins.

**Applying deflection to topography**

gFlex returns the *total* instantaneous deflection, not an increment.
In a time-stepping loop, apply only the *change* to topography:

.. code-block:: python

   flex.qs = qs
   flex.run()
   w_new = flex.w
   topo += w_new - w_prev
   w_prev = w_new.copy()

**QGIS Processing provider**

For GIS-based workflows where load and elastic-thickness data are already
rasters, `processing_gflex <https://github.com/awickert/processing_gflex>`_
exposes gFlex as a no-code algorithm in the QGIS Processing Toolbox and
Graphical Modeler.  It supports all 2-D solution methods (FD, FFT, SAS),
variable or scalar :math:`T_e`, all boundary conditions, and in-plane
stresses.  Installation::

   pip install "gflex>=2.0.0"   # auto-installed by the plugin on first use

Requires QGIS ≥ 3.16.  Usable from the Toolbox, Graphical Modeler, or
headlessly via ``processing.run()``.

BMI interface
-------------

:class:`~gflex.BmiGflex` exposes the CSDMS Basic Model Interface, enabling
gFlex to be coupled with other models in the CSDMS framework.  It requires
the optional ``bmipy`` dependency (``pip install gflex[bmi]``).

**BMI variables**

Grid 0 — spatial flexure grid:

.. list-table::
   :header-rows: 1
   :widths: 40 10 10 40

   * - Name
     - Direction
     - Units
     - Description
   * - ``load__normal_component_of_stress``
     - input
     - Pa
     - Surface-normal load stress :math:`q_s`.  Material-agnostic: convert
       ice, water, sediment, etc. to Pa before calling :meth:`set_value`.
   * - ``lithosphere__elastic_thickness``
     - input
     - m
     - Elastic thickness :math:`T_e`.  Usually set once at initialisation;
       updating it between :meth:`update` calls invalidates the cached LU
       factorisation so the next solve rebuilds the stiffness matrix.
   * - ``lithosphere__vertical_displacement``
     - output
     - m
     - Deflection :math:`w` (downward negative).

Grid 1 — scalar physical constants (single-element arrays):

.. list-table::
   :header-rows: 1
   :widths: 40 10 10 40

   * - Name
     - Direction
     - Units
     - Description
   * - ``lithosphere__young_modulus``
     - input
     - Pa
     - Young's modulus :math:`E`.
   * - ``lithosphere__poisson_ratio``
     - input
     - 1
     - Poisson's ratio :math:`\nu` (dimensionless).
   * - ``mantle__mass-per-volume_density``
     - input
     - kg m⁻³
     - Mantle density :math:`\rho_m`.
   * - ``infill_material__mass-per-volume_density``
     - input
     - kg m⁻³
     - Infill material density :math:`\rho_\text{fill}`.  The only constant
       with a runtime-update use case: a basin transitioning from subaerial
       (:math:`\rho_\text{fill}=0`) to subaqueous
       (:math:`\rho_\text{fill}=1030`) during a simulation.
   * - ``planet_surface__gravitational_acceleration``
     - input
     - m s⁻²
     - Gravitational acceleration :math:`g`.

Updating any scalar constant via :meth:`set_value` propagates to the solver
immediately and invalidates the cached LU factorisation, so the next
:meth:`update` rebuilds the stiffness matrix automatically.

**Coupling example**

.. code-block:: python

   from gflex import BmiGflex
   import numpy as np

   bmi = BmiGflex()
   bmi.initialize("my_config.yaml")

   n = bmi.get_grid_size(0)
   load = np.zeros(n)

   for step in range(n_steps):
       load[:] = rho_ice * g * ice_thickness.ravel()
       bmi.set_value("load__normal_component_of_stress", load)
       bmi.update()
       w = np.empty(n)
       bmi.get_value("lithosphere__vertical_displacement", w)
       # … apply w to topography …

   bmi.finalize()

.. autoclass:: gflex.BmiGflex
   :members: initialize, update, finalize, get_value, set_value

Landlab component
-----------------

The Landlab component (``landlab.components.gFlex``) exposes gFlex within
the `Landlab <https://landlab.readthedocs.io>`_ Earth-surface modelling
framework.  It uses the same underlying solver and the same CSDMS Standard
Name fields as the BMI, but follows Landlab conventions: construction via
``__init__`` and time-stepping via ``run_one_step()``.

**Comparison with the BMI**

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * -
     - gFlex BMI
     - Landlab component
   * - Lifecycle
     - ``initialize`` / ``update`` / ``finalize``
     - ``__init__`` / ``run_one_step``
   * - Load input
     - ``set_value("load__normal_component_of_stress", q)``
     - ``grid.at_node["load__normal_component_of_stress"][:] = q``
   * - :math:`T_e` update
     - ``set_value("lithosphere__elastic_thickness", te)``
     - ``grid.at_node["lithosphere__elastic_thickness"][:] = te``
   * - Deflection output
     - ``get_value("lithosphere__vertical_displacement", w)``
     - ``grid.at_node["lithosphere__vertical_displacement"]``

**Field names**

.. list-table::
   :header-rows: 1
   :widths: 40 10 10 40

   * - CSDMS Standard Name
     - Direction
     - Units
     - Notes
   * - ``load__normal_component_of_stress``
     - input
     - Pa
     - Required.  Maps to gFlex internal :math:`q_s`.
   * - ``lithosphere__elastic_thickness``
     - input
     - m
     - Optional.  Re-read on every ``run_one_step()`` call if present,
       enabling runtime :math:`T_e` updates without re-initialisation.
   * - ``lithosphere__vertical_displacement``
     - output
     - m
     - Total deflection :math:`w` (downward negative).

``rho_fill`` defaults to ``0.0`` (air; no infill) in the Landlab component.
Set it explicitly at construction for marine (``rho_fill=1030``) or
sediment-filled (``rho_fill=2000``–``2700``) basins.

**Coupling example**

.. code-block:: python

   import numpy as np
   from landlab import RasterModelGrid
   from landlab.components import gFlex

   mg = RasterModelGrid((100, 100), xy_spacing=5000.0)
   mg.add_zeros("load__normal_component_of_stress", at="node")
   mg.add_zeros("topographic__elevation", at="node")

   gf = gFlex(mg, Youngs_modulus=65e9, Poissons_ratio=0.25,
              rho_mantle=3300, rho_fill=0, elastic_thickness=35e3)

   w_prev = np.zeros(mg.number_of_nodes)

   for step in range(n_steps):
       mg.at_node["load__normal_component_of_stress"][:] = (
           rho_ice * g * ice_thickness.ravel()
       )
       gf.run_one_step()
       w = mg.at_node["lithosphere__vertical_displacement"]
       mg.at_node["topographic__elevation"] += w - w_prev
       w_prev = w.copy()

**Installation**

.. code-block:: bash

   pip install landlab

Requires a Landlab release that includes the gFlex v2 component
(``landlab/landlab#2420``).
