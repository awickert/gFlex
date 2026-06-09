Configuration Files
===================

gFlex can be driven by a YAML configuration file rather than (or alongside)
the programmatic Python API.  Pass the file path to the ``gflex`` CLI or
to an :class:`~gflex.F1D` / :class:`~gflex.F2D` constructor.  The file
must have a ``.yaml`` or ``.yml`` extension.

----

Parameters
----------

``mode`` section
~~~~~~~~~~~~~~~~

``dimension``
    ``1`` or ``2``.  Selects a 1-D (profile) or 2-D (map-view) flexural
    solution.

``method``
    Solution method:

    * ``fd`` — Finite Difference.  Supports spatially variable elastic
      thickness.  Requires a grid (``dx``, and ``dy`` in 2-D).
    * ``fft`` — Spectral (Fast Fourier Transform).  Requires
      scalar (uniform) :math:`T_e`.  Spectrally accurate and fast.  When
      all boundary conditions are ``periodic`` the domain tiles exactly;
      for any other boundary condition the load is zero-padded by
      :math:`4\alpha` on each side, approximating the ``no_outside_loads``
      condition (the default when BCs are unset).  Setting only *some* BCs to
      ``periodic`` triggers a ``UserWarning`` and falls back to zero-padding.
      In-plane stresses (:math:`\sigma_{xx}`, :math:`\sigma_{yy}`,
      :math:`\sigma_{xy}`) are supported.
    * ``sas`` — Superposition of Analytical Solutions.  Constant elastic
      thickness only; fast and analytically exact.
    * ``sas_ng`` — SAS on an unstructured point set (NG = "no grid").
      Load and output locations are arbitrary (x, q0) or (x, y, q0)
      columns; see ``loads`` below.

----

``parameter`` section
~~~~~~~~~~~~~~~~~~~~~

``youngs_modulus``
    Young's modulus :math:`E` [Pa].  Typical lithospheric value: 65 GPa
    (``6.5e10``).

``poissons_ratio``
    Poisson's ratio :math:`\nu` [dimensionless].  Typical value: 0.25.

``gravitational_acceleration``
    Gravitational acceleration :math:`g` [m s⁻²].  Earth standard: 9.8.

``mantle_density``
    Density of the mantle :math:`\rho_m` [kg m⁻³].  Typical value: 3300.

``infill_material_density``
    Density of the material that fills (or vacates) the flexural depression
    :math:`\rho_\text{fill}` [kg m⁻³].  Common values:

    * ``0`` — air (no infill)
    * ``1030`` — seawater
    * ``2000``–``2700`` — sediment

    If the infill density varies spatially (e.g., at a subsiding shoreline
    that progressively floods), iterate externally: flex, update the
    inundation mask, re-flex, repeat.

----

``input`` section
~~~~~~~~~~~~~~~~~

``loads``
    Path to the load file.

    * *Gridded methods (FD, SAS)*: a space-delimited array of surface
      normal stresses [Pa] (:math:`\rho g h`).  Grid cell area
      (:math:`\Delta x \times \Delta y`) is applied internally to convert
      stress to force.
    * *SAS_NG*: a space-delimited file with columns ``(x, q0)`` in 1-D or
      ``(x, y, q0)`` in 2-D, where q0 is a point force [N].

    Paths are resolved relative to the directory containing the
    configuration file.

``elastic_thickness``
    Elastic thickness [m].  Either a scalar value or a path to a
    space-delimited array.  Arrays are required for FD solutions with
    spatially variable *Te*.  Use :func:`~gflex.pad_domain` to extend a
    variable-*Te* grid with a smooth boundary buffer before running (works
    for both 1-D and 2-D; scalar *Te* is also supported).

``xw``, ``yw``
    *(SAS_NG only)*  Vectors of x (and y for 2-D) coordinates at which to
    evaluate deflection.  If omitted, deflection is evaluated at the load
    points.

----

``output`` section
~~~~~~~~~~~~~~~~~~

``deflection_out``
    Path for writing deflection output as a space-delimited ASCII file.
    Leave blank to suppress file output.

``plot``
    Controls inline plotting after the run:

    * ``q`` — plot the applied load.
    * ``w`` — plot the deflection.
    * ``both`` — deflection and load in separate subplots.
    * ``combo`` — *(1-D only)* deflection with the load overlaid.

    Any other value (or blank) suppresses plotting.

----

``numerical`` section
~~~~~~~~~~~~~~~~~~~~~

``grid_spacing_x``
    Grid cell size in the x-direction [m].

``boundary_condition_west``, ``boundary_condition_east``
    Boundary conditions on the west and east edges.

    For FD: ``zero_displacement_zero_slope`` (alias ``clamped``),
    ``zero_displacement_zero_moment`` (alias ``pinned``),
    ``zero_moment_zero_shear`` (alias ``free``),
    ``zero_slope_zero_shear`` (alias ``mirror``),
    ``no_outside_loads`` (alias ``infinite``), or ``periodic``.
    When ``no_outside_loads`` is set on an FD edge, the solver
    automatically pads the domain by one flexural wavelength, applies a
    clamped outer boundary, solves, and crops ``w`` back to the original
    domain — the padding is invisible to the caller.
    See :doc:`boundary_conditions` for the physical meaning of each and
    guidance on choosing.

    For SAS / SAS_NG: ``no_outside_loads`` is always assumed (the
    analytical solution is inherently infinite-plate); leave blank or set
    explicitly.

.. note::

   **In-plane stresses** (:math:`\sigma_{xx}`, :math:`\sigma_{yy}`,
   :math:`\sigma_{xy}` [Pa]) cannot be set from a configuration file.
   They must be assigned programmatically before calling
   :meth:`~gflex.F2D.initialize`::

       flex.sigma_xx = 1e6   # east–west compression [Pa]
       flex.sigma_yy = 0.
       flex.sigma_xy = 0.

   All three default to zero if not set.  They are supported by ``fd``
   and ``fft`` in 2-D, and by ``fd`` and ``fft`` in 1-D (``sigma_xx``
   only).  Setting them with ``sas`` or ``sas_ng`` raises a warning and
   has no effect.  See :doc:`theory` for the governing equations.

----

``numerical2D`` section
~~~~~~~~~~~~~~~~~~~~~~~

``grid_spacing_y``
    Grid cell size in the y-direction [m].

``boundary_condition_north``, ``boundary_condition_south``
    Same options as ``boundary_condition_west`` / ``boundary_condition_east``.

``latlon``
    ``true`` / ``false``.  Interpret input coordinates as geographic
    latitude and longitude.  Default: ``false``.

``planetary_radius``
    Planetary radius [m].  Required when ``latlon = true``.  Earth:
    6 371 000 m.

----

``verbosity`` section
~~~~~~~~~~~~~~~~~~~~~

``verbose``
    ``true`` / ``false``.  Print progress messages during the run.
    Default: ``true``.

``debug``
    ``true`` / ``false``.  Print internal arrays and solver diagnostics.
    Default: ``false``.

``quiet``
    ``true`` / ``false``.  Suppress all output.  Overrides ``verbose`` and
    ``debug``.  Default: ``false``.

----

Complete examples
-----------------

1-D finite-difference example (YAML)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

    # All units are SI.

    mode:
      dimension: 1
      method: fd

    parameter:
      youngs_modulus: 6.5e10
      poissons_ratio: 0.25
      gravitational_acceleration: 9.8
      mantle_density: 3300
      infill_material_density: 0

    input:
      loads: q0_sample/1D/central_block.txt
      elastic_thickness: Te_sample/1D/8km_20km_ramp.txt

    output:
      deflection_out: ""
      plot: combo          # overlay deflection and load (1-D only)

    numerical:
      grid_spacing_x: 6000
      boundary_condition_west: periodic
      boundary_condition_east: periodic

    verbosity:
      verbose: false
      debug: false
      quiet: false

2-D finite-difference example (YAML)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

    # All units are SI.

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
      loads: q0_sample/2D/diag.txt
      elastic_thickness: Te_sample/2D/fault_24-30.txt

    output:
      deflection_out: ""
      plot: both

    numerical:
      grid_spacing_x: 4000
      boundary_condition_west: zero_moment_zero_shear
      boundary_condition_east: zero_displacement_zero_slope

    numerical2D:
      grid_spacing_y: 4000
      boundary_condition_north: mirror
      boundary_condition_south: zero_slope_zero_shear

    verbosity:
      verbose: false
      debug: false
      quiet: false

