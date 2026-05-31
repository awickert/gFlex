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

    * ``FD`` — Finite Difference.  Supports spatially variable elastic
      thickness.  Requires a grid (``dx``, and ``dy`` in 2-D).
    * ``FFT`` — Spectral (Fast Fourier Transform).  Requires
      scalar (uniform) :math:`T_e`.  Spectrally accurate and fast.  When
      all boundary conditions are ``periodic`` the domain tiles exactly;
      for any other boundary condition the load is zero-padded by
      :math:`4\alpha` on each side, approximating the ``no_outside_loads``
      condition.  In-plane stresses (:math:`\sigma_{xx}`, :math:`\sigma_{yy}`,
      :math:`\sigma_{xy}`) are supported.
    * ``SAS`` — Superposition of Analytical Solutions.  Constant elastic
      thickness only; fast and analytically exact.
    * ``SAS_NG`` — SAS on an unstructured point set (NG = "no grid").
      Load and output locations are arbitrary (x, q0) or (x, y, q0)
      columns; see ``Loads`` below.

----

``parameter`` section
~~~~~~~~~~~~~~~~~~~~~

``YoungsModulus``
    Young's modulus :math:`E` [Pa].  Typical lithospheric value: 65 GPa
    (``6.5e10``).

``PoissonsRatio``
    Poisson's ratio :math:`\nu` [dimensionless].  Typical value: 0.25.

``GravAccel``
    Gravitational acceleration :math:`g` [m s⁻²].  Earth standard: 9.8.

``MantleDensity``
    Density of the mantle :math:`\rho_m` [kg m⁻³].  Typical value: 3300.

``InfillMaterialDensity``
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

``Loads``
    Path to the load file.

    * *Gridded methods (FD, SAS)*: a space-delimited array of surface
      normal stresses [Pa] (:math:`\rho g h`).  Grid cell area
      (:math:`\Delta x \times \Delta y`) is applied internally to convert
      stress to force.
    * *SAS_NG*: a space-delimited file with columns ``(x, q0)`` in 1-D or
      ``(x, y, q0)`` in 2-D, where q0 is a point force [N].

    Paths are resolved relative to the directory containing the
    configuration file.

``ElasticThickness``
    Elastic thickness [m].  Either a scalar value or a path to a
    space-delimited array.  Arrays are required for FD solutions with
    spatially variable *Te*.  Use :func:`~gflex.smooth_pad_Te` and
    :func:`~gflex.pad_domain` (2-D) or :func:`~gflex.smooth_pad_Te_1d` and
    :func:`~gflex.pad_domain_1d` (1-D) to extend a variable-*Te* grid with a
    smooth boundary buffer before running.

``xw``, ``yw``
    *(SAS_NG only)*  Vectors of x (and y for 2-D) coordinates at which to
    evaluate deflection.  If omitted, deflection is evaluated at the load
    points.

----

``output`` section
~~~~~~~~~~~~~~~~~~

``DeflectionOut``
    Path for writing deflection output as a space-delimited ASCII file.
    Leave blank to suppress file output.

``Plot``
    Controls inline plotting after the run:

    * ``q`` — plot the applied load.
    * ``w`` — plot the deflection.
    * ``both`` — deflection and load in separate subplots.
    * ``combo`` — *(1-D only)* deflection with the load overlaid.

    Any other value (or blank) suppresses plotting.

----

``numerical`` section
~~~~~~~~~~~~~~~~~~~~~

``GridSpacing_x``
    Grid cell size in the x-direction [m].

``BoundaryCondition_West``, ``BoundaryCondition_East``
    Boundary conditions on the west and east edges.

    For FD: ``zero_displacement_zero_slope``, ``zero_displacement_zero_moment``,
    ``zero_moment_zero_shear``, ``zero_slope_zero_shear``, ``mirror``, or ``periodic``.
    See :doc:`boundary_conditions` for the physical meaning of each and
    guidance on choosing.

    For SAS / SAS_NG: ``no_outside_loads`` (assumed if left blank).

.. note::

   **In-plane stresses** (:math:`\sigma_{xx}`, :math:`\sigma_{yy}`,
   :math:`\sigma_{xy}` [Pa]) cannot be set from a configuration file.
   They must be assigned programmatically before calling
   :meth:`~gflex.F2D.initialize`::

       flex.sigma_xx = 1e6   # east–west compression [Pa]
       flex.sigma_yy = 0.
       flex.sigma_xy = 0.

   All three default to zero if not set.  They are supported by ``FD``
   and ``FFT`` in 2-D, and by ``FD`` and ``FFT`` in 1-D (``sigma_xx``
   only).  Setting them with ``SAS`` or ``SAS_NG`` raises a warning and
   has no effect.  See :doc:`theory` for the governing equations.

----

``numerical2D`` section
~~~~~~~~~~~~~~~~~~~~~~~

``GridSpacing_y``
    Grid cell size in the y-direction [m].

``BoundaryCondition_North``, ``BoundaryCondition_South``
    Same options as ``BoundaryCondition_West`` / ``BoundaryCondition_East``.

``latlon``
    ``true`` / ``false``.  Interpret input coordinates as geographic
    latitude and longitude.  Default: ``false``.

``PlanetaryRadius``
    Planetary radius [m].  Required when ``latlon = true``.  Earth:
    6 371 000 m.

----

``verbosity`` section
~~~~~~~~~~~~~~~~~~~~~

``Verbose``
    ``true`` / ``false``.  Print progress messages during the run.
    Default: ``true``.

``Debug``
    ``true`` / ``false``.  Print internal arrays and solver diagnostics.
    Default: ``false``.

``Quiet``
    ``true`` / ``false``.  Suppress all output.  Overrides ``Verbose`` and
    ``Debug``.  Default: ``false``.

----

Complete examples
-----------------

1-D finite-difference example (YAML)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

    # All units are SI.

    mode:
      dimension: 1
      method: FD

    parameter:
      YoungsModulus: 6.5e10
      PoissonsRatio: 0.25
      GravAccel: 9.8
      MantleDensity: 3300
      InfillMaterialDensity: 0

    input:
      Loads: q0_sample/1D/central_block.txt
      ElasticThickness: Te_sample/1D/8km_20km_ramp.txt

    output:
      DeflectionOut: ""
      Plot: combo          # overlay deflection and load (1-D only)

    numerical:
      GridSpacing_x: 6000
      BoundaryCondition_West: periodic
      BoundaryCondition_East: periodic

    verbosity:
      Verbose: false
      Debug: false
      Quiet: false

2-D finite-difference example (YAML)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

    # All units are SI.

    mode:
      dimension: 2
      method: FD

    parameter:
      YoungsModulus: 6.5e10
      PoissonsRatio: 0.25
      GravAccel: 9.8
      MantleDensity: 3300
      InfillMaterialDensity: 0

    input:
      Loads: q0_sample/2D/diag.txt
      ElasticThickness: Te_sample/2D/fault_24-30.txt

    output:
      DeflectionOut: ""
      Plot: both

    numerical:
      GridSpacing_x: 4000
      BoundaryCondition_West: zero_moment_zero_shear
      BoundaryCondition_East: zero_displacement_zero_slope

    numerical2D:
      GridSpacing_y: 4000
      BoundaryCondition_North: mirror
      BoundaryCondition_South: zero_slope_zero_shear

    verbosity:
      Verbose: false
      Debug: false
      Quiet: false

