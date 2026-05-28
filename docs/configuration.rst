Configuration Files
===================

gFlex can be driven by a configuration file rather than (or alongside)
the programmatic Python API.  Pass the file path to the ``gflex`` CLI or
to an :class:`~gflex.F1D` / :class:`~gflex.F2D` constructor.

Two formats are supported.  The file extension determines which parser is
used:

* **YAML** (``.yaml`` or ``.yml``) — recommended for new configurations.
* **INI** (any other extension) — the original legacy format.

Both formats use the same section names and parameter keys.

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
    * ``SAS`` — Superposition of Analytical Solutions.  Constant elastic
      thickness only; fast and highly accurate.
    * ``SAS_NG`` — SAS on an unstructured point set (NG = "no grid").
      Load and output locations are arbitrary (x, q0) or (x, y, q0)
      columns; see ``Loads`` below.

``PlateSolutionType``
    *(2-D only)*  Plate bending equation variant:

    * ``vWC1994`` — van Wees & Cloetingh (1994); recommended.
    * ``G2009`` — Govers et al. (2009); less robust near boundaries.

----

``parameter`` section
~~~~~~~~~~~~~~~~~~~~~

``YoungsModulus``
    Young's modulus [Pa].  Typical lithospheric value: 65 GPa (``6.5e10``).

``PoissonsRatio``
    Poisson's ratio [dimensionless].  Typical value: 0.25.

``GravAccel``
    Gravitational acceleration [m s⁻²].  Earth standard: 9.8.

``MantleDensity``
    Density of the mantle [kg m⁻³].  Typical value: 3300.

``InfillMaterialDensity``
    Density of the material that fills (or vacates) the flexural depression
    [kg m⁻³].  Common values:

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
      normal stresses [Pa] (ρ × g × h).  Grid cell area (dx × dy) is
      applied internally to convert stress to force.
    * *SAS_NG*: a space-delimited file with columns ``(x, q0)`` in 1-D or
      ``(x, y, q0)`` in 2-D, where q0 is a point force [N].

    Paths are resolved relative to the directory containing the
    configuration file.

``ElasticThickness``
    Elastic thickness [m].  Either a scalar value or a path to a
    space-delimited array.  Arrays are required for FD solutions with
    spatially variable *Te*.  Use :func:`~gflex.smooth_pad_Te` and
    :func:`~gflex.pad_domain` to extend a variable-*Te* grid with a smooth
    boundary buffer before running.

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

    * ``q0`` — plot the applied load.
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

    For FD solutions:

    * ``0Slope0Shear`` — free end; zero slope and shear force.
    * ``0Moment0Shear`` — free end; zero bending moment and shear force.
      This is the natural (Neumann) boundary condition.
    * ``0Displacement0Slope`` — clamped; zero displacement and slope.
    * ``Mirror`` — mirror the domain at this edge.
    * ``Periodic`` — wrap-around (useful for convergence testing).

    For SAS / SAS_NG: ``NoOutsideLoads`` (assumed if left blank).

    Flexural solutions can be sensitive to boundary conditions; choose
    carefully and consider using :func:`~gflex.pad_domain` to push
    boundaries away from the region of interest.

``Solver``
    Linear system solver for FD:

    * ``direct`` — sparse direct solver; recommended for most grids.
    * ``iterative`` — lower peak memory on very large grids; slower.

``ConvergenceTolerance``
    Maximum allowable change between successive iterative solver steps [m].
    Only used when ``Solver = iterative``.  Set to ``0`` to run a fixed
    number of iterations without a tolerance check.

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
      BoundaryCondition_West: Periodic
      BoundaryCondition_East: Periodic
      Solver: direct
      ConvergenceTolerance: 0.001

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
      PlateSolutionType: vWC1994   # van Wees & Cloetingh (1994); recommended

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
      BoundaryCondition_West: 0Moment0Shear
      BoundaryCondition_East: 0Displacement0Slope
      Solver: direct
      ConvergenceTolerance: 1.0e-3

    numerical2D:
      GridSpacing_y: 4000
      BoundaryCondition_North: Mirror
      BoundaryCondition_South: 0Slope0Shear

    verbosity:
      Verbose: false
      Debug: false
      Quiet: false

INI format (legacy)
~~~~~~~~~~~~~~~~~~~

The INI format uses the same section names and parameter keys as YAML, with
``[section]`` headers and ``key=value`` pairs.  Comments begin with ``;``.
The file extension is not restricted — the absence of a ``.yaml`` / ``.yml``
suffix causes gFlex to treat the file as INI.

Equivalent 1-D example in INI:

.. code-block:: ini

    ; All units are SI.

    [mode]
    dimension=1
    method=FD

    [parameter]
    YoungsModulus=6.5E10
    PoissonsRatio=0.25
    GravAccel=9.8
    MantleDensity=3300
    InfillMaterialDensity=0

    [input]
    Loads=q0_sample/1D/central_block.txt
    ElasticThickness=Te_sample/1D/8km_20km_ramp.txt

    [output]
    DeflectionOut=
    Plot=combo

    [numerical]
    GridSpacing_x=6000
    BoundaryCondition_West=Periodic
    BoundaryCondition_East=Periodic
    Solver=direct
    ConvergenceTolerance=0.001

    [verbosity]
    Verbose=false
    Debug=false
    Quiet=false
