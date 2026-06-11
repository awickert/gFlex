GRASS GIS Interfaces
====================

gFlex ships two GRASS GIS add-on modules that expose the full solver
suite from within a GRASS session: **r.flexure** for gridded (raster)
loads and **v.flexure** for scattered point loads.  Both require
**gFlex ≥ 2.0.0**.

Source code and issue tracker:
`<https://github.com/awickert/grass-addons-gflex>`_

----

Installation
------------

Install gFlex first (the GRASS modules call into it)::

   pip install "gflex>=2.0.0"

Then install the add-ons from the GRASS session::

   g.extension extension=r.flexure
   g.extension extension=v.flexure

----

r.flexure — raster loads
-------------------------

*r.flexure* computes the elastic deflection of the lithosphere in
response to a raster load field.  It supports spatially variable
elastic thickness (FD only) and all boundary conditions.

**Methods** (set with ``method=``):

``fd`` — Finite Difference
    Supports spatially variable elastic thickness via the van Wees &
    Cloetingh (1994) stencil.  Fastest for large grids.  All boundary
    conditions and in-plane stresses are supported.

``fft`` — Fast Fourier Transform
    Requires scalar (uniform) :math:`T_e`.  Fastest option for large
    uniform-:math:`T_e` grids; spectrally accurate.  Supports in-plane
    stresses.  Each axis is handled independently: setting both edges
    of a pair to ``periodic`` makes that axis exactly periodic; any
    other setting zero-pads that axis to approximate
    ``no_outside_loads``.

``sas`` — Superposition of Analytical Solutions
    Requires scalar :math:`T_e`.  Analytically exact; best for
    smaller grids or when accuracy near individual loads matters.
    Boundary conditions are not applicable (the solver inherently
    assumes an infinite plate).

**Key parameters**

``input``
    Raster of surface-normal load stress [Pa],
    i.e. :math:`\rho_\text{load}\,g\,h`.  Computed by the user before
    running the module.

``te``
    Elastic thickness [m or km, set via ``te_units``].  A scalar or
    (for FD) a raster of spatially variable :math:`T_e`.

``g``
    Gravitational acceleration [m s⁻²].

``rho_m``, ``rho_fill``
    Mantle and infill densities [kg m⁻³].

``output``
    Raster of lithospheric deflection [m].  Negative values indicate
    downward displacement.

``northbc``, ``southbc``, ``eastbc``, ``westbc``
    Boundary conditions.  SAS ignores these; FD and FFT accept the
    values below.

**Boundary conditions** (FD and FFT)

+------------------+-------------------------------------+-----------------------------------------------------+
| Name             | gFlex canonical                     | Physical meaning                                    |
+==================+=====================================+=====================================================+
| ``infinite``     | ``no_outside_loads``                | Domain padded by one flexural wavelength; default.  |
+------------------+-------------------------------------+-----------------------------------------------------+
| ``clamped``      | ``zero_displacement_zero_slope``    | Plate fixed and horizontal at edge.                 |
+------------------+-------------------------------------+-----------------------------------------------------+
| ``pinned``       | ``zero_displacement_zero_moment``   | Plate hinged (no displacement, free rotation).      |
+------------------+-------------------------------------+-----------------------------------------------------+
| ``free``         | ``zero_moment_zero_shear``          | Free broken-plate end (passive margin, trench).     |
+------------------+-------------------------------------+-----------------------------------------------------+
| ``mirror``       | ``zero_slope_zero_shear``           | Load and :math:`T_e` reflected across boundary.     |
+------------------+-------------------------------------+-----------------------------------------------------+
| ``periodic``     | ``periodic``                        | Opposite edges connect; must be set in pairs (FD).  |
+------------------+-------------------------------------+-----------------------------------------------------+

**In-plane stresses** (``sigma_xx``, ``sigma_yy``, ``sigma_xy`` [Pa])
add tectonic membrane-stress terms to the governing equation
(FD and FFT only; default 0).

**Latitude/longitude** coordinates are supported via the ``-l`` flag.
The module computes a single :math:`dx` at the midpoint latitude of
the domain; this approximation degrades near the poles.

Set the computational region with *g.region* before running.

**Example — asymmetric volcanic load on heterogeneous lithosphere**

A circular volcanic edifice on a lithosphere whose elastic thickness
rises from 15 km in the west to 35 km in the east across a sigmoid
transition.  The thinner western lithosphere deflects more steeply
over a shorter wavelength; the stiffer eastern side produces a broader
depression with a more prominent forebulge.

.. code-block:: bash

   # Domain: 750 × 750 km at 5 km resolution (projected CRS required)
   g.region n=750000 s=0 e=750000 w=0 res=5000

   # Sigmoid Te: 15 km (west) → 35 km (east)
   r.mapcalc expression="te_sigmoid = 15000 + 20000 / (1 + exp(-(x()-375000) / 20000))"

   # Circular volcanic load (30 MPa) within 61.3 km of the domain centre
   r.mapcalc expression="load_volcano = if(sqrt((x()-375000)^2+(y()-375000)^2) <= 61300, 30000000.0, 0)"

   # FD flexural deflection; default infinite BCs approximate an unbounded plate
   r.flexure method=fd \
       input=load_volcano \
       te=te_sigmoid te_units=m \
       output=w_volcano

   d.rast w_volcano
   d.legend w_volcano

For a broken-plate scenario (e.g. a passive margin or oceanic trench),
set the relevant edge to ``free``::

   r.flexure method=fd input=load_volcano te=te_sigmoid te_units=m \
       output=w_volcano_free_south southbc=free

----

v.flexure — point loads
------------------------

*v.flexure* computes the elastic deflection of the lithosphere due to
a set of point loads stored as a GRASS vector map.  It uses the
Superposition of Analytical Solutions for Non-Gridded points
(SAS_NG) method: each point load is a concentrated force, and the
far-field deflection is the superposition of Kelvin–Bessel (kei)
functions.

Elastic thickness must be **scalar** for *v.flexure*; use *r.flexure*
for spatially variable :math:`T_e`.

**Key parameters**

``input``
    Vector points map.  Each point carries a load in units of force
    [N].  For a distributed load discretized as points, multiply the
    load stress [Pa] by the tributary area [m²] to get [N].

``column``
    Name of the attribute column containing the point force values [N].

``te``
    Scalar elastic thickness [m or km, set via ``te_units``].

``output``
    Raster map of deflections [m] at the spacing and extent of the
    current GRASS computational region.

``w_points`` *(optional)*
    Existing vector points map at which deflection is also evaluated
    (GPS stations, boreholes, tide gauges, etc.).  Deflection values
    are written into the column specified by ``w_column`` (default:
    ``w``).  The raster grid and the w_points locations are solved
    in a single gFlex call — no extra computational cost.

In latitude/longitude coordinates, *v.flexure* automatically computes
great-circle distances between load and output points.

**Example — three seamounts along a volcanic chain**

Three point loads represent seamounts spaced 500 km apart.  Each
force (~5 × 10¹⁶ N) corresponds roughly to a seamount with a 25 km
base radius, 3 km height, and basalt density.  Because the flexural
parameter :math:`\alpha \approx 41\,\text{km}` for
:math:`T_e = 25\,\text{km}` is much smaller than the 500 km spacing,
the three depressions are essentially independent and their deflections
add linearly everywhere.

.. code-block:: bash

   # Domain: 2000 × 1000 km at 20 km resolution (projected CRS required)
   g.region n=1000000 s=0 e=2000000 w=0 res=20000

   # Three seamounts along the chain
   echo "500000 500000
   1000000 500000
   1500000 500000" | v.in.ascii format=point input=- output=seamount_chain

   # Add load column
   v.db.addtable map=seamount_chain
   v.db.addcolumn map=seamount_chain columns="force double precision"
   v.db.update map=seamount_chain column=force value=5e16

   # SAS_NG deflection
   v.flexure input=seamount_chain column=force \
       te=25 te_units=km \
       output=w_chain_rast

   # Optionally evaluate at tide-gauge or GPS stations as well
   v.flexure input=seamount_chain column=force \
       te=25 te_units=km \
       output=w_chain_rast \
       w_points=tide_gauges w_column=w_flexure

   d.rast w_chain_rast
   d.legend w_chain_rast

----

References
----------

Wickert, A. D. (2016), Open-source modular solutions for flexural
isostasy: gFlex v1.0, *Geoscientific Model Development*, *9*\(3), 997–1017,
`doi:10.5194/gmd-9-997-2016 <https://doi.org/10.5194/gmd-9-997-2016>`_.

van Wees, J. D., and S. Cloetingh (1994), A finite-difference technique
to incorporate spatial variations in rigidity and planar faults into
3-D models for lithospheric flexure, *Geophysical Journal International*,
*117*\(1), 179–195,
`doi:10.1111/j.1365-246X.1994.tb03311.x <https://doi.org/10.1111/j.1365-246X.1994.tb03311.x>`_.
