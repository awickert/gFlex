Boundary Conditions
===================

Boundary conditions specify what happens at the edges of the modelled domain.
gFlex supports six named conditions for the finite-difference (FD) solver;
each imposes constraints on which plate mechanical quantities — deflection,
slope, bending moment, and shear force — vanish at that edge.

The spectral (FFT) and analytical-superposition (SAS / SAS_NG) methods do
not use these named conditions.  FFT zero-pads the domain by
:math:`4\alpha` on each side (approximating ``NoOutsideLoads`` except when
all edges are set to ``Periodic``), and SAS / SAS_NG always assume
``NoOutsideLoads``.

Flexural solutions can be sensitive to boundary conditions.  When in doubt,
use :func:`~gflex.pad_domain` (2-D) or :func:`~gflex.pad_domain_1d` (1-D)
to push the boundaries far from the region of interest, and choose
``0Moment0Shear`` (free end) to minimise their influence.

----

Quantity names
--------------

The boundary condition names encode which plate mechanical quantities vanish
at that edge.  The table below maps each name component to its physical
meaning and derivative order.

.. list-table::
   :widths: 22 12 42 12
   :header-rows: 1

   * - Quantity
     - Symbol
     - Definition (uniform :math:`D`, 1-D)
     - :math:`w`-derivative order
   * - Deflection
     - :math:`w`
     - :math:`w`
     - 0th
   * - Slope
     - :math:`S`
     - :math:`\frac{\mathrm{d}w}{\mathrm{d}x}`
     - 1st
   * - Bending moment
     - :math:`M`
     - :math:`D\,\frac{\mathrm{d}^2 w}{\mathrm{d}x^2}`
     - 2nd
   * - Shear force
     - :math:`V`
     - :math:`D\,\frac{\mathrm{d}^3 w}{\mathrm{d}x^3}`
     - 3rd

----

Conditions
----------

.. list-table::
   :widths: 24 18 17 14 27
   :header-rows: 1

   * - gFlex name
     - Zero at boundary
     - Structural mechanics
     - Geophysical
     - Description
   * - ``0Displacement0Slope``
     - :math:`w`,\ :math:`S`
     - clamped end
     - —
     - No deflection, no rotation
   * - ``0Displacement0Moment``
     - :math:`w`,\ :math:`M`
     - simply supported
     - —
     - No deflection, free to rotate
   * - ``0Moment0Shear``
     - :math:`M`,\ :math:`V`
     - free end
     - broken plate
     - Free, unsupported plate end
   * - ``0Slope0Shear``
     - :math:`S`,\ :math:`V`
     - guided end
     - —
     - Plate is level at edge, free to deflect; no shear
   * - ``Mirror``
     - :math:`S`,\ :math:`V`
     - —
     - —
     - Even reflection; model half of a symmetric system
   * - ``Periodic``
     - —
     - —
     - —
     - Domain tiles infinitely in both directions

.. list-table::
   :widths: 28 42 30
   :header-rows: 1

   * - gFlex name
     - Geological description / use
     - Examples
   * - ``0Displacement0Slope``
     - Plate is fully fixed at the boundary — no deflection and no rotation.
       Functions as a conservative far-field constraint when the domain edge
       is far from the load and the plate is expected to be undisturbed.
     - Far boundary in a synthetic model
   * - ``0Displacement0Moment``
     - Plate is pinned to its reference elevation and free to rotate — the
       natural condition for sine-series (Discrete Sine Transform) solutions.
       Pinning both ends in deflection but freeing them in moment yields an
       odd-periodic extension that supports spectral computation without edge
       artefacts.
     - —
   * - ``0Moment0Shear``
     - Natural free-edge condition: no moment or shear transmitted.  The
       far-field condition for a plate decaying to its reference level, and
       the foundation for broken-plate flexure when paired with an
       edge-applied load that supplies (:math:`M_0`, :math:`V_0`) through
       the loading vector.
     - Far-field boundary of an interior-loaded domain; passive/rifted
       margin; broken-plate flexure (T&S) with edge load at the boundary
       node; subduction trench/outer rise (slab pull as edge-applied
       vertical load)
   * - ``0Slope0Shear``
     - Plate is horizontal at the boundary and free to deflect; no shear
       transmitted.  Primarily a mathematical convenience.
     - —
   * - ``Mirror``
     - Symmetry plane: the system is identical on both sides of the
       boundary.  Model only half of a symmetric domain, halving
       computation.  Naturally compatible with cosine-series (Discrete
       Cosine Transform) solutions.
     - One flank of a mountain range or orogenic belt; one side of a
       continental ice sheet; half of a foreland basin profile
   * - ``Periodic``
     - Domain tiles infinitely in both directions; the solution wraps
       around.  Native to FFT-based spectral solutions.
     - Regularly spaced seamount or volcanic chain; broad-scale FFT
       calculations

The following figure from Wickert (2016) shows schematics of five of the six
conditions (``0Displacement0Moment`` was added after publication):

.. figure:: _static/fig4_bc_schematics.png
   :width: 55%
   :align: center
   :alt: Schematics of the five finite-difference boundary condition types

   Schematics of five FD boundary condition types (a–e) from Wickert (2016),
   Fig. 4 (`CC BY 3.0 <https://creativecommons.org/licenses/by/3.0/>`_).
   ``0Displacement0Moment`` (added post-publication) is not shown.

0Displacement0Slope
~~~~~~~~~~~~~~~~~~~

Zero displacement and zero slope: the plate is fully clamped at the
boundary — no deflection and no rotation.

.. figure:: _static/bc_diagram_0Displacement0Slope.svg
   :width: 80%
   :align: center
   :alt: Diagram of the 0Displacement0Slope (clamped end) boundary condition

   *Clamped end* — zero deflection and zero slope at the boundary.

0Displacement0Moment
~~~~~~~~~~~~~~~~~~~~

Zero displacement and zero bending moment: the classical simply-supported
(pinned) plate end.  The plate is held at zero deflection but is free to
rotate, so no moment is transmitted.  Implemented as a Dirichlet condition
(:math:`w = 0`) at the boundary node and an odd-reflection ghost
(:math:`w_\text{ghost} = -w_\text{interior}`) at the first interior node
to enforce zero curvature.

``Mirror`` and ``0Displacement0Moment`` are reflection boundary conditions
but encode opposite parities.  ``Mirror`` uses an *even* reflection
(:math:`w_\text{ghost} = +w_\text{interior}`): the symmetry plane lies
between the last real node and its ghost, the plate is horizontal at the
boundary, and the deflection there is generally non-zero — making it the
correct choice for modelling one half of a symmetric system.
``0Displacement0Moment`` uses an *odd* reflection
(:math:`w_\text{ghost} = -w_\text{interior}`): the boundary node is the
fixed point of the reflection, so :math:`w = 0` there by definition, and
the plate is free to rotate — the simply-supported end.  Sine modes
satisfy ``0Displacement0Moment``; cosine modes satisfy ``Mirror``.

.. figure:: _static/bc_diagram_0Displacement0Moment.svg
   :width: 80%
   :align: center
   :alt: Diagram of the 0Displacement0Moment (simply supported) boundary condition

   *Simply supported* — zero deflection, free to rotate; no bending moment transmitted.

0Moment0Shear
~~~~~~~~~~~~~

The natural condition at a free edge: no bending moment and no shear
force are transmitted across the boundary (Wickert, 2016, Table 1).  It
is the far-field condition into which a flexed plate decays.  Combined
with an edge-applied load — a vertical point force supplies
:math:`V_0`, a closely-spaced couple supplies :math:`M_0` — it produces
the classical broken-plate response of Turcotte and Schubert: the load
carries the inhomogeneity through the loading vector while the BC matrix
stays homogeneous.

.. figure:: _static/bc_diagram_0Moment0Shear.svg
   :width: 80%
   :align: center
   :alt: Diagram of the 0Moment0Shear (free end) boundary condition

   *Free end* — no bending moment and no shear force; the plate ends freely ("broken plate").

0Slope0Shear
~~~~~~~~~~~~

Zero slope and zero shear force: the plate is level at the boundary but
free to deflect there, with no shear transmitted.  Wickert (2016) calls
this "free displacement of a horizontally clamped boundary."  It is
superficially similar to ``Mirror`` (both enforce zero slope), but uses a
different finite-difference stencil and produces noticeably different
solutions; prefer ``Mirror`` for symmetry problems.

.. figure:: _static/bc_diagram_0Slope0Shear.svg
   :width: 80%
   :align: center
   :alt: Diagram of the 0Slope0Shear (guided end) boundary condition

   *Guided end* — zero slope, no shear force; plate is level and free to deflect.

Mirror
~~~~~~

Even reflection at the boundary: model only half of a symmetric system
(e.g., one flank of a mountain range or ice sheet).  See
`0Displacement0Moment`_ above for the contrast with the odd-reflection
simply-supported condition.

.. figure:: _static/bc_diagram_Mirror.svg
   :width: 80%
   :align: center
   :alt: Diagram of the Mirror (symmetry plane) boundary condition

   *Symmetry plane* — even reflection; use when the system is symmetric about the boundary.

Periodic
~~~~~~~~

Wrap-around: the domain tiles infinitely in both directions.

.. figure:: _static/bc_diagram_Periodic.svg
   :width: 80%
   :align: center
   :alt: Diagram of the Periodic boundary condition

   *Periodic* — the domain wraps around; east and west edges are connected.
