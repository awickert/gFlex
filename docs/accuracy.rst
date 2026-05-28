Numerical Accuracy
==================

2-D finite-difference solver: second-order convergence
-------------------------------------------------------

The 2-D finite-difference solver (``Method = 'FD'``, ``PlateSolutionType =
'vWC1994'``) is second-order accurate in space: halving the grid spacing *dx*
reduces the numerical error by a factor of approximately four (O(dx²)).

This is verified two ways in gFlex:

1. **Method of Manufactured Solutions (MMS)** — a formal unit test
   (:func:`tests.test_2D_FD.test_2d_fd_convergence_order`) that runs with every
   CI build.  A cosine deflection field is chosen as the exact solution; the
   corresponding load is derived analytically and fed to the solver.  Three grid
   refinements (N = 30, 60, 120 cells per side) on a periodic domain confirm a
   convergence rate > 1.8 at all refinement levels.

2. **Kelvin-function benchmark** — a grid-convergence study comparing the FD
   solver to the analytical Kelvin-function solution for an infinite plate under
   a point load.  Results are summarised in the table below.

Kelvin-function benchmark setup
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 30 70
   :header-rows: 0

   * - Young's modulus *E*
     - 65 GPa
   * - Elastic thickness *Te*
     - 10 km
   * - Poisson's ratio *ν*
     - 0.25
   * - Mantle density *ρ*\ :sub:`m`
     - 3300 kg m⁻³
   * - Infill density *ρ*\ :sub:`fill`
     - 0 kg m⁻³ (air)
   * - *g*
     - 9.81 m s⁻²
   * - Flexural parameter *α*
     - ≈ 21 km
   * - Domain
     - 600 km × 600 km, point load at centre
   * - Boundary conditions
     - ``0Moment0Shear`` on all four sides
   * - Grid spacings tested
     - *dx* = 10 000, 5 000, 2 500, 1 250 m

Convergence orders (coarse → fine)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The table gives the local convergence order between successive grid refinements
at several distances from the load, expressed as multiples of the flexural
parameter *α*.

.. list-table::
   :header-rows: 1
   :widths: 15 20 20 20

   * - *r* / *α*
     - dx 10→5 km
     - dx 5→2.5 km
     - dx 2.5→1.25 km
   * - 0.97
     - 2.20
     - 2.07
     - 2.01
   * - 1.46
     - 2.68
     - 2.18
     - 2.04
   * - 1.95
     - 1.07
     - 1.88
     - 1.97
   * - 2.92
     - −0.07
     - 1.77
     - 1.95

At the finest spacing (dx = 1 250 m, dx/α ≈ 0.06) relative errors are
0.02–0.1 % everywhere except very close to the singularity of the point load.

The sub-second-order rates at *r/α* = 1.95 and 2.92 on the coarsest grids
reflect pre-asymptotic effects near the forebulge zero-crossing, where the
solution changes sign and the absolute error temporarily dominates the relative
error.  All rates converge to ≈ 2 at finer resolution, confirming asymptotic
second-order behaviour.

Practical guidance
~~~~~~~~~~~~~~~~~~

For most geoscience applications, **dx/α ≤ 0.1** (ten cells per flexural
parameter) is sufficient to keep errors below 1 %.  At **dx/α ≈ 0.5** (two
cells per flexural parameter) errors can be several percent, particularly near
the load centre and the forebulge.  Use :func:`gflex.flexural_wavelengths` to
estimate *α* before choosing a grid spacing.
