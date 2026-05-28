Numerical Accuracy
==================

.. figure:: _static/fig3_fd_vs_sas.png
   :width: 90%
   :align: center
   :alt: Comparison of FD and analytical (SAS) deflection solutions in 1-D and 2-D

   Comparison of numerical (FD) and analytical (SAS) solutions in one dimension
   **(a)** and two dimensions **(c)**, and their differences **(b, d)**, for a
   100 km central line load / circular load.  The 1–2 m offset in **(b)** is due
   primarily to the ``NoOutsideLoads`` BC of the analytical solution versus the
   ``0Displacement0Slope`` BC of the FD solution; the cross-shaped residual in
   **(d)** reflects boundary effects along the longer diagonal boundaries.
   Reproduced from Wickert (2016), Fig. 3;
   `CC BY 3.0 <https://creativecommons.org/licenses/by/3.0/>`_.

----

2-D finite-difference solver: second-order convergence
-------------------------------------------------------

The 2-D finite-difference solver (``Method = 'FD'``, ``PlateSolutionType =
'vWC1994'``) is second-order accurate in space: halving the grid spacing
:math:`\Delta x` reduces the numerical error by a factor of approximately
four (:math:`\mathcal{O}(\Delta x^2)`).

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

   * - Young's modulus :math:`E`
     - 65 GPa
   * - Elastic thickness :math:`T_e`
     - 10 km
   * - Poisson's ratio :math:`\nu`
     - 0.25
   * - Mantle density :math:`\rho_m`
     - 3300 kg m⁻³
   * - Infill density :math:`\rho_\text{fill}`
     - 0 kg m⁻³ (air)
   * - :math:`g`
     - 9.81 m s⁻²
   * - Flexural parameter :math:`\alpha`
     - ≈ 21 km
   * - Domain
     - 600 km × 600 km, point load at centre
   * - Boundary conditions
     - ``0Moment0Shear`` on all four sides
   * - Grid spacings tested
     - :math:`\Delta x` = 10 000, 5 000, 2 500, 1 250 m

Convergence orders (coarse → fine)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The table gives the local convergence order between successive grid refinements
at several distances from the load, expressed as multiples of the flexural
parameter *α*.

.. list-table::
   :header-rows: 1
   :widths: 15 20 20 20

   * - :math:`r/\alpha`
     - :math:`\Delta x` 10→5 km
     - :math:`\Delta x` 5→2.5 km
     - :math:`\Delta x` 2.5→1.25 km
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

At the finest spacing (:math:`\Delta x = 1250` m,
:math:`\Delta x/\alpha \approx 0.06`) relative errors are 0.02–0.1 %
everywhere except very close to the singularity of the point load.

The sub-second-order rates at :math:`r/\alpha` = 1.95 and 2.92 on the
coarsest grids reflect pre-asymptotic effects near the forebulge
zero-crossing, where the solution changes sign and the absolute error
temporarily dominates the relative error.  All rates converge to ≈ 2 at
finer resolution, confirming asymptotic second-order behaviour.

Practical guidance
~~~~~~~~~~~~~~~~~~

For most geoscience applications, :math:`\Delta x / \alpha \leq 0.1` (ten
cells per flexural parameter) is sufficient to keep errors below 1 %.  At
:math:`\Delta x / \alpha \approx 0.5` (two cells per flexural parameter)
errors can be several percent, particularly near the load centre and the
forebulge.  Use :func:`gflex.flexural_wavelengths` to estimate :math:`\alpha`
before choosing a grid spacing.
