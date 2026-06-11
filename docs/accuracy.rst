Numerical Accuracy
==================

Each solution method in gFlex has a distinct accuracy profile.

Analytical solutions (``sas`` and ``sas_ng``)
----------------------------------------------

The superposition-of-analytical-solutions methods are **exact** for a plate
of constant elastic thickness under the ``no_outside_loads`` boundary condition
(zero deflection at infinity).  In one dimension the Green's function is the
exponential sinusoid of Eq. (3) in Wickert (2016); in two dimensions it is the
zeroth-order Kelvin function :math:`\mathrm{kei}`, evaluated via
:func:`scipy.special.kei`.  Floating-point errors in these special-function
evaluations are at machine precision and negligible in practice.

The only meaningful source of error is a **domain that is too small**: if the
load produces non-trivial deflection at the domain boundary, the
``no_outside_loads`` assumption is violated and the solution is physically
wrong — not because of numerics, but because the boundary condition does not
match the situation.  Use :func:`gflex.flexural_wavelengths` to estimate the
flexural parameter :math:`\alpha` and ensure the domain extends several
:math:`\alpha` beyond the loaded region.  Unlike the finite difference solver,
there is no grid-spacing error; output sampling density affects only the
display, not the solution.

FFT spectral solver
-------------------

The FFT solver is **spectrally accurate** for uniform elastic thickness: for a
smooth load field the solution error decays faster than any finite power of the
grid spacing, limited in practice only by floating-point arithmetic.  For
periodic boundary conditions (``periodic`` on all sides) the solution is exact
to machine precision.  For all other boundary conditions the load is
zero-padded by :math:`4\alpha` on each side before the transform, which
approximates the ``no_outside_loads`` condition; the padding introduces a small
error near the domain edges that decreases as the pad width increases relative
to the load's flexural footprint.  For typical geoscience applications the
:math:`4\alpha` default padding is more than sufficient.

Because the FFT assumes a single, spatially uniform flexural rigidity
:math:`D`, it cannot represent variable-:math:`T_e` problems; those require
the finite difference solver.

Finite difference solver (``fd``)
---------------------------------

.. figure:: _static/fig3_fd_vs_sas.png
   :width: 90%
   :align: center
   :alt: Comparison of FD and analytical (SAS) deflection solutions in 1-D and 2-D

   Comparison of numerical (FD) and analytical (SAS) solutions in one dimension
   **(a)** and two dimensions **(c)**, and their differences **(b, d)**, for a
   100 km central line load / circular load.  The 1–2 m offset in **(b)** is due
   primarily to the ``no_outside_loads`` BC of the analytical solution versus the
   ``zero_displacement_zero_slope`` BC of the FD solution; the cross-shaped residual in
   **(d)** reflects boundary effects along the longer diagonal boundaries.
   Reproduced from Wickert (2016), Fig. 3;
   `CC BY 3.0 <https://creativecommons.org/licenses/by/3.0/>`_.

----

2-D finite-difference solver: second-order convergence
-------------------------------------------------------

The 2-D finite-difference solver (``method = 'fd'``) is second-order accurate in space: halving the grid spacing
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
     - 9.8 m s⁻²
   * - Flexural parameter :math:`\alpha`
     - ≈ 21 km
   * - Domain
     - 600 km × 600 km, point load at centre
   * - Boundary conditions
     - ``zero_moment_zero_shear`` on all four sides
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

----

1-D ``zero_displacement_zero_slope`` boundary condition: MMS verification
--------------------------------------------------------------------------

The ``zero_displacement_zero_slope`` (clamped) boundary condition enforces
:math:`w = 0` and :math:`dw/dx = 0` at the plate edge.  The manufactured
solution

.. math::

   w_\mathrm{exact}(\xi) = -W_0\,\xi^2\,(1-\xi)^2, \quad \xi = x/L,

satisfies both conditions at :math:`\xi = 0` and :math:`\xi = 1` exactly.
Its fourth derivative is constant, so the manufactured load

.. math::

   q_s(\xi) = \frac{24\,D\,W_0}{L^4}
              + \Delta\rho\,g\,W_0\,\xi^2(1-\xi)^2

includes a spatially varying elastic-foundation term — making this a
nontrivial test of the full governing equation, not just the biharmonic
operator alone.  The error metric is the :math:`L^\infty` relative error:

.. math::

   e = \frac{\max|w_\mathrm{num} - w_\mathrm{exact}|}
            {\max|w_\mathrm{exact}|}

Physical parameters: :math:`T_e = 30` km, :math:`E = 65` GPa,
:math:`\nu = 0.25`, :math:`\rho_m = 3300` kg m⁻³, :math:`\rho_\mathrm{fill}
= 0`, :math:`g = 9.8` m s⁻², :math:`L = 600` km.

Results
~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 10 12 16

   * - :math:`n_x`
     - :math:`\Delta x` [km]
     - :math:`L^\infty` error
   * - 26
     - 24.0
     - 1.84 × 10⁻³
   * - 51
     - 12.0
     - 4.55 × 10⁻⁴
   * - 101
     - 6.0
     - 1.14 × 10⁻⁴
   * - 201
     - 3.0
     - 2.85 × 10⁻⁵
   * - 401
     - 1.5
     - 7.12 × 10⁻⁶
   * - 801
     - 0.75
     - 1.78 × 10⁻⁶

Convergence slope (finest three points): :math:`\mathcal{O}(\Delta x^{2.00})`.
The test function
:func:`tests.test_bc_mms.TestClampedBC1D.test_convergence_order`
runs with every CI build.

----

1-D and 2-D ``zero_moment_zero_shear`` boundary condition: MMS verification
----------------------------------------------------------------------------

The manufactured solution

.. math::

   w_\mathrm{exact}(\xi) = -W_0\,\xi^4\,(1-\xi)^4, \quad \xi = x/L,

satisfies all four free-end boundary conditions
(:math:`w''=w'''=0`) at both ends exactly.  Its fourth derivative is

.. math::

   \frac{d^4 w}{dx^4} = -\frac{W_0}{L^4}
   \bigl(24 - 480\xi + 2160\xi^2 - 3360\xi^3 + 1680\xi^4\bigr),

giving the manufactured load

.. math::

   q_s(\xi) = \frac{D\,W_0}{L^4}
              \bigl(24 - 480\xi + 2160\xi^2 - 3360\xi^3 + 1680\xi^4\bigr)
              + \Delta\rho\,g\,W_0\,\xi^4(1-\xi)^4.

The 2-D extension uses the separable solution
:math:`w = -W_0\,f(\xi)\,f(\eta)` with :math:`f(t) = t^4(1-t)^4`,
which satisfies the free-end condition on all four sides.

Physical parameters: :math:`T_e = 30` km, :math:`E = 65` GPa,
:math:`\nu = 0.25`, :math:`\rho_m = 3300` kg m⁻³, :math:`\rho_\mathrm{fill}
= 0`, :math:`g = 9.8` m s⁻², :math:`L = 600` km,
:math:`W_0 = 25600` m (giving :math:`|w_\mathrm{exact}|_\mathrm{max} = 100`
m).

Results (1-D)
~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 10 12 16

   * - :math:`n_x`
     - :math:`\Delta x` [km]
     - :math:`L^\infty` error
   * - 26
     - 24.0
     - 1.72 × 10⁻²
   * - 51
     - 12.0
     - 4.40 × 10⁻³
   * - 101
     - 6.0
     - 1.11 × 10⁻³
   * - 201
     - 3.0
     - 2.77 × 10⁻⁴
   * - 401
     - 1.5
     - 6.93 × 10⁻⁵
   * - 801
     - 0.75
     - 1.73 × 10⁻⁵

Convergence slope (finest three points): :math:`\mathcal{O}(\Delta x^{2.00})`.

Results (2-D)
~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 10 12 16

   * - :math:`n_x = n_y`
     - :math:`\Delta x` [km]
     - :math:`L^\infty` error
   * - 26
     - 24.0
     - 1.76 × 10⁻²
   * - 51
     - 12.0
     - 4.53 × 10⁻³
   * - 101
     - 6.0
     - 1.15 × 10⁻³
   * - 201
     - 3.0
     - 2.85 × 10⁻⁴

Convergence slope (finest two points): :math:`\mathcal{O}(\Delta x^{2.01})`.
The test functions
:func:`tests.test_bc_mms.TestFreeEndBC1D.test_convergence_order` and
:func:`tests.test_bc_mms.TestFreeEndBC2D.test_convergence_order` run with
every CI build.

----

2-D ``zero_displacement_zero_slope`` boundary condition: MMS verification
--------------------------------------------------------------------------

The exact solution

.. math::

   w_\mathrm{exact}(\xi, \eta)
     = -W_0\, g(\xi)\, g(\eta),
   \quad g(t) = t^2(1-t)^2,
   \quad \xi = x/L,\ \eta = y/L,

satisfies all four clamped boundary conditions (:math:`w = 0`,
:math:`\partial w/\partial n = 0`) exactly.  Because :math:`g''''(t) = 24`
(constant), the manufactured load

.. math::

   q_s(\xi,\eta) =
     \frac{D\,W_0}{L^4}
     \bigl[24\,g(\eta) + 2\,g''(\xi)\,g''(\eta) + 24\,g(\xi)\bigr]
     + \Delta\rho\,g\,W_0\,g(\xi)\,g(\eta)

includes a spatially varying elastic-foundation term.  The error metric is
the :math:`L^\infty` relative error.

Physical parameters: :math:`T_e = 30` km, :math:`E = 65` GPa,
:math:`\nu = 0.25`, :math:`\rho_m = 3300` kg m⁻³, :math:`\rho_\mathrm{fill}
= 0`, :math:`g = 9.8` m s⁻², :math:`L = 600` km, :math:`W_0 = 1600` m
(max :math:`|w_\mathrm{exact}| = 6.25` m).

Results
~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 10 12 16

   * - :math:`n_x = n_y`
     - :math:`\Delta x` [km]
     - :math:`L^\infty` error
   * - 26
     - 24.0
     - 1.61 × 10⁻³
   * - 51
     - 12.0
     - 4.31 × 10⁻⁴
   * - 101
     - 6.0
     - 1.10 × 10⁻⁴
   * - 201
     - 3.0
     - 2.75 × 10⁻⁵

Convergence slope (finest two points): :math:`\mathcal{O}(\Delta x^{1.99})`.
The test function
:func:`tests.test_bc_mms.TestClampedBC2D.test_convergence_order`
runs with every CI build.

----

``zero_slope_zero_shear`` (mirror) boundary condition: MMS verification
------------------------------------------------------------------------

The ``zero_slope_zero_shear`` (mirror) boundary condition enforces
:math:`dw/dx = 0` and :math:`d^3w/dx^3 = 0` at the plate edge.  A
cosine satisfies both conditions identically, giving the manufactured
solution

.. math::

   w_\mathrm{exact}(\xi) = -W_0\,\cos(\pi\xi), \quad \xi = x/L,

and the manufactured load

.. math::

   q_s(\xi) = \Bigl[D\!\left(\frac{\pi}{L}\right)^4 + \Delta\rho\,g\Bigr]
              W_0\cos(\pi\xi).

The 2-D extension uses the separable solution
:math:`w = -W_0\cos(\pi\xi)\cos(\pi\eta)`, for which the load is

.. math::

   q_s(\xi,\eta) = \Bigl[4D\!\left(\frac{\pi}{L}\right)^4 + \Delta\rho\,g\Bigr]
                   W_0\cos(\pi\xi)\cos(\pi\eta).

Physical parameters: :math:`T_e = 30` km, :math:`E = 65` GPa,
:math:`\nu = 0.25`, :math:`\rho_m = 3300` kg m⁻³, :math:`\rho_\mathrm{fill}
= 0`, :math:`g = 9.8` m s⁻², :math:`L = 600` km, :math:`W_0 = 1` m.

Results (1-D)
~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 10 12 16

   * - :math:`n_x`
     - :math:`\Delta x` [km]
     - :math:`L^\infty` error
   * - 26
     - 24.0
     - 9.50 × 10⁻⁶
   * - 51
     - 12.0
     - 2.38 × 10⁻⁶
   * - 101
     - 6.0
     - 5.94 × 10⁻⁷
   * - 201
     - 3.0
     - 1.49 × 10⁻⁷
   * - 401
     - 1.5
     - 3.81 × 10⁻⁸

Convergence slopes (coarse → fine): :math:`\mathcal{O}(\Delta x^{2.06})`,
:math:`\mathcal{O}(\Delta x^{2.03})`, :math:`\mathcal{O}(\Delta x^{2.01})`,
:math:`\mathcal{O}(\Delta x^{1.97})`.

Results (2-D)
~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 10 12 16

   * - :math:`n_x = n_y`
     - :math:`\Delta x` [km]
     - :math:`L^\infty` error
   * - 26
     - 24.0
     - 3.76 × 10⁻⁵
   * - 51
     - 12.0
     - 9.40 × 10⁻⁶
   * - 101
     - 6.0
     - 2.35 × 10⁻⁶
   * - 201
     - 3.0
     - 5.88 × 10⁻⁷

Convergence slopes (coarse → fine): :math:`\mathcal{O}(\Delta x^{2.06})`,
:math:`\mathcal{O}(\Delta x^{2.03})`, :math:`\mathcal{O}(\Delta x^{2.01})`.

The test functions :func:`tests.test_bc_mms.TestMirrorBC1D.test_convergence_order`
and :func:`tests.test_bc_mms.TestMirrorBC2D.test_convergence_order` run with
every CI build and assert slopes > 1.8 at all refinement levels.

----

Variable elastic thickness: 1-D MMS verification
-------------------------------------------------

The finite-difference solver supports spatially variable elastic thickness
:math:`T_e(x)`.  To verify convergence, a linearly varying rigidity

.. math::

   D(\xi) = D_0\,(1 + a\xi), \quad a = 0.5, \quad \xi = x/L,

is constructed by setting

.. math::

   T_e(\xi) = T_{e,\mathrm{ref}}\,(1 + a\xi)^{1/3},

so that :math:`D(\xi) = E T_e^3/[12(1-\nu^2)]` ranges from :math:`D_0`
at :math:`\xi = 0` to :math:`1.5\,D_0` at :math:`\xi = 1`.  The exact
solution is the clamped shape function

.. math::

   w_\mathrm{exact}(\xi) = -W_0\,\xi^2(1-\xi)^2.

Expanding the governing equation
:math:`\partial^2/\partial x^2\!\left[D\,\partial^2 w/\partial x^2\right]
+ \Delta\rho g\,w = -q_s` for linear :math:`D` (so :math:`D'' = 0`) and
noting that :math:`w'''' = 24/L^4` (constant) gives the manufactured load

.. math::

   q_s(\xi) = \frac{D_0}{L^4}\bigl[24(1-a) + 72a\xi\bigr]
              + \Delta\rho g\,\xi^2(1-\xi)^2.

Physical parameters: :math:`T_{e,\mathrm{ref}} = 30` km, :math:`E = 65` GPa,
:math:`\nu = 0.25`, :math:`\rho_m = 3300` kg m⁻³, :math:`\rho_\mathrm{fill}
= 0`, :math:`g = 9.8` m s⁻², :math:`L = 600` km, :math:`W_0 = 1` m,
clamped BCs on both ends.

Results
~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 10 12 16

   * - :math:`n_x`
     - :math:`\Delta x` [km]
     - :math:`L^\infty` error
   * - 26
     - 24.0
     - 2.07 × 10⁻³
   * - 51
     - 12.0
     - 5.10 × 10⁻⁴
   * - 101
     - 6.0
     - 1.27 × 10⁻⁴
   * - 201
     - 3.0
     - 3.17 × 10⁻⁵
   * - 401
     - 1.5
     - 7.92 × 10⁻⁶
   * - 801
     - 0.75
     - 1.98 × 10⁻⁶

Convergence slopes (finest four points): :math:`\mathcal{O}(\Delta x^{2.01})`,
:math:`\mathcal{O}(\Delta x^{2.01})`, :math:`\mathcal{O}(\Delta x^{2.00})`.

The test function :func:`tests.test_bc_mms.TestVariableTe1D.test_convergence_order`
runs with every CI build.

----

Variable elastic thickness: 2-D MMS verification
-------------------------------------------------

The 2-D variable-:math:`T_e` test uses a rigidity that varies in both
directions:

.. math::

   D(\xi,\eta) = D_0\,(1 + a\xi)(1 + b\eta),
   \quad a = 0.5,\quad b = 0.3,\quad \xi = x/L,\ \eta = y/L,

constructed by setting

.. math::

   T_e(\xi,\eta) = T_{e,\mathrm{ref}}\,\bigl[(1+a\xi)(1+b\eta)\bigr]^{1/3}.

The exact solution is the same separable clamped shape as the constant-:math:`T_e`
2-D case:

.. math::

   w_\mathrm{exact}(\xi,\eta) = -g(\xi)\,g(\eta),
   \quad g(t) = t^2(1-t)^2.

Because :math:`D` varies in both :math:`x` and :math:`y`, the governing
equation (van Wees & Cloetingh 1994) contains cross-derivative terms that are
absent in the 1-D case.  For bilinear :math:`D` — where
:math:`\partial^2 D/\partial x^2 = \partial^2 D/\partial y^2 = 0` — the
equation reduces to

.. math::

   \frac{\partial^2}{\partial x^2}\!\left[D\frac{\partial^2 w}{\partial x^2}\right]
   + 2\frac{\partial^2}{\partial x\,\partial y}\!\left[D\frac{\partial^2 w}{\partial x\,\partial y}\right]
   + \frac{\partial^2}{\partial y^2}\!\left[D\frac{\partial^2 w}{\partial y^2}\right]
   + 2(\partial_x D)(\partial_{xyy}w)
   + 2(\partial_y D)(\partial_{xxy}w)
   + 2(1-\nu)\frac{\partial^2 D}{\partial x\,\partial y}\frac{\partial^2 w}{\partial x\,\partial y}
   + \Delta\rho\,g\,w = -q_s,

where the last three terms arise from spatial variation of :math:`D` and are
absent when :math:`D` is uniform.  The final term uses coefficient
:math:`2(1-\nu)` rather than :math:`2`, reflecting the full van Wees &
Cloetingh form (see :doc:`theory_and_numerics`).  A non-zero cross-derivative
term exercises the off-diagonal entries of the variable-:math:`T_e`
finite-difference stencil and provides a stronger verification than the 1-D
strip test alone.  The manufactured load :math:`q_s` is derived analytically
from the equation above with :math:`w = w_\mathrm{exact}`.

Physical parameters: :math:`T_{e,\mathrm{ref}} = 30` km, :math:`E = 65` GPa,
:math:`\nu = 0.25`, :math:`\rho_m = 3300` kg m⁻³, :math:`\rho_\mathrm{fill}
= 0`, :math:`g = 9.8` m s⁻², :math:`L = 600` km, clamped BCs on all four sides.

Results
~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 10 12 16

   * - :math:`n_x = n_y`
     - :math:`\Delta x` [km]
     - :math:`L^\infty` error
   * - 26
     - 24.0
     - 1.86 × 10⁻³
   * - 51
     - 12.0
     - 4.98 × 10⁻⁴
   * - 101
     - 6.0
     - 1.26 × 10⁻⁴
   * - 201
     - 3.0
     - 3.16 × 10⁻⁵

Convergence slopes (coarse → fine): :math:`\mathcal{O}(\Delta x^{1.90})`,
:math:`\mathcal{O}(\Delta x^{1.98})`, :math:`\mathcal{O}(\Delta x^{2.00})`.

The test function :func:`tests.test_bc_mms.TestVariableTe2D.test_convergence_order`
runs with every CI build.
