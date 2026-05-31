Theory
======

.. note::

   The derivations and notation follow `Wickert (2016)
   <https://doi.org/10.5194/gmd-9-997-2016>`_ for the methods present at
   v1.0 (SAS, SAS_NG, FD without in-plane stresses).  In-plane stress
   support and the FFT spectral solver, added in later versions, are
   synthesized here following the same conventions.  Please refer to
   Wickert (2016) for the full derivations, figures, and authoritative
   treatment of the v1.0 methods.  This page is AI-generated.

Flexure of the lithosphere is the process by which loads bend the elastic
outer shell of Earth or other planets (Watts, 2001; Watters and McGovern,
2006).  The sources of these loads are wide-ranging, encompassing volcanic
islands and seamounts, mountain-belt-forming thrust sheets, sedimentary
basins, continental ice sheets, lakes, seas and oceans, erosional unloading,
and mantle plumes — among others.  By bending over distances of several tens
to hundreds of kilometers, the lithosphere low-pass filters a discontinuous
surface loading field into a smoothed solid-Earth response.

Governing equations
-------------------

gFlex solves the thin elastic plate (Kirchhoff–Love) equation for flexural
isostasy.  In compact tensor form the two-dimensional governing equation is

.. math::

   \nabla^2(D\,\nabla^2 w) - T_e\,\boldsymbol{\sigma} : \nabla\nabla w + \Delta\rho\,g\,w = q,

where :math:`\boldsymbol{\sigma} : \nabla\nabla w` is the double contraction
(Frobenius inner product) of the in-plane stress tensor with the Hessian of
:math:`w`.  For spatially uniform :math:`D`,
:math:`\nabla^2(D\,\nabla^2 w) = D\nabla^4 w`; expanding the stress term and
writing the one- and two-dimensional forms for uniform :math:`D`
(Wickert, 2016, Eqs. 1–2, extended to include in-plane stresses):

.. math::

   D \frac{\mathrm{d}^4 w}{\mathrm{d}x^4}
   - \sigma_{xx}\,T_e\,\frac{\mathrm{d}^2 w}{\mathrm{d}x^2}
   + \Delta\rho\,g\,w = q,

.. math::

   D \nabla^4 w
   - \sigma_{xx}\,T_e\,\frac{\partial^2 w}{\partial x^2}
   - \sigma_{yy}\,T_e\,\frac{\partial^2 w}{\partial y^2}
   - 2\sigma_{xy}\,T_e\,\frac{\partial^2 w}{\partial x \partial y}
   + \Delta\rho\,g\,w = q.

Here, :math:`w` [m] is vertical deflection of the plate (:math:`w < 0` is
downward into the mantle), :math:`q` [Pa] is the applied surface normal
stress, :math:`\Delta\rho = \rho_m - \rho_f` [kg m⁻³] is the density of the
mantle minus the density of the infilling material, :math:`g` [m s⁻²] is
gravitational acceleration, and :math:`D` [N m] is the flexural rigidity
(uniform for the analytical solutions, spatially variable for the finite
difference solver).  The terms :math:`\sigma_{xx}`, :math:`\sigma_{yy}` [Pa]
are the normal components of the in-plane stress tensor and
:math:`\sigma_{xy}` [Pa] is its shear component; each is multiplied by the
elastic thickness :math:`T_e` [m] to convert a depth-averaged stress to a
force per unit width acting on the plate edge.  In one dimension only
:math:`\sigma_{xx}` enters; :math:`\sigma_{yy}` and :math:`\sigma_{xy}` are
absent.  All three stress terms default to zero, recovering the purely
biharmonic equations of Wickert (2016); they are supported by the
finite difference and FFT solvers but not by the analytical (SAS/SAS\_NG)
solutions.

The :math:`\Delta\rho\,g\,w` term represents the feedback by which flexural
subsidence can lead a depression to be filled by material — for example,
seawater or sediment — which leads to additional flexural subsidence.  If the
infilling material is not uniform in density or spatial extent, one may solve
for the flexural response with :math:`\rho_f = \rho_\text{air} \approx 0`,
add loads based on conditions that match a given inundation or deposition
rule, and then re-calculate flexure iteratively until convergence is achieved.

For finite difference solutions, where :math:`D` may vary spatially,
expanding :math:`\nabla^2(D\,\nabla^2 w)` by the chain rule gives

.. math::

   \nabla^2(D\,\nabla^2 w)
   = D\nabla^4 w
   + 2\,\nabla D \cdot \nabla(\nabla^2 w)
   + \nabla^2 D \cdot \nabla^2 w,

showing that spatial variation in :math:`D` contributes correction terms
proportional to the gradient and Laplacian of :math:`D`; these vanish for
uniform :math:`D`.  In one dimension the compact form is

.. math::

   \frac{\partial^2}{\partial x^2}\!\left(D\,\frac{\partial^2 w}{\partial x^2}\right)
   - \sigma_{xx}\,T_e\,\frac{\partial^2 w}{\partial x^2}
   + \Delta\rho\,g\,w = q,

which expands by the product rule to

.. math::

   D \frac{\partial^4 w}{\partial x^4}
   + 2 \frac{\partial D}{\partial x} \frac{\partial^3 w}{\partial x^3}
   + \frac{\partial^2 D}{\partial x^2} \frac{\partial^2 w}{\partial x^2}
   - \sigma_{xx}\,T_e\,\frac{\partial^2 w}{\partial x^2}
   + \Delta\rho\,g\,w = q.

The two-dimensional variable-:math:`D` form (van Wees and Cloetingh, 1994)
includes the chain-rule terms above and adds a further cross-derivative
coupling between the second derivatives of :math:`D` and the Hessian of
:math:`w` that arises from the plate bending constitutive relations:

.. math::

   \nabla^2(D\,\nabla^2 w)
   - (1 - \nu)\!\left[
     \frac{\partial^2 D}{\partial x^2}\frac{\partial^2 w}{\partial y^2}
     - 2\frac{\partial^2 D}{\partial x\,\partial y}\frac{\partial^2 w}{\partial x\,\partial y}
     + \frac{\partial^2 D}{\partial y^2}\frac{\partial^2 w}{\partial x^2}
   \right]
   - \sigma_{xx}\,T_e\,\frac{\partial^2 w}{\partial x^2}
   - \sigma_{yy}\,T_e\,\frac{\partial^2 w}{\partial y^2}
   - 2\sigma_{xy}\,T_e\,\frac{\partial^2 w}{\partial x\,\partial y}
   + \Delta\rho\,g\,w = q,

where the bracketed term vanishes when :math:`D` is uniform, recovering
:math:`D\nabla^4 w`.  Both the 1-D and 2-D variable-:math:`D` equations are
discretized using a second-order centered finite difference approximation,
reducing the problem to the sparse linear matrix equation

.. math::

   \mathbf{A}\mathbf{W} = \mathbf{Q},

where :math:`\mathbf{A}` is a sparse matrix of finite difference operators,
:math:`\mathbf{W}` is a vector of deflections, and :math:`\mathbf{Q}` is a
vector of imposed loads.

Flexural rigidity and elastic thickness
---------------------------------------

The scalar flexural rigidity :math:`D` [N m] is the key parameter that
controls the flexural response, and is a function of :math:`T_e`, :math:`E`,
and :math:`\nu` (Turcotte and Schubert, 2002):

.. math::

   D = \frac{E\,T_e^3}{12\!\left(1 - \nu^2\right)}.

Here, :math:`E` [Pa] is Young's modulus, :math:`T_e` [m] is the elastic
thickness, and :math:`\nu` is Poisson's ratio.  Because :math:`D` scales as
:math:`T_e^3`, changes in the effective elastic thickness of the lithosphere,
cubed, are more significant than changes in Poisson's ratio, squared, or
Young's modulus.  gFlex contains the additional simplifying assumption that
:math:`E` and :math:`\nu` are uniform constants, permitting variations in
scalar flexural rigidity to map to variations in effective elastic thickness
via this equation.

Flexural parameter
------------------

For a plate of constant :math:`T_e`, the *flexural parameter* :math:`\alpha`
[m] provides the characteristic length scale of the deflection.  Following
Vening Meinesz (1931, after Hertz, 1884):

.. math::

   \alpha_\text{1D} = \left[\frac{4D}{\Delta\rho\,g}\right]^{1/4}, \qquad
   \alpha_\text{2D} = \left[\frac{D}{\Delta\rho\,g}\right]^{1/4}.

The significance of the flexural parameter is that the flexural wavelength
:math:`\lambda_\alpha` is related to it as :math:`\lambda_\alpha = 2\pi\alpha`.
The distance from a point load to the first flexural bulge (forebulge) that it
creates around its local depression, for example, is a flexural half-wavelength,
:math:`\pi\alpha`.  This nature of plate bending as an exponentially decaying
periodic function can be seen most easily in the one-dimensional analytical
(constant :math:`T_e`) solution, Eq. (3) below.

Use :func:`gflex.flexural_wavelengths` to compute :math:`\alpha`,
:math:`\lambda_\alpha`, and the first zero-crossing for a given set of elastic
parameters before choosing a grid spacing or domain size.

Solution methods
----------------

gFlex solves for lithospheric flexure in two major ways.  First, it can
produce analytical solutions to flexural isostasy generated by superposition
of local solutions to point loads in the spatial domain.  Second, it can
compute finite difference solutions for both constant and arbitrarily varying
lithospheric elastic thickness structures.  These solutions are formulated for
both one-dimensional (line load, assumed to extend infinitely in an orientation
orthogonal to the line along which the equation is solved) and two-dimensional
(point load) cases (Wickert, 2016).

Superposition of Analytical Solutions (``SAS``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first solution type takes advantage of the linear nature of the analytical
solution for flexure of a plate of constant thickness and elastic properties
when subjected to a point or line load.  Because the flexure equation is
linear, these solutions may be superposed (i.e., summed) in space to compute
the full flexural response to any arbitrary load.

In one dimension, the response at position :math:`x` to a line load
:math:`q` [Pa] at position :math:`x_i` is (Wickert, 2016, Eq. 3)

.. math::

   w_i = q\,\frac{\alpha_\text{1D}^3}{8D}\,
         e^{-|x - x_i| / \alpha_\text{1D}}
         \!\left[
           \cos\!\frac{|x - x_i|}{\alpha_\text{1D}}
         + \sin\!\frac{|x - x_i|}{\alpha_\text{1D}}
         \right].

In two dimensions, the response at :math:`(x, y)` to a point load at
:math:`(x_i, y_j)` involves the zeroth-order Kelvin function
:math:`\mathrm{kei}` (Brotchie and Silvester, 1969; Abramowitz and Stegun,
1972):

.. math::

   w_{i,j} = q\,\frac{\alpha_\text{2D}^2}{2\pi D}\,
              \mathrm{kei}\!\left(
                \frac{\sqrt{(x - x_i)^2 + (y - y_j)^2}}{\,\alpha_\text{2D}}
              \right).

The implicit boundary condition for all analytical solutions is
``NoOutsideLoads``: the plate is undeflected at infinity
(:math:`w_\infty = 0`).

Superposition of Analytical Solutions, No Grid (``SAS_NG``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The same analytical superposition as ``SAS``, but load and output locations
are supplied as an arbitrary (x, q₀) or (x, y, q₀) point set rather than a
regular grid.  This solution type is termed SAS_NG — superposition of
analytical solutions, no grid — because it lacks the grid uniformity that
permits a solution template to be used, and so its computational time is not
optimized in the same way (Wickert, 2016, Sect. 2.5).

Finite Difference (``FD``)
~~~~~~~~~~~~~~~~~~~~~~~~~~

Finite difference solutions employ the variable-:math:`D` expansions of the
governing equations on a regular Cartesian grid and, following van Wees and
Cloetingh (1994), permit computations with spatially varying flexural rigidity.
The grid spacings :math:`\Delta x` and :math:`\Delta y` may differ from one
another, but each must be constant.  The resulting sparse linear system is
solved directly using a sparse LU factorization or, at the user's choice,
iteratively.  Six boundary conditions are available; see :doc:`configuration`
for details and physical interpretations.

Fast Fourier Transform (``FFT``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The FFT spectral solver exploits the fact that, for a
uniform flexural rigidity :math:`D`, the thin plate equation is a convolution
in the spatial domain — and therefore diagonal in the wavenumber domain.
Taking the two-dimensional Fourier transform of the governing equation with
in-plane stresses yields an algebraic expression for the deflection spectrum:

.. math::

   \hat{W}(k_x, k_y) =
   \frac{-\hat{Q}(k_x,k_y)}{%
     D\!\left(k_x^2 + k_y^2\right)^{\!2}
     + \sigma_{xx}\,T_e\,k_x^2
     + \sigma_{yy}\,T_e\,k_y^2
     + 2\sigma_{xy}\,T_e\,k_x k_y
     + \Delta\rho\,g},

where :math:`k_x` and :math:`k_y` are angular wavenumbers [m⁻¹].  The
deflection field :math:`w` is recovered by the inverse FFT.

Two boundary condition modes are available.  When all four sides are set to
``Periodic``, the load array is transformed as-is and the solution is
inherently periodic.  For all other boundary conditions, the load is
zero-padded by four flexural parameters (:math:`4\alpha`) on each side before
the transform and trimmed back afterward, producing a result equivalent to the
``NoOutsideLoads`` condition of the analytical solutions.

Because the transfer function assumes a single, spatially uniform value of
:math:`D`, the FFT solver requires scalar (constant) :math:`T_e`.  Spatially
variable elastic thickness requires the finite difference solver.
