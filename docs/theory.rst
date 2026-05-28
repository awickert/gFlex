Theory
======

gFlex solves the thin elastic plate equation governing lithospheric flexural
isostasy: the downward (or upward) bending of Earth's rigid outer shell in
response to an applied surface load.

Governing equation
------------------

For a two-dimensional plate, the deflection :math:`w(x, y)` [m] satisfies

.. math::

   \nabla^2 \!\left(D\, \nabla^2 w\right)
   - \frac{\partial^2}{\partial x^2}\!\left(\sigma_{xx}\,\frac{\partial w}{\partial x}\right)
   - \frac{\partial^2}{\partial y^2}\!\left(\sigma_{yy}\,\frac{\partial w}{\partial y}\right)
   - 2\frac{\partial^2}{\partial x\,\partial y}\!\left(\sigma_{xy}\,\frac{\partial w}{\partial x}\right)
   + \Delta\rho\,g\,w = q,

where :math:`q(x,y)` [Pa] is the applied surface normal stress,
:math:`\Delta\rho = \rho_m - \rho_\text{fill}` [kg m⁻³] is the density
contrast between mantle and infill material, :math:`g` [m s⁻²] is
gravitational acceleration, and :math:`\sigma_{xx}`, :math:`\sigma_{yy}`,
:math:`\sigma_{xy}` [Pa] are depth-integrated tectonic stresses in the plate
(Wickert, 2016, after van Wees & Cloetingh, 1994).

In the common case of no in-plane tectonic stresses and constant flexural
rigidity, this reduces to the isotropic biharmonic equation

.. math::

   D\,\nabla^4 w + \Delta\rho\,g\,w = q.

The one-dimensional (profile) form is

.. math::

   D\,\frac{\mathrm{d}^4 w}{\mathrm{d}x^4} + \Delta\rho\,g\,w = q(x).

Flexural rigidity
-----------------

The flexural rigidity :math:`D` [N m] depends on the elastic properties of
the plate:

.. math::

   D = \frac{E\,T_e^3}{12\!\left(1 - \nu^2\right)},

where :math:`E` [Pa] is Young's modulus, :math:`T_e` [m] is the elastic
thickness, and :math:`\nu` is Poisson's ratio.  :math:`D` may vary spatially
when :math:`T_e` varies, as is common in the real lithosphere.

Flexural parameter
------------------

For a constant-rigidity plate the characteristic length scale of the
deflection — the *flexural parameter* :math:`\alpha` — is

.. math::

   \alpha_\text{1D} = \left(\frac{4D}{\Delta\rho\,g}\right)^{1/4}, \qquad
   \alpha_\text{2D} = \left(\frac{D}{\Delta\rho\,g}\right)^{1/4}.

The corresponding *flexural wavelength* is
:math:`\lambda = 2\pi\alpha`, and the first zero-crossing of the deflection
(the forebulge) falls near :math:`\pi\alpha`.

Use :func:`gflex.flexural_wavelengths` to compute :math:`\alpha`,
:math:`\lambda`, and the first zero-crossing for a given set of elastic
parameters.

Solution methods
----------------

gFlex offers three methods, selectable via ``Method`` in the configuration or
Python API.

Finite Difference (``FD``)
~~~~~~~~~~~~~~~~~~~~~~~~~~

The finite-difference method discretises the governing PDE on a regular
Cartesian grid and solves the resulting linear system.  Spatially variable
:math:`T_e` (and hence variable :math:`D`) is fully supported.  Five boundary
conditions are available (see :doc:`configuration`).  The direct sparse
solver is recommended; an iterative solver is available for very large grids.

Superposition of Analytical Solutions (``SAS``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For a plate with *constant* :math:`T_e`, the deflection at any point is the
superposition of responses from individual load elements.

In one dimension, the response to a line load :math:`F` [N m⁻¹] at distance
:math:`r` is

.. math::

   w(r) = \frac{\alpha^3}{8D}\,F\,
          e^{-r/\alpha}\!\left(\cos\frac{r}{\alpha} + \sin\frac{r}{\alpha}\right).

In two dimensions, the response to a point load :math:`F` [N] at distance
:math:`r` involves the Kelvin function :math:`\mathrm{kei}`:

.. math::

   w(r) = \frac{\alpha^2}{2\pi D}\,F\;\mathrm{kei}\!\left(\frac{r}{\alpha}\right).

Because SAS bypasses sparse-matrix assembly, it is fast and avoids
boundary-condition artefacts.  The implicit boundary condition is
``NoOutsideLoads``: the plate is undeflected at infinity.

Superposition of Analytical Solutions, No Grid (``SAS_NG``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The same analytical superposition as ``SAS``, but load and output locations
are supplied as arbitrary (x, q₀) or (x, y, q₀) point sets rather than a
regular grid.  Useful when load data are sparsely or irregularly distributed.

Sign convention
---------------

Positive :math:`w` is *upward*; negative :math:`w` is *downward* (into the
mantle).  A surface load applied downward (:math:`q > 0`) produces
:math:`w < 0`.  The output array ``flex.w`` follows this convention.
