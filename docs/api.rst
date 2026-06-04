API Reference
=============

Primary solvers
---------------

Both solvers follow the same lifecycle::

   flex.initialize()
   flex.run()
   w = flex.w          # read deflection before finalize clears it
   flex.output()       # optional: save to file or display plots
   flex.finalize()     # releases w, qs, and the coefficient matrix

.. warning::

   :meth:`finalize` deletes ``w``, ``qs``, and the cached coefficient
   matrix.  Read ``w`` (and call :meth:`~gflex.base.Flexure.output` if
   needed) **before** calling ``finalize``.  Accessing ``w`` afterwards
   raises :exc:`AttributeError`.

.. autoclass:: gflex.F2D
   :members: initialize, run, finalize

.. autoclass:: gflex.F1D
   :members: initialize, run, finalize

Output
------

.. automethod:: gflex.base.Flexure.output

In-plane stresses
-----------------

In-plane stresses are set as attributes directly on the solver instance
before calling :meth:`initialize`.  They are not available as configuration
file keys.

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Attribute
     - Solvers
     - Description
   * - ``sigma_xx``
     - FD, FFT (1-D and 2-D)
     - Normal stress in the x-direction :math:`\sigma_{xx}` [Pa].
       Default ``0``.
   * - ``sigma_yy``
     - FD, FFT (2-D only)
     - Normal stress in the y-direction :math:`\sigma_{yy}` [Pa].
       Default ``0``.
   * - ``sigma_xy``
     - FD, FFT (2-D only)
     - Shear stress :math:`\sigma_{xy}` [Pa].  Default ``0``.

All three default to zero if not assigned; setting any of them with
``SAS`` or ``SAS_NG`` raises a :exc:`RuntimeWarning` and has no effect.
See :doc:`theory` for the governing equations that include these terms.

Domain-padding utilities
------------------------

These functions help when running :class:`~gflex.F1D` or :class:`~gflex.F2D`
with a spatially variable elastic thickness grid.  A smooth padding zone
reduces spurious deflections at the domain boundary caused by sharp rigidity
gradients, and ensures that the flexural forebulge can develop freely before
reaching the boundary.

**2-D (F2D)**

.. autofunction:: gflex.pad_domain

.. autofunction:: gflex.smooth_pad_Te

.. autofunction:: gflex.recommended_pad_width

**1-D (F1D)**

.. autofunction:: gflex.pad_domain_1d

.. autofunction:: gflex.smooth_pad_Te_1d

.. autofunction:: gflex.recommended_pad_width_1d

FD boundary-condition warnings
------------------------------

When running :class:`~gflex.F1D` or :class:`~gflex.F2D` with the finite-difference
solver, gFlex issues :exc:`UserWarning` messages for two categories of potentially
problematic boundary conditions.

**BC-type warnings** fire whenever a side carries a BC whose physical
interpretation deserves verification:

* ``'zero_moment_zero_shear'`` (alias ``'free'``) — assumes a free broken plate end (zero moment and shear
  force).  Physically appropriate for rifted or passive continental margins,
  subduction trenches with an applied edge load (slab pull), and broken-plate
  flexure (Turcotte & Schubert).  Often applied uncritically elsewhere in the
  literature — verify that one of these settings applies.

**Proximity warnings** fire for ``'zero_displacement_zero_slope'`` (alias ``'clamped'``) boundaries when the
nearest loaded cell is within one flexural wavelength
(:math:`\lambda = 2\pi\alpha`, where
:math:`\alpha = (4D / \Delta\rho g)^{1/4}`)
of the boundary.  Within this distance the flexural forebulge — which peaks at
:math:`\approx \pi\alpha` from the load — will be suppressed by the zero-displacement
condition, contaminating the solution.  The warning message reports the distance as
a fraction of the local flexural wavelength and directs you to the domain-padding
utilities.

**Warning deduplication in model-coupling loops**

Python's default warning filter shows each unique warning *once per call site*
per interpreter session.  In a time-stepping or iterative-coupling loop such as::

   for load in loads:
       flex.qs = load
       flex.run()          # warning fires on first iteration, silenced thereafter

the proximity warning fires on the first iteration and is not repeated, even if
the load subsequently moves closer to the boundary.  To re-enable the warning on
every call::

   import warnings
   warnings.filterwarnings("always", category=UserWarning, module="gflex")

**Suppressing warnings you have verified**

Once you have confirmed that a boundary condition is appropriate for your setup,
suppress the corresponding warning by message text::

   import warnings
   warnings.filterwarnings("ignore", message=".*zero_moment_zero_shear.*")

To suppress all gFlex warnings::

   warnings.filterwarnings("ignore", module="gflex")

Or use a context manager for a single run::

   with warnings.catch_warnings():
       warnings.simplefilter("ignore")
       flex.run()

LU factorization cache
----------------------

For coupling workflows that call :meth:`~gflex.F1D.run` (or
:meth:`~gflex.F2D.run`) repeatedly with the same grid, elastic thickness,
and boundary conditions, the sparse-LU factorization of the coefficient
matrix can be cached to avoid re-factorizing on every call.  Set the
attribute **before** calling :meth:`~gflex.base.Flexure.initialize`:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Value
     - Behaviour
   * - ``False`` *(default)*
     - No caching.  The matrix is factorized on every :meth:`run` call.
   * - ``True``
     - Cache with hash check.  The factorization is reused when a
       hash of the coefficient matrix matches the stored hash; it is
       recomputed when any of Te, dx/dy, BCs, or physical parameters change.
   * - ``"no_check"``
     - Cache without hash check.  The stored factorization is reused on
       every call without computing a hash.  Gives the maximum performance
       benefit.  The coefficient matrix is freed from memory immediately
       after factorization; only the LU factorization is retained.  Smart
       invalidation (see below) still applies: reassigning a
       matrix-determining input clears the cache and triggers a rebuild on
       the next call.

**Smart invalidation**

Reassigning any matrix-determining attribute — ``T_e``, ``E``, ``nu``,
``g``, ``rho_m``, ``rho_fill``, ``dx``, ``dy``, boundary conditions, or
in-plane stresses — automatically clears the cached coefficient matrix and
LU factorization.  No explicit cache management is needed between solves
when only ``qs`` changes.

.. note::

   Smart invalidation is triggered by *assignment* (``flex.T_e = new_array``),
   not by in-place mutation of a NumPy array (``flex.T_e[5] = 40e3``).  If
   you mutate an array in place, reassign it afterwards to ensure the cache
   is correctly invalidated::

      flex.T_e[5] = 40e3
      flex.T_e = flex.T_e   # trigger invalidation

Example (coupling loop)::

   flex.cache_factorization = True
   flex.initialize()

   for load in load_sequence:
       flex.qs = load
       flex.run()
       w = flex.w
       # ... process w ...

   flex.finalize()

The cache is cleared by :meth:`~gflex.base.Flexure.finalize`.

Flexural wavelength
-------------------

.. autofunction:: gflex.flexural_wavelengths

BMI interface
-------------

:class:`~gflex.BmiGflex` exposes the CSDMS Basic Model Interface, enabling
gFlex to be coupled with other models in the CSDMS framework.  It requires
the optional ``bmipy`` dependency (``pip install gflex[bmi]``).

.. autoclass:: gflex.BmiGflex
   :members: initialize, update, finalize, get_value, set_value
