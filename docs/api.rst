API Reference
=============

Primary solvers
---------------

Both solvers follow the same three-step workflow::

   flex.initialize()
   flex.run()
   flex.finalize()

After :meth:`finalize`, the deflection is available as ``flex.w``.  Call
:meth:`~gflex.base.Flexure.output` to save results to file or display plots.

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

These functions help when running :class:`~gflex.F2D` with a spatially
variable elastic thickness grid.  A smooth padding zone reduces spurious
deflections at the domain boundary caused by sharp rigidity gradients.

.. autofunction:: gflex.pad_domain

.. autofunction:: gflex.smooth_pad_Te

.. autofunction:: gflex.recommended_pad_width

FD boundary-condition warnings
------------------------------

When running :class:`~gflex.F1D` or :class:`~gflex.F2D` with the finite-difference
solver, gFlex issues :exc:`UserWarning` messages for two categories of potentially
problematic boundary conditions.

**BC-type warnings** fire whenever a side carries a BC whose physical
interpretation deserves verification:

* ``'0Moment0Shear'`` — assumes a free broken plate end (zero moment and shear
  force).  This is physically appropriate for a rifted margin or spreading ridge
  where the plate really is broken, but is often applied uncritically elsewhere
  in the literature.
* ``'0Slope0Shear'`` — requires the plate to be simultaneously horizontal and
  shear-free at the boundary.  No clear geological analog is known for this
  combination in a non-trivial deflection setting.

**Proximity warnings** fire for ``'0Displacement0Slope'`` boundaries when the
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
   warnings.filterwarnings("ignore", message=".*0Moment0Shear.*")

To suppress all gFlex warnings::

   warnings.filterwarnings("ignore", module="gflex")

Or use a context manager for a single run::

   with warnings.catch_warnings():
       warnings.simplefilter("ignore")
       flex.run()

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
