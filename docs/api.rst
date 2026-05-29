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
