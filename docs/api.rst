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
