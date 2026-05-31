Numerical Methods
=================

gFlex solves the flexure equation by three strategies: superposition of
analytical solutions (SAS / SAS_NG), a spectral FFT solver, and a
finite-difference (FD) solver.  See :doc:`theory` for the governing
equations and a comparison of the strategies.

This page documents the FD interior stencils — the fixed-coefficient
templates applied at every grid point away from the domain boundary.
Stencils that enforce the domain boundary conditions are described in
:doc:`boundary_conditions`.

----

Finite-difference interior stencils
------------------------------------

The FD solver replaces each spatial derivative in the governing PDE
with a second-order centred finite difference (Fornberg, 1988, Table 1)
to assemble a sparse linear system (Wickert, 2016, Eq. 11)

.. math::

   \mathbf{A}\,\mathbf{w} = \mathbf{q},

where :math:`\mathbf{A}` is the coefficient (operator) matrix,
:math:`\mathbf{w}` is the vector of unknown deflections, and
:math:`\mathbf{q}` is the imposed load.

**1-D governing equation**

The one-dimensional flexural isostasy equation for a plate of variable
rigidity :math:`D(x)` (Wickert, 2016, Eq. A19) is

.. math::

   \frac{d^2}{dx^2}\!\left(D\frac{d^2 w}{dx^2}\right) + \Delta\rho\,g\,w = q.

Expanding the outer derivative (Wickert, 2016, Eq. 9) and adding the
in-plane (end-load) normal stress :math:`\sigma_{xx}`:

.. math::

   D\frac{d^4 w}{dx^4}
   + 2\frac{dD}{dx}\frac{d^3 w}{dx^3}
   + \frac{d^2 D}{dx^2}\frac{d^2 w}{dx^2}
   - \sigma_{xx}\,T_e\,\frac{d^2 w}{dx^2}
   + \Delta\rho\,g\,w = q.

For uniform :math:`D`, applying the second-order central-difference
approximation of the fourth derivative at node :math:`i`,

.. math::

   \frac{d^4 w}{dx^4}\bigg|_i
   \approx
   \frac{w_{i-2} - 4\,w_{i-1} + 6\,w_i - 4\,w_{i+1} + w_{i+2}}{\Delta x^4},

gives the stencil equation at each interior node:

.. math::

   \frac{D}{\Delta x^4}
   \bigl(w_{i-2} - 4\,w_{i-1} + 6\,w_i - 4\,w_{i+1} + w_{i+2}\bigr)
   - \frac{\sigma_{xx}\,T_e}{\Delta x^2}
   \bigl(w_{i-1} - 2\,w_i + w_{i+1}\bigr)
   + \Delta\rho\,g\,w_i = q_i.

With :math:`\sigma_{xx} = 0` this is the biharmonic stencil with
integer weights :math:`(+1,\,-4,\,+6,\,-4,\,+1)` scaled by
:math:`D/\Delta x^4`, with :math:`\Delta\rho\,g` absorbed into the
centre coefficient.  For variable :math:`D`, the :math:`D`-gradient
terms additionally couple the stencil weights across neighbouring
nodes.

**1-D interior stencil** — five-point stencil spanning two nodes on
each side of the evaluation point:

.. figure:: _static/stencil_1d_interior.svg
   :width: 90%
   :align: center
   :alt: 1-D interior finite-difference stencil

   Five-point 1-D stencil:
   :math:`\tfrac{D}{\Delta x^4}(+1,\,-4,\,+6,\,-4,\,+1)` applied at
   each interior node (uniform :math:`D`).

**2-D governing equation**

The two-dimensional equation follows van Wees and Cloetingh (1994) and,
in compact form (Wickert, 2016, Eq. A20), is

.. math::

   \hat{\boldsymbol{\kappa}}^T\,\mathbf{D}\,\boldsymbol{\kappa}
   + \Delta\rho\,g\,w = q,

where :math:`\hat{\boldsymbol{\kappa}}` is the curvature operator and
:math:`\mathbf{D}` is the plate rigidity matrix (Wickert, 2016,
Eqs. A9–A12).  Expanded for variable :math:`D` (Wickert, 2016,
Eq. 10):

.. math::

   \begin{aligned}
   &D\frac{\partial^4 w}{\partial x^4}
    + D\frac{\partial^4 w}{\partial y^4}
    + 2D\frac{\partial^4 w}{\partial x^2\partial y^2} \\
   &\quad
    + 2\frac{\partial D}{\partial x}\frac{\partial^3 w}{\partial x^3}
    + \frac{\partial^2 D}{\partial x^2}\frac{\partial^2 w}{\partial x^2}
    + 2\frac{\partial D}{\partial y}\frac{\partial^3 w}{\partial y^3}
    + \frac{\partial^2 D}{\partial y^2}\frac{\partial^2 w}{\partial y^2} \\
   &\quad
    + \nu\frac{\partial^2 D}{\partial y^2}\frac{\partial^2 w}{\partial x^2}
    + \nu\frac{\partial^2 D}{\partial x^2}\frac{\partial^2 w}{\partial y^2}
    + 2\frac{\partial D}{\partial x}\frac{\partial^3 w}{\partial x\partial y^2}
    + 2\frac{\partial D}{\partial y}\frac{\partial^3 w}{\partial x^2\partial y} \\
   &\quad
    + 2(1-\nu)\frac{\partial^2 D}{\partial x\partial y}
      \frac{\partial^2 w}{\partial x\partial y}
    + \Delta\rho\,g\,w = q.
   \end{aligned}

In-plane stresses :math:`\sigma_{xx}`, :math:`\sigma_{yy}`, and
:math:`\sigma_{xy}` add to the left-hand side:

.. math::

   -\sigma_{xx}\,T_e\,\frac{\partial^2 w}{\partial x^2}
   - \sigma_{yy}\,T_e\,\frac{\partial^2 w}{\partial y^2}
   - 2\,\sigma_{xy}\,T_e\,\frac{\partial^2 w}{\partial x\partial y}.

For uniform :math:`D` on a square grid
(:math:`\Delta x = \Delta y = h`), the equation reduces to
:math:`D\nabla^4 w + \Delta\rho\,g\,w = q`, and the thirteen-point
stencil equation at interior node :math:`(i,j)` is

.. math::

   \frac{D}{h^4}\bigl(
     w_{i-2,j} + w_{i+2,j} + w_{i,j-2} + w_{i,j+2}
   - 8\,w_{i-1,j} - 8\,w_{i+1,j}
   - 8\,w_{i,j-1} - 8\,w_{i,j+1} \\
   + 2\,w_{i-1,j-1} + 2\,w_{i+1,j-1}
   + 2\,w_{i-1,j+1} + 2\,w_{i+1,j+1}
   + 20\,w_{i,j}\bigr)
   + \Delta\rho\,g\,w_{i,j} = q_{i,j}.

The integer weights arise from expanding :math:`h^4\nabla^4`: centre
:math:`+20`, four nearest axis neighbours :math:`-8`, four
next-nearest axis neighbours :math:`+1`, and four diagonal
neighbours :math:`+2`.

**2-D interior stencil** — thirteen-point stencil on a square
:math:`(\Delta x = \Delta y = h)` grid: centre coefficient +20, the
four nearest neighbours −8, the four next-nearest along each axis +1,
and the four diagonal neighbours +2.

.. figure:: _static/stencil_2d_interior.svg
   :width: 70%
   :align: center
   :alt: 2-D interior finite-difference stencil

   Thirteen-point 2-D stencil for the biharmonic operator on a square
   grid (uniform :math:`D`).
