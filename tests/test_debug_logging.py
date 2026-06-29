"""Regression tests: debug logging must not crash with a scalar elastic thickness.

`fd_solve()` logs the Te shape under `debug=True`.  Because logging evaluates
its arguments before checking the level, `self.T_e.shape` raised
``AttributeError: 'float' object has no attribute 'shape'`` for a scalar Te
regardless of the configured log level.  Both solvers now use
``np.shape(self.T_e)``, which returns ``()`` for a scalar.
"""

import numpy as np

from gflex.f1d import F1D
from gflex.f2d import F2D

E, NU, TE, RHO_M, RHO_F, G = 65.0e9, 0.25, 30.0e3, 3300.0, 0.0, 9.8
DX = 30.0e3


def _common(flex):
    flex.quiet = True
    flex.verbose = False
    flex.debug = True          # exercises the Te-shape debug log in fd_solve
    flex.method = "fd"
    flex.g, flex.E, flex.nu = G, E, NU
    flex.rho_m, flex.rho_fill = RHO_M, RHO_F
    flex.T_e = TE              # scalar: the case that previously crashed


def test_1d_debug_scalar_te_no_crash():
    flex = F1D()
    _common(flex)
    flex.dx = DX
    flex.qs = np.zeros(40)
    flex.qs[20] = 1.0e8
    flex.bc_west = flex.bc_east = "zero_displacement_zero_slope"
    flex.initialize()
    flex.run()
    assert np.all(np.isfinite(flex.w))


def test_2d_debug_scalar_te_no_crash():
    flex = F2D()
    _common(flex)
    flex.dx = flex.dy = DX
    flex.qs = np.zeros((40, 40))
    flex.qs[20, 20] = 1.0e8
    flex.bc_west = flex.bc_east = "zero_displacement_zero_slope"
    flex.bc_north = flex.bc_south = "zero_displacement_zero_slope"
    flex.initialize()
    flex.run()
    assert np.all(np.isfinite(flex.w))
