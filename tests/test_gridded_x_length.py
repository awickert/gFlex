"""The x-coordinate array must have exactly one entry per grid node.

``np.arange(0, N*dx, dx)`` can return N+1 entries for some float ``dx``
(floating-point endpoint rounding), desynchronising ``_x_local`` from ``w``.
``np.arange(N)*dx`` is exact.
"""

import warnings

import numpy as np

from gflex.f1d import F1D

# dx / n known to overrun under the old arange(0, n*dx, dx) idiom.
DX_OVERRUN = 333.3
NX_OVERRUN = 13


def test_gridded_x_length_matches_grid_for_float_dx():
    f = F1D()
    f.dx = DX_OVERRUN
    f.qs = np.zeros(NX_OVERRUN)
    f.gridded_x()
    assert len(f._x_local) == NX_OVERRUN


def test_x_local_matches_w_after_no_outside_loads_run():
    f = F1D()
    f.quiet = True
    f.verbose = False
    f.debug = False
    f.method = "fd"
    f.g, f.E, f.nu, f.rho_m, f.rho_fill = 9.8, 65e9, 0.25, 3300.0, 0.0
    f.dx = DX_OVERRUN
    f.T_e = 1.0e3
    qs = np.zeros(NX_OVERRUN)
    qs[NX_OVERRUN // 2] = 1.0e6
    f.qs = qs
    f.bc_west = f.bc_east = "no_outside_loads"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        f.initialize()
        f.run()
    assert len(f._x_local) == len(f.w) == NX_OVERRUN
