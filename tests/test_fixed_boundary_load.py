"""A load on a fixed-displacement boundary node is reacted by the support.

For displacement-fixing BCs (clamped ``zero_displacement_zero_slope``, pinned
``zero_displacement_zero_moment``, and dict ``{"displacement": ...}``), the
boundary node's deflection is prescribed.  A load applied there must be taken
up by the support, not deflect the node, so ``w`` keeps its prescribed value.
Free / mirror edges do not fix displacement and are unaffected.
"""

import warnings

import numpy as np
import pytest

from gflex.f1d import F1D
from gflex.f2d import F2D

E, NU, RHO_M, RHO_F, G = 65.0e9, 0.25, 3300.0, 0.0, 9.8
QS_BIG = 5.0e8   # load placed directly on the boundary node


def _run_1d(bc_west, te=30.0e3, qs0=QS_BIG, record=False):
    f = F1D()
    f.quiet = True
    f.verbose = False
    f.debug = False
    f.method = "fd"
    f.g, f.E, f.nu, f.rho_m, f.rho_fill = G, E, NU, RHO_M, RHO_F
    f.dx = 30.0e3
    f.T_e = te
    qs = np.zeros(40)
    qs[0] = qs0
    f.qs = qs
    f.bc_west = bc_west
    f.bc_east = "zero_displacement_zero_slope"
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        f.initialize()
        f.run()
    w0 = f.w[0]
    if record:
        fired = any("fixed-displacement" in str(c.message) for c in caught)
        return w0, fired
    return w0


# ---------------------------------------------------------------------------
# Fixed-displacement boundaries hold their value under a boundary load
# ---------------------------------------------------------------------------

def test_1d_clamped_holds_zero_under_boundary_load():
    assert abs(_run_1d("zero_displacement_zero_slope")) < 1e-9


def test_1d_pinned_holds_zero_under_boundary_load():
    assert abs(_run_1d("zero_displacement_zero_moment")) < 1e-9


def test_1d_dict_displacement_holds_value_under_boundary_load():
    w0 = _run_1d({"displacement": 100.0, "slope": 0.0})
    assert abs(w0 - 100.0) < 1e-9


def test_2d_clamped_holds_zero_under_edge_load():
    n = 31
    f = F2D()
    f.quiet = True
    f.verbose = False
    f.debug = False
    f.method = "fd"
    f.solver = "direct"
    f.g, f.E, f.nu, f.rho_m, f.rho_fill = G, E, NU, RHO_M, RHO_F
    f.dx = f.dy = 30.0e3
    f.T_e = 30.0e3
    qs = np.zeros((n, n))
    qs[15, 0] = QS_BIG          # load on a west-edge node
    f.qs = qs
    f.bc_west = f.bc_east = f.bc_north = f.bc_south = "zero_displacement_zero_slope"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        f.initialize()
        f.run()
    assert abs(f.w[15, 0]) < 1e-9


# ---------------------------------------------------------------------------
# Non-displacement edges are unaffected
# ---------------------------------------------------------------------------

def test_1d_free_end_boundary_load_not_cancelled():
    """A free end does not fix displacement: a boundary load deflects it."""
    w0, fired = _run_1d("zero_moment_zero_shear", record=True)
    assert abs(w0) > 1.0          # genuinely deflects
    assert not fired              # no fixed-node warning


# ---------------------------------------------------------------------------
# Warning behaviour
# ---------------------------------------------------------------------------

def test_warning_fires_on_load_at_fixed_boundary():
    _, fired = _run_1d("zero_displacement_zero_slope", record=True)
    assert fired


def test_no_warning_when_no_load_on_fixed_boundary():
    _, fired = _run_1d("zero_displacement_zero_slope", qs0=0.0, record=True)
    assert not fired
