"""Read-only T_e and cache-correctness regression tests.

The coefficient-matrix / LU cache is keyed to ``T_e``; an in-place edit of a
``T_e`` array would silently leave the cache stale and the next solve would
return deflections for the *old* Te.  To make that impossible, a ``T_e`` array
is stored read-only, so in-place edits raise.  The solver works on a separate
writeable copy (``_te``), and ``self.T_e`` keeps its shape across a run.
"""

import warnings

import numpy as np
import pytest

from gflex.f1d import F1D
from gflex.f2d import F2D

E, NU, RHO_M, RHO_F, G = 65.0e9, 0.25, 3300.0, 0.0, 9.8


def _make_1d(te):
    f = F1D()
    f.quiet = True
    f.verbose = False
    f.debug = False
    f.method = "fd"
    f.g, f.E, f.nu, f.rho_m, f.rho_fill = G, E, NU, RHO_M, RHO_F
    f.dx = 30.0e3
    f.qs = np.zeros(40)
    f.qs[20] = 1.0e8
    f.bc_west = f.bc_east = "zero_displacement_zero_slope"
    f.T_e = te
    return f


def _make_2d(te):
    f = F2D()
    f.quiet = True
    f.verbose = False
    f.debug = False
    f.method = "fd"
    f.solver = "direct"
    f.g, f.E, f.nu, f.rho_m, f.rho_fill = G, E, NU, RHO_M, RHO_F
    f.dx = f.dy = 30.0e3
    f.qs = np.zeros((31, 31))
    f.qs[15, 15] = 1.0e9
    f.bc_west = f.bc_east = f.bc_north = f.bc_south = "zero_displacement_zero_slope"
    f.T_e = te
    return f


def _run(f):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        f.initialize()
        f.run()
    return f.w.copy()


# ---------------------------------------------------------------------------
# In-place mutation of a T_e array is forbidden (would desync the cache)
# ---------------------------------------------------------------------------

def test_1d_inplace_te_edit_raises():
    f = _make_1d(30.0e3 * np.ones(40))
    with pytest.raises(ValueError):
        f.T_e[5] = 40.0e3


def test_2d_inplace_te_edit_raises():
    f = _make_2d(30.0e3 * np.ones((31, 31)))
    with pytest.raises(ValueError):
        f.T_e[5, 5] = 40.0e3


def test_inplace_scaling_raises():
    f = _make_1d(30.0e3 * np.ones(40))
    with pytest.raises(ValueError):
        f.T_e *= 2.0


@pytest.mark.parametrize("attr", ["sigma_xx", "sigma_yy", "sigma_xy"])
def test_inplace_sigma_edit_raises(attr):
    """In-plane stress arrays are matrix-determining; in-place edits must raise."""
    f = _make_1d(30.0e3 * np.ones(40))
    setattr(f, attr, 1.0e6 * np.ones(40))
    with pytest.raises(ValueError):
        getattr(f, attr)[0] = 2.0e6


# ---------------------------------------------------------------------------
# Reassignment changes Te correctly (cache invalidated; not stale)
# ---------------------------------------------------------------------------

def test_1d_reassignment_invalidates_cache():
    f = _make_1d(20.0e3 * np.ones(40))
    w_a = _run(f)
    f.T_e = 60.0e3 * np.ones(40)        # reassign -> setter invalidates cache
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        f.run()
    w_b = f.w.copy()
    assert not np.allclose(w_a, w_b), "reassigning Te did not change the solution"
    w_fresh = _run(_make_1d(60.0e3 * np.ones(40)))
    np.testing.assert_allclose(w_b, w_fresh, rtol=1e-12, atol=0)


# ---------------------------------------------------------------------------
# self.T_e stays the pristine user grid across a run (no padding leak)
# ---------------------------------------------------------------------------

def test_2d_te_shape_stable_across_run():
    te = 30.0e3 * np.ones((31, 31))
    f = _make_2d(te)
    _run(f)
    assert f.T_e.shape == (31, 31), f"T_e shape changed to {f.T_e.shape} after run"


def test_1d_te_unchanged_across_run():
    te_vals = 25.0e3 * (1.0 + 0.4 * np.linspace(0, 1, 40))
    f = _make_1d(te_vals)
    _run(f)
    np.testing.assert_array_equal(f.T_e, te_vals)
