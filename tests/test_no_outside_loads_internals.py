#! /usr/bin/env python
"""Internal-state and coupling-loop tests for the 'no_outside_loads' FD feature.

Checks:
  - self.nx is restored to the inner-domain size after a 1-D FD run
  - self.nx / self.ny are restored after a 2-D FD run
  - A second run on the same object (coupling loop) gives a consistent result
  - _auto_pad_Te_2d matches smooth_pad_Te for the symmetric all-sides case
  - 'no_outside_loads' BC with method='fft' uses zero-padding (same as unset)
"""

import numpy as np
import pytest

from gflex.f1d import F1D
from gflex.f2d import F2D, _auto_pad_Te_2d, smooth_pad_Te

E      = 65e9
nu     = 0.25
rho_m  = 3300.0
rho_f  = 0.0
g      = 9.8
Te     = 35e3
dx     = 5000.0


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_1d(qs, bc_w, bc_e, Te_val=None, finalize=True):
    flex = F1D()
    flex.quiet  = True
    flex.method = "fd"
    flex.solver = "direct"
    flex.g = g; flex.E = E; flex.nu = nu
    flex.rho_m = rho_m; flex.rho_fill = rho_f
    flex.T_e = Te_val if Te_val is not None else Te
    flex.qs  = qs.copy()
    flex.dx  = dx
    flex.bc_west = bc_w
    flex.bc_east = bc_e
    flex.initialize()
    flex.run()
    if finalize:
        w = flex.w.copy()
        flex.finalize()
        return w
    return flex


def _make_2d(qs, bc_w, bc_e, bc_n, bc_s, Te_val=None, finalize=True):
    flex = F2D()
    flex.quiet  = True
    flex.method = "fd"
    flex.solver = "direct"
    flex.g = g; flex.E = E; flex.nu = nu
    flex.rho_m = rho_m; flex.rho_fill = rho_f
    flex.T_e = Te_val if Te_val is not None else Te
    flex.qs  = qs.copy()
    flex.dx  = flex.dy = dx
    flex.bc_west  = bc_w; flex.bc_east  = bc_e
    flex.bc_north = bc_n; flex.bc_south = bc_s
    flex.initialize()
    flex.run()
    if finalize:
        w = flex.w.copy()
        flex.finalize()
        return w
    return flex


# ---------------------------------------------------------------------------
# nx / ny restoration after auto-padding
# ---------------------------------------------------------------------------

def test_1d_nx_matches_w_after_run():
    """After a no_outside_loads FD run, self.nx must equal len(self.w)."""
    qs = np.zeros(100)
    qs[50] = 1e6
    flex = _make_1d(qs, "no_outside_loads", "no_outside_loads", finalize=False)
    assert flex.nx == len(flex.w), (
        f"self.nx={flex.nx} but len(self.w)={len(flex.w)}"
    )
    assert flex.nx == len(qs)
    flex.finalize()


def test_1d_nx_matches_w_west_only():
    """West-only padding: self.nx still equals len(self.w)."""
    qs = np.zeros(80)
    qs[40] = 1e6
    flex = _make_1d(qs, "no_outside_loads", "zero_displacement_zero_slope", finalize=False)
    assert flex.nx == len(flex.w) == len(qs)
    flex.finalize()


def test_2d_nx_ny_match_w_after_run():
    """After a no_outside_loads FD run, self.nx/ny must match self.w.shape."""
    qs = np.zeros((50, 60))
    qs[25, 30] = 1e6
    flex = _make_2d(qs, "no_outside_loads", "no_outside_loads",
                    "no_outside_loads", "no_outside_loads", finalize=False)
    assert flex.w.shape == qs.shape
    assert flex.ny == qs.shape[0], f"self.ny={flex.ny}, expected {qs.shape[0]}"
    assert flex.nx == qs.shape[1], f"self.nx={flex.nx}, expected {qs.shape[1]}"
    flex.finalize()


def test_2d_nx_ny_match_partial_padding():
    """Partial (west+east) no_outside_loads: self.nx/ny still correct."""
    qs = np.zeros((50, 60))
    qs[25, 30] = 1e6
    flex = _make_2d(qs, "no_outside_loads", "no_outside_loads",
                    "mirror", "mirror", finalize=False)
    assert flex.w.shape == qs.shape
    assert flex.ny == qs.shape[0]
    assert flex.nx == qs.shape[1]
    flex.finalize()


# ---------------------------------------------------------------------------
# Coupling-loop: two runs on the same object
# ---------------------------------------------------------------------------

def test_1d_coupling_loop_second_run_matches():
    """Two consecutive 1-D FD runs give identical deflection for identical loads."""
    qs = np.zeros(120)
    qs[60] = 1e6
    flex = _make_1d(qs, "no_outside_loads", "no_outside_loads", finalize=False)
    w1 = flex.w.copy()
    # Second run: change qs slightly then restore; verify shape stays correct
    flex.qs = qs.copy()
    flex.run()
    w2 = flex.w.copy()
    assert w2.shape == w1.shape
    assert np.allclose(w1, w2, rtol=1e-10)
    assert flex.nx == len(qs)
    flex.finalize()


def test_2d_coupling_loop_second_run_matches():
    """Two consecutive 2-D FD runs give identical deflection for identical loads."""
    qs = np.zeros((60, 60))
    qs[30, 30] = 1e6
    flex = _make_2d(qs, "no_outside_loads", "no_outside_loads",
                    "no_outside_loads", "no_outside_loads", finalize=False)
    w1 = flex.w.copy()
    flex.qs = qs.copy()
    flex.run()
    w2 = flex.w.copy()
    assert w2.shape == w1.shape
    assert np.allclose(w1, w2, rtol=1e-10)
    assert flex.ny == qs.shape[0]
    assert flex.nx == qs.shape[1]
    flex.finalize()


# ---------------------------------------------------------------------------
# _auto_pad_Te_2d vs smooth_pad_Te (symmetric case)
# ---------------------------------------------------------------------------

def test_auto_pad_te_2d_matches_smooth_pad_te():
    """_auto_pad_Te_2d with equal sides matches smooth_pad_Te."""
    rng = np.random.default_rng(42)
    Te_inner = 25e3 + 20e3 * rng.random((20, 20))
    p = 8
    ref = smooth_pad_Te(Te_inner, p)
    got = _auto_pad_Te_2d(Te_inner, p, p, p, p)
    assert got.shape == ref.shape
    np.testing.assert_allclose(got, ref, rtol=1e-12)


def test_auto_pad_te_2d_shape_asymmetric():
    """_auto_pad_Te_2d produces the correct shape for asymmetric pads."""
    Te_inner = np.ones((10, 12)) * 30e3
    pn, ps, pw, pe = 3, 5, 0, 4
    out = _auto_pad_Te_2d(Te_inner, pn, ps, pw, pe)
    assert out.shape == (10 + pn + ps, 12 + pw + pe)


def test_auto_pad_te_2d_inner_domain_preserved():
    """_auto_pad_Te_2d leaves the inner domain unchanged."""
    Te_inner = np.arange(1, 101, dtype=float).reshape(10, 10)
    out = _auto_pad_Te_2d(Te_inner, 4, 4, 4, 4)
    np.testing.assert_array_equal(out[4:14, 4:14], Te_inner)


def test_auto_pad_te_2d_outer_corners_near_te_out():
    """Outermost corner cells should be close to Te_out (fully outside)."""
    Te_inner = np.full((8, 8), 40e3)
    Te_out = 30e3
    out = _auto_pad_Te_2d(Te_inner, 6, 6, 6, 6, Te_out=Te_out)
    # Very corners should be Te_out (k=0 in both loops → max(0,0)=0 → Te_out)
    assert abs(out[0, 0] - Te_out) < 1.0
    assert abs(out[0, -1] - Te_out) < 1.0
    assert abs(out[-1, 0] - Te_out) < 1.0
    assert abs(out[-1, -1] - Te_out) < 1.0


# ---------------------------------------------------------------------------
# FFT: 'no_outside_loads' BC is equivalent to unset (both zero-pad)
# ---------------------------------------------------------------------------

def _run_fft_2d(qs, bc_w, bc_e, bc_n, bc_s):
    flex = F2D()
    flex.quiet  = True
    flex.method = "fft"
    flex.solver = "direct"
    flex.g = g; flex.E = E; flex.nu = nu
    flex.rho_m = rho_m; flex.rho_fill = rho_f
    flex.T_e   = Te
    flex.qs    = qs.copy()
    flex.dx    = flex.dy = dx
    flex.bc_west  = bc_w; flex.bc_east  = bc_e
    flex.bc_north = bc_n; flex.bc_south = bc_s
    flex.initialize()
    flex.run()
    w = flex.w.copy()
    flex.finalize()
    return w


def test_fft_2d_no_outside_loads_same_as_unset():
    """'no_outside_loads' FFT BC produces the same result as unset BCs."""
    qs = np.zeros((80, 80))
    qs[40, 40] = 1e6
    w_unset = _run_fft_2d(qs, "", "", "", "")
    w_nol   = _run_fft_2d(qs, "no_outside_loads", "no_outside_loads",
                          "no_outside_loads", "no_outside_loads")
    np.testing.assert_allclose(w_nol, w_unset, rtol=1e-12)


def test_fft_2d_no_outside_loads_mixed_no_warn():
    """'no_outside_loads' on all four sides with FFT: no partial-periodic warning."""
    import warnings
    qs = np.zeros((60, 60))
    qs[30, 30] = 1e6
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        _run_fft_2d(qs, "no_outside_loads", "no_outside_loads",
                    "no_outside_loads", "no_outside_loads")
    periodic_warns = [str(w.message) for w in caught if "no_outside_loads" in str(w.message)]
    assert not periodic_warns, f"unexpected warnings: {periodic_warns}"


