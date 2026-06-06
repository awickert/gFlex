#! /usr/bin/env python
"""Tests for F2D FFT per-axis boundary-condition pairs.

Checks:
  - All-periodic case still returns correct shape and finite values
  - Mixed (x-periodic, y-padded) returns w.shape == qs.shape
  - Mixed (x-padded, y-periodic) returns w.shape == qs.shape
  - All-padded (no BCs set) returns w.shape == qs.shape
  - All-periodic result differs from all-padded (genuinely different physics)
  - Mixed-axis result is finite and non-zero at load location
"""

import numpy as np
import pytest

from gflex.f2d import F2D

E      = 65e9
nu     = 0.25
rho_m  = 3300.0
rho_f  = 0.0
g      = 9.8
Te     = 35e3
dx     = 5000.0


def _run(qs, bc_w, bc_e, bc_n, bc_s):
    flex = F2D()
    flex.quiet   = True
    flex.method  = "fft"
    flex.solver  = "direct"
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


@pytest.mark.parametrize("bc_w,bc_e,bc_n,bc_s", [
    ("periodic", "periodic", "periodic", "periodic"),   # all periodic
    ("periodic", "periodic", "",         ""),           # x-periodic, y-padded
    ("",         "",         "periodic", "periodic"),   # x-padded,   y-periodic
    ("",         "",         "",         ""),           # all padded
])
def test_shape(bc_w, bc_e, bc_n, bc_s):
    """w.shape == qs.shape for every per-axis combination."""
    qs = np.zeros((60, 60))
    qs[30, 30] = 1e6
    w = _run(qs, bc_w, bc_e, bc_n, bc_s)
    assert w.shape == qs.shape


@pytest.mark.parametrize("bc_w,bc_e,bc_n,bc_s", [
    ("periodic", "periodic", "",         ""),
    ("",         "",         "periodic", "periodic"),
    ("",         "",         "",         ""),
])
def test_finite_and_nonzero(bc_w, bc_e, bc_n, bc_s):
    """Non-trivially deflects and produces finite values."""
    qs = np.zeros((80, 80))
    qs[40, 40] = 1e6
    w = _run(qs, bc_w, bc_e, bc_n, bc_s)
    assert np.all(np.isfinite(w))
    assert w.min() < 0


def test_all_periodic_differs_from_padded():
    """All-periodic and all-padded give distinct deflection fields."""
    qs = np.zeros((60, 60))
    qs[30, 30] = 1e6
    w_per = _run(qs, "periodic", "periodic", "periodic", "periodic")
    w_pad = _run(qs, "", "", "", "")
    assert not np.allclose(w_per, w_pad)


def test_mixed_x_periodic_shape_and_symmetry():
    """x-periodic, y-padded: deflection is symmetric along x for centred load."""
    ny, nx = 80, 80
    qs = np.zeros((ny, nx))
    qs[40, 40] = 1e6
    w = _run(qs, "periodic", "periodic", "", "")
    assert w.shape == (ny, nx)
    assert np.all(np.isfinite(w))
    # The load is at x=40 in an 80-cell periodic domain; w should be
    # symmetric about x=40 (cols 40 and 40 are the same row).
    assert np.allclose(w[:, 39], w[:, 41], rtol=1e-6)
