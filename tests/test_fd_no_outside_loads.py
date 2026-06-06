#! /usr/bin/env python
"""Tests for FD 'no_outside_loads' auto-padding in F1D and F2D.

Checks:
  - Output w.shape matches input qs.shape (padding is invisible to the caller)
  - FD result agrees with the SAS reference within ~2% for a centred load
  - Partial padding (one side 'no_outside_loads', other side explicit BC) works
  - Mixed BCs (some sides 'no_outside_loads', some 'mirror') produce correct shape
  - Variable-Te run completes without errors
"""

import numpy as np
import pytest

from gflex.f1d import F1D
from gflex.f2d import F2D

# ---------------------------------------------------------------------------
# Shared physical parameters
# ---------------------------------------------------------------------------
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

def _run_1d(qs, bc_w, bc_e, Te_val=None):
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
    w = flex.w.copy()
    flex.finalize()
    return w


def _run_2d(qs, bc_w, bc_e, bc_n, bc_s, Te_val=None):
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
    w = flex.w.copy()
    flex.finalize()
    return w


def _run_1d_sas(qs):
    flex = F1D()
    flex.quiet  = True
    flex.method = "sas"
    flex.solver = "direct"
    flex.g = g; flex.E = E; flex.nu = nu
    flex.rho_m = rho_m; flex.rho_fill = rho_f
    flex.T_e = Te
    flex.qs  = qs.copy()
    flex.dx  = dx
    flex.bc_west = flex.bc_east = "no_outside_loads"
    flex.initialize()
    flex.run()
    w = flex.w.copy()
    flex.finalize()
    return w


# ---------------------------------------------------------------------------
# 1-D: output shape
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("bc_w,bc_e", [
    ("no_outside_loads", "no_outside_loads"),
    ("no_outside_loads", "mirror"),
    ("mirror",           "no_outside_loads"),
])
def test_1d_output_shape(bc_w, bc_e):
    """w.shape must equal qs.shape regardless of which sides are auto-padded."""
    qs = np.zeros(120)
    qs[60] = 1e6
    w = _run_1d(qs, bc_w, bc_e)
    assert w.shape == qs.shape


# ---------------------------------------------------------------------------
# 1-D: accuracy vs SAS reference
# ---------------------------------------------------------------------------

def test_1d_both_sides_vs_sas():
    """FD 'no_outside_loads' on both sides matches SAS within 2% at the peak."""
    qs = np.zeros(300)
    qs[150] = 1e6
    w_fd  = _run_1d(qs, "no_outside_loads", "no_outside_loads")
    w_sas = _run_1d_sas(qs)
    rel_err = abs(w_fd.min() - w_sas.min()) / abs(w_sas.min())
    assert rel_err < 0.02, f"relative peak error {rel_err:.4f} exceeds 2 %"


# ---------------------------------------------------------------------------
# 1-D: partial padding
# ---------------------------------------------------------------------------

def test_1d_west_only_shape_and_finite():
    """West-only 'no_outside_loads' gives correct shape and finite values."""
    qs = np.zeros(100)
    qs[50] = 1e6
    w = _run_1d(qs, "no_outside_loads", "zero_displacement_zero_slope")
    assert w.shape == qs.shape
    assert np.all(np.isfinite(w))


def test_1d_east_only_shape_and_finite():
    qs = np.zeros(100)
    qs[50] = 1e6
    w = _run_1d(qs, "zero_displacement_zero_slope", "no_outside_loads")
    assert w.shape == qs.shape
    assert np.all(np.isfinite(w))


# ---------------------------------------------------------------------------
# 1-D: variable Te
# ---------------------------------------------------------------------------

def test_1d_variable_te_no_outside_loads():
    """Variable Te with 'no_outside_loads' on both sides runs and returns correct shape."""
    nx = 100
    qs = np.zeros(nx)
    qs[50] = 1e6
    Te_arr = np.linspace(20e3, 50e3, nx)
    w = _run_1d(qs, "no_outside_loads", "no_outside_loads", Te_val=Te_arr)
    assert w.shape == qs.shape
    assert np.all(np.isfinite(w))
    assert w.min() < 0  # loads depress the plate


# ---------------------------------------------------------------------------
# 2-D: output shape
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("bc_w,bc_e,bc_n,bc_s", [
    ("no_outside_loads", "no_outside_loads", "no_outside_loads", "no_outside_loads"),
    ("no_outside_loads", "mirror",           "mirror",           "mirror"),
    ("mirror",           "no_outside_loads", "mirror",           "mirror"),
    ("mirror",           "mirror",           "no_outside_loads", "mirror"),
    ("mirror",           "mirror",           "mirror",           "no_outside_loads"),
    ("no_outside_loads", "no_outside_loads", "mirror",           "mirror"),
])
def test_2d_output_shape(bc_w, bc_e, bc_n, bc_s):
    """2-D w.shape must equal qs.shape for any auto-padded combination."""
    qs = np.zeros((60, 60))
    qs[30, 30] = 1e6
    w = _run_2d(qs, bc_w, bc_e, bc_n, bc_s)
    assert w.shape == qs.shape


# ---------------------------------------------------------------------------
# 2-D: accuracy vs SAS reference
# ---------------------------------------------------------------------------

def test_2d_all_sides_vs_sas():
    """FD 'no_outside_loads' on all four sides: peak deflection within 5% of SAS."""
    qs = np.zeros((200, 200))
    qs[100, 100] = 1e6
    w_fd = _run_2d(
        qs,
        "no_outside_loads", "no_outside_loads",
        "no_outside_loads", "no_outside_loads",
    )
    # SAS reference (via F2D sas)
    flex_sas = F2D()
    flex_sas.quiet  = True
    flex_sas.method = "sas"
    flex_sas.g = g; flex_sas.E = E; flex_sas.nu = nu
    flex_sas.rho_m = rho_m; flex_sas.rho_fill = rho_f
    flex_sas.T_e = Te
    flex_sas.qs  = qs.copy()
    flex_sas.dx  = flex_sas.dy = dx
    flex_sas.bc_west = flex_sas.bc_east = "no_outside_loads"
    flex_sas.bc_north = flex_sas.bc_south = "no_outside_loads"
    flex_sas.initialize()
    flex_sas.run()
    w_sas = flex_sas.w.copy()
    flex_sas.finalize()

    rel_err = abs(w_fd.min() - w_sas.min()) / abs(w_sas.min())
    assert rel_err < 0.05, f"relative peak error {rel_err:.4f} exceeds 5 %"


# ---------------------------------------------------------------------------
# 2-D: variable Te
# ---------------------------------------------------------------------------

def test_2d_variable_te_no_outside_loads():
    """Variable 2-D Te with all-sides 'no_outside_loads' runs and returns correct shape."""
    ny, nx = 50, 50
    qs = np.zeros((ny, nx))
    qs[25, 25] = 1e6
    y_idx, x_idx = np.mgrid[0:ny, 0:nx]
    Te_arr = 25e3 + 20e3 * (x_idx / nx)
    w = _run_2d(
        qs,
        "no_outside_loads", "no_outside_loads",
        "no_outside_loads", "no_outside_loads",
        Te_val=Te_arr,
    )
    assert w.shape == qs.shape
    assert np.all(np.isfinite(w))
    assert w.min() < 0
