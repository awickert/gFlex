"""Smoke tests for gFlex output() and plotting().

All tests use a monkeypatched plt.show so no windows open.  Each test closes
all figures in its teardown to prevent figure accumulation.

Valid plot_choice strings are 'q', 'w', 'both', and (1-D only) 'combo'.
"""

import numpy as np
import pytest
from matplotlib import pyplot as plt

from gflex.f1d import F1D
from gflex.f2d import F2D

# ---------------------------------------------------------------------------
# Shared elastic / grid parameters
# ---------------------------------------------------------------------------

_g       = 9.8
_E       = 65e9
_nu      = 0.25
_rho_m   = 3300.0
_rho_f   = 0.0
_Te      = 35e3
_dx = _dy = 5000.0


# ---------------------------------------------------------------------------
# Helper: build and run minimal flex objects
# ---------------------------------------------------------------------------

def _1d_fd():
    N = 60
    qs = np.zeros(N); qs[25:35] = 1e6
    flex = F1D()
    flex.quiet = True;  flex.method = "fd";  flex.solver = "direct"
    flex.g = _g;  flex.E = _E;  flex.nu = _nu
    flex.rho_m = _rho_m;  flex.rho_fill = _rho_f
    flex.te = _Te * np.ones(N);  flex.qs = qs;  flex.dx = _dx
    flex.bc_west = flex.bc_east = "periodic"
    flex.initialize();  flex.run()
    return flex


def _1d_sas():
    N = 60
    qs = np.zeros(N); qs[25:35] = 1e6
    flex = F1D()
    flex.quiet = True;  flex.method = "sas"
    flex.g = _g;  flex.E = _E;  flex.nu = _nu
    flex.rho_m = _rho_m;  flex.rho_fill = _rho_f
    flex.te = _Te;  flex.qs = qs;  flex.dx = _dx
    flex.initialize();  flex.run()
    return flex


def _1d_fft():
    N = 64
    qs = np.zeros(N); qs[28:36] = 1e6
    flex = F1D()
    flex.quiet = True;  flex.method = "fft"
    flex.g = _g;  flex.E = _E;  flex.nu = _nu
    flex.rho_m = _rho_m;  flex.rho_fill = _rho_f
    flex.te = _Te;  flex.qs = qs;  flex.dx = _dx
    flex.bc_west = flex.bc_east = "periodic"
    flex.initialize();  flex.run()
    return flex


def _1d_sas_ng():
    N = 60
    x = np.arange(N) * _dx
    q = np.zeros(N); q[25:35] = 1e6 * _dx   # point forces [N/m]
    flex = F1D()
    flex.quiet = True;  flex.method = "sas_ng"
    flex.g = _g;  flex.E = _E;  flex.nu = _nu
    flex.rho_m = _rho_m;  flex.rho_fill = _rho_f
    flex.te = _Te
    flex.x = x;  flex.q = q;  flex.xw = x.copy()
    flex.initialize();  flex.run()
    return flex


def _2d_fd():
    N = 20
    qs = np.zeros((N, N)); qs[7:13, 7:13] = 1e6
    flex = F2D()
    flex.quiet = True;  flex.method = "fd";  flex.solver = "direct"
    flex.g = _g;  flex.E = _E;  flex.nu = _nu
    flex.rho_m = _rho_m;  flex.rho_fill = _rho_f
    flex.te = _Te * np.ones((N, N));  flex.qs = qs
    flex.dx = flex.dy = _dx
    flex.bc_west = flex.bc_east = flex.bc_north = flex.bc_south = "periodic"
    flex.initialize();  flex.run()
    return flex


def _2d_sas():
    N = 20
    qs = np.zeros((N, N)); qs[7:13, 7:13] = 1e6
    flex = F2D()
    flex.quiet = True;  flex.method = "sas"
    flex.g = _g;  flex.E = _E;  flex.nu = _nu
    flex.rho_m = _rho_m;  flex.rho_fill = _rho_f
    flex.te = _Te;  flex.qs = qs;  flex.dx = flex.dy = _dx
    flex.initialize();  flex.run()
    return flex


def _2d_fft():
    N = 20
    qs = np.zeros((N, N)); qs[7:13, 7:13] = 1e6
    flex = F2D()
    flex.quiet = True;  flex.method = "fft"
    flex.g = _g;  flex.E = _E;  flex.nu = _nu
    flex.rho_m = _rho_m;  flex.rho_fill = _rho_f
    flex.te = _Te;  flex.qs = qs;  flex.dx = flex.dy = _dx   # scalar Te required
    flex.bc_west = flex.bc_east = flex.bc_north = flex.bc_south = "periodic"
    flex.initialize();  flex.run()
    return flex


def _2d_sas_ng():
    N = 10
    x = np.arange(N, dtype=float) * _dx
    y = np.arange(N, dtype=float) * _dy
    X, Y = np.meshgrid(x, y)
    xf, yf = X.ravel(), Y.ravel()
    q = np.zeros(N * N)
    mask = (np.abs(xf - x[N // 2]) < 2 * _dx) & (np.abs(yf - y[N // 2]) < 2 * _dy)
    q[mask] = 1e6 * _dx * _dy    # point forces [N]
    flex = F2D()
    flex.quiet = True;  flex.method = "sas_ng"
    flex.g = _g;  flex.E = _E;  flex.nu = _nu
    flex.rho_m = _rho_m;  flex.rho_fill = _rho_f
    flex.te = _Te
    flex.x = xf;  flex.y = yf;  flex.q = q
    flex.xw = xf.copy();  flex.yw = yf.copy()
    flex.initialize();  flex.run()
    return flex


# ---------------------------------------------------------------------------
# Autouse fixture: suppress plt.show and close figures after every test
# ---------------------------------------------------------------------------

@pytest.fixture(autouse=True)
def _no_show_close(monkeypatch):
    """Prevent plt.show from opening windows; close all figures on teardown."""
    monkeypatch.setattr("matplotlib.pyplot.show", lambda: None)
    yield
    plt.close("all")


# ---------------------------------------------------------------------------
# 1-D plots: q, w, both — all four methods
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("make", [_1d_fd, _1d_sas, _1d_fft, _1d_sas_ng],
                          ids=["FD", "SAS", "FFT", "SAS_NG"])
@pytest.mark.parametrize("choice", ["q", "w", "both"])
def test_1d_plot(make, choice):
    """1-D plotting smoke test: no exception for any valid plot_choice × solver."""
    flex = make()
    flex.plot_choice = choice
    flex.output()


# ---------------------------------------------------------------------------
# 1-D combo plot — all four methods
# (SAS_NG only plots deflection in this mode but should not raise)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("make", [_1d_fd, _1d_sas, _1d_fft, _1d_sas_ng],
                          ids=["FD", "SAS", "FFT", "SAS_NG"])
def test_1d_plot_combo(make):
    """1-D combo plot smoke test."""
    flex = make()
    flex.plot_choice = "combo"
    flex.output()


# ---------------------------------------------------------------------------
# 2-D plots: q, w, both — all four methods
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("make", [_2d_fd, _2d_sas, _2d_fft, _2d_sas_ng],
                          ids=["FD", "SAS", "FFT", "SAS_NG"])
@pytest.mark.parametrize("choice", ["q", "w", "both"])
def test_2d_plot(make, choice):
    """2-D plotting smoke test: no exception for any valid plot_choice × solver."""
    flex = make()
    flex.plot_choice = choice
    flex.output()


# ---------------------------------------------------------------------------
# output() with no plot_choice and no w_out_file — should be a silent no-op
# ---------------------------------------------------------------------------

def test_output_no_plotchoice():
    """output() with no plot_choice set is a silent no-op that creates no figures."""
    flex = _1d_fd()
    flex.output()
    assert plt.get_fignums() == []


# ---------------------------------------------------------------------------
# Invalid plot_choice — prints a message but does not raise
# ---------------------------------------------------------------------------

def test_invalid_plotchoice_silent():
    """An unrecognised plot_choice string does not raise an exception."""
    flex = _1d_fd()
    flex.plot_choice = "not_a_valid_choice"
    flex.output()


# ---------------------------------------------------------------------------
# File output
# ---------------------------------------------------------------------------

def test_output_file_1d_ascii(tmp_path):
    """w_out_file writes a 1-D deflection array as ASCII."""
    flex = _1d_fd()
    out = tmp_path / "deflection.txt"
    flex.w_out_file = str(out)
    flex.output()
    assert out.exists()
    loaded = np.loadtxt(str(out))
    assert loaded.shape == flex.w.shape


def test_output_file_1d_npy(tmp_path):
    """w_out_file ending in .npy writes a binary NumPy array that round-trips."""
    flex = _1d_fd()
    out = tmp_path / "deflection.npy"
    flex.w_out_file = str(out)
    flex.output()
    assert out.exists()
    np.testing.assert_array_equal(np.load(str(out)), flex.w)


def test_output_file_2d_ascii(tmp_path):
    """w_out_file writes a 2-D deflection grid as ASCII."""
    flex = _2d_fd()
    out = tmp_path / "deflection_2d.txt"
    flex.w_out_file = str(out)
    flex.output()
    assert out.exists()
    loaded = np.loadtxt(str(out))
    assert loaded.shape == flex.w.shape
