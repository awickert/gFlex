#! /usr/bin/env python
"""Tests for the 1-D finite-difference and FFT solvers, including end-load (sigma_xx) support."""

import numpy as np
import pytest

from gflex.f1d import F1D


def _run_flex_1d(Te, qs, dx, method="FD", bc_w="0Moment0Shear",
                 bc_e="0Moment0Shear", sigma_xx=None):
    """Helper: run a 1-D flexure calculation and return the flex object."""
    flex = F1D()
    flex.Quiet = True
    flex.Method = method
    flex.Solver = "direct"
    flex.g = 9.8
    flex.E = 65e9
    flex.nu = 0.25
    flex.rho_m = 3300.0
    flex.rho_fill = 0.0
    flex.Te = Te
    flex.qs = qs.copy()
    flex.dx = dx
    flex.BC_W = bc_w
    flex.BC_E = bc_e
    if sigma_xx is not None:
        flex.sigma_xx = sigma_xx
    flex.initialize()
    flex.run()
    flex.finalize()
    return flex


# ---------------------------------------------------------------------------
# Basic smoke test
# ---------------------------------------------------------------------------

def test_1d_fd_runs():
    """1-D FD solver produces a non-trivial, finite deflection."""
    N = 200
    dx = 4000.0
    qs = np.zeros(N)
    qs[90:110] = 1e6
    flex = _run_flex_1d(Te=30e3, qs=qs, dx=dx)
    assert flex.w.ndim == 1
    assert flex.w.size == N
    assert not np.any(np.isnan(flex.w))
    assert flex.w.min() < 0


# ---------------------------------------------------------------------------
# End-load (sigma_xx) correctness tests
# ---------------------------------------------------------------------------

def test_1d_sigma_xx_zero_matches_default():
    """sigma_xx=0 gives identical results to not setting sigma_xx."""
    N = 150
    dx = 4000.0
    qs = np.zeros(N)
    qs[60:90] = 1e6
    flex_default = _run_flex_1d(Te=30e3, qs=qs, dx=dx,
                                bc_w="Periodic", bc_e="Periodic")
    flex_zero = _run_flex_1d(Te=30e3, qs=qs, dx=dx,
                             bc_w="Periodic", bc_e="Periodic",
                             sigma_xx=0.0)
    np.testing.assert_array_equal(flex_default.w, flex_zero.w)


def test_1d_end_load_fft_vs_fd():
    """FFT and FD solutions agree for sigma_xx with periodic BCs.

    The FFT method applies the exact spectral transfer function
    W(k) = -Q(k) / (D k⁴ + σ_xx T_e k² + Δρg), making it the analytical
    reference for the FD end-load implementation.  Agreement to within FD
    truncation error (~rtol 1e-4 for this grid) confirms sign, scaling, and
    placement of the σ_xx term in the FD stencil.
    """
    dx = 4000.0
    N = 256
    L = N * dx
    Te = 30e3
    sigma_xx = 2e8  # tensile [Pa]

    n_waves = 3
    k = 2.0 * np.pi * n_waves / L
    x = (np.arange(N) + 0.5) * dx
    qs = np.sin(k * x) * 1e6

    flex_fft = _run_flex_1d(Te=Te, qs=qs, dx=dx, method="FFT",
                             bc_w="Periodic", bc_e="Periodic",
                             sigma_xx=sigma_xx)
    flex_fd = _run_flex_1d(Te=Te, qs=qs, dx=dx, method="FD",
                            bc_w="Periodic", bc_e="Periodic",
                            sigma_xx=sigma_xx)

    # rtol=5e-4: FD discretisation error in k⁴ is O((kdx)²) ≈ 0.05 % for
    # this grid; FFT is the exact spectral answer.
    np.testing.assert_allclose(flex_fd.w, flex_fft.w, rtol=5e-4)


def test_1d_end_load_monotonicity():
    """Tensile sigma_xx reduces deflection; compressive increases it."""
    N = 200
    dx = 4000.0
    qs = np.zeros(N)
    qs[90:110] = 1e6

    flex_0 = _run_flex_1d(Te=30e3, qs=qs, dx=dx)
    flex_t = _run_flex_1d(Te=30e3, qs=qs, dx=dx, sigma_xx=+1e8)   # tensile
    flex_c = _run_flex_1d(Te=30e3, qs=qs, dx=dx, sigma_xx=-1e8)   # compressive

    assert flex_t.w.min() > flex_0.w.min(), "tensile end load should reduce subsidence"
    assert flex_c.w.min() < flex_0.w.min(), "compressive end load should increase subsidence"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
