#! /usr/bin/env python
"""Tests for the 1-D FFT solver, including end-load (sigma_xx) support."""

import numpy as np
import pytest

from gflex.f1d import F1D


# ---------------------------------------------------------------------------
# Shared parameters
# ---------------------------------------------------------------------------

E    = 65e9
nu   = 0.25
Te   = 30e3
rho_m    = 3300.0
rho_fill = 0.0
g    = 9.8
dx   = 4000.0

D    = E * Te**3 / (12.0 * (1.0 - nu**2))
drho = rho_m - rho_fill


def _run_flex_1d(qs, method="fft", bc_w="", bc_e="", sigma_xx=None):
    """Run a 1-D flexure calculation with the shared parameter set."""
    flex = F1D()
    flex.quiet = True
    flex.method = method
    flex.solver = "direct"
    flex.g = g
    flex.E = E
    flex.nu = nu
    flex.rho_m = rho_m
    flex.rho_fill = rho_fill
    flex.T_e = Te
    flex.qs = qs.copy()
    flex.dx = dx
    flex.bc_west = bc_w
    flex.bc_east = bc_e
    if sigma_xx is not None:
        flex.sigma_xx = sigma_xx
    flex.initialize()
    flex.run()
    return flex


# ---------------------------------------------------------------------------
# periodic BC: exact match against analytical transfer function
# ---------------------------------------------------------------------------

def test_fft_periodic_no_end_load_exact():
    """FFT (periodic) without sigma_xx matches the exact spectral formula."""
    N = 256
    L = N * dx
    n_waves = 3
    k = 2.0 * np.pi * n_waves / L
    x = (np.arange(N) + 0.5) * dx
    qs = np.sin(k * x) * 1e6

    flex = _run_flex_1d(qs, bc_w="periodic", bc_e="periodic")

    w_exact = -1e6 / (D * k**4 + drho * g) * np.sin(k * x)
    np.testing.assert_allclose(flex.w, w_exact, rtol=1e-10)


def test_fft_periodic_end_load_exact():
    """FFT (periodic) with sigma_xx matches the exact spectral formula."""
    N = 256
    L = N * dx
    n_waves = 3
    k = 2.0 * np.pi * n_waves / L
    x = (np.arange(N) + 0.5) * dx
    qs = np.sin(k * x) * 1e6
    sigma_xx = 2e8  # tensile [Pa]

    flex = _run_flex_1d(qs, bc_w="periodic", bc_e="periodic", sigma_xx=sigma_xx)

    w_exact = -1e6 / (D * k**4 + sigma_xx * Te * k**2 + drho * g) * np.sin(k * x)
    np.testing.assert_allclose(flex.w, w_exact, rtol=1e-10)


# ---------------------------------------------------------------------------
# Zero-padded mode: internal consistency
# ---------------------------------------------------------------------------

def test_fft_padded_matches_manual_padding():
    """Auto-padded FFT (non-periodic BCs) matches FFT with periodic on a manually padded domain.

    The FFT solver pads internally when BCs are not periodic.  This test
    verifies that internal padding produces the same result as manually
    replicating the padding and running with periodic BCs — i.e., the two
    paths through the code are equivalent.
    """
    N = 100
    qs = np.zeros(N)
    qs[40:60] = 1e6
    sigma_xx = 1e8

    # Internal auto-padding (blank BCs → zero-padded)
    flex_auto = _run_flex_1d(qs, sigma_xx=sigma_xx)

    # Manual padding to match what the code does internally
    alpha = (4.0 * D / (drho * g)) ** 0.25
    pad = int(np.ceil(4.0 * alpha / dx))
    qs_padded = np.pad(qs, pad, mode="constant")

    flex_manual = _run_flex_1d(qs_padded, bc_w="periodic", bc_e="periodic",
                                sigma_xx=sigma_xx)

    np.testing.assert_allclose(flex_auto.w, flex_manual.w[pad : pad + N],
                                rtol=1e-10)


# ---------------------------------------------------------------------------
# Zero-padded mode: physics check against FD
# ---------------------------------------------------------------------------

def test_fft_padded_end_load_vs_fd():
    """Zero-padded FFT with sigma_xx agrees with FD/periodic on the same padded domain.

    The auto-padded FFT and FD with periodic BCs on the manually padded
    domain solve the same discrete problem, so they must agree to within FD
    truncation error (~rtol 5e-4).  This confirms sigma_xx is applied
    correctly in the zero-padded (non-periodic) FFT path.

    Note: comparing to FD with physical BCs such as zero_moment_zero_shear on a
    finite domain is not appropriate here — those BCs impose a plate end,
    which is a genuinely different physical assumption from the infinite
    plate (no_outside_loads) that the zero-padded FFT represents.
    """
    N = 200
    qs = np.zeros(N)
    qs[90:110] = 1e6
    sigma_xx = 1e8

    # FFT: non-periodic BCs trigger internal zero-padding by 4α
    flex_fft = _run_flex_1d(qs, sigma_xx=sigma_xx)

    # FD/periodic on the same manually padded domain — equivalent problem
    alpha = (4.0 * D / (drho * g)) ** 0.25
    pad = int(np.ceil(4.0 * alpha / dx))
    qs_large = np.pad(qs, pad, mode="constant")

    flex_fd = _run_flex_1d(qs_large, method="fd",
                            bc_w="periodic", bc_e="periodic",
                            sigma_xx=sigma_xx)

    # atol=0.01 m: the rectangular load excites many wavenumbers, each with
    # O((k·dx)²) FD truncation error.  Their sum gives ~6 mm peak absolute
    # error — tiny relative to the ~15 m peak deflection (~0.04%) but large
    # relative to the near-zero forebulge (~0.02 m), so rtol alone fails.
    np.testing.assert_allclose(flex_fft.w, flex_fd.w[pad : pad + N], atol=0.01)


# ---------------------------------------------------------------------------
# Monotonicity: sigma_xx sign effect
# ---------------------------------------------------------------------------

def test_fft_end_load_monotonicity():
    """Tensile sigma_xx reduces deflection; compressive increases it."""
    N = 200
    qs = np.zeros(N)
    qs[90:110] = 1e6

    flex_0 = _run_flex_1d(qs)
    flex_t = _run_flex_1d(qs, sigma_xx=+1e8)   # tensile
    flex_c = _run_flex_1d(qs, sigma_xx=-1e8)   # compressive

    assert flex_t.w.min() > flex_0.w.min(), "tensile end load should reduce subsidence"
    assert flex_c.w.min() < flex_0.w.min(), "compressive end load should increase subsidence"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
