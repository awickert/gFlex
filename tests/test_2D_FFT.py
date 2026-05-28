#! /usr/bin/env python
"""Tests for the 2-D FFT solver."""

import numpy as np
import pytest

from gflex.f2d import F2D


# ---------------------------------------------------------------------------
# Shared parameters
# ---------------------------------------------------------------------------

E        = 65e9
nu       = 0.25
Te       = 30e3
rho_m    = 3300.0
rho_fill = 0.0
g        = 9.8
dx       = 4000.0
dy       = 4000.0

D    = E * Te**3 / (12.0 * (1.0 - nu**2))
drho = rho_m - rho_fill


def _run(qs, method="FFT", bc_w="", bc_e="", bc_n="", bc_s="",
         dx_=dx, dy_=dy, sigma_xx=None, sigma_yy=None):
    flex = F2D()
    flex.Quiet = True
    flex.Method = method
    flex.g = g
    flex.E = E
    flex.nu = nu
    flex.rho_m = rho_m
    flex.rho_fill = rho_fill
    flex.Te = Te
    flex.qs = qs.copy()
    flex.dx = dx_
    flex.dy = dy_
    flex.BC_W = bc_w
    flex.BC_E = bc_e
    flex.BC_N = bc_n
    flex.BC_S = bc_s
    if sigma_xx is not None:
        flex.sigma_xx = sigma_xx
    if sigma_yy is not None:
        flex.sigma_yy = sigma_yy
    flex.initialize()
    flex.run()
    flex.finalize()
    return flex


# ---------------------------------------------------------------------------
# Periodic: exact match against analytical transfer function
# ---------------------------------------------------------------------------

def test_fft_2d_periodic_exact():
    """2-D FFT (Periodic) matches the exact spectral formula for a cosine load.

    For a separable load q(x,y) = q0·cos(kx·x)·cos(ky·y) the exact solution is
    w(x,y) = -q0 / (D(kx²+ky²)² + Δρg) · cos(kx·x)·cos(ky·y).
    """
    Nx, Ny = 64, 64
    Lx, Ly = Nx * dx, Ny * dy
    nx_waves, ny_waves = 2, 3
    kx = 2.0 * np.pi * nx_waves / Lx
    ky = 2.0 * np.pi * ny_waves / Ly
    x = (np.arange(Nx) + 0.5) * dx
    y = (np.arange(Ny) + 0.5) * dy
    X, Y = np.meshgrid(x, y)
    q0 = 1e6
    qs = q0 * np.cos(kx * X) * np.cos(ky * Y)

    flex = _run(qs, bc_w="Periodic", bc_e="Periodic", bc_n="Periodic", bc_s="Periodic")

    w_exact = -q0 / (D * (kx**2 + ky**2)**2 + drho * g) * np.cos(kx * X) * np.cos(ky * Y)
    np.testing.assert_allclose(flex.w, w_exact, rtol=1e-10)


# ---------------------------------------------------------------------------
# Periodic with sigma_xx: exact match
# ---------------------------------------------------------------------------

def test_fft_2d_periodic_sigma_xx_exact():
    """2-D FFT (Periodic) with sigma_xx matches the exact spectral formula."""
    Nx, Ny = 64, 64
    Lx, Ly = Nx * dx, Ny * dy
    nx_waves, ny_waves = 2, 3
    kx = 2.0 * np.pi * nx_waves / Lx
    ky = 2.0 * np.pi * ny_waves / Ly
    x = (np.arange(Nx) + 0.5) * dx
    y = (np.arange(Ny) + 0.5) * dy
    X, Y = np.meshgrid(x, y)
    q0 = 1e6
    sigma_xx = 2e8
    qs = q0 * np.cos(kx * X) * np.cos(ky * Y)

    flex = _run(qs, bc_w="Periodic", bc_e="Periodic", bc_n="Periodic", bc_s="Periodic",
                sigma_xx=sigma_xx)

    K2 = kx**2 + ky**2
    w_exact = -q0 / (D * K2**2 + sigma_xx * Te * kx**2 + drho * g) * np.cos(kx * X) * np.cos(ky * Y)
    np.testing.assert_allclose(flex.w, w_exact, rtol=1e-10)


# ---------------------------------------------------------------------------
# Zero-padded: agrees with SAS (both use NoOutsideLoads assumption)
# ---------------------------------------------------------------------------

def test_fft_2d_padded_vs_sas():
    """Zero-padded 2-D FFT agrees with SAS for a central square load.

    Both methods assume no loads outside the domain (infinite plate).
    Agreement to rtol=1e-3 confirms the padding is sufficient and the
    transfer function is correct.
    """
    N = 80
    qs = np.zeros((N, N))
    qs[35:45, 35:45] = 1e6

    flex_fft = _run(qs)
    flex_sas = _run(qs, method="SAS")

    # Compare interior, away from any near-boundary differences
    m = 10
    np.testing.assert_allclose(
        flex_fft.w[m:N-m, m:N-m],
        flex_sas.w[m:N-m, m:N-m],
        rtol=1e-3,
    )


# ---------------------------------------------------------------------------
# Non-square grid (dx ≠ dy): zero-padded FFT vs SAS
# ---------------------------------------------------------------------------

def test_fft_2d_nonsquare_grid():
    """2-D FFT with dx ≠ dy agrees with SAS in the interior."""
    Ny, Nx = 60, 80
    qs = np.zeros((Ny, Nx))
    qs[25:35, 35:45] = 1e6

    flex_fft = _run(qs, dx_=4000., dy_=6000.)
    flex_sas = _run(qs, method="SAS", dx_=4000., dy_=6000.)

    m = 8
    np.testing.assert_allclose(
        flex_fft.w[m:Ny-m, m:Nx-m],
        flex_sas.w[m:Ny-m, m:Nx-m],
        rtol=2e-3,
    )


# ---------------------------------------------------------------------------
# Monotonicity: sigma_xx and sigma_yy sign effects
# ---------------------------------------------------------------------------

def test_fft_2d_sigma_xx_monotonicity():
    """Tensile sigma_xx reduces deflection; compressive increases it."""
    N = 60
    qs = np.zeros((N, N))
    qs[25:35, 25:35] = 1e6

    flex_0 = _run(qs)
    flex_t = _run(qs, sigma_xx=+1e8)
    flex_c = _run(qs, sigma_xx=-1e8)

    assert flex_t.w.min() > flex_0.w.min(), "tensile sigma_xx should reduce subsidence"
    assert flex_c.w.min() < flex_0.w.min(), "compressive sigma_xx should increase subsidence"


def test_fft_2d_sigma_yy_monotonicity():
    """Tensile sigma_yy reduces deflection; compressive increases it."""
    N = 60
    qs = np.zeros((N, N))
    qs[25:35, 25:35] = 1e6

    flex_0 = _run(qs)
    flex_t = _run(qs, sigma_yy=+1e8)
    flex_c = _run(qs, sigma_yy=-1e8)

    assert flex_t.w.min() > flex_0.w.min(), "tensile sigma_yy should reduce subsidence"
    assert flex_c.w.min() < flex_0.w.min(), "compressive sigma_yy should increase subsidence"


# ---------------------------------------------------------------------------
# Zero-padded: internal auto-padding equals explicit manual padding
# ---------------------------------------------------------------------------

def test_fft_2d_padded_matches_manual_padding():
    """2-D auto-padded FFT matches FFT/Periodic on the manually padded domain.

    The FFT solver pads by 4α in both x and y when BCs are not Periodic
    (where α = (4D/Δρg)^0.25 is the 1-D flexural parameter used for the
    padding estimate).  This test verifies that the internal padding produces
    bit-for-bit the same result as running FFT/Periodic on an explicitly
    constructed padded domain — i.e., the two code paths are equivalent.
    """
    N  = 60
    qs = np.zeros((N, N))
    qs[25:35, 25:35] = 1e6

    flex_auto = _run(qs)   # non-Periodic BCs → internal zero-padding

    # Replicate the padding formula used inside f2d.FFT
    D_val = E * Te**3 / (12.0 * (1.0 - nu**2))
    alpha = (4.0 * D_val / (drho * g)) ** 0.25   # 1-D alpha for padding estimate
    pad_x = int(np.ceil(4.0 * alpha / dx))
    pad_y = int(np.ceil(4.0 * alpha / dy))
    qs_padded = np.pad(qs, ((pad_y, pad_y), (pad_x, pad_x)), mode="constant")

    flex_manual = _run(qs_padded, bc_w="Periodic", bc_e="Periodic",
                                  bc_n="Periodic", bc_s="Periodic")

    np.testing.assert_allclose(
        flex_auto.w,
        flex_manual.w[pad_y : pad_y + N, pad_x : pad_x + N],
        rtol=1e-10,
    )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
