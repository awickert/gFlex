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


def _run(qs, method="fft", bc_w="", bc_e="", bc_n="", bc_s="",
         dx_=dx, dy_=dy, sigma_xx=None, sigma_yy=None, sigma_xy=None):
    flex = F2D()
    flex.quiet = True
    flex.method = method
    flex.g = g
    flex.E = E
    flex.nu = nu
    flex.rho_m = rho_m
    flex.rho_fill = rho_fill
    flex.te = Te
    flex.qs = qs.copy()
    flex.dx = dx_
    flex.dy = dy_
    flex.bc_west = bc_w
    flex.bc_east = bc_e
    flex.bc_north = bc_n
    flex.bc_south = bc_s
    if sigma_xx is not None:
        flex.sigma_xx = sigma_xx
    if sigma_yy is not None:
        flex.sigma_yy = sigma_yy
    if sigma_xy is not None:
        flex.sigma_xy = sigma_xy
    flex.initialize()
    flex.run()
    return flex


# ---------------------------------------------------------------------------
# periodic: exact match against analytical transfer function
# ---------------------------------------------------------------------------

def test_fft_2d_periodic_exact():
    """2-D FFT (periodic) matches the exact spectral formula for a cosine load.

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

    flex = _run(qs, bc_w="periodic", bc_e="periodic", bc_n="periodic", bc_s="periodic")

    w_exact = -q0 / (D * (kx**2 + ky**2)**2 + drho * g) * np.cos(kx * X) * np.cos(ky * Y)
    np.testing.assert_allclose(flex.w, w_exact, rtol=1e-10)


# ---------------------------------------------------------------------------
# periodic with non-zero rho_fill: exact match
# ---------------------------------------------------------------------------

def test_fft_2d_periodic_rho_fill_exact():
    """2-D FFT (periodic) with rho_fill ≠ 0 matches the exact spectral formula.

    The Winkler restoring term uses Δρ = rho_m − rho_fill, not bare rho_m.
    With rho_fill = 1030 kg/m³ (seawater) the effective restoring modulus
    is weaker, so deflections are deeper.  Agreement to rtol = 1e-10
    against the formula with the correct Δρ confirms the subtraction is
    applied throughout the solver.
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
    rho_fill_water = 1030.0
    drho_water = rho_m - rho_fill_water
    qs = q0 * np.cos(kx * X) * np.cos(ky * Y)

    flex = F2D()
    flex.quiet = True
    flex.method = "fft"
    flex.g = g
    flex.E = E
    flex.nu = nu
    flex.rho_m = rho_m
    flex.rho_fill = rho_fill_water
    flex.te = Te
    flex.qs = qs.copy()
    flex.dx = dx
    flex.dy = dy
    flex.bc_west = flex.bc_east = flex.bc_north = flex.bc_south = "periodic"
    flex.initialize()
    flex.run()
    w = flex.w
    flex.finalize()

    K2 = kx**2 + ky**2
    w_exact = -q0 / (D * K2**2 + drho_water * g) * np.cos(kx * X) * np.cos(ky * Y)
    np.testing.assert_allclose(w, w_exact, rtol=1e-10)

    # Deflection must be deeper (more negative) than with air fill (rho_fill=0)
    w_air = _run(qs, bc_w="periodic", bc_e="periodic", bc_n="periodic", bc_s="periodic").w
    assert w.min() < w_air.min(), (
        "water fill (weaker Winkler restoring) should produce deeper deflection than air fill"
    )


# ---------------------------------------------------------------------------
# periodic with sigma_xx: exact match
# ---------------------------------------------------------------------------

def test_fft_2d_periodic_sigma_xx_exact():
    """2-D FFT (periodic) with sigma_xx matches the exact spectral formula."""
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

    flex = _run(qs, bc_w="periodic", bc_e="periodic", bc_n="periodic", bc_s="periodic",
                sigma_xx=sigma_xx)

    K2 = kx**2 + ky**2
    w_exact = -q0 / (D * K2**2 + sigma_xx * Te * kx**2 + drho * g) * np.cos(kx * X) * np.cos(ky * Y)
    np.testing.assert_allclose(flex.w, w_exact, rtol=1e-10)


# ---------------------------------------------------------------------------
# Zero-padded: agrees with SAS (both use no_outside_loads assumption)
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
    flex_sas = _run(qs, method="sas")

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
    flex_sas = _run(qs, method="sas", dx_=4000., dy_=6000.)

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


def test_fft_2d_periodic_sigma_xy_exact():
    """2-D FFT (periodic) with sigma_xy matches the exact spectral formula for both diagonals.

    A separable cos(kx·X)·cos(ky·Y) load has rfft2 peaks at both (kx,+ky)
    and (kx,−ky), which see different denominators when σ_xy ≠ 0 and cannot
    be reduced to a simple closed form.  The two orthogonal diagonal waves
    each sit at a single rfft2 bin and avoid this:

      cos(kx·X + ky·Y): peak at (Kx=+kx, Ky=+ky) → denom gets +2σ_xy·Te·kx·ky
      cos(kx·X − ky·Y): peak at (Kx=+kx, Ky=−ky) → denom gets −2σ_xy·Te·kx·ky

    Testing both diagonals confirms the sign of 2σ_xy·Te·Kx·Ky is directional
    (not symmetric), which a sign error in the stencil or transfer function
    would break for one diagonal but not the other.
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
    sigma_xy = 2e8
    K2 = kx**2 + ky**2

    # Diagonal (+): Kx=+kx, Ky=+ky → denom includes +2σ_xy·Te·kx·ky
    qs_plus = q0 * np.cos(kx * X + ky * Y)
    flex_plus = _run(qs_plus, bc_w="periodic", bc_e="periodic",
                     bc_n="periodic", bc_s="periodic",
                     sigma_xx=0.0, sigma_yy=0.0, sigma_xy=sigma_xy)
    denom_plus  = D * K2**2 + 2.0 * sigma_xy * Te * kx * ky + drho * g
    w_exact_plus = -q0 / denom_plus * np.cos(kx * X + ky * Y)
    np.testing.assert_allclose(flex_plus.w, w_exact_plus, rtol=1e-10)

    # Diagonal (−): Kx=+kx, Ky=−ky → denom includes −2σ_xy·Te·kx·ky
    qs_minus = q0 * np.cos(kx * X - ky * Y)
    flex_minus = _run(qs_minus, bc_w="periodic", bc_e="periodic",
                      bc_n="periodic", bc_s="periodic",
                      sigma_xx=0.0, sigma_yy=0.0, sigma_xy=sigma_xy)
    denom_minus  = D * K2**2 - 2.0 * sigma_xy * Te * kx * ky + drho * g
    w_exact_minus = -q0 / denom_minus * np.cos(kx * X - ky * Y)
    np.testing.assert_allclose(flex_minus.w, w_exact_minus, rtol=1e-10)

    # The two diagonals must give different deflection amplitudes when σ_xy ≠ 0
    assert not np.allclose(flex_plus.w, flex_minus.w), (
        "orthogonal diagonals should differ when sigma_xy != 0"
    )


# ---------------------------------------------------------------------------
# Zero-padded: internal auto-padding equals explicit manual padding
# ---------------------------------------------------------------------------

def test_fft_2d_padded_matches_manual_padding():
    """2-D auto-padded FFT matches FFT/periodic on the manually padded domain.

    The FFT solver pads by 4α in both x and y when BCs are not periodic
    (where α = (4D/Δρg)^0.25 is the 1-D flexural parameter used for the
    padding estimate).  This test verifies that the internal padding produces
    bit-for-bit the same result as running FFT/periodic on an explicitly
    constructed padded domain — i.e., the two code paths are equivalent.
    """
    N  = 60
    qs = np.zeros((N, N))
    qs[25:35, 25:35] = 1e6

    flex_auto = _run(qs)   # non-periodic BCs → internal zero-padding

    # Replicate the padding formula used inside f2d.FFT
    D_val = E * Te**3 / (12.0 * (1.0 - nu**2))
    alpha = (4.0 * D_val / (drho * g)) ** 0.25   # 1-D alpha for padding estimate
    pad_x = int(np.ceil(4.0 * alpha / dx))
    pad_y = int(np.ceil(4.0 * alpha / dy))
    qs_padded = np.pad(qs, ((pad_y, pad_y), (pad_x, pad_x)), mode="constant")

    flex_manual = _run(qs_padded, bc_w="periodic", bc_e="periodic",
                                  bc_n="periodic", bc_s="periodic")

    np.testing.assert_allclose(
        flex_auto.w,
        flex_manual.w[pad_y : pad_y + N, pad_x : pad_x + N],
        rtol=1e-10,
    )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
