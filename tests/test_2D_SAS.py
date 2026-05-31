#! /usr/bin/env python
"""Tests for the 2-D SAS (Superposition of Analytical Solutions) solver.

Green's-function validation
----------------------------
A single concentrated load P = qs · dx · dy reduces SAS to one
Kelvin-function evaluation.  The exact 2-D infinite-plate deflection is::

    w(r) = P · (α²/2πD) · kei(r/α)

where α = (D / Δρg)^0.25 is the 2-D flexural parameter and kei is the
Kelvin function of the second kind of zero order.  kei(r/α) < 0 for r ≥ 0,
so a positive (downward) load produces negative (downward) deflection.
Because SAS implements this formula directly, the comparison is exact to
floating-point precision and validates α, the coefficient α²/2πD, and
the kei normalisation together.

Interior agreement with FD
---------------------------
On a large domain with a block load at the centre, SAS (infinite plate)
and FD (finite plate, 0Moment0Shear BC) agree in the interior once the
load sits many flexural wavelengths from the boundary.  BC corrections
decay exponentially with distance from the edge.
"""

import numpy as np
import pytest
from scipy.special import kei

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

# 2-D flexural parameter and Kelvin-function Green's-function coefficient
alpha = (D / (drho * g)) ** 0.25          # ≈ 46.9 km
coeff = alpha**2 / (2.0 * np.pi * D)      # α²/2πD


def _run(qs, method="SAS", bc_w="0Moment0Shear", bc_e="0Moment0Shear",
         bc_n="0Moment0Shear", bc_s="0Moment0Shear"):
    """Run a 2-D flexure calculation and return the flex object."""
    flex = F2D()
    flex.quiet   = True
    flex.method  = method
    flex.solver  = "direct"
    flex.g        = g
    flex.E        = E
    flex.nu       = nu
    flex.rho_m    = rho_m
    flex.rho_fill = rho_fill
    flex.te       = Te
    flex.qs       = qs.copy()
    flex.dx       = dx
    flex.dy       = dy
    flex.bc_west     = bc_w
    flex.bc_east     = bc_e
    flex.bc_north     = bc_n
    flex.bc_south     = bc_s
    flex.initialize()
    flex.run()
    flex.finalize()
    return flex


# ---------------------------------------------------------------------------
# Analytical Kelvin-function validation
# ---------------------------------------------------------------------------

def test_2d_sas_kelvin_green_function_exact():
    """2-D SAS matches the analytical Kelvin-function Green's function exactly.

    For a single concentrated load P at (x₀, y₀), the exact deflection is:

        w(r) = P · (α²/2πD) · kei(r/α)

    where r = √((x−x₀)²+(y−y₀)²).  A single load cell means SAS performs
    exactly one kei evaluation with no superposition approximation.  Agreement
    to rtol = 1e-10 confirms α, the coefficient α²/2πD, and the kei
    normalisation are all correct.  kei(0) = −π/4 is well-defined, so the
    load cell itself is included in the comparison.
    """
    N        = 100
    j0, i0  = N // 2, N // 2
    qs       = np.zeros((N, N))
    qs[j0, i0] = 1e6
    P = qs[j0, i0] * dx * dy   # total load [N]

    flex = _run(qs)

    x = np.arange(N) * dx
    y = np.arange(N) * dy
    X, Y = np.meshgrid(x, y)
    r = np.sqrt((X - x[i0])**2 + (Y - y[j0])**2)
    w_exact = P * coeff * kei(r / alpha)

    np.testing.assert_allclose(flex.w, w_exact, rtol=1e-10)


# ---------------------------------------------------------------------------
# Interior agreement with FD
# ---------------------------------------------------------------------------

def test_2d_sas_vs_fd_large_domain():
    """2-D SAS and FD agree to within 0.5 % in the interior of a large domain.

    A 10×10-cell (40 km × 40 km) block load at the centre of a 200×200 grid
    sits ~8.5 α from each boundary (α ≈ 46.9 km ≈ 11.7 cells).  The
    0Moment0Shear BC correction at the 70-cell margin decays as exp(−6 α)
    ≈ 0.25 %, and FD truncation for the dominant block-load modes adds
    another ~0.1 %, so rtol = 5e-3 is conservative.
    """
    N  = 200
    qs = np.zeros((N, N))
    qs[95:105, 95:105] = 1e6

    flex_sas = _run(qs)
    flex_fd  = _run(qs, method="FD")

    m = 70
    np.testing.assert_allclose(
        flex_sas.w[m : N - m, m : N - m],
        flex_fd.w[m : N - m, m : N - m],
        rtol=5e-3,
    )


# ---------------------------------------------------------------------------
# SAS_NG: non-gridded (arbitrary point-load) solver
# ---------------------------------------------------------------------------

def _run_sas_ng(x, y, q, xw, yw):
    """Run a 2-D SAS_NG (non-gridded) flexure calculation."""
    flex = F2D()
    flex.quiet   = True
    flex.method  = "SAS_NG"
    flex.g        = g
    flex.E        = E
    flex.nu       = nu
    flex.rho_m    = rho_m
    flex.rho_fill = rho_fill
    flex.te       = Te
    flex.x        = x.copy()
    flex.y        = y.copy()
    flex.q        = q.copy()
    flex.xw       = xw.copy()
    flex.yw       = yw.copy()
    flex.initialize()
    flex.run()
    flex.finalize()
    return flex


def test_2d_sas_ng_green_function_exact():
    """2-D SAS_NG matches the analytical Kelvin-function Green's function exactly.

    SAS_NG takes total load q [N] directly (no dx·dy multiplication), so a
    single point load P at (x₀, y₀) gives:

        w(r) = P · (α²/2πD) · kei(r/α)

    Agreement to rtol = 1e-10 confirms α, the coefficient α²/2πD, and the
    kei normalisation are all correct for the non-gridded 2-D code path.
    """
    x_load = np.array([50.0 * dx])   # single load at (200 km, 200 km)
    y_load = np.array([50.0 * dy])
    P      = 1e6 * dx * dy            # total load [N]
    q      = np.array([P])

    N  = 100
    xw = np.arange(N, dtype=float) * dx
    yw = np.arange(N, dtype=float) * dy
    Xw, Yw = np.meshgrid(xw, yw)
    xw_flat = Xw.ravel()
    yw_flat = Yw.ravel()

    flex = _run_sas_ng(x_load, y_load, q, xw_flat, yw_flat)

    r       = np.sqrt((xw_flat - x_load[0])**2 + (yw_flat - y_load[0])**2)
    w_exact = P * coeff * kei(r / alpha)
    np.testing.assert_allclose(flex.w, w_exact, rtol=1e-10)


def test_2d_sas_ng_equals_sas():
    """2-D SAS_NG equals gridded SAS when point loads are placed at cell centres.

    SAS_NG with q[k] = qs[j,i]·dx·dy placed at grid-cell positions (i·dx, j·dy)
    performs the same superposition as gridded SAS (which multiplies qs by dx·dy
    internally), so the two outputs must be identical to floating-point precision.
    """
    N  = 60
    qs = np.zeros((N, N))
    qs[25:35, 25:35] = 1e6   # block load

    flex_sas = _run(qs)   # gridded SAS

    i_grid, j_grid = np.meshgrid(np.arange(N), np.arange(N))
    x_centres = (i_grid * dx).ravel()
    y_centres = (j_grid * dy).ravel()
    mask      = qs.ravel() != 0
    x_load    = x_centres[mask]
    y_load    = y_centres[mask]
    q_load    = qs.ravel()[mask] * dx * dy   # Pa → N

    flex_ng = _run_sas_ng(x_load, y_load, q_load, x_centres, y_centres)

    np.testing.assert_allclose(
        flex_ng.w.reshape(N, N), flex_sas.w, rtol=1e-10
    )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
