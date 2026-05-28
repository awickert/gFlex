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
    flex.Quiet   = True
    flex.Method  = method
    flex.Solver  = "direct"
    flex.g        = g
    flex.E        = E
    flex.nu       = nu
    flex.rho_m    = rho_m
    flex.rho_fill = rho_fill
    flex.Te       = Te
    flex.qs       = qs.copy()
    flex.dx       = dx
    flex.dy       = dy
    flex.BC_W     = bc_w
    flex.BC_E     = bc_e
    flex.BC_N     = bc_n
    flex.BC_S     = bc_s
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


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
