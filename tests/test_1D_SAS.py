#! /usr/bin/env python
"""Tests for the 1-D SAS (Superposition of Analytical Solutions) solver.

Green's-function validation
----------------------------
A single concentrated load P = qs · dx reduces SAS to one Green's-function
evaluation.  The exact 1-D infinite-plate deflection is::

    w(r) = −P · (α³/8D) · exp(−r/α) · (cos(r/α) + sin(r/α))

where α = (4D / Δρg)^0.25 is the 1-D flexural parameter and
r = |x − x₀|.  Because SAS implements this formula directly, the
comparison is exact to floating-point precision and validates α, the
coefficient α³/8D, and the exponential-cosine-sine formula together.

Interior agreement with FD
---------------------------
On a large domain with a block load at the centre, SAS (infinite plate)
and FD (finite plate, 0Moment0Shear BC) agree in the interior once the
load sits many flexural wavelengths from the boundary.  BC corrections
decay exponentially with distance from the edge.
"""

import numpy as np
import pytest

from gflex.f1d import F1D


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

D    = E * Te**3 / (12.0 * (1.0 - nu**2))
drho = rho_m - rho_fill

# 1-D flexural parameter and Green's-function coefficient
alpha = (4.0 * D / (drho * g)) ** 0.25   # ≈ 66.3 km
coeff = alpha**3 / (8.0 * D)


def _run(qs, method="SAS", bc_w="0Moment0Shear", bc_e="0Moment0Shear"):
    """Run a 1-D flexure calculation and return the flex object."""
    flex = F1D()
    flex.Quiet  = True
    flex.Method = method
    flex.Solver = "direct"
    flex.g       = g
    flex.E       = E
    flex.nu      = nu
    flex.rho_m   = rho_m
    flex.rho_fill = rho_fill
    flex.Te      = Te
    flex.qs      = qs.copy()
    flex.dx      = dx
    flex.BC_W    = bc_w
    flex.BC_E    = bc_e
    flex.initialize()
    flex.run()
    flex.finalize()
    return flex


# ---------------------------------------------------------------------------
# Analytical Green's-function validation
# ---------------------------------------------------------------------------

def test_1d_sas_green_function_exact():
    """1-D SAS matches the analytical exponential Green's function exactly.

    For a single concentrated load P at x₀, the exact deflection is:

        w(r) = −P · (α³/8D) · exp(−r/α) · (cos(r/α) + sin(r/α))

    A single load cell means SAS performs exactly one Green's-function
    evaluation with no superposition approximation.  Agreement to
    rtol = 1e-10 confirms α, the coefficient α³/8D, and the
    exponential-cosine-sine formula are all correct.
    """
    N  = 200
    i0 = N // 2
    qs = np.zeros(N)
    qs[i0] = 1e6
    P = qs[i0] * dx   # total load [N/m]

    flex = _run(qs)

    x = np.arange(N) * dx
    r = np.abs(x - x[i0])
    w_exact = (
        -P * coeff
        * np.exp(-r / alpha)
        * (np.cos(r / alpha) + np.sin(r / alpha))
    )

    np.testing.assert_allclose(flex.w, w_exact, rtol=1e-10)


# ---------------------------------------------------------------------------
# Interior agreement with FD
# ---------------------------------------------------------------------------

def test_1d_sas_vs_fd_large_domain():
    """1-D SAS and FD agree to within 10 mm absolute in the interior.

    A 20-cell (80 km) block load at the centre of a 400-cell domain sits
    ~12 α from each boundary (α ≈ 66 km ≈ 16.6 cells).  The 0Moment0Shear
    BC correction at the 100-cell margin decays as exp(−6 α) ≈ 0.24 %.
    FD truncation error for the block-load modes gives a peak absolute error
    of ~6.6 mm — about 0.09 % of the ~7.5 m peak deflection — but the
    relative error blows up near the forebulge zero-crossing, so atol
    is used instead of rtol, following the same reasoning as
    test_fft_padded_end_load_vs_fd in test_1D_FFT.py.
    """
    N  = 400
    qs = np.zeros(N)
    qs[190:210] = 1e6

    flex_sas = _run(qs)
    flex_fd  = _run(qs, method="FD")

    m = 100
    np.testing.assert_allclose(
        flex_sas.w[m : N - m],
        flex_fd.w[m : N - m],
        atol=0.01,
    )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
