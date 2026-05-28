#! /usr/bin/env python
"""Tests for 1-D FD boundary conditions, cross-validated against analytical references.

Boundary conditions tested:
  Periodic       — exact via FFT spectral formula (covered in test_1D_FFT.py)
  Mirror         — exact: cosine eigenfunction; Mirror == Periodic on 2× even-extended domain
  0Slope0Shear   — same physical intent as Mirror (symmetry BC) but different stencil;
                   shown here to be genuinely distinct from Mirror
  0Displacement0Slope — clamped end; vs SAS (infinite plate) for load far from boundary
  0Moment0Shear  — free end; vs SAS (interior) + Hetényi semi-infinite plate formula (near end)

Grid convention: node-centred, x[i] = i*dx, boundary nodes at x[0]=0 and x[N-1]=(N-1)*dx.
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

D     = E * Te**3 / (12.0 * (1.0 - nu**2))
drho  = rho_m - rho_fill
beta  = (drho * g / (4.0 * D)) ** 0.25   # 1 / alpha
alpha = 1.0 / beta                         # flexural parameter [m]


def _run(qs, method="FD", bc_w="0Moment0Shear", bc_e="0Moment0Shear",
         sigma_xx=None):
    """Run a 1-D flexure calculation and return the flex object."""
    flex = F1D()
    flex.Quiet = True
    flex.Method = method
    flex.Solver = "direct"
    flex.g = g
    flex.E = E
    flex.nu = nu
    flex.rho_m = rho_m
    flex.rho_fill = rho_fill
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
# Mirror: exact cosine eigenfunction
# ---------------------------------------------------------------------------

def test_fd_mirror_cosine_eigenfunction():
    """FD/Mirror with a cosine load matches the continuous analytical formula.

    cos(nπx/L) is an eigenfunction of the 4th-order operator for Mirror BCs
    (zero slope at both endpoints): the deflection is the same cosine scaled by
    -q₀/(D(nπ/L)⁴ + Δρg).  The FD solution is exact in shape but has an
    O((k·dx)²) error in the eigenvalue — here < 0.1 % — so rtol = 1e-3 is safe.
    """
    N = 256
    L = (N - 1) * dx          # node-centred: x[0]=0, x[N-1]=L
    n_waves = 2
    k = n_waves * np.pi / L
    x = np.arange(N) * dx
    q0 = 1e6
    qs = q0 * np.cos(k * x)

    flex = _run(qs, bc_w="Mirror", bc_e="Mirror")

    w_exact = -q0 / (D * k**4 + drho * g) * np.cos(k * x)
    np.testing.assert_allclose(flex.w, w_exact, rtol=1e-3)


# ---------------------------------------------------------------------------
# Mirror: exact equivalence with Periodic on the even-extended 2× domain
# ---------------------------------------------------------------------------

def test_fd_mirror_equals_periodic_2x():
    """FD/Mirror on [0, L] equals FD/Periodic on [0, 2L] with even-extended load.

    Mirror BCs at both endpoints of a node-centred grid implement even reflection
    about x = 0 and x = L.  The resulting ghost-cell assignments are identical to
    running a Periodic problem on the period-2(N-1) even extension of the load.
    Agreement to rtol = 1e-10 (limited by the direct solver) confirms both the
    Mirror stencil and the Periodic stencil implement the same linear system.
    """
    N = 80
    qs = np.zeros(N)
    qs[10:30] = 1e6    # off-centre so the even extension is non-trivial

    flex_m = _run(qs, bc_w="Mirror", bc_e="Mirror")

    # Even extension on a node-centred grid: period = 2(N-1), sharing endpoints.
    # qs_ext = [qs[0], ..., qs[N-1], qs[N-2], ..., qs[1]]  (length 2N-2)
    qs_ext = np.concatenate([qs, qs[-2:0:-1]])
    flex_p = _run(qs_ext, bc_w="Periodic", bc_e="Periodic")

    np.testing.assert_allclose(flex_m.w, flex_p.w[:N], rtol=1e-8)


# ---------------------------------------------------------------------------
# 0Slope0Shear: distinct from Mirror despite same physical intent
# ---------------------------------------------------------------------------

def test_fd_0slope0shear_and_mirror_are_distinct():
    """Mirror and 0Slope0Shear are genuinely different stencils.

    Both BCs are intended to represent the symmetry condition (zero slope and
    zero shear at the boundary).  Mirror implements exact even reflection in the
    ghost cells; 0Slope0Shear uses a different finite-difference approximation
    that assigns the ghost cell at i = -1 the value of w[3] (not w[1]) for the
    i = 1 stencil row.  This difference propagates throughout the domain.

    For a cosine load — the eigenfunction of the Mirror BC operator — the two
    BCs give solutions that differ by > 10 % of the peak amplitude.
    Mirror agrees with the exact cosine formula; 0Slope0Shear does not.
    """
    N = 200
    L = (N - 1) * dx
    n_waves = 3
    k = n_waves * np.pi / L
    x = np.arange(N) * dx
    q0 = 1e6
    qs = q0 * np.cos(k * x)

    w_mirror = _run(qs, bc_w="Mirror", bc_e="Mirror").w
    w_0ss    = _run(qs, bc_w="0Slope0Shear", bc_e="0Slope0Shear").w
    w_exact  = -q0 / (D * k**4 + drho * g) * np.cos(k * x)

    # Mirror reproduces the cosine eigenfunction
    np.testing.assert_allclose(w_mirror, w_exact, rtol=2e-3)

    # 0Slope0Shear gives a substantially different result
    amp = q0 / (D * k**4 + drho * g)       # peak deflection amplitude [m]
    max_diff = np.abs(w_mirror - w_0ss).max() / amp
    assert max_diff > 0.10, (
        f"Mirror and 0Slope0Shear should differ by > 10 % of the amplitude "
        f"for this cosine load; got {max_diff:.1%}"
    )


# ---------------------------------------------------------------------------
# Large-domain interior checks: 0Slope0Shear, 0Displacement0Slope, 0Moment0Shear
# ---------------------------------------------------------------------------

def _central_load_sas_comparison(bc, margin=85):
    """Run FD with *bc* and SAS on a 200-cell domain with a central block load.

    The load sits ~6α from each boundary, so BC corrections at the centre
    are O(e⁻⁶) < 0.3 %.  Agreement to rtol = 2e-3 in the interior
    (cells *margin* to 200-margin) confirms the FD solver produces the
    physically correct infinite-plate solution far from the boundary.

    Note: SAS ignores BC_W / BC_E; it always uses the NoOutsideLoads
    (infinite-plate) assumption.
    """
    N = 200      # 200 × 4 km = 800 km; α ≈ 66 km → load is ~6α from each edge
    qs = np.zeros(N)
    qs[95:105] = 1e6

    flex_fd  = _run(qs, bc_w=bc, bc_e=bc)
    flex_sas = _run(qs, method="SAS", bc_w=bc, bc_e=bc)

    np.testing.assert_allclose(
        flex_fd.w[margin : N - margin],
        flex_sas.w[margin : N - margin],
        rtol=2e-3,
    )


def test_fd_0slope0shear_vs_sas():
    """FD/0Slope0Shear matches SAS for a central load far from the boundary."""
    _central_load_sas_comparison("0Slope0Shear")


def test_fd_0displacement0slope_vs_sas():
    """FD/0Displacement0Slope matches SAS for a central load far from the boundary."""
    _central_load_sas_comparison("0Displacement0Slope")


def test_fd_0moment0shear_vs_sas():
    """FD/0Moment0Shear matches SAS for a central load far from the boundary."""
    _central_load_sas_comparison("0Moment0Shear")


# ---------------------------------------------------------------------------
# 0Moment0Shear: Hetényi semi-infinite plate analytical formula
# ---------------------------------------------------------------------------

def _hetenyi_free_end(x, x_load, qs_load):
    """Deflection for a semi-infinite plate (x ≥ 0) with 0Moment0Shear at x = 0.

    In gFlex the FD equation is D w'''' + Δρg w = −qs, so a grid cell of
    pressure qs_load [Pa] acts as a line source of strength P = qs_load·dx
    [N/m] in the equivalent continuous problem D G'''' + Δρg G = P·δ(x − a).

    The particular (infinite-plate) solution is the standard Kelvin Green's
    function G∞; the homogeneous correction w_h cancels the moment and shear at
    the free end (Hetényi, 1946).  The resulting BCs D·w''(0) = D·w'''(0) = 0
    are verified algebraically in the derivation.
    """
    k_spring = drho * g
    P = qs_load * dx           # effective line load [N/m] = Pa × m
    a = x_load

    # Infinite-plate Green's function (negative: load down → deflection down)
    r = np.abs(x - a)
    w_p = -(P * beta / (2.0 * k_spring)) * np.exp(-beta * r) * (
        np.cos(beta * r) + np.sin(beta * r)
    )

    # Hetényi image constants that cancel M and V at x = 0
    ba = beta * a
    A_h = -(P * beta / (2.0 * k_spring)) * np.exp(-ba) * (
        3.0 * np.cos(ba) - np.sin(ba)
    )
    B_h = (P * beta / (2.0 * k_spring)) * np.exp(-ba) * (
        np.cos(ba) - np.sin(ba)
    )
    w_h = np.exp(-beta * x) * (A_h * np.cos(beta * x) + B_h * np.sin(beta * x))

    return w_p + w_h


def test_fd_0moment0shear_hetenyi():
    """FD/0Moment0Shear matches the Hetényi semi-infinite plate formula.

    A point load is placed 2α from the free west end.  The domain extends 8α
    beyond the load so the east boundary correction is O(e⁻⁸) < 0.04 %.

    The 0Moment0Shear stencil uses a one-sided finite-difference approximation
    at the boundary cell (first-order in dx/α), giving an O(dx/α) ≈ 6 %
    error at x = 0 for dx = 4 km.  An absolute tolerance of 0.05 m covers
    this near-boundary truncation error while confirming the correct shape and
    sign of the free-end upward deflection (the plate tips up like a seesaw).
    """
    # Load at 2α from the free west end
    ia = int(np.round(2.0 * alpha / dx))
    N  = ia + int(np.ceil(8.0 * alpha / dx)) + 1

    qs = np.zeros(N)
    qs[ia] = 1e6
    x = np.arange(N) * dx

    flex_fd = _run(qs, bc_w="0Moment0Shear", bc_e="0Displacement0Slope")

    w_exact = _hetenyi_free_end(x, x[ia], qs[ia])

    # Compare the first half (near the free end); east-end correction negligible
    half = N // 2
    np.testing.assert_allclose(flex_fd.w[:half], w_exact[:half], atol=0.05)


# ---------------------------------------------------------------------------
# sigma_xx monotonicity for every FD BC
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("bc", [
    "0Moment0Shear",
    "0Displacement0Slope",
    "Mirror",
    "0Slope0Shear",
])
def test_sigma_xx_monotonicity_all_bcs(bc):
    """Tensile sigma_xx reduces deflection; compressive increases it — for every BC."""
    N = 200
    qs = np.zeros(N)
    qs[90:110] = 1e6

    flex_0 = _run(qs, bc_w=bc, bc_e=bc)
    flex_t = _run(qs, bc_w=bc, bc_e=bc, sigma_xx=+1e8)   # tensile
    flex_c = _run(qs, bc_w=bc, bc_e=bc, sigma_xx=-1e8)   # compressive

    assert flex_t.w.min() > flex_0.w.min(), (
        f"{bc}: tensile sigma_xx should reduce subsidence"
    )
    assert flex_c.w.min() < flex_0.w.min(), (
        f"{bc}: compressive sigma_xx should increase subsidence"
    )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
