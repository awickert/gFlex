#!/usr/bin/env python
"""Analytical tests for inhomogeneous (prescribed-value) 1-D boundary conditions.

STATUS: These tests FAIL until issues #51 / #52 are implemented.
        They define the target behaviour; they will pass once the
        prescribed-BC machinery is in place.

Semi-infinite plate derivation
-------------------------------
For a plate on a Winkler foundation with no distributed load,
D·d⁴w/dx⁴ + Δρg·w = 0, the solution decaying as x → ∞ is

    w(x) = e^{-λx} [C₁ cos(λx) + C₂ sin(λx)],   λ = (Δρg / 4D)^{1/4}

Derivatives at x = 0:
    w(0)    = C₁
    w'(0)   = λ(C₂ − C₁)                     [slope]
    w''(0)  = −2λ²C₂  →  M₀ = D·w''(0)       [bending moment]
    w'''(0) = 2λ³(C₁+C₂)  →  V₀ = D·w'''(0)  [shear force]

Four canonical cases
--------------------
A. Prescribed moment M₀, zero shear:
     C₁ = M₀/(2Dλ²),  C₂ = −M₀/(2Dλ²)
     w(x) = [M₀/(2Dλ²)] e^{-λx} [cos(λx) − sin(λx)]

B. Prescribed shear V₀, zero moment:
     C₁ = V₀/(2Dλ³),  C₂ = 0
     w(x) = [V₀/(2Dλ³)] e^{-λx} cos(λx)

C. Prescribed displacement w₀, zero slope (clamped, offset):
     C₁ = C₂ = w₀
     w(x) = w₀ e^{-λx} [cos(λx) + sin(λx)]

D. Prescribed slope θ₀, zero displacement:
     C₁ = 0,  C₂ = θ₀/λ
     w(x) = (θ₀/λ) e^{-λx} sin(λx)
"""

import warnings

import numpy as np

from gflex.f1d import F1D

# ---------------------------------------------------------------------------
# Physical parameters  (consistent with the rest of the test suite)
# ---------------------------------------------------------------------------

E        = 65.0e9     # Young's modulus, Pa
NU       = 0.25
TE       = 30.0e3     # elastic thickness, m
RHO_M    = 3300.0     # mantle density, kg m⁻³
RHO_F    = 0.0        # infill density (air)
G        = 9.8        # m s⁻²

D      = E * TE**3 / (12.0 * (1.0 - NU**2))
DRHOG  = (RHO_M - RHO_F) * G
LAM    = (DRHOG / (4.0 * D)) ** 0.25   # λ = 1/α
ALPHA  = 1.0 / LAM                      # flexural parameter, m

# Domain long enough that the far-end BC contributes < 0.1 % truncation error.
# At x = 10α: e^{-10} ≈ 4.5×10⁻⁵.
L_DOMAIN = 10.0 * ALPHA   # m
NX       = 401             # Δx ≈ α/40 → O(Δx²) error ≈ 0.06 %

REL_TOL = 0.01   # 1 % L-inf relative error tolerance


# ---------------------------------------------------------------------------
# Exact semi-infinite solutions
# ---------------------------------------------------------------------------

def _w_exact(x, C1, C2):
    lx = LAM * x
    return np.exp(-lx) * (C1 * np.cos(lx) + C2 * np.sin(lx))


def w_prescribed_moment(x, M0):
    """Case A: M(0) = M₀, V(0) = 0."""
    c = M0 / (2.0 * D * LAM**2)
    return _w_exact(x, c, -c)


def w_prescribed_shear(x, V0):
    """Case B: M(0) = 0, V(0) = V₀."""
    c = V0 / (2.0 * D * LAM**3)
    return _w_exact(x, c, 0.0)


def w_prescribed_displacement(x, w0):
    """Case C: w(0) = w₀, w'(0) = 0."""
    return _w_exact(x, w0, w0)


def w_prescribed_slope(x, theta0):
    """Case D: w(0) = 0, w'(0) = θ₀."""
    return _w_exact(x, 0.0, theta0 / LAM)


# ---------------------------------------------------------------------------
# gFlex runner
# ---------------------------------------------------------------------------

def _run(bc_west, bc_east, nx=NX, L=L_DOMAIN, qs=None):
    dx = L / (nx - 1)
    x  = np.arange(nx) * dx

    flex = F1D()
    flex.quiet    = True
    flex.verbose  = False
    flex.debug    = False
    flex.method   = "fd"
    flex.g        = G
    flex.E        = E
    flex.nu       = NU
    flex.rho_m    = RHO_M
    flex.rho_fill = RHO_F
    flex.T_e      = TE
    flex.dx       = dx
    # Deflection driven by BCs (zero load) unless a load is supplied.
    flex.qs       = np.zeros(nx) if qs is None else np.asarray(qs, dtype=float)
    flex.bc_west  = bc_west
    flex.bc_east  = bc_east

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        flex.initialize()
        flex.run()
        w = flex.w
        flex.finalize()

    return x, w


def _check(w_num, w_ex, label):
    scale = np.max(np.abs(w_ex))
    err   = np.max(np.abs(w_num - w_ex)) / scale
    assert err < REL_TOL, (
        f"{label}: L-inf relative error {err:.3e} exceeds {REL_TOL:.0%}"
    )


# ---------------------------------------------------------------------------
# Case A — prescribed moment, zero shear
# ---------------------------------------------------------------------------

class TestPrescribedMoment:
    """M(0) = M₀, V(0) = 0; zero-moment/zero-shear (free) at far end."""

    M0 = 1.0e12   # N·m / m

    def test_deflection_profile(self):
        x, w_num = _run(
            bc_west={"moment": self.M0, "shear": 0.0},
            bc_east="zero_moment_zero_shear",
        )
        _check(w_num, w_prescribed_moment(x, self.M0), "Case A profile")

    def test_west_boundary_displacement(self):
        """Boundary node should have the analytically determined displacement."""
        x, w_num = _run(
            bc_west={"moment": self.M0, "shear": 0.0},
            bc_east="zero_moment_zero_shear",
        )
        w0_exact = self.M0 / (2.0 * D * LAM**2)   # C₁
        assert abs(w_num[0] - w0_exact) / abs(w0_exact) < REL_TOL


# ---------------------------------------------------------------------------
# Case B — zero moment, prescribed shear
# ---------------------------------------------------------------------------

class TestPrescribedShear:
    """M(0) = 0, V(0) = V₀; free at far end."""

    V0 = 1.0e8   # N / m

    def test_deflection_profile(self):
        x, w_num = _run(
            bc_west={"moment": 0.0, "shear": self.V0},
            bc_east="zero_moment_zero_shear",
        )
        _check(w_num, w_prescribed_shear(x, self.V0), "Case B profile")

    def test_west_boundary_displacement(self):
        x, w_num = _run(
            bc_west={"moment": 0.0, "shear": self.V0},
            bc_east="zero_moment_zero_shear",
        )
        w0_exact = self.V0 / (2.0 * D * LAM**3)   # C₁
        assert abs(w_num[0] - w0_exact) / abs(w0_exact) < REL_TOL


# ---------------------------------------------------------------------------
# Case C — prescribed displacement, zero slope  (clamped, offset)
# ---------------------------------------------------------------------------

class TestPrescribedDisplacement:
    """w(0) = w₀, w'(0) = 0; clamped at far end."""

    W0 = 100.0   # m

    def test_deflection_profile(self):
        x, w_num = _run(
            bc_west={"displacement": self.W0, "slope": 0.0},
            bc_east="zero_displacement_zero_slope",
        )
        _check(w_num, w_prescribed_displacement(x, self.W0), "Case C profile")

    def test_west_boundary_displacement(self):
        _, w_num = _run(
            bc_west={"displacement": self.W0, "slope": 0.0},
            bc_east="zero_displacement_zero_slope",
        )
        assert abs(w_num[0] - self.W0) / self.W0 < REL_TOL


# ---------------------------------------------------------------------------
# Case D — zero displacement, prescribed slope
# ---------------------------------------------------------------------------

class TestPrescribedSlope:
    """w(0) = 0, w'(0) = θ₀; free at far end."""

    THETA0 = 1.0e-3   # rad

    def test_deflection_profile(self):
        x, w_num = _run(
            bc_west={"displacement": 0.0, "slope": self.THETA0},
            bc_east="zero_moment_zero_shear",
        )
        _check(w_num, w_prescribed_slope(x, self.THETA0), "Case D profile")


# ---------------------------------------------------------------------------
# Case E — prescribed displacement + slope at BOTH ends (east-path regression)
# ---------------------------------------------------------------------------
#
# Manufactured solution  w(x) = c·x·(x − L):
#   * w'''' = 0          → the FD bending operator reproduces it exactly, so
#                          the only error source is boundary-condition handling;
#   * w(0) = w(L) = 0    → the load q = −Δρg·w vanishes at the boundary nodes,
#                          so the prescribed displacement is uncontaminated;
#   * w'(0) = −cL, w'(L) = +cL  → both slope-BC paths are exercised, with the
#                          west path acting as a control on the east path.
#
# This pins the east-boundary slope RHS sign.  The earlier sign produced a
# ~28 % L-inf error localized just inside the east boundary; the corrected
# code reproduces w to solver precision.  Cases A–D above only prescribe at
# the west end, so the east inhomogeneous path was previously untested.

class TestPrescribedDisplacementSlopeBothEnds:
    """Quadratic MMS reproduced to machine precision (east-slope regression)."""

    C = 1.0e-9   # curvature scale, m⁻¹

    def test_quadratic_mms_exact(self):
        nx  = 201
        dx  = L_DOMAIN / (nx - 1)
        x   = np.arange(nx) * dx
        w_m = self.C * x * (x - L_DOMAIN)        # manufactured deflection
        wp  = self.C * (2.0 * x - L_DOMAIN)      # slope w'(x)
        qs  = -DRHOG * w_m                       # interior: Δρg·w = −qs

        _, w_num = _run(
            bc_west={"displacement": float(w_m[0]),  "slope": float(wp[0])},
            bc_east={"displacement": float(w_m[-1]), "slope": float(wp[-1])},
            nx=nx, qs=qs,
        )

        scale = np.max(np.abs(w_m))
        err   = np.max(np.abs(w_num - w_m)) / scale
        assert err < 1.0e-8, (
            f"quadratic MMS L-inf relative error {err:.3e} exceeds 1e-8 "
            "(east-slope BC RHS sign regression?)"
        )
