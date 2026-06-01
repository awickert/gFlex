#!/usr/bin/env python
"""Analytical tests for inhomogeneous (prescribed-value) 2-D boundary conditions.

STATUS: These tests FAIL until 2-D inhomogeneous BC machinery is implemented.
        They define the target behaviour; they will pass once the
        prescribed-BC machinery is in place.

Mirror-structure counterpart to test_1D_inhomogeneous_BC.py
-----------------------------------------------------------
The domain is long in x (10α) and narrow in y (NY=11 cells).  The north and
south edges carry ``bc_north="mirror"`` / ``bc_south="mirror"``, which forces
dw/dy = 0 there.  With a uniform load (qs = 0) and y-independent BCs the 2-D
plate equation reduces to its 1-D form, so the numerical solution should be
y-uniform and the centre row should match the 1-D semi-infinite exact solution.

Semi-infinite plate derivation  (same as 1-D file)
---------------------------------------------------
D·d⁴w/dx⁴ + Δρg·w = 0, solution decaying as x → ∞:

    w(x) = e^{-λx} [C₁ cos(λx) + C₂ sin(λx)],   λ = (Δρg / 4D)^{1/4}

Four canonical cases (A–D) are identical to the 1-D tests.
"""

import warnings

import numpy as np

from gflex.f2d import F2D

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

# Domain long enough that far-end BC contributes < 0.1 % truncation error.
L_DOMAIN = 10.0 * ALPHA   # x-extent, m
NX       = 401             # Δx ≈ α/40
NY       = 11              # narrow y-extent; solution is y-uniform

REL_TOL = 0.01   # 1 % L-inf relative error tolerance


# ---------------------------------------------------------------------------
# Exact semi-infinite solutions  (same as 1-D)
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
# gFlex 2-D runner
# ---------------------------------------------------------------------------

def _run(bc_west, bc_east, nx=NX, ny=NY, L=L_DOMAIN):
    dx = L / (nx - 1)
    x  = np.arange(nx) * dx

    flex = F2D()
    flex.quiet    = True
    flex.method   = "fd"
    flex.solver   = "direct"
    flex.g        = G
    flex.E        = E
    flex.nu       = NU
    flex.rho_m    = RHO_M
    flex.rho_fill = RHO_F
    flex.te       = TE
    flex.dx       = dx
    flex.dy       = dx                       # square cells
    flex.qs       = np.zeros((ny, nx))       # deflection driven by BCs only
    flex.bc_west  = bc_west
    flex.bc_east  = bc_east
    flex.bc_north = "mirror"                 # enforce y-uniformity
    flex.bc_south = "mirror"

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        flex.initialize()
        flex.run()
        flex.finalize()

    return x, flex.w


def _check(w_num, w_ex, label):
    """Compare centre row of 2-D solution to 1-D exact; L-inf relative error."""
    centre = w_num[NY // 2, :]
    scale  = np.max(np.abs(w_ex))
    err    = np.max(np.abs(centre - w_ex)) / scale
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
        """Centre-row boundary node should match the analytical displacement."""
        x, w_num = _run(
            bc_west={"moment": self.M0, "shear": 0.0},
            bc_east="zero_moment_zero_shear",
        )
        w0_exact = self.M0 / (2.0 * D * LAM**2)   # C₁
        assert abs(w_num[NY // 2, 0] - w0_exact) / abs(w0_exact) < REL_TOL


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
        assert abs(w_num[NY // 2, 0] - w0_exact) / abs(w0_exact) < REL_TOL


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
        assert abs(w_num[NY // 2, 0] - self.W0) / self.W0 < REL_TOL


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
# Variable-Te Dirichlet exactness tests
# ---------------------------------------------------------------------------

_VTE_NX = 51
_VTE_NY = 11
_VTE_W0 = 75.0    # prescribed displacement, m

_DISP_BC = {"displacement": _VTE_W0, "slope": 0.0}
_CLAMP   = "zero_displacement_zero_slope"
_MIRROR  = "mirror"


def _make_variable_te(ny, nx):
    """Smooth 2-D Te array: ±30 % variation in both directions."""
    eta = np.linspace(0, 1, ny)[:, np.newaxis]
    xi  = np.linspace(0, 1, nx)[np.newaxis, :]
    return TE * (1.0 + 0.3 * np.sin(np.pi * eta) * np.cos(np.pi * xi))


def _run_vte(bc_west, bc_east, bc_north, bc_south,
             nx=_VTE_NX, ny=_VTE_NY):
    dx = L_DOMAIN / (nx - 1)
    flex = F2D()
    flex.quiet    = True
    flex.method   = "fd"
    flex.solver   = "direct"
    flex.g        = G
    flex.E        = E
    flex.nu       = NU
    flex.rho_m    = RHO_M
    flex.rho_fill = RHO_F
    flex.te       = _make_variable_te(ny, nx)
    flex.dx       = dx
    flex.dy       = dx
    flex.qs       = np.zeros((ny, nx))
    flex.bc_west  = bc_west
    flex.bc_east  = bc_east
    flex.bc_north = bc_north
    flex.bc_south = bc_south
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        flex.initialize()
        flex.run()
        flex.finalize()
    return flex.w


class TestVariableTeDirichletExact:
    """Prescribed displacement is enforced exactly under spatially varying Te.

    The RHS correction for displacement/slope BCs is a pure Dirichlet
    decoupling: correction[edge] = cj0i0_coeff[edge] * w0.  The coefficient
    depends on D, but w = w0 holds to floating-point precision regardless of
    how D varies.  Tests cover each single edge, both pairs of opposite edges,
    all four adjacent-edge (corner) pairs, and all four edges simultaneously.
    """

    # --- single edges ---

    def test_west_exact(self):
        w = _run_vte(_DISP_BC, _CLAMP, _MIRROR, _MIRROR)
        assert np.all(w[:, 0] == _VTE_W0), "west edge not exactly W0"

    def test_east_exact(self):
        w = _run_vte(_CLAMP, _DISP_BC, _MIRROR, _MIRROR)
        assert np.all(w[:, -1] == _VTE_W0), "east edge not exactly W0"

    def test_north_exact(self):
        w = _run_vte(_CLAMP, _CLAMP, _DISP_BC, _MIRROR)
        assert np.all(w[0, :] == _VTE_W0), "north edge not exactly W0"

    def test_south_exact(self):
        w = _run_vte(_CLAMP, _CLAMP, _MIRROR, _DISP_BC)
        assert np.all(w[-1, :] == _VTE_W0), "south edge not exactly W0"

    # --- opposite-edge pairs ---

    def test_west_east_exact(self):
        w = _run_vte(_DISP_BC, _DISP_BC, _MIRROR, _MIRROR)
        assert np.all(w[:, 0]  == _VTE_W0), "west edge not exactly W0"
        assert np.all(w[:, -1] == _VTE_W0), "east edge not exactly W0"

    def test_north_south_exact(self):
        w = _run_vte(_CLAMP, _CLAMP, _DISP_BC, _DISP_BC)
        assert np.all(w[0, :]  == _VTE_W0), "north edge not exactly W0"
        assert np.all(w[-1, :] == _VTE_W0), "south edge not exactly W0"

    # --- adjacent-edge (corner) pairs ---

    def test_northwest_exact(self):
        w = _run_vte(_DISP_BC, _CLAMP, _DISP_BC, _MIRROR)
        assert np.all(w[:, 0] == _VTE_W0), "west edge not exactly W0"
        assert np.all(w[0, :] == _VTE_W0), "north edge not exactly W0"

    def test_northeast_exact(self):
        w = _run_vte(_CLAMP, _DISP_BC, _DISP_BC, _MIRROR)
        assert np.all(w[:, -1] == _VTE_W0), "east edge not exactly W0"
        assert np.all(w[0, :]  == _VTE_W0), "north edge not exactly W0"

    def test_southwest_exact(self):
        w = _run_vte(_DISP_BC, _CLAMP, _MIRROR, _DISP_BC)
        assert np.all(w[:, 0]  == _VTE_W0), "west edge not exactly W0"
        assert np.all(w[-1, :] == _VTE_W0), "south edge not exactly W0"

    def test_southeast_exact(self):
        w = _run_vte(_CLAMP, _DISP_BC, _MIRROR, _DISP_BC)
        assert np.all(w[:, -1] == _VTE_W0), "east edge not exactly W0"
        assert np.all(w[-1, :] == _VTE_W0), "south edge not exactly W0"

    # --- all four edges ---

    def test_all_four_exact(self):
        w = _run_vte(_DISP_BC, _DISP_BC, _DISP_BC, _DISP_BC)
        assert np.all(w[:, 0]  == _VTE_W0), "west edge not exactly W0"
        assert np.all(w[:, -1] == _VTE_W0), "east edge not exactly W0"
        assert np.all(w[0, :]  == _VTE_W0), "north edge not exactly W0"
        assert np.all(w[-1, :] == _VTE_W0), "south edge not exactly W0"


# ---------------------------------------------------------------------------
# North / South domain helpers
# ---------------------------------------------------------------------------
# Domain long in y (NY_NS cells) and narrow in x (NX_NS cells);
# mirror on east/west enforces x-uniformity so the 2-D problem reduces
# to the 1-D semi-infinite plate equation in the y direction.

_NS_NX = 11
_NS_NY = 401   # Δy ≈ α/40, matching the x-direction resolution


def _run_ns(bc_north, bc_south, nx=_NS_NX, ny=_NS_NY, L=L_DOMAIN):
    dy = L / (ny - 1)
    flex = F2D()
    flex.quiet    = True
    flex.method   = "fd"
    flex.solver   = "direct"
    flex.g        = G
    flex.E        = E
    flex.nu       = NU
    flex.rho_m    = RHO_M
    flex.rho_fill = RHO_F
    flex.te       = TE
    flex.dx       = dy          # square cells
    flex.dy       = dy
    flex.qs       = np.zeros((ny, nx))
    flex.bc_west  = "mirror"
    flex.bc_east  = "mirror"
    flex.bc_north = bc_north
    flex.bc_south = bc_south
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        flex.initialize()
        flex.run()
        flex.finalize()
    y = np.arange(ny) * dy
    return y, flex.w


def _check_col(w_num, w_ex, label):
    """Compare centre column of 2-D solution to 1-D exact; L-inf relative error."""
    centre = w_num[:, _NS_NX // 2]
    scale  = np.max(np.abs(w_ex))
    err    = np.max(np.abs(centre - w_ex)) / scale
    assert err < REL_TOL, (
        f"{label}: L-inf relative error {err:.3e} exceeds {REL_TOL:.0%}"
    )


# ---------------------------------------------------------------------------
# Square-domain runner for Neumann superposition tests
# ---------------------------------------------------------------------------

_SQ_N = 31   # small square domain — solution need not decay to machine zero


def _run_sq(bc_west, bc_east, bc_north, bc_south, n=_SQ_N):
    L  = 3.0 * ALPHA          # domain just large enough for meaningful deflection
    dx = L / (n - 1)
    flex = F2D()
    flex.quiet    = True
    flex.method   = "fd"
    flex.solver   = "direct"
    flex.g        = G
    flex.E        = E
    flex.nu       = NU
    flex.rho_m    = RHO_M
    flex.rho_fill = RHO_F
    flex.te       = TE
    flex.dx       = dx
    flex.dy       = dx
    flex.qs       = np.zeros((n, n))
    flex.bc_west  = bc_west
    flex.bc_east  = bc_east
    flex.bc_north = bc_north
    flex.bc_south = bc_south
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        flex.initialize()
        flex.run()
        flex.finalize()
    return flex.w


# ---------------------------------------------------------------------------
# East edge analytical tests  (Cases A – D)
# ---------------------------------------------------------------------------
# Exact solutions use ξ = L_DOMAIN − x (distance from east boundary).
#
# Sign note — shear and slope are ODD under x → L − x because
# d³/dx³ and d/dx each pick up a sign flip:
#   east V₀ → w(ξ) = −V₀/(2Dλ³) e^{−λξ} cos(λξ) = −w_prescribed_shear(ξ, V₀)
#   east θ₀ → w(ξ) = −(θ₀/λ) e^{−λξ} sin(λξ)    = −w_prescribed_slope(ξ, θ₀)
# Moment and displacement are EVEN and carry no sign change.
# ---------------------------------------------------------------------------

class TestEastPrescribedMoment:
    """M(L) = M₀, V(L) = 0; zero-moment/zero-shear (free) at west end."""

    M0 = 1.0e12   # N·m / m

    def test_deflection_profile(self):
        x, w_num = _run(
            bc_west="zero_moment_zero_shear",
            bc_east={"moment": self.M0, "shear": 0.0},
        )
        _check(w_num, w_prescribed_moment(L_DOMAIN - x, self.M0), "East Case A profile")

    def test_east_boundary_displacement(self):
        x, w_num = _run(
            bc_west="zero_moment_zero_shear",
            bc_east={"moment": self.M0, "shear": 0.0},
        )
        w0_exact = self.M0 / (2.0 * D * LAM**2)
        assert abs(w_num[NY // 2, -1] - w0_exact) / abs(w0_exact) < REL_TOL


class TestEastPrescribedShear:
    """M(L) = 0, V(L) = V₀; free at west end.

    d³w/dx³ at east reverses sign under ξ = L − x, so V₀ > 0 at east
    produces negative deflection: w(ξ) = −w_prescribed_shear(ξ, V₀).
    """

    V0 = 1.0e8   # N / m

    def test_deflection_profile(self):
        x, w_num = _run(
            bc_west="zero_moment_zero_shear",
            bc_east={"moment": 0.0, "shear": self.V0},
        )
        _check(w_num, -w_prescribed_shear(L_DOMAIN - x, self.V0), "East Case B profile")

    def test_east_boundary_displacement(self):
        x, w_num = _run(
            bc_west="zero_moment_zero_shear",
            bc_east={"moment": 0.0, "shear": self.V0},
        )
        w0_exact = -self.V0 / (2.0 * D * LAM**3)   # negative: see class docstring
        assert abs(w_num[NY // 2, -1] - w0_exact) / abs(w0_exact) < REL_TOL


class TestEastPrescribedDisplacement:
    """w(L) = W₀, dw/dx(L) = 0; clamped (zero displacement/slope) at west end."""

    W0 = 100.0   # m

    def test_deflection_profile(self):
        x, w_num = _run(
            bc_west="zero_displacement_zero_slope",
            bc_east={"displacement": self.W0, "slope": 0.0},
        )
        _check(w_num, w_prescribed_displacement(L_DOMAIN - x, self.W0), "East Case C profile")

    def test_east_boundary_displacement(self):
        _, w_num = _run(
            bc_west="zero_displacement_zero_slope",
            bc_east={"displacement": self.W0, "slope": 0.0},
        )
        assert abs(w_num[NY // 2, -1] - self.W0) / self.W0 < REL_TOL


class TestEastPrescribedSlope:
    """w(L) = 0, dw/dx(L) = +θ₀; free at west end.

    dw/dx = −dw/dξ at east, so θ₀ > 0 gives negative deflection:
    w(ξ) = −w_prescribed_slope(ξ, θ₀).
    """

    THETA0 = 1.0e-3   # rad

    def test_deflection_profile(self):
        x, w_num = _run(
            bc_west="zero_moment_zero_shear",
            bc_east={"displacement": 0.0, "slope": self.THETA0},
        )
        _check(w_num, -w_prescribed_slope(L_DOMAIN - x, self.THETA0), "East Case D profile")


# ---------------------------------------------------------------------------
# North edge analytical tests  (Cases A – C)
# ---------------------------------------------------------------------------
# North outward normal is −y, the same sense as west (−x), so ALL four
# sign conventions carry over unchanged (y replaces x, ny replaces nx).
# ---------------------------------------------------------------------------

class TestNorthPrescribedMoment:
    """M(y=0) = M₀, V(y=0) = 0; free at south end."""

    M0 = 1.0e12

    def test_deflection_profile(self):
        y, w_num = _run_ns(
            bc_north={"moment": self.M0, "shear": 0.0},
            bc_south="zero_moment_zero_shear",
        )
        _check_col(w_num, w_prescribed_moment(y, self.M0), "North Case A profile")

    def test_north_boundary_displacement(self):
        y, w_num = _run_ns(
            bc_north={"moment": self.M0, "shear": 0.0},
            bc_south="zero_moment_zero_shear",
        )
        w0_exact = self.M0 / (2.0 * D * LAM**2)
        assert abs(w_num[0, _NS_NX // 2] - w0_exact) / abs(w0_exact) < REL_TOL


class TestNorthPrescribedShear:
    """M(y=0) = 0, V(y=0) = V₀; free at south end."""

    V0 = 1.0e8

    def test_deflection_profile(self):
        y, w_num = _run_ns(
            bc_north={"moment": 0.0, "shear": self.V0},
            bc_south="zero_moment_zero_shear",
        )
        _check_col(w_num, w_prescribed_shear(y, self.V0), "North Case B profile")


class TestNorthPrescribedDisplacement:
    """w(y=0) = W₀, dw/dy(y=0) = 0; clamped at south end."""

    W0 = 100.0

    def test_deflection_profile(self):
        y, w_num = _run_ns(
            bc_north={"displacement": self.W0, "slope": 0.0},
            bc_south="zero_displacement_zero_slope",
        )
        _check_col(w_num, w_prescribed_displacement(y, self.W0), "North Case C profile")

    def test_north_boundary_displacement(self):
        y, w_num = _run_ns(
            bc_north={"displacement": self.W0, "slope": 0.0},
            bc_south="zero_displacement_zero_slope",
        )
        assert abs(w_num[0, _NS_NX // 2] - self.W0) / self.W0 < REL_TOL


class TestNorthPrescribedSlope:
    """w(y=0) = 0, dw/dy(y=0) = +θ₀; free at south end.

    Outward normal at north is −y (same sense as west −x), so the slope
    sign convention carries over unchanged: θ₀ > 0 gives positive deflection.
    w(y) = w_prescribed_slope(y, θ₀).
    """

    THETA0 = 1.0e-3   # rad

    def test_deflection_profile(self):
        y, w_num = _run_ns(
            bc_north={"displacement": 0.0, "slope": self.THETA0},
            bc_south="zero_moment_zero_shear",
        )
        _check_col(w_num, w_prescribed_slope(y, self.THETA0), "North Case D profile")


# ---------------------------------------------------------------------------
# South edge analytical tests  (Cases A – D)
# ---------------------------------------------------------------------------
# South outward normal is +y, the same sense as east (+x), so shear and
# slope are ODD under ξ = L_DOMAIN − y and pick up a sign flip.
# Moment and displacement are EVEN and carry no sign change.
# ---------------------------------------------------------------------------

class TestSouthPrescribedMoment:
    """M(y=L) = M₀, V(y=L) = 0; free at north end."""

    M0 = 1.0e12

    def test_deflection_profile(self):
        y, w_num = _run_ns(
            bc_north="zero_moment_zero_shear",
            bc_south={"moment": self.M0, "shear": 0.0},
        )
        _check_col(w_num, w_prescribed_moment(L_DOMAIN - y, self.M0), "South Case A profile")

    def test_south_boundary_displacement(self):
        y, w_num = _run_ns(
            bc_north="zero_moment_zero_shear",
            bc_south={"moment": self.M0, "shear": 0.0},
        )
        w0_exact = self.M0 / (2.0 * D * LAM**2)
        assert abs(w_num[-1, _NS_NX // 2] - w0_exact) / abs(w0_exact) < REL_TOL


class TestSouthPrescribedShear:
    """M(y=L) = 0, V(y=L) = V₀; free at north end.

    d³w/dy³ at south reverses sign under ξ = L − y, so V₀ > 0 gives
    negative deflection: w(ξ) = −w_prescribed_shear(ξ, V₀).
    """

    V0 = 1.0e8

    def test_deflection_profile(self):
        y, w_num = _run_ns(
            bc_north="zero_moment_zero_shear",
            bc_south={"moment": 0.0, "shear": self.V0},
        )
        _check_col(w_num, -w_prescribed_shear(L_DOMAIN - y, self.V0), "South Case B profile")

class TestSouthPrescribedDisplacement:
    """w(y=L) = W₀, dw/dy(y=L) = 0; clamped at north end."""

    W0 = 100.0

    def test_deflection_profile(self):
        y, w_num = _run_ns(
            bc_north="zero_displacement_zero_slope",
            bc_south={"displacement": self.W0, "slope": 0.0},
        )
        _check_col(w_num, w_prescribed_displacement(L_DOMAIN - y, self.W0), "South Case C profile")

    def test_south_boundary_displacement(self):
        y, w_num = _run_ns(
            bc_north="zero_displacement_zero_slope",
            bc_south={"displacement": self.W0, "slope": 0.0},
        )
        assert abs(w_num[-1, _NS_NX // 2] - self.W0) / self.W0 < REL_TOL


class TestSouthPrescribedSlope:
    """w(y=L) = 0, dw/dy(y=L) = +θ₀; free at north end.

    dw/dy at south reverses sign under ξ = L − y, so θ₀ > 0 gives
    negative deflection: w(ξ) = −w_prescribed_slope(ξ, θ₀).
    """

    THETA0 = 1.0e-3   # rad

    def test_deflection_profile(self):
        y, w_num = _run_ns(
            bc_north="zero_moment_zero_shear",
            bc_south={"displacement": 0.0, "slope": self.THETA0},
        )
        _check_col(w_num, -w_prescribed_slope(L_DOMAIN - y, self.THETA0), "South Case D profile")



# ---------------------------------------------------------------------------
# Mixed Neumann + Dirichlet: one end prescribed moment (ZMZS type),
# the other prescribed displacement (ZDSZS type).
# ---------------------------------------------------------------------------
# For L = 10α each solution decays to < e^{-10} ≈ 4.5 e-5 of its peak
# before reaching the opposite end, so the superposition approximation
# w ≈ w_near + w_far holds to < 0.01 % — well within the 1 % tolerance.
#
# Corner interactions exercised:
#   west-moment + north/south mirror  → ZMZS+mirror correction (NW/SW)
#   east-disp   + north/south mirror  → ZDSZS+mirror: all inf-flagged, no-op
#   west-disp   + north/south mirror  → ZDSZS+mirror: all inf-flagged, no-op
#   east-moment + north/south mirror  → ZMZS+mirror correction (NE/SE)
# ---------------------------------------------------------------------------

class TestMixedNeumannDirichlet:
    """Neumann (moment/shear) at one end, Dirichlet (displacement/slope) at the other.

    W0 is chosen so both ends produce equal peak deflection, ensuring neither
    contribution dominates the L-inf error comparison.
    """

    M0 = 1.0e12                        # prescribed moment (Neumann end)
    W0 = M0 / (2.0 * D * LAM**2)      # prescribed displacement (Dirichlet end)

    # --- east / west ---

    def test_west_moment_east_displacement(self):
        """West: moment (ZMZS); East: displacement (ZDSZS)."""
        x, w_num = _run(
            bc_west={"moment": self.M0, "shear": 0.0},
            bc_east={"displacement": self.W0, "slope": 0.0},
        )
        w_ex = (w_prescribed_moment(x, self.M0)
                + w_prescribed_displacement(L_DOMAIN - x, self.W0))
        _check(w_num, w_ex, "Mixed west-moment east-disp")

    def test_west_displacement_east_moment(self):
        """West: displacement (ZDSZS); East: moment (ZMZS)."""
        x, w_num = _run(
            bc_west={"displacement": self.W0, "slope": 0.0},
            bc_east={"moment": self.M0, "shear": 0.0},
        )
        w_ex = (w_prescribed_displacement(x, self.W0)
                + w_prescribed_moment(L_DOMAIN - x, self.M0))
        _check(w_num, w_ex, "Mixed west-disp east-moment")

    # --- north / south ---

    def test_north_moment_south_displacement(self):
        """North: moment (ZMZS); South: displacement (ZDSZS)."""
        y, w_num = _run_ns(
            bc_north={"moment": self.M0, "shear": 0.0},
            bc_south={"displacement": self.W0, "slope": 0.0},
        )
        w_ex = (w_prescribed_moment(y, self.M0)
                + w_prescribed_displacement(L_DOMAIN - y, self.W0))
        _check_col(w_num, w_ex, "Mixed north-moment south-disp")

    def test_north_displacement_south_moment(self):
        """North: displacement (ZDSZS); South: moment (ZMZS)."""
        y, w_num = _run_ns(
            bc_north={"displacement": self.W0, "slope": 0.0},
            bc_south={"moment": self.M0, "shear": 0.0},
        )
        w_ex = (w_prescribed_displacement(y, self.W0)
                + w_prescribed_moment(L_DOMAIN - y, self.M0))
        _check_col(w_num, w_ex, "Mixed north-disp south-moment")


# ---------------------------------------------------------------------------
# Neumann superposition — linearity across all edge combinations
# ---------------------------------------------------------------------------
# With qs = 0 and identical coefficient matrices (all edges normalise to
# zero_moment_zero_shear), the RHS correction is additive and the solution
# satisfies w(A + B) = w(A) + w(B) to floating-point precision.
# A failure here indicates a bug in the multi-edge Neumann corner treatment
# (cross-derivative stencil terms at adjacent prescribed edges).
# ---------------------------------------------------------------------------

_FREE_NS = "zero_moment_zero_shear"
_M0_SQ   = 1.0e12


class TestNeumannSuperposition:
    """w(A+B) == w(A) + w(B) to floating-point precision for moment BCs."""

    # --- corner pairs ---

    def test_west_north(self):
        w_w = _run_sq({"moment": _M0_SQ, "shear": 0.0}, _FREE_NS, _FREE_NS, _FREE_NS)
        w_n = _run_sq(_FREE_NS, _FREE_NS, {"moment": _M0_SQ, "shear": 0.0}, _FREE_NS)
        w   = _run_sq({"moment": _M0_SQ, "shear": 0.0}, _FREE_NS,
                       {"moment": _M0_SQ, "shear": 0.0}, _FREE_NS)
        np.testing.assert_allclose(w, w_w + w_n, rtol=1e-10, atol=0)

    def test_west_south(self):
        w_w = _run_sq({"moment": _M0_SQ, "shear": 0.0}, _FREE_NS, _FREE_NS, _FREE_NS)
        w_s = _run_sq(_FREE_NS, _FREE_NS, _FREE_NS, {"moment": _M0_SQ, "shear": 0.0})
        w   = _run_sq({"moment": _M0_SQ, "shear": 0.0}, _FREE_NS,
                       _FREE_NS, {"moment": _M0_SQ, "shear": 0.0})
        np.testing.assert_allclose(w, w_w + w_s, rtol=1e-10, atol=0)

    def test_east_north(self):
        w_e = _run_sq(_FREE_NS, {"moment": _M0_SQ, "shear": 0.0}, _FREE_NS, _FREE_NS)
        w_n = _run_sq(_FREE_NS, _FREE_NS, {"moment": _M0_SQ, "shear": 0.0}, _FREE_NS)
        w   = _run_sq(_FREE_NS, {"moment": _M0_SQ, "shear": 0.0},
                       {"moment": _M0_SQ, "shear": 0.0}, _FREE_NS)
        np.testing.assert_allclose(w, w_e + w_n, rtol=1e-10, atol=0)

    def test_east_south(self):
        w_e = _run_sq(_FREE_NS, {"moment": _M0_SQ, "shear": 0.0}, _FREE_NS, _FREE_NS)
        w_s = _run_sq(_FREE_NS, _FREE_NS, _FREE_NS, {"moment": _M0_SQ, "shear": 0.0})
        w   = _run_sq(_FREE_NS, {"moment": _M0_SQ, "shear": 0.0},
                       _FREE_NS, {"moment": _M0_SQ, "shear": 0.0})
        np.testing.assert_allclose(w, w_e + w_s, rtol=1e-10, atol=0)

    # --- opposite-edge pairs ---

    def test_west_east(self):
        w_w = _run_sq({"moment": _M0_SQ, "shear": 0.0}, _FREE_NS, _FREE_NS, _FREE_NS)
        w_e = _run_sq(_FREE_NS, {"moment": _M0_SQ, "shear": 0.0}, _FREE_NS, _FREE_NS)
        w   = _run_sq({"moment": _M0_SQ, "shear": 0.0}, {"moment": _M0_SQ, "shear": 0.0},
                       _FREE_NS, _FREE_NS)
        np.testing.assert_allclose(w, w_w + w_e, rtol=1e-10, atol=0)

    def test_north_south(self):
        w_n = _run_sq(_FREE_NS, _FREE_NS, {"moment": _M0_SQ, "shear": 0.0}, _FREE_NS)
        w_s = _run_sq(_FREE_NS, _FREE_NS, _FREE_NS, {"moment": _M0_SQ, "shear": 0.0})
        w   = _run_sq(_FREE_NS, _FREE_NS,
                       {"moment": _M0_SQ, "shear": 0.0}, {"moment": _M0_SQ, "shear": 0.0})
        np.testing.assert_allclose(w, w_n + w_s, rtol=1e-10, atol=0)

    # --- all four ---

    def test_all_four(self):
        w_w = _run_sq({"moment": _M0_SQ, "shear": 0.0}, _FREE_NS, _FREE_NS, _FREE_NS)
        w_e = _run_sq(_FREE_NS, {"moment": _M0_SQ, "shear": 0.0}, _FREE_NS, _FREE_NS)
        w_n = _run_sq(_FREE_NS, _FREE_NS, {"moment": _M0_SQ, "shear": 0.0}, _FREE_NS)
        w_s = _run_sq(_FREE_NS, _FREE_NS, _FREE_NS, {"moment": _M0_SQ, "shear": 0.0})
        w   = _run_sq(
            {"moment": _M0_SQ, "shear": 0.0}, {"moment": _M0_SQ, "shear": 0.0},
            {"moment": _M0_SQ, "shear": 0.0}, {"moment": _M0_SQ, "shear": 0.0},
        )
        np.testing.assert_allclose(w, w_w + w_e + w_n + w_s, rtol=1e-10, atol=0)
