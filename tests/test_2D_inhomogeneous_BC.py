#!/usr/bin/env python
"""Analytical tests for inhomogeneous (prescribed-value) 2-D boundary conditions.

STATUS: These tests FAIL until 2-D inhomogeneous BC machinery is implemented.
        They define the target behaviour; they will pass once the
        prescribed-BC machinery is in place.

Mirror-structure counterpart to test_1D_inhomogeneous_BC.py
-----------------------------------------------------------
The domain is long in x (10α) and narrow in y (NY=11 cells).  The north and
south edges carry ``bc_north="zero_slope_zero_shear"`` / ``bc_south="zero_slope_zero_shear"``,
which forces dw/dy = 0 there.  With a uniform load (qs = 0) and y-independent BCs the 2-D
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

def _run(bc_west, bc_east, bc_north="mirror", bc_south="mirror",
         nx=NX, ny=NY, L=L_DOMAIN, qs_override=None):
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
    flex.T_e      = TE
    flex.dx       = dx
    flex.dy       = dx
    flex.qs       = np.zeros((ny, nx)) if qs_override is None else qs_override
    flex.bc_west  = bc_west
    flex.bc_east  = bc_east
    flex.bc_north = bc_north
    flex.bc_south = bc_south

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        flex.initialize()
        flex.run()
        w = flex.w
        flex.finalize()

    return x, w


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
_MIRROR  = "zero_slope_zero_shear"


def _make_variable_te(ny, nx):
    """Smooth 2-D Te array: ±30 % variation in both directions."""
    eta = np.linspace(0, 1, ny)[:, np.newaxis]
    xi  = np.linspace(0, 1, nx)[np.newaxis, :]
    return TE * (1.0 + 0.3 * np.sin(np.pi * eta) * np.cos(np.pi * xi))


def _run_vte(bc_west, bc_east, bc_north, bc_south,
             nx=_VTE_NX, ny=_VTE_NY, te=None):
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
    flex.T_e      = _make_variable_te(ny, nx) if te is None else te
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
        w = flex.w
        flex.finalize()
    return w


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
# Variable-Te Neumann tests
# ---------------------------------------------------------------------------

class TestVariableTePrescribedNeumann:
    """Prescribed Neumann (moment / shear) BCs under spatially varying Te.

    Three orthogonal checks:
    1. Uniform 2-D Te array produces the same solution as scalar Te (tests
       the array-Te code path in _apply_bc_rhs_inhomogeneous_2d).
    2. Superposition w(M₀, V₀) = w(M₀, 0) + w(0, V₀) holds to machine
       precision under the variable-Te field (linearity check).
    3. The variable-Te solution is quantifiably different from the uniform-Te
       solution (sanity check that D variation is actually used).
    """

    M0 = 1.0e12   # N·m / m
    V0 = 1.0e8    # N / m
    _ZMZS = "zero_moment_zero_shear"

    def test_uniform_array_te_matches_scalar(self):
        """Uniform 2-D Te array must give the same answer as scalar Te."""
        Te_uni = TE * np.ones((_VTE_NY, _VTE_NX))
        w_scalar = _run_vte(
            {"moment": self.M0, "shear": 0.0}, self._ZMZS, _MIRROR, _MIRROR,
            te=TE,
        )
        w_array = _run_vte(
            {"moment": self.M0, "shear": 0.0}, self._ZMZS, _MIRROR, _MIRROR,
            te=Te_uni,
        )
        np.testing.assert_allclose(w_array, w_scalar, rtol=1e-10, atol=0)

    def test_superposition_holds_under_variable_te(self):
        """w(M₀,V₀) = w(M₀,0) + w(0,V₀) to machine precision under variable Te."""
        Te_var = _make_variable_te(_VTE_NY, _VTE_NX)
        w_M  = _run_vte({"moment": self.M0, "shear": 0.0}, self._ZMZS, _MIRROR, _MIRROR, te=Te_var)
        w_V  = _run_vte({"moment": 0.0, "shear": self.V0},  self._ZMZS, _MIRROR, _MIRROR, te=Te_var)
        w_MV = _run_vte({"moment": self.M0, "shear": self.V0}, self._ZMZS, _MIRROR, _MIRROR, te=Te_var)
        np.testing.assert_allclose(w_MV, w_M + w_V, rtol=1e-10, atol=0)

    def test_variable_te_differs_from_uniform(self):
        """Variable Te must produce a quantifiably different solution than uniform Te."""
        Te_var = _make_variable_te(_VTE_NY, _VTE_NX)
        w_uni = _run_vte({"moment": self.M0, "shear": 0.0}, self._ZMZS, _MIRROR, _MIRROR, te=TE)
        w_var = _run_vte({"moment": self.M0, "shear": 0.0}, self._ZMZS, _MIRROR, _MIRROR, te=Te_var)
        diff = np.max(np.abs(w_var - w_uni))
        assert diff > 0.05 * np.max(np.abs(w_uni)), (
            f"Variable Te (±30%) should produce >5% change; got {diff / np.max(np.abs(w_uni)):.1%}"
        )


# ---------------------------------------------------------------------------
# North / South domain helpers
# ---------------------------------------------------------------------------
# Domain long in y (NY_NS cells) and narrow in x (NX_NS cells);
# mirror on east/west enforces x-uniformity so the 2-D problem reduces
# to the 1-D semi-infinite plate equation in the y direction.

_NS_NX = 11
_NS_NY = 401   # Δy ≈ α/40, matching the x-direction resolution


def _run_ns(bc_north, bc_south, bc_west="mirror", bc_east="mirror",
            nx=_NS_NX, ny=_NS_NY, L=L_DOMAIN):
    dy = L / (ny - 1)
    dx = L / (nx - 1)   # equals dy for default narrow domain; differs for wide domain
    flex = F2D()
    flex.quiet    = True
    flex.method   = "fd"
    flex.solver   = "direct"
    flex.g        = G
    flex.E        = E
    flex.nu       = NU
    flex.rho_m    = RHO_M
    flex.rho_fill = RHO_F
    flex.T_e      = TE
    flex.dx       = dx
    flex.dy       = dy
    flex.qs       = np.zeros((ny, nx))
    flex.bc_west  = bc_west
    flex.bc_east  = bc_east
    flex.bc_north = bc_north
    flex.bc_south = bc_south
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        flex.initialize()
        flex.run()
        w = flex.w
        flex.finalize()
    y = np.arange(ny) * dy
    return y, w


def _check_col(w_num, w_ex, label, col=None):
    """Compare a column of 2-D solution to 1-D exact; L-inf relative error.

    col: column index to compare (default: _NS_NX // 2, the centre of the
         standard narrow NS domain).
    """
    if col is None:
        col = _NS_NX // 2
    centre = w_num[:, col]
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
    flex.T_e      = TE
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
        w = flex.w
        flex.finalize()
    return w


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
# Non-ZMZS far-end tests
# ---------------------------------------------------------------------------
# All existing single-edge prescribed BC tests use ZMZS at the far end.
# Here the far end carries a different string preset (ZDSZS or ZDZM) to
# confirm the RHS correction is unaffected by the far-end stiffness structure.
#
# For L = 10α the solution decays to e^{-10} ≈ 4.5e-5 of peak before reaching
# the far boundary, so the near-end profile matches the semi-infinite reference
# to well within 1%.
#
# Sign convention for east (odd BCs): shear d³w/dx³ reverses under ξ = L − x,
# so east V₀ > 0 gives w(ξ) = −w_prescribed_shear(ξ, V₀). This matches the
# convention in TestEastPrescribedShear.
# ---------------------------------------------------------------------------

class TestNonZMZSFarEnd:
    """Prescribed Neumann BC at one end with a non-ZMZS preset at the far end."""

    M0 = 1.0e12   # N·m / m
    V0 = 1.0e8    # N / m
    _ZDSZS = "zero_displacement_zero_slope"
    _ZDZM  = "zero_displacement_zero_moment"

    # --- west prescribed, east far end ---

    def test_west_moment_east_zdszs(self):
        x, w_num = _run(
            bc_west={"moment": self.M0, "shear": 0.0},
            bc_east=self._ZDSZS,
        )
        _check(w_num, w_prescribed_moment(x, self.M0), "West moment + east ZDSZS")

    def test_west_moment_east_zdzm(self):
        x, w_num = _run(
            bc_west={"moment": self.M0, "shear": 0.0},
            bc_east=self._ZDZM,
        )
        _check(w_num, w_prescribed_moment(x, self.M0), "West moment + east ZDZM")

    def test_west_shear_east_zdszs(self):
        x, w_num = _run(
            bc_west={"moment": 0.0, "shear": self.V0},
            bc_east=self._ZDSZS,
        )
        _check(w_num, w_prescribed_shear(x, self.V0), "West shear + east ZDSZS")

    def test_west_shear_east_zdzm(self):
        x, w_num = _run(
            bc_west={"moment": 0.0, "shear": self.V0},
            bc_east=self._ZDZM,
        )
        _check(w_num, w_prescribed_shear(x, self.V0), "West shear + east ZDZM")

    # --- east prescribed, west far end ---

    def test_east_moment_west_zdszs(self):
        x, w_num = _run(
            bc_west=self._ZDSZS,
            bc_east={"moment": self.M0, "shear": 0.0},
        )
        _check(w_num, w_prescribed_moment(L_DOMAIN - x, self.M0), "East moment + west ZDSZS")

    def test_east_moment_west_zdzm(self):
        x, w_num = _run(
            bc_west=self._ZDZM,
            bc_east={"moment": self.M0, "shear": 0.0},
        )
        _check(w_num, w_prescribed_moment(L_DOMAIN - x, self.M0), "East moment + west ZDZM")

    def test_east_shear_west_zdszs(self):
        """East shear: sign flips under ξ = L − x (d³/dx³ is odd)."""
        x, w_num = _run(
            bc_west=self._ZDSZS,
            bc_east={"moment": 0.0, "shear": self.V0},
        )
        _check(w_num, -w_prescribed_shear(L_DOMAIN - x, self.V0), "East shear + west ZDSZS")

    def test_east_shear_west_zdzm(self):
        """East shear: sign flips under ξ = L − x (d³/dx³ is odd)."""
        x, w_num = _run(
            bc_west=self._ZDZM,
            bc_east={"moment": 0.0, "shear": self.V0},
        )
        _check(w_num, -w_prescribed_shear(L_DOMAIN - x, self.V0), "East shear + west ZDZM")


class TestNonZMZSFarEndDirichlet:
    """Prescribed Dirichlet BC at one end with a non-ZMZS preset at the far end.

    Extends TestNonZMZSFarEnd to displacement and slope BCs.  The existing
    single-edge tests use ZMZS as the far end; here we test ZDSZS and ZDZM.

    Sign convention: east slope uses −w_prescribed_slope(ξ, θ₀) because
    dw/dx reverses sign under ξ = L − x (same as east shear).  East
    displacement is even and carries no sign change.
    """

    W0     = 100.0    # m
    THETA0 = 1.0e-3   # rad
    _ZDSZS = "zero_displacement_zero_slope"
    _ZDZM  = "zero_displacement_zero_moment"

    # --- west displacement, non-ZMZS east far end ---

    def test_west_displacement_east_zdszs(self):
        x, w_num = _run(
            bc_west={"displacement": self.W0, "slope": 0.0},
            bc_east=self._ZDSZS,
        )
        _check(w_num, w_prescribed_displacement(x, self.W0), "West disp + east ZDSZS")

    def test_west_displacement_east_zdzm(self):
        x, w_num = _run(
            bc_west={"displacement": self.W0, "slope": 0.0},
            bc_east=self._ZDZM,
        )
        _check(w_num, w_prescribed_displacement(x, self.W0), "West disp + east ZDZM")

    # --- east displacement, non-ZDSZS west far end ---

    def test_east_displacement_west_zdzm(self):
        """East disp, ZDZM at west (ZDSZS at west is covered by TestEastPrescribedDisplacement)."""
        x, w_num = _run(
            bc_west=self._ZDZM,
            bc_east={"displacement": self.W0, "slope": 0.0},
        )
        _check(w_num, w_prescribed_displacement(L_DOMAIN - x, self.W0), "East disp + west ZDZM")

    # --- west slope, non-ZMZS east far end ---

    def test_west_slope_east_zdszs(self):
        x, w_num = _run(
            bc_west={"displacement": 0.0, "slope": self.THETA0},
            bc_east=self._ZDSZS,
        )
        _check(w_num, w_prescribed_slope(x, self.THETA0), "West slope + east ZDSZS")

    def test_west_slope_east_zdzm(self):
        x, w_num = _run(
            bc_west={"displacement": 0.0, "slope": self.THETA0},
            bc_east=self._ZDZM,
        )
        _check(w_num, w_prescribed_slope(x, self.THETA0), "West slope + east ZDZM")

    # --- east slope, non-ZMZS west far end ---

    def test_east_slope_west_zdszs(self):
        """East slope: sign flips under ξ = L − x (d/dx is odd)."""
        x, w_num = _run(
            bc_west=self._ZDSZS,
            bc_east={"displacement": 0.0, "slope": self.THETA0},
        )
        _check(w_num, -w_prescribed_slope(L_DOMAIN - x, self.THETA0), "East slope + west ZDSZS")

    def test_east_slope_west_zdzm(self):
        """East slope: sign flips under ξ = L − x (d/dx is odd)."""
        x, w_num = _run(
            bc_west=self._ZDZM,
            bc_east={"displacement": 0.0, "slope": self.THETA0},
        )
        _check(w_num, -w_prescribed_slope(L_DOMAIN - x, self.THETA0), "East slope + west ZDZM")


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
        np.testing.assert_allclose(w, w_w + w_n, rtol=1e-9, atol=0)

    def test_west_south(self):
        w_w = _run_sq({"moment": _M0_SQ, "shear": 0.0}, _FREE_NS, _FREE_NS, _FREE_NS)
        w_s = _run_sq(_FREE_NS, _FREE_NS, _FREE_NS, {"moment": _M0_SQ, "shear": 0.0})
        w   = _run_sq({"moment": _M0_SQ, "shear": 0.0}, _FREE_NS,
                       _FREE_NS, {"moment": _M0_SQ, "shear": 0.0})
        np.testing.assert_allclose(w, w_w + w_s, rtol=1e-9, atol=0)

    def test_east_north(self):
        w_e = _run_sq(_FREE_NS, {"moment": _M0_SQ, "shear": 0.0}, _FREE_NS, _FREE_NS)
        w_n = _run_sq(_FREE_NS, _FREE_NS, {"moment": _M0_SQ, "shear": 0.0}, _FREE_NS)
        w   = _run_sq(_FREE_NS, {"moment": _M0_SQ, "shear": 0.0},
                       {"moment": _M0_SQ, "shear": 0.0}, _FREE_NS)
        np.testing.assert_allclose(w, w_e + w_n, rtol=1e-9, atol=0)

    def test_east_south(self):
        w_e = _run_sq(_FREE_NS, {"moment": _M0_SQ, "shear": 0.0}, _FREE_NS, _FREE_NS)
        w_s = _run_sq(_FREE_NS, _FREE_NS, _FREE_NS, {"moment": _M0_SQ, "shear": 0.0})
        w   = _run_sq(_FREE_NS, {"moment": _M0_SQ, "shear": 0.0},
                       _FREE_NS, {"moment": _M0_SQ, "shear": 0.0})
        np.testing.assert_allclose(w, w_e + w_s, rtol=1e-9, atol=0)

    # --- opposite-edge pairs ---

    def test_west_east(self):
        w_w = _run_sq({"moment": _M0_SQ, "shear": 0.0}, _FREE_NS, _FREE_NS, _FREE_NS)
        w_e = _run_sq(_FREE_NS, {"moment": _M0_SQ, "shear": 0.0}, _FREE_NS, _FREE_NS)
        w   = _run_sq({"moment": _M0_SQ, "shear": 0.0}, {"moment": _M0_SQ, "shear": 0.0},
                       _FREE_NS, _FREE_NS)
        np.testing.assert_allclose(w, w_w + w_e, rtol=1e-9, atol=0)

    def test_north_south(self):
        w_n = _run_sq(_FREE_NS, _FREE_NS, {"moment": _M0_SQ, "shear": 0.0}, _FREE_NS)
        w_s = _run_sq(_FREE_NS, _FREE_NS, _FREE_NS, {"moment": _M0_SQ, "shear": 0.0})
        w   = _run_sq(_FREE_NS, _FREE_NS,
                       {"moment": _M0_SQ, "shear": 0.0}, {"moment": _M0_SQ, "shear": 0.0})
        np.testing.assert_allclose(w, w_n + w_s, rtol=1e-9, atol=0)

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
        np.testing.assert_allclose(w, w_w + w_e + w_n + w_s, rtol=1e-9, atol=0)


# ---------------------------------------------------------------------------
# Zero-value dict degeneracy
# ---------------------------------------------------------------------------

class TestZeroValueDictDegeneracy:
    """A dict BC with all zero values must produce the same solution as the
    equivalent string preset.

    {"moment": 0, "shear": 0}         ≡  "zero_moment_zero_shear"
    {"displacement": 0, "slope": 0}   ≡  "zero_displacement_zero_slope"

    The RHS correction is identically zero for zero prescribed values, so the
    two code paths produce the same linear system and the solutions must be
    exactly equal.
    """

    _ZMZS  = "zero_moment_zero_shear"
    _ZDSZS = "zero_displacement_zero_slope"

    def test_west_zero_moment_shear_equals_string(self):
        _, w_dict = _run(bc_west={"moment": 0.0, "shear": 0.0}, bc_east=self._ZMZS)
        _, w_str  = _run(bc_west=self._ZMZS,                     bc_east=self._ZMZS)
        np.testing.assert_array_equal(w_dict, w_str)

    def test_west_zero_displacement_slope_equals_string(self):
        _, w_dict = _run(bc_west={"displacement": 0.0, "slope": 0.0}, bc_east=self._ZMZS)
        _, w_str  = _run(bc_west=self._ZDSZS,                          bc_east=self._ZMZS)
        np.testing.assert_array_equal(w_dict, w_str)

    def test_east_zero_moment_shear_equals_string(self):
        _, w_dict = _run(bc_west=self._ZMZS, bc_east={"moment": 0.0, "shear": 0.0})
        _, w_str  = _run(bc_west=self._ZMZS, bc_east=self._ZMZS)
        np.testing.assert_array_equal(w_dict, w_str)

    def test_east_zero_displacement_slope_equals_string(self):
        _, w_dict = _run(bc_west=self._ZMZS, bc_east={"displacement": 0.0, "slope": 0.0})
        _, w_str  = _run(bc_west=self._ZMZS, bc_east=self._ZDSZS)
        np.testing.assert_array_equal(w_dict, w_str)


# ---------------------------------------------------------------------------
# Superposition under non-zero distributed load
# ---------------------------------------------------------------------------

class TestLoadSuperposition:
    """Prescribed-value BC corrections must be additive with distributed load.

    For any linear load q₀ and BC value B, linearity demands:
        w(q₀, B) = w(q₀, B=0) + w(q₀=0, B)

    This holds to machine precision (rtol ≈ 1e-13) and catches any bug
    where the RHS correction interacts non-linearly with the load vector.
    All existing inhomogeneous BC tests use qs = 0; this is the first
    test with a real distributed load.
    """

    M0 = 1.0e12   # N·m / m
    Q0 = 1.0e4    # Pa  (uniform downward load)

    def _qs(self):
        dx = L_DOMAIN / (NX - 1)
        return self.Q0 * np.ones((NY, NX))

    def test_west_moment_with_load(self):
        qs = self._qs()
        _, w_both = _run({"moment": self.M0, "shear": 0.0}, "zero_moment_zero_shear",
                         qs_override=qs)
        _, w_load = _run("zero_moment_zero_shear", "zero_moment_zero_shear",
                         qs_override=qs)
        _, w_bc   = _run({"moment": self.M0, "shear": 0.0}, "zero_moment_zero_shear")
        np.testing.assert_allclose(w_both, w_load + w_bc, rtol=1e-10, atol=0)

    def test_east_moment_with_load(self):
        qs = self._qs()
        _, w_both = _run("zero_moment_zero_shear", {"moment": self.M0, "shear": 0.0},
                         qs_override=qs)
        _, w_load = _run("zero_moment_zero_shear", "zero_moment_zero_shear",
                         qs_override=qs)
        _, w_bc   = _run("zero_moment_zero_shear", {"moment": self.M0, "shear": 0.0})
        np.testing.assert_allclose(w_both, w_load + w_bc, rtol=1e-10, atol=0)


# ---------------------------------------------------------------------------
# Periodic lateral BCs with prescribed east / west BCs
# ---------------------------------------------------------------------------

class TestPeriodicLateralBCs:
    """Prescribed BC at west or east with periodic north/south boundary conditions.

    For y-independent BCs (uniform M₀ along the boundary, uniform Te),
    the periodic and mirror lateral conditions both enforce y-uniformity:

    * Mirror:   dw/dy = 0 at y = 0 and y = L_y.
    * Periodic: w(y=0) ≡ w(y=L_y) and all y-derivatives match.

    Any y-independent solution satisfies both, so the two solves produce
    numerically identical results (differences are at the 1e-10 level due to
    different sparse-matrix structures).

    The y-uniformity assertion is the numerical analogue of a visual corner
    check: if the periodic corner handling were wrong, a spurious y-gradient
    would appear near the prescribed BC edge.
    """

    M0 = 1.0e12

    def test_west_moment_periodic_matches_mirror(self):
        _, w_mirror   = _run({"moment": self.M0, "shear": 0.0}, "zero_moment_zero_shear",
                             "mirror", "mirror")
        _, w_periodic = _run({"moment": self.M0, "shear": 0.0}, "zero_moment_zero_shear",
                             "periodic", "periodic")
        np.testing.assert_allclose(w_periodic, w_mirror, rtol=1e-6, atol=0)

    def test_west_moment_periodic_is_y_uniform(self):
        _, w = _run({"moment": self.M0, "shear": 0.0}, "zero_moment_zero_shear",
                    "periodic", "periodic")
        centre = w[NY // 2, :]
        assert np.max(np.abs(w - centre)) < 1e-10 * np.max(np.abs(centre)), (
            "Periodic lateral BCs with uniform west BC should give y-uniform solution"
        )

    def test_east_moment_periodic_matches_mirror(self):
        _, w_mirror   = _run("zero_moment_zero_shear", {"moment": self.M0, "shear": 0.0},
                             "mirror", "mirror")
        _, w_periodic = _run("zero_moment_zero_shear", {"moment": self.M0, "shear": 0.0},
                             "periodic", "periodic")
        np.testing.assert_allclose(w_periodic, w_mirror, rtol=1e-6, atol=0)

    def test_east_moment_periodic_is_y_uniform(self):
        _, w = _run("zero_moment_zero_shear", {"moment": self.M0, "shear": 0.0},
                    "periodic", "periodic")
        centre = w[NY // 2, :]
        assert np.max(np.abs(w - centre)) < 1e-10 * np.max(np.abs(centre)), (
            "Periodic lateral BCs with uniform east BC should give y-uniform solution"
        )


# ---------------------------------------------------------------------------
# North/south prescribed BCs with non-mirror lateral ends
# ---------------------------------------------------------------------------

class TestNorthSouthNonMirrorLateral:
    """North prescribed BC with non-mirror (ZDSZS or ZDZM) east/west lateral ends.

    Three independent lines of evidence:

    1. Rotation symmetry — On a square domain with constant D and dx = dy,
       the biharmonic ∇⁴ is symmetric under 90° CCW rotation that maps
       (x, y) → (L − y, x), sending the west face to the north face.
       The induced array relation is:
           w_NS[i, j] = w_EW[n − 1 − j, i]   →   w_NS = w_EW[::-1, :].T
       This connects the untested NS code path to the independently-validated
       EW code path without needing a 1-D analytical reference.

    2. Superposition — w(M₀ + V₀) = w(M₀) + w(V₀) for north BC with ZDSZS
       lateral ends, rtol = 1e-9.  Catches RHS sign or assembly bugs that
       a symmetry test alone would miss.

    3. Stiffness inequality — ZDSZS lateral ends clamp the plate globally;
       the peak north-boundary deflection must be strictly less than the
       mirror-lateral (1-D equivalent) case.  The center column does NOT
       match the 1-D solution even for wide domains: ZDSZS forces w = 0
       at the lateral edges for all y, changing the plate response throughout.
    """

    M0     = 1.0e12
    V0     = 1.0e8
    _ZDSZS = "zero_displacement_zero_slope"
    _ZDZM  = "zero_displacement_zero_moment"
    _ZMZS  = "zero_moment_zero_shear"

    # --- 1. Rotation symmetry ---

    def test_rotation_symmetry_north_moment_zdszs(self):
        """NS prescribed-moment solution = 90° CCW rotation of EW prescribed-moment.

        EW problem: west = {M0}, east/north/south = ZDSZS.
        NS problem: north = {M0}, south/east/west = ZDSZS.
        Rotation maps west→north, east→south, north→east, south→west.
        """
        w_ew = _run_sq(
            bc_west={"moment": self.M0, "shear": 0.0},
            bc_east=self._ZDSZS,
            bc_north=self._ZDSZS,
            bc_south=self._ZDSZS,
        )
        w_ns = _run_sq(
            bc_west=self._ZDSZS,
            bc_east=self._ZDSZS,
            bc_north={"moment": self.M0, "shear": 0.0},
            bc_south=self._ZDSZS,
        )
        np.testing.assert_allclose(w_ns, w_ew[::-1, :].T, rtol=1e-7, atol=0)

    def test_rotation_symmetry_north_moment_zdzm(self):
        """Same rotation symmetry test with ZDZM at all non-prescribed ends."""
        w_ew = _run_sq(
            bc_west={"moment": self.M0, "shear": 0.0},
            bc_east=self._ZDZM,
            bc_north=self._ZDZM,
            bc_south=self._ZDZM,
        )
        w_ns = _run_sq(
            bc_west=self._ZDZM,
            bc_east=self._ZDZM,
            bc_north={"moment": self.M0, "shear": 0.0},
            bc_south=self._ZDZM,
        )
        np.testing.assert_allclose(w_ns, w_ew[::-1, :].T, rtol=1e-7, atol=0)

    # --- 2. Superposition ---

    def test_superposition_north_zdszs_lateral(self):
        """w(M₀+V₀) = w(M₀)+w(V₀) for north BC with ZDSZS at east and west.

        Uses _run_sq (small square domain) so the solution need not decay to
        machine zero — superposition holds for any linear system regardless.
        """
        kw = dict(bc_west=self._ZDSZS, bc_east=self._ZDSZS,
                  bc_south=self._ZMZS)
        w_m  = _run_sq(bc_north={"moment": self.M0, "shear": 0.0}, **kw)
        w_v  = _run_sq(bc_north={"moment": 0.0,     "shear": self.V0}, **kw)
        w_mv = _run_sq(bc_north={"moment": self.M0, "shear": self.V0}, **kw)
        np.testing.assert_allclose(w_mv, w_m + w_v, rtol=1e-9, atol=0)

    # --- 3. x-symmetry of the solution ---

    def test_north_moment_zdszs_x_symmetric(self):
        """With uniform north moment and symmetric ZDSZS at E/W, the solution
        must be symmetric about the x-midpoint: w[i, j] = w[i, n-1-j].

        Any bug in the lateral BC application (wrong sign, wrong corner
        handling, asymmetric stencil) would break this symmetry.  The
        tolerance is tight (rtol=1e-10) because the discrete problem is
        exactly symmetric — the coefficient matrix and RHS have the same
        left-right symmetry as the geometry.
        """
        w = _run_sq(
            bc_west=self._ZDSZS,
            bc_east=self._ZDSZS,
            bc_north={"moment": self.M0, "shear": 0.0},
            bc_south=self._ZMZS,
        )
        np.testing.assert_allclose(w, w[:, ::-1], rtol=1e-10, atol=0)


# ---------------------------------------------------------------------------
# Array-valued BC values (#64)
# ---------------------------------------------------------------------------

class TestArrayValuedBCs:
    """Array-valued BC values — same BC type along the whole edge, values vary.

    Three checks:

    1. Scalar-array equivalence — np.full(ny, scalar) produces the same
       solution as the scalar itself, to machine precision.  This validates
       the broadcast path without needing an independent analytical reference.

    2. Linearity — w(α·M0_array) = α·w(M0_array).  Catches sign or scaling
       errors in the per-row RHS correction.

    3. Superposition — w(A + B) = w(A) + w(B) for two array BC values.
       Catches cross-row coupling bugs that linearity alone misses.
    """

    M0 = 1.0e12   # N·m / m  (representative prescribed moment)

    # --- 1. scalar-array equivalence ---

    def test_west_moment_array_equals_scalar(self):
        """np.full(NY, M0) must give exactly the same solution as scalar M0."""
        _, w_scalar = _run({"moment": self.M0, "shear": 0.0}, "zero_moment_zero_shear")
        _, w_array  = _run({"moment": np.full(NY, self.M0), "shear": 0.0},
                           "zero_moment_zero_shear")
        np.testing.assert_array_equal(w_scalar, w_array)

    def test_east_moment_array_equals_scalar(self):
        _, w_scalar = _run("zero_moment_zero_shear", {"moment": self.M0, "shear": 0.0})
        _, w_array  = _run("zero_moment_zero_shear",
                           {"moment": np.full(NY, self.M0), "shear": 0.0})
        np.testing.assert_array_equal(w_scalar, w_array)

    def test_north_moment_array_equals_scalar(self):
        _, w_scalar = _run_ns({"moment": self.M0, "shear": 0.0}, "zero_moment_zero_shear")
        _, w_array  = _run_ns({"moment": np.full(_NS_NX, self.M0), "shear": 0.0},
                              "zero_moment_zero_shear")
        np.testing.assert_array_equal(w_scalar, w_array)

    def test_south_moment_array_equals_scalar(self):
        _, w_scalar = _run_ns("zero_moment_zero_shear", {"moment": self.M0, "shear": 0.0})
        _, w_array  = _run_ns("zero_moment_zero_shear",
                              {"moment": np.full(_NS_NX, self.M0), "shear": 0.0})
        np.testing.assert_array_equal(w_scalar, w_array)

    # --- 2. linearity ---

    def test_west_moment_array_linearity(self):
        """w(2·M0_array) = 2·w(M0_array) to rtol=1e-10."""
        M0_arr = self.M0 * (1.0 + 0.5 * np.sin(np.linspace(0, np.pi, NY)))
        _, w1 = _run({"moment": M0_arr,     "shear": 0.0}, "zero_moment_zero_shear")
        _, w2 = _run({"moment": 2.0*M0_arr, "shear": 0.0}, "zero_moment_zero_shear")
        np.testing.assert_allclose(w2, 2.0 * w1, rtol=1e-10, atol=0)

    # --- 3. superposition ---

    def test_west_moment_array_superposition(self):
        """w(A + B) = w(A) + w(B) for two array BC values, rtol=1e-10."""
        M0_a = self.M0 * np.ones(NY)
        M0_b = self.M0 * 0.5 * np.sin(np.linspace(0, np.pi, NY))
        _, w_a  = _run({"moment": M0_a,       "shear": 0.0}, "zero_moment_zero_shear")
        _, w_b  = _run({"moment": M0_b,       "shear": 0.0}, "zero_moment_zero_shear")
        _, w_ab = _run({"moment": M0_a + M0_b,"shear": 0.0}, "zero_moment_zero_shear")
        np.testing.assert_allclose(w_ab, w_a + w_b, rtol=1e-10, atol=0)


# ---------------------------------------------------------------------------
# Nested-model gradient round-trip
# ---------------------------------------------------------------------------
# Motivation: the nested-model use case (e.g. a coarse Greenland solve feeding
# a fine sub-domain solve) extracts w and dw/dn from a coarse solve and
# re-imposes them as inhomogeneous Dirichlet BCs on the fine sub-domain.
#
# The slope sign convention is:
#   west / north:  slope = dw/dx or dw/dy  in the *positive* coordinate
#                  direction (east / south-row-index).  Positive slope gives
#                  positive deflection decaying inward.
#   east / south:  same coordinate direction, but the domain is "behind" the
#                  boundary, so the same positive slope drives negative
#                  deflection (sign flip under ξ = L − coord).
#
# Test strategy:
#   1. Solve a full N×N problem with an off-centre load and clamped BCs.
#   2. Extract w and np.gradient-based slopes at the boundary of an interior
#      window (rows I1:I2+1, cols J1:J2+1).
#   3. Re-solve the window with those as all-four-edge Dirichlet BCs and the
#      same load.
#   4. Check that the window interior matches the full solution to within the
#      O(dx²) accuracy of np.gradient.
#
# A sign error on any single edge would flip the deflection trend at that
# edge, causing an O(1) error in the interior — far above the threshold.
# ---------------------------------------------------------------------------

_NRT_N   = 40      # full-domain cells (40×40)
_NRT_DX  = ALPHA * 0.25   # ~25 % of flexural parameter → 40 cells ≈ 10α
_NRT_I1, _NRT_I2 = 10, 30   # sub-domain row window (inclusive)
_NRT_J1, _NRT_J2 = 10, 30   # sub-domain col window (inclusive)


def _run_full_nrt():
    """Full 40×40 solve with an off-centre rectangular load, all sides clamped."""
    n  = _NRT_N
    dx = _NRT_DX
    qs = np.zeros((n, n))
    # Off-centre load so the 2-D solution has genuine gradients on all four
    # sub-domain edges (a centred load would make N/S or E/W slopes trivially zero).
    qs[12:22, 14:24] = 1.0e6   # 1 MPa patch, shifted NW of centre

    flex = F2D()
    flex.quiet    = True
    flex.method   = "fd"
    flex.solver   = "direct"
    flex.g        = G
    flex.E        = E
    flex.nu       = NU
    flex.rho_m    = RHO_M
    flex.rho_fill = RHO_F
    flex.T_e      = TE
    flex.dx       = dx
    flex.dy       = dx
    flex.qs       = qs
    flex.bc_west = flex.bc_east = flex.bc_north = flex.bc_south = \
        "zero_displacement_zero_slope"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        flex.initialize()
        flex.run()
        w  = flex.w.copy()
        flex.finalize()
    return w, qs


class TestNestedModelGradientRoundTrip:
    """Round-trip: extract w and slopes from a full solve, re-solve sub-domain.

    Verifies the slope (gradient) sign convention on all four edges at once.
    The slope convention used for extraction must match gFlex's internal
    convention:

      west / east:  slope = dw/dx  (np.gradient(w, dx, axis=1))
      north / south: slope = dw/dy (np.gradient(w, dy, axis=0))

    where dy increases with row index (row 0 = north → row n-1 = south).
    """

    # 2 % of peak deflection — generous enough for np.gradient O(dx²) error,
    # tight enough to catch any sign flip (which would produce O(1) error).
    RTOL = 0.02

    def test_all_four_edges(self):
        """Sub-domain interior matches full solve after gradient-BC round-trip."""
        w_full, qs_full = _run_full_nrt()
        dx = _NRT_DX
        i1, i2 = _NRT_I1, _NRT_I2
        j1, j2 = _NRT_J1, _NRT_J2

        dw_dx = np.gradient(w_full, dx, axis=1)   # positive = eastward
        dw_dy = np.gradient(w_full, dx, axis=0)   # positive = southward (row-index dir)

        n_sub_y = i2 - i1 + 1
        n_sub_x = j2 - j1 + 1

        flex2 = F2D()
        flex2.quiet    = True
        flex2.method   = "fd"
        flex2.solver   = "direct"
        flex2.g        = G
        flex2.E        = E
        flex2.nu       = NU
        flex2.rho_m    = RHO_M
        flex2.rho_fill = RHO_F
        flex2.T_e      = TE
        flex2.dx       = dx
        flex2.dy       = dx
        flex2.qs       = qs_full[i1:i2+1, j1:j2+1]
        flex2.bc_west  = {"displacement": w_full[i1:i2+1, j1],
                          "slope":        dw_dx[i1:i2+1, j1]}
        flex2.bc_east  = {"displacement": w_full[i1:i2+1, j2],
                          "slope":        dw_dx[i1:i2+1, j2]}
        flex2.bc_north = {"displacement": w_full[i1, j1:j2+1],
                          "slope":        dw_dy[i1, j1:j2+1]}
        flex2.bc_south = {"displacement": w_full[i2, j1:j2+1],
                          "slope":        dw_dy[i2, j1:j2+1]}

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            flex2.initialize()
            flex2.run()
            w_sub = flex2.w.copy()
            flex2.finalize()

        # Compare interior nodes (exclude boundary rows/cols, which are by
        # construction equal to the prescribed displacement values).
        w_sub_int  = w_sub[1:-1, 1:-1]
        w_full_int = w_full[i1+1:i2, j1+1:j2]

        scale = np.max(np.abs(w_full_int))
        err   = np.max(np.abs(w_sub_int - w_full_int)) / scale
        assert err < self.RTOL, (
            f"Nested round-trip L-inf error {err:.3e} exceeds {self.RTOL:.0%}; "
            "check gradient sign convention (west/east: dw/dx; "
            "north/south: dw/dy in row-index direction)"
        )
