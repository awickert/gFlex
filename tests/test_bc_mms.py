#!/usr/bin/env python3
"""BC-specific MMS convergence and Dirichlet-exactness tests.

Covers boundary conditions that were either recently corrected or newly
introduced, for which the existing periodic-MMS tests in test_1D_FD and
test_2D_FD provide no signal, plus variable-Te FD convergence.

BCs tested
----------
zero_displacement_zero_slope  (clamped)  — corrected 2026; convergence + Dirichlet exactness
zero_moment_zero_shear        (free end) — corrected 2026; convergence
zero_displacement_zero_moment (pinned)   — new in 2026; convergence
zero_slope_zero_shear         (mirror)   — convergence; MMS uses w=cos(πx/L)

Variable Te
-----------
variable_te_1d  — linearly-varying D; clamped BCs; MMS load derived analytically
"""

import warnings

import numpy as np

from gflex.f1d import F1D
from gflex.f2d import F2D

# ---------------------------------------------------------------------------
# Shared physical parameters
# ---------------------------------------------------------------------------

E        = 65.0e9
NU       = 0.25
TE       = 30.0e3
RHO_M    = 3300.0
RHO_F    = 0.0
G        = 9.8
L        = 600.0e3   # domain length (both sides)

D      = E * TE**3 / (12.0 * (1.0 - NU**2))
DRHO_G = (RHO_M - RHO_F) * G

# Grid sizes: small enough for CI, enough levels for a reliable slope estimate
NX_1D = [50, 100, 200]
N_2D  = [25, 50, 100]


# ---------------------------------------------------------------------------
# MMS shape functions
# ---------------------------------------------------------------------------

def _g(t):
    """g(t) = t²(1-t)²  — satisfies w=0, dw/dt=0 at t=0,1 (clamped)."""
    return t**2 * (1 - t)**2

def _g2(t):
    """g''(t) = 2 - 12t + 12t²."""
    return 2 - 12*t + 12*t**2

def _f(t):
    """f(t) = t⁴(1-t)⁴  — satisfies w=w'=w''=w'''=0 at t=0,1 (free end)."""
    return t**4 * (1 - t)**4

def _f2(t):
    """f''(t)."""
    return 12*t**2 - 80*t**3 + 180*t**4 - 168*t**5 + 56*t**6

def _f4(t):
    """f''''(t)."""
    return 24 - 480*t + 2160*t**2 - 3360*t**3 + 1680*t**4

def _m(t):
    """m(t) = cos(π t) — satisfies dw/dt=0, d³w/dt³=0 at t=0,1 (mirror BC)."""
    return np.cos(np.pi * t)


# ---------------------------------------------------------------------------
# 1-D MMS load builders
# ---------------------------------------------------------------------------

def _mms_1d_clamped(nx):
    dx = L / (nx - 1)
    xi = np.arange(nx) / (nx - 1)
    w_ex = -_g(xi)                         # W0 = 1
    qs   = D * 24.0 / L**4 + DRHO_G * _g(xi)
    return dx, w_ex, qs

def _mms_1d_free_end(nx):
    dx = L / (nx - 1)
    xi = np.arange(nx) / (nx - 1)
    w_ex = -_f(xi)
    qs   = D / L**4 * _f4(xi) + DRHO_G * _f(xi)
    return dx, w_ex, qs

def _mms_1d_pinned(nx):
    dx = L / (nx - 1)
    xi = np.arange(nx) / (nx - 1)
    w_ex = -np.sin(np.pi * xi)
    qs   = (D * (np.pi / L)**4 + DRHO_G) * np.sin(np.pi * xi)
    return dx, w_ex, qs

def _mms_1d_mirror(nx):
    """Mirror BC MMS: w_exact = -cos(π ξ), satisfies dw/dx=d³w/dx³=0 at both ends."""
    dx   = L / (nx - 1)
    xi   = np.arange(nx) / (nx - 1)
    k    = np.pi / L
    w_ex = -_m(xi)
    qs   = (D * k**4 + DRHO_G) * _m(xi)
    return dx, w_ex, qs

# Fractional D variation for the variable-Te MMS: D(ξ) = D*(1 + _VAR_A*ξ).
_VAR_A = 0.5

def _mms_1d_variable_te(nx):
    """Variable-Te MMS (clamped BCs): D(ξ) = D*(1+0.5ξ), w_exact = -ξ²(1-ξ)².

    The manufactured load is derived analytically from
    d²/dx²[D(x) d²w/dx²] + Δρg·w = −q, with D linear in ξ = x/L.  For
    D(ξ) = D₀(1+aξ) and w = -g(ξ):

        q(ξ) = (D₀/L⁴)·(24(1−a) + 72aξ) + Δρg·g(ξ)
    """
    dx   = L / (nx - 1)
    xi   = np.arange(nx) / (nx - 1)
    a    = _VAR_A
    te   = TE * (1.0 + a * xi) ** (1.0/3.0)   # T_e(ξ): D(ξ) = D*(1+a*ξ)
    w_ex = -_g(xi)
    qs   = (D / L**4) * (24*(1-a) + 72*a*xi) + DRHO_G * _g(xi)
    return dx, te, w_ex, qs


# ---------------------------------------------------------------------------
# 2-D MMS load builders  (square domain, rows = y/eta, cols = x/xi)
# ---------------------------------------------------------------------------

def _mms_2d_clamped(n):
    dx  = L / (n - 1)
    eta = np.arange(n)[:, np.newaxis] / (n - 1)   # (n,1)
    xi  = np.arange(n)[np.newaxis, :] / (n - 1)   # (1,n)
    w_ex = -_g(eta) * _g(xi)
    qs   = (D / L**4 * (24*_g(eta) + 2*_g2(eta)*_g2(xi) + 24*_g(xi))
            + DRHO_G * _g(eta) * _g(xi))
    return dx, w_ex, qs

def _mms_2d_free_end(n):
    dx  = L / (n - 1)
    eta = np.arange(n)[:, np.newaxis] / (n - 1)
    xi  = np.arange(n)[np.newaxis, :] / (n - 1)
    w_ex = -_f(eta) * _f(xi)
    qs   = (D / L**4 * (_f4(eta)*_f(xi) + 2*_f2(eta)*_f2(xi) + _f(eta)*_f4(xi))
            + DRHO_G * _f(eta) * _f(xi))
    return dx, w_ex, qs

def _mms_2d_pinned(n):
    dx  = L / (n - 1)
    eta = np.arange(n)[:, np.newaxis] / (n - 1)
    xi  = np.arange(n)[np.newaxis, :] / (n - 1)
    w_ex = -np.sin(np.pi * eta) * np.sin(np.pi * xi)
    qs   = (4*D*(np.pi/L)**4 + DRHO_G) * np.sin(np.pi * eta) * np.sin(np.pi * xi)
    return dx, w_ex, qs

def _mms_2d_mirror(n):
    """Mirror BC MMS 2-D: w = -cos(πη)cos(πξ), all four sides are mirror."""
    dx  = L / (n - 1)
    eta = np.arange(n)[:, np.newaxis] / (n - 1)
    xi  = np.arange(n)[np.newaxis, :] / (n - 1)
    k   = np.pi / L
    w_ex = -_m(eta) * _m(xi)
    qs   = (4*D * k**4 + DRHO_G) * _m(eta) * _m(xi)
    return dx, w_ex, qs


# ---------------------------------------------------------------------------
# gFlex runners
# ---------------------------------------------------------------------------

def _run_1d(dx, qs, bc):
    s = F1D()
    s.quiet    = True
    s.verbose  = False
    s.debug    = False
    s.method   = "fd"
    s.g        = G
    s.E        = E
    s.nu       = NU
    s.rho_m    = RHO_M
    s.rho_fill = RHO_F
    s.T_e      = TE
    s.dx       = dx
    s.qs       = qs
    s.bc_west  = bc
    s.bc_east  = bc
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s.initialize()
        s.run()
        w = s.w
        s.finalize()
    return w


def _run_1d_te(dx, te, qs, bc):
    """Like _run_1d but accepts a Te array for variable-rigidity tests."""
    s = F1D()
    s.quiet    = True
    s.verbose  = False
    s.debug    = False
    s.method   = "fd"
    s.g        = G
    s.E        = E
    s.nu       = NU
    s.rho_m    = RHO_M
    s.rho_fill = RHO_F
    s.T_e      = te
    s.dx       = dx
    s.qs       = qs
    s.bc_west  = bc
    s.bc_east  = bc
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s.initialize()
        s.run()
        w = s.w
        s.finalize()
    return w


def _run_2d(dx, qs, bc):
    s = F2D()
    s.quiet    = True
    s.method   = "fd"
    s.solver   = "direct"
    s.g        = G
    s.E        = E
    s.nu       = NU
    s.rho_m    = RHO_M
    s.rho_fill = RHO_F
    s.T_e      = TE
    s.dx       = dx
    s.dy       = dx
    s.qs       = qs
    s.bc_west  = bc
    s.bc_east  = bc
    s.bc_north = bc
    s.bc_south = bc
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s.initialize()
        s.run()
        w = s.w
        s.finalize()
    return w


# ---------------------------------------------------------------------------
# Convergence check helper
# ---------------------------------------------------------------------------

def _convergence_slopes(dxs, errs):
    return [
        np.log(errs[i] / errs[i+1]) / np.log(dxs[i] / dxs[i+1])
        for i in range(len(dxs) - 1)
    ]


def _assert_second_order(dxs, errs, bc_name, dim):
    slopes = _convergence_slopes(dxs, errs)
    for slope in slopes:
        assert slope > 1.8, (
            f"{bc_name} {dim}: expected O(dx²) convergence (slope ≥ 1.8), "
            f"got {slope:.2f}  (errors: {errs})"
        )


# ===========================================================================
# zero_displacement_zero_slope  (clamped)
# ===========================================================================

class TestClampedBC1D:

    def test_convergence_order(self):
        """1-D clamped BC achieves O(dx²) MMS convergence."""
        dxs, errs = [], []
        for nx in NX_1D:
            dx, w_ex, qs = _mms_1d_clamped(nx)
            w = _run_1d(dx, qs, "zero_displacement_zero_slope")
            errs.append(np.max(np.abs(w - w_ex)) / np.max(np.abs(w_ex)))
            dxs.append(dx)
        _assert_second_order(dxs, errs, "zero_displacement_zero_slope", "1-D")

    def test_dirichlet_exact(self):
        """Interior-only point load: w[0] and w[-1] are exactly zero."""
        nx = 101
        dx = L / (nx - 1)
        qs = np.zeros(nx)
        qs[nx // 2] = 1e6
        w = _run_1d(dx, qs, "zero_displacement_zero_slope")
        assert w[0]  == 0.0, f"w[0]  = {w[0]:.3e}; expected exactly 0"
        assert w[-1] == 0.0, f"w[-1] = {w[-1]:.3e}; expected exactly 0"


class TestClampedBC2D:

    def test_convergence_order(self):
        """2-D clamped BC achieves O(dx²) MMS convergence."""
        dxs, errs = [], []
        for n in N_2D:
            dx, w_ex, qs = _mms_2d_clamped(n)
            w = _run_2d(dx, qs, "zero_displacement_zero_slope")
            errs.append(np.max(np.abs(w - w_ex)) / np.max(np.abs(w_ex)))
            dxs.append(dx)
        _assert_second_order(dxs, errs, "zero_displacement_zero_slope", "2-D")

    def test_dirichlet_exact(self):
        """Interior-only point load: all four boundary edges are exactly zero."""
        n = 51
        dx = L / (n - 1)
        qs = np.zeros((n, n))
        qs[n // 2, n // 2] = 1e6
        w = _run_2d(dx, qs, "zero_displacement_zero_slope")
        assert np.all(w[:, 0]  == 0.0), "west edge not exactly zero"
        assert np.all(w[:, -1] == 0.0), "east edge not exactly zero"
        assert np.all(w[0, :]  == 0.0), "north edge not exactly zero"
        assert np.all(w[-1, :] == 0.0), "south edge not exactly zero"


# ===========================================================================
# zero_moment_zero_shear  (free end)
# ===========================================================================

class TestFreeEndBC1D:

    def test_convergence_order(self):
        """1-D free-end BC achieves O(dx²) MMS convergence."""
        dxs, errs = [], []
        for nx in NX_1D:
            dx, w_ex, qs = _mms_1d_free_end(nx)
            w = _run_1d(dx, qs, "zero_moment_zero_shear")
            errs.append(np.max(np.abs(w - w_ex)) / np.max(np.abs(w_ex)))
            dxs.append(dx)
        _assert_second_order(dxs, errs, "zero_moment_zero_shear", "1-D")


class TestFreeEndBC2D:

    def test_convergence_order(self):
        """2-D free-end BC achieves O(dx²) MMS convergence."""
        dxs, errs = [], []
        for n in N_2D:
            dx, w_ex, qs = _mms_2d_free_end(n)
            w = _run_2d(dx, qs, "zero_moment_zero_shear")
            errs.append(np.max(np.abs(w - w_ex)) / np.max(np.abs(w_ex)))
            dxs.append(dx)
        _assert_second_order(dxs, errs, "zero_moment_zero_shear", "2-D")


# ===========================================================================
# zero_displacement_zero_moment  (pinned / simply-supported)
# ===========================================================================

class TestPinnedBC1D:

    def test_convergence_order(self):
        """1-D pinned BC achieves O(dx²) MMS convergence."""
        dxs, errs = [], []
        for nx in NX_1D:
            dx, w_ex, qs = _mms_1d_pinned(nx)
            w = _run_1d(dx, qs, "zero_displacement_zero_moment")
            errs.append(np.max(np.abs(w - w_ex)) / np.max(np.abs(w_ex)))
            dxs.append(dx)
        _assert_second_order(dxs, errs, "zero_displacement_zero_moment", "1-D")


class TestPinnedBC2D:

    def test_convergence_order(self):
        """2-D pinned BC achieves O(dx²) MMS convergence."""
        dxs, errs = [], []
        for n in N_2D:
            dx, w_ex, qs = _mms_2d_pinned(n)
            w = _run_2d(dx, qs, "zero_displacement_zero_moment")
            errs.append(np.max(np.abs(w - w_ex)) / np.max(np.abs(w_ex)))
            dxs.append(dx)
        _assert_second_order(dxs, errs, "zero_displacement_zero_moment", "2-D")


# ===========================================================================
# zero_slope_zero_shear  (mirror / symmetric)
# ===========================================================================

class TestMirrorBC1D:

    def test_convergence_order(self):
        """1-D mirror BC achieves O(dx²) MMS convergence."""
        dxs, errs = [], []
        for nx in NX_1D:
            dx, w_ex, qs = _mms_1d_mirror(nx)
            w = _run_1d(dx, qs, "zero_slope_zero_shear")
            errs.append(np.max(np.abs(w - w_ex)) / np.max(np.abs(w_ex)))
            dxs.append(dx)
        _assert_second_order(dxs, errs, "zero_slope_zero_shear", "1-D")


class TestMirrorBC2D:

    def test_convergence_order(self):
        """2-D mirror BC achieves O(dx²) MMS convergence."""
        dxs, errs = [], []
        for n in N_2D:
            dx, w_ex, qs = _mms_2d_mirror(n)
            w = _run_2d(dx, qs, "zero_slope_zero_shear")
            errs.append(np.max(np.abs(w - w_ex)) / np.max(np.abs(w_ex)))
            dxs.append(dx)
        _assert_second_order(dxs, errs, "zero_slope_zero_shear", "2-D")


# ===========================================================================
# Variable elastic thickness (FD only)
# ===========================================================================

class TestVariableTe1D:

    def test_convergence_order(self):
        """Variable-Te 1-D FD achieves O(dx²) MMS convergence.

        Uses linearly-varying D(ξ) = D₀(1+0.5ξ) with clamped BCs.  The
        manufactured load is derived analytically, so any deviation from
        second-order convergence indicates a bug in the variable-Te stencil.
        """
        dxs, errs = [], []
        for nx in NX_1D:
            dx, te, w_ex, qs = _mms_1d_variable_te(nx)
            w = _run_1d_te(dx, te, qs, "zero_displacement_zero_slope")
            errs.append(np.max(np.abs(w - w_ex)) / np.max(np.abs(w_ex)))
            dxs.append(dx)
        _assert_second_order(dxs, errs, "variable-Te", "1-D")
