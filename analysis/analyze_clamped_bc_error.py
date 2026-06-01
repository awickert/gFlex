#!/usr/bin/env python3
"""
Error analysis: zero_displacement_zero_slope ghost-node fix (2026)

Quantifies the error of the original (pre-fix) and corrected ghost-node
implementations using a Method of Manufactured Solutions (MMS) exact solution
for a 1-D and 2-D clamped plate on an elastic foundation.

Manufactured solution (1-D)
---------------------------
  xi = x / L,   L = (N-1) * dx

  w_exact(xi) = -W0 * xi**2 * (1 - xi)**2

satisfies both BCs at each end exactly:
  w(0) = 0,  w(L) = 0          (zero displacement)
  w'(0) = 0, w'(L) = 0         (zero slope)

The 4th derivative is constant:
  d4w/dx4 = 24 * W0 / L**4

so the load that produces this solution under the full flexure equation
  D * d4w/dx4 + drho*g * w = -qs
is:
  qs(xi) = 24*D*W0/L**4 + drho*g * W0 * xi**2 * (1 - xi)**2

Both terms are spatially varying (the second through the foundation restoring
force), making this a nontrivial test of the full governing equation.

Manufactured solution (2-D)
---------------------------
  xi = x/L,  eta = y/L   (square domain)

  w_exact(xi, eta) = -W0 * g(xi) * g(eta),   g(t) = t**2 * (1-t)**2

satisfies all BCs on all four edges: w = dw/dn = 0 at xi,eta in {0,1}.
The corresponding load follows from the 2-D biharmonic:

  qs[i,j] = D*W0/L**4 * [24*g(eta[i]) + 2*g''(xi[j])*g''(eta[i]) + 24*g(xi[j])]
             + drho*g * W0 * g(xi[j]) * g(eta[i])

where g''(t) = 2 - 12t + 12t**2

Error metric
------------
  e_rel = max|w_num - w_exact| / max|w_exact|   (L-inf relative error)

The old implementation left the boundary row coupled to interior nodes and
did not apply the even-reflection ghost at the first interior node; this script
shows how the relative error and its convergence rate differ from the corrected
version.
"""

import numpy as np
import matplotlib.pyplot as plt

from gflex.f1d import F1D
from gflex.f2d import F2D


# ---------------------------------------------------------------------------
# Reconstruct the old (pre-fix) BC behaviour as a subclass
# ---------------------------------------------------------------------------

class F1D_OldClampedBC(F1D):
    """F1D with the pre-2026 ghost-node treatment for zero_displacement_zero_slope.

    The original code dropped ghost-node stencil terms entirely rather than
    folding them via even reflection, and did not decouple the boundary row
    from interior nodes.
    """

    def _bc_zero_displacement_zero_slope(self):
        if self.bc_west == "zero_displacement_zero_slope":
            i = 0
            self.l2[i] = np.nan
            self.l1[i] = np.nan
            # r1[0] and r2[0] NOT zeroed: boundary row left coupled to interior
            i = 1
            self.l2[i] = np.nan
            # c0[1] NOT modified: ghost reflection w[-1]=w[1] not applied
        if self.bc_east == "zero_displacement_zero_slope":
            i = -2
            # c0[-2] NOT modified: ghost reflection w[N]=w[N-2] not applied
            self.r2[i] = np.nan
            i = -1
            # l2[-1] and l1[-1] NOT zeroed: boundary row left coupled to interior
            self.r1[i] = np.nan
            self.r2[i] = np.nan


# ---------------------------------------------------------------------------
# MMS helpers
# ---------------------------------------------------------------------------

def mms_exact(xi, W0):
    """Exact deflection: -W0 * xi**2 * (1-xi)**2."""
    return -W0 * xi**2 * (1 - xi)**2


def mms_load(xi, W0, D, L, drho_g):
    """Manufactured load qs that produces mms_exact under the full flexure equation."""
    return 24.0 * D * W0 / L**4 + drho_g * W0 * xi**2 * (1 - xi)**2


def run_mms(solver_class, nx, te, E, nu, rho_m, rho_fill, g, W0, L):
    """Solve the MMS problem and return (w_numerical, w_exact, L-inf relative error)."""
    dx = L / (nx - 1)
    xi = np.arange(nx) / (nx - 1)

    D = E * te**3 / (12.0 * (1.0 - nu**2))
    drho_g = (rho_m - rho_fill) * g

    w_ex = mms_exact(xi, W0)
    qs = mms_load(xi, W0, D, L, drho_g)

    s = solver_class()
    s.dx = dx
    s.te = te
    s.E = E
    s.nu = nu
    s.rho_m = rho_m
    s.rho_fill = rho_fill
    s.g = g
    s.qs = qs
    s.bc_west = "zero_displacement_zero_slope"
    s.bc_east = "zero_displacement_zero_slope"
    s.method = "fd"
    s.quiet = True
    s.verbose = False
    s.debug = False
    s.initialize()
    s.run()
    s.finalize()

    err = np.max(np.abs(s.w - w_ex)) / np.max(np.abs(w_ex))
    return s.w, w_ex, err, dx, xi


# ---------------------------------------------------------------------------
# Physical parameters
# ---------------------------------------------------------------------------

te     = 30.0e3     # elastic thickness, m
E      = 65.0e9     # Young's modulus, Pa
nu     = 0.25
rho_m  = 3300.0     # mantle density, kg m-3
rho_fill = 0.0      # infill (air)
g      = 9.81       # m s-2
L      = 600.0e3    # plate length, m
W0     = 1600.0     # MMS amplitude, m  →  max deflection ~100 m

D      = E * te**3 / (12.0 * (1.0 - nu**2))
drho_g = (rho_m - rho_fill) * g
alpha  = (4.0 * D / drho_g) ** 0.25

# ---------------------------------------------------------------------------
# 1. Single-resolution comparison
# ---------------------------------------------------------------------------

nx_ref = 201

w_new, w_ex, err_new, dx_ref, xi_ref = run_mms(
    F1D, nx_ref, te, E, nu, rho_m, rho_fill, g, W0, L
)
w_old, _,    err_old, _,      _      = run_mms(
    F1D_OldClampedBC, nx_ref, te, E, nu, rho_m, rho_fill, g, W0, L
)

print("=" * 62)
print("  MMS error analysis: zero_displacement_zero_slope fix")
print("=" * 62)
print(f"  L = {L/1e3:.0f} km,  te = {te/1e3:.0f} km,  "
      f"alpha = {alpha/1e3:.0f} km")
print(f"  W0 = {W0:.0f} m  ->  |w_exact|_max = {np.max(np.abs(w_ex)):.1f} m")
print(f"  nx = {nx_ref},  dx = {dx_ref/1e3:.2f} km")
print()
print(f"  {'Implementation':<20}  {'L-inf rel. error':>18}")
print(f"  {'-'*20}  {'-'*18}")
print(f"  {'Old (pre-fix)':<20}  {err_old:>18.4e}")
print(f"  {'New (corrected)':<20}  {err_new:>18.4e}")
print(f"  {'Improvement':<20}  {err_old/err_new:>17.1f}x")
print()
print(f"  Boundary values (should be 0):")
print(f"    Old:  w[0]={w_old[0]:+.4e} m,  w[-1]={w_old[-1]:+.4e} m")
print(f"    New:  w[0]={w_new[0]:+.4e} m,  w[-1]={w_new[-1]:+.4e} m")
print()

# ---------------------------------------------------------------------------
# 2. Convergence study
# ---------------------------------------------------------------------------

nx_vals = [26, 51, 101, 201, 401, 801]

errs_new = []
errs_old = []
dx_vals  = []

for nx in nx_vals:
    _, _, e_new, dx, _ = run_mms(F1D,               nx, te, E, nu, rho_m, rho_fill, g, W0, L)
    _, _, e_old, _,  _ = run_mms(F1D_OldClampedBC,  nx, te, E, nu, rho_m, rho_fill, g, W0, L)
    errs_new.append(e_new)
    errs_old.append(e_old)
    dx_vals.append(dx)
    print(f"  nx={nx:4d}  dx={dx/1e3:6.2f} km  |  "
          f"old={e_old:.3e}  new={e_new:.3e}")

dx_arr      = np.array(dx_vals)
errs_new_arr = np.array(errs_new)
errs_old_arr = np.array(errs_old)

# Fit convergence slopes over the finest half of the resolution range
n_fit = len(nx_vals) // 2
slope_new = np.polyfit(np.log(dx_arr[-n_fit:]), np.log(errs_new_arr[-n_fit:]), 1)[0]
slope_old = np.polyfit(np.log(dx_arr[-n_fit:]), np.log(errs_old_arr[-n_fit:]), 1)[0]

print()
print(f"  Convergence order (finest {n_fit} points):")
print(f"    Old: O(dx^{slope_old:.2f})")
print(f"    New: O(dx^{slope_new:.2f})")
print()

# ---------------------------------------------------------------------------
# 3. Plots
# ---------------------------------------------------------------------------

fig, axes = plt.subplots(1, 3, figsize=(14, 4))
x_km = xi_ref * L / 1e3

# Panel 1: deflection profiles
ax = axes[0]
ax.plot(x_km, w_ex,    "k-",  lw=2,   label="Exact (MMS)")
ax.plot(x_km, w_new,   "b--", lw=1.5, label="New (corrected)")
ax.plot(x_km, w_old,   "r:",  lw=1.5, label="Old (pre-fix)")
ax.set_xlabel("x [km]")
ax.set_ylabel("w [m]")
ax.set_title("Deflection profiles")
ax.legend(fontsize=8)

# Panel 2: residuals at reference resolution
ax = axes[1]
ax.plot(x_km, w_new - w_ex, "b-",  lw=1.5, label=f"New  (err={err_new:.2e})")
ax.plot(x_km, w_old - w_ex, "r--", lw=1.5, label=f"Old  (err={err_old:.2e})")
ax.axhline(0, color="k", lw=0.5)
ax.set_xlabel("x [km]")
ax.set_ylabel("w_num - w_exact [m]")
ax.set_title(f"Residuals at dx = {dx_ref/1e3:.1f} km")
ax.legend(fontsize=8)

# Panel 3: convergence
ax = axes[2]
ax.loglog(dx_arr / 1e3, errs_new_arr, "bs-",  label=f"New  (O(dx^{slope_new:.1f}))")
ax.loglog(dx_arr / 1e3, errs_old_arr, "r^--", label=f"Old  (O(dx^{slope_old:.1f}))")
# Reference lines
ref_x = np.array([dx_arr.min(), dx_arr.max()]) / 1e3
ax.loglog(ref_x, 1e-5 * (ref_x / ref_x[0]) ** 2, "k:",  lw=0.8, label="O(dx²)")
ax.loglog(ref_x, 3e-4 * (ref_x / ref_x[0]) ** 1, "k--", lw=0.8, label="O(dx¹)")
ax.set_xlabel("dx [km]")
ax.set_ylabel("L-inf relative error")
ax.set_title("Convergence with grid refinement")
ax.legend(fontsize=8)

fig.suptitle(
    "MMS error analysis: zero_displacement_zero_slope\n"
    f"te={te/1e3:.0f} km,  L={L/1e3:.0f} km,  "
    f"alpha={alpha/1e3:.0f} km,  W0={W0:.0f} m",
    fontsize=10,
)
fig.tight_layout()
fig.savefig("analysis/results/clamped_bc_error.png", dpi=150, bbox_inches="tight")
print("  Figure saved to analysis/results/clamped_bc_error.png")


# ===========================================================================
# 2-D: F2D_OldClampedBC subclass replicating pre-fix behaviour
# ===========================================================================

class F2D_OldClampedBC(F2D):
    """F2D with the pre-2026 ghost-node treatment for zero_displacement_zero_slope.

    The original code left the boundary row/column coupled to interior nodes
    (off-diagonal stencil entries were not set to inf) and did not apply the
    even-reflection ghost at the first interior row/column.  The north/south
    blocks additionally contained a "[!= inf] = nan" pattern that was harmless
    in practice (those entries land outside the assembled matrix range).

    Strategy: call super()._apply_bc_flexure() (the corrected code) then undo
    the two corrections — (a) restore in-domain off-diagonal coupling at the
    boundary row/column, (b) subtract the ghost coefficient from cj0i0 at the
    first interior row/column.
    """

    _WEST_RESTORE = [
        "cj0i_2", "cj0i_1", "cj0i1", "cj0i2",
        "cj1i_1", "cj1i0", "cj1i1", "cj2i0",
    ]
    _EAST_RESTORE = [
        "cj_2i0", "cj_1i_1", "cj_1i0", "cj_1i1",
        "cj0i_2", "cj0i_1", "cj0i1", "cj0i2",
    ]
    _NS_RESTORE = [
        "cj_2i0", "cj_1i_1", "cj_1i0", "cj_1i1",
        "cj0i_2", "cj0i_1", "cj0i1", "cj0i2",
        "cj1i_1", "cj1i0", "cj1i1", "cj2i0",
    ]

    def _apply_bc_flexure(self):
        all_arrays = self._WEST_RESTORE + ["cj_2i0", "cj_1i_1", "cj_1i0", "cj_1i1"]
        all_arrays = list(dict.fromkeys(all_arrays + self._NS_RESTORE))
        saved = {a: getattr(self, a).copy() for a in all_arrays}

        super()._apply_bc_flexure()

        if self.bc_west == "zero_displacement_zero_slope":
            for attr in self._WEST_RESTORE:
                getattr(self, attr)[:, 0] = saved[attr][:, 0]
            self.cj0i0[:, 1] -= self.cj_2i0_coeff_ij[:, 1]

        if self.bc_east == "zero_displacement_zero_slope":
            for attr in self._EAST_RESTORE:
                getattr(self, attr)[:, -1] = saved[attr][:, -1]
            self.cj0i0[:, -2] -= self.cj2i0_coeff_ij[:, -2]

        if self.bc_north == "zero_displacement_zero_slope":
            for attr in self._NS_RESTORE:
                getattr(self, attr)[0, :] = saved[attr][0, :]
            self.cj0i0[1, :] -= self.cj0i_2_coeff_ij[1, :]

        if self.bc_south == "zero_displacement_zero_slope":
            for attr in self._NS_RESTORE:
                getattr(self, attr)[-1, :] = saved[attr][-1, :]
            self.cj0i0[-2, :] -= self.cj0i2_coeff_ij[-2, :]


# ---------------------------------------------------------------------------
# 2-D MMS helpers:  g(t) = t^2*(1-t)^2,  g''(t) = 2 - 12t + 12t^2
# ---------------------------------------------------------------------------

def _g(t):
    return t**2 * (1 - t)**2


def _g2(t):
    """Second derivative of g."""
    return 2 - 12*t + 12*t**2


def run_mms_2d(solver_class, ny, nx, te, E, nu, rho_m, rho_fill, g_acc, W0, L):
    """Solve the 2-D MMS problem; return (w_num, w_exact, L-inf rel error, dx)."""
    dx = L / (nx - 1)
    dy = L / (ny - 1)
    xi  = np.arange(nx) / (nx - 1)   # x-direction
    eta = np.arange(ny) / (ny - 1)   # y-direction

    D      = E * te**3 / (12.0 * (1.0 - nu**2))
    drho_g = (rho_m - rho_fill) * g_acc

    # w_ex[i, j] = -W0 * g(eta[i]) * g(xi[j]),  shape (ny, nx)
    gxi  = _g(xi)[np.newaxis, :]    # (1, nx)
    geta = _g(eta)[:, np.newaxis]   # (ny, 1)
    w_ex = -W0 * geta * gxi

    g2xi  = _g2(xi)[np.newaxis, :]
    g2eta = _g2(eta)[:, np.newaxis]
    qs = (D * W0 / L**4 * (24 * geta + 2 * g2xi * g2eta + 24 * gxi)
          + drho_g * W0 * gxi * geta)

    s = solver_class()
    s.quiet    = True
    s.method   = "fd"
    s.solver   = "direct"
    s.g        = g_acc
    s.E        = E
    s.nu       = nu
    s.rho_m    = rho_m
    s.rho_fill = rho_fill
    s.te       = te
    s.qs       = qs
    s.dx       = dx
    s.dy       = dy
    s.bc_west  = "zero_displacement_zero_slope"
    s.bc_east  = "zero_displacement_zero_slope"
    s.bc_north = "zero_displacement_zero_slope"
    s.bc_south = "zero_displacement_zero_slope"
    s.initialize()
    s.run()
    s.finalize()

    err = np.max(np.abs(s.w - w_ex)) / np.max(np.abs(w_ex))
    return s.w, w_ex, err, dx


# ===========================================================================
# 2-D analysis
# ===========================================================================

print("=" * 62)
print("  MMS error analysis: zero_displacement_zero_slope fix — 2-D")
print("=" * 62)

n_vals_2d = [26, 51, 101, 201]

errs_new_2d = []
errs_old_2d = []
dx_vals_2d  = []

print(f"  {'n':>4}  {'dx [km]':>8}  {'old error':>12}  {'new error':>12}")
print(f"  {'-'*4}  {'-'*8}  {'-'*12}  {'-'*12}")

for n in n_vals_2d:
    _, _, e_new, dx = run_mms_2d(F2D,               n, n, te, E, nu, rho_m, rho_fill, g, W0, L)
    _, _, e_old, _  = run_mms_2d(F2D_OldClampedBC,  n, n, te, E, nu, rho_m, rho_fill, g, W0, L)
    errs_new_2d.append(e_new)
    errs_old_2d.append(e_old)
    dx_vals_2d.append(dx)
    print(f"  {n:>4}  {dx/1e3:>8.2f}  {e_old:>12.3e}  {e_new:>12.3e}")

dx_arr_2d       = np.array(dx_vals_2d)
errs_new_arr_2d = np.array(errs_new_2d)
errs_old_arr_2d = np.array(errs_old_2d)

n_fit_2d = len(n_vals_2d) // 2
slope_new_2d = np.polyfit(np.log(dx_arr_2d[-n_fit_2d:]), np.log(errs_new_arr_2d[-n_fit_2d:]), 1)[0]
slope_old_2d = np.polyfit(np.log(dx_arr_2d[-n_fit_2d:]), np.log(errs_old_arr_2d[-n_fit_2d:]), 1)[0]

print()
print(f"  Convergence order (finest {n_fit_2d} points):")
print(f"    Old: O(dx^{slope_old_2d:.2f})")
print(f"    New: O(dx^{slope_new_2d:.2f})")
print()

# Reference run for figure
n_ref_2d = 101
w_new_2d, w_ex_2d, err_new_2d_ref, dx_ref_2d = run_mms_2d(
    F2D,              n_ref_2d, n_ref_2d, te, E, nu, rho_m, rho_fill, g, W0, L
)
w_old_2d, _,       err_old_2d_ref, _         = run_mms_2d(
    F2D_OldClampedBC, n_ref_2d, n_ref_2d, te, E, nu, rho_m, rho_fill, g, W0, L
)

xi_ref_2d = np.arange(n_ref_2d) / (n_ref_2d - 1)
x_km_2d   = xi_ref_2d * L / 1e3

# 2-D figure
fig2, axes2 = plt.subplots(1, 3, figsize=(14, 4))
mid = n_ref_2d // 2

ax = axes2[0]
ax.plot(x_km_2d, w_ex_2d[mid, :],          "k-",  lw=2,   label="Exact (MMS)")
ax.plot(x_km_2d, w_new_2d[mid, :],         "b--", lw=1.5, label="New (corrected)")
ax.plot(x_km_2d, w_old_2d[mid, :],         "r:",  lw=1.5, label="Old (pre-fix)")
ax.set_xlabel("x [km]")
ax.set_ylabel("w [m]  (centre row)")
ax.set_title("Deflection profiles")
ax.legend(fontsize=8)

ax = axes2[1]
ax.plot(x_km_2d, (w_new_2d - w_ex_2d)[mid, :], "b-",  lw=1.5,
        label=f"New  (err={err_new_2d_ref:.2e})")
ax.plot(x_km_2d, (w_old_2d - w_ex_2d)[mid, :], "r--", lw=1.5,
        label=f"Old  (err={err_old_2d_ref:.2e})")
ax.axhline(0, color="k", lw=0.5)
ax.set_xlabel("x [km]")
ax.set_ylabel("w_num − w_exact [m]  (centre row)")
ax.set_title(f"Residuals at dx = {dx_ref_2d/1e3:.1f} km")
ax.legend(fontsize=8)

ax = axes2[2]
ax.loglog(dx_arr_2d/1e3, errs_new_arr_2d, "bs-",  label=f"New  (O(dx^{slope_new_2d:.2f}))")
ax.loglog(dx_arr_2d/1e3, errs_old_arr_2d, "r^--", label=f"Old  (O(dx^{slope_old_2d:.2f}))")
ref_x2 = np.array([dx_arr_2d.min(), dx_arr_2d.max()]) / 1e3
ax.loglog(ref_x2, 1e-4 * (ref_x2 / ref_x2[0])**2, "k:",  lw=0.8, label="O(dx²)")
ax.loglog(ref_x2, 1e-2 * (ref_x2 / ref_x2[0])**1, "k--", lw=0.8, label="O(dx¹)")
ax.set_xlabel("dx [km]")
ax.set_ylabel("L-inf relative error")
ax.set_title("Convergence with grid refinement")
ax.legend(fontsize=8)

fig2.suptitle(
    "MMS error analysis: zero_displacement_zero_slope  (2-D)\n"
    f"te={te/1e3:.0f} km,  L={L/1e3:.0f} km,  "
    f"alpha={alpha/1e3:.0f} km,  W0={W0:.0f} m",
    fontsize=10,
)
fig2.tight_layout()
fig2.savefig("analysis/results/clamped_bc_error_2d.png", dpi=150, bbox_inches="tight")
print("  Figure saved to analysis/results/clamped_bc_error_2d.png")
