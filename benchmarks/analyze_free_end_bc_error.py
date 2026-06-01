#!/usr/bin/env python3
"""
Error analysis: zero_moment_zero_shear ghost-node fix (2026)

Quantifies the error of the original (pre-fix) and corrected ghost-node
implementations using a Method of Manufactured Solutions (MMS) exact solution
for a 1-D and 2-D free-end plate on an elastic foundation.

Manufactured solution (1-D)
---------------------------
  xi = x / L,   L = (N-1) * dx

  w_exact(xi) = -W0 * xi**4 * (1 - xi)**4

satisfies all four boundary conditions at each end exactly:
  w''(0)  = 0,  w''(L)  = 0          (zero moment)
  w'''(0) = 0,  w'''(L) = 0          (zero shear)

The 4th derivative is:
  d4w/dx4 = -W0/L**4 * (24 - 480*xi + 2160*xi**2 - 3360*xi**3 + 1680*xi**4)

so the load that produces this solution under the full flexure equation
  D * d4w/dx4 + drho*g * w = -qs
is:
  qs(xi) = D*W0/L**4 * (24 - 480*xi + 2160*xi**2 - 3360*xi**3 + 1680*xi**4)
           + drho*g * W0 * xi**4 * (1 - xi)**4

Manufactured solution (2-D)
---------------------------
  xi = x/L,  eta = y/L   (square domain)

  w_exact(xi, eta) = -W0 * f(xi) * f(eta),   f(t) = t**4 * (1-t)**4

satisfies all BCs on all four edges: w'' = w''' = 0 at xi,eta in {0,1}.
The corresponding load follows from the 2-D biharmonic:

  qs[i,j] = D*W0/L**4 * [f''''(xi[j])*f(eta[i]) + 2*f''(xi[j])*f''(eta[i])
                           + f(xi[j])*f''''(eta[i])]
             + drho*g * W0 * f(xi[j]) * f(eta[i])

where f''(t) = 12t**2 - 80t**3 + 180t**4 - 168t**5 + 56t**6
      f''''(t) = 24 - 480t + 2160t**2 - 3360t**3 + 1680t**4

Error metric
------------
  e_rel = max|w_num - w_exact| / max|w_exact|   (L-inf relative error)

The old implementation used the shear-condition ghost at the first interior
node (staggered one cell inward from the boundary) instead of the moment-
condition ghost at the boundary.  This script quantifies the difference.
"""

import numpy as np
import matplotlib.pyplot as plt

from gflex.f1d import F1D
from gflex.f2d import F2D


# ---------------------------------------------------------------------------
# Polynomial helpers for the MMS function f(t) = t^4*(1-t)^4
# ---------------------------------------------------------------------------

def _f(t):
    return t**4 * (1 - t)**4


def _f2(t):
    """Second derivative of f."""
    return 12*t**2 - 80*t**3 + 180*t**4 - 168*t**5 + 56*t**6


def _f4(t):
    """Fourth derivative of f."""
    return 24 - 480*t + 2160*t**2 - 3360*t**3 + 1680*t**4


# ---------------------------------------------------------------------------
# 1-D subclass replicating the old (pre-fix) shear-ghost treatment
# ---------------------------------------------------------------------------

class F1D_OldFreeEndBC(F1D):
    """F1D with the pre-2026 ghost-node treatment for zero_moment_zero_shear.

    The original code eliminated the ghost at the first interior node using
    the shear condition (d³w/dx³=0) staggered to that node, not the moment
    condition at the boundary.
    """

    def _bc_zero_moment_zero_shear(self):
        if self._bc_west_norm == "zero_moment_zero_shear":
            # Boundary node (i=0): identical to corrected code
            i = 0
            self.l2[i] += np.nan
            self.l1[i] += np.nan
            self.c0[i] += 4 * self.l2_coeff_i[i] + 2 * self.l1_coeff_i[i]
            self.r1[i] += -4 * self.l2_coeff_i[i] - self.l1_coeff_i[i]
            self.r2[i] += self.l2_coeff_i[i]
            # First interior node (i=1): OLD shear ghost
            # d³w/dx³=0 at x=dx → w[-1] = 2w[0] - 2w[2] + w[3]
            i = 1
            self.l2[i] += np.nan
            self.l1[i] += 2 * self.l2_coeff_i[i]
            # c0 unchanged  (new code does -= l2_coeff; old did nothing)
            self.r1[i] += -2 * self.l2_coeff_i[i]
            self.r2[i] += self.l2_coeff_i[i]

        if self._bc_east_norm == "zero_moment_zero_shear":
            # First interior node (i=N-2): OLD shear ghost
            # d³w/dx³=0 at x=(N-2)*dx → w[N] = 2w[N-1] - 2w[N-3] + w[N-4]
            i = -2
            self.l2[i] += self.r2_coeff_i[i]
            self.l1[i] += -2 * self.r2_coeff_i[i]
            # c0 unchanged  (new code does -= r2_coeff; old did nothing)
            self.r1[i] += 2 * self.r2_coeff_i[i]
            self.r2[i] += np.nan
            # Boundary node (i=N-1): identical to corrected code
            i = -1
            self.l2[i] += self.r2_coeff_i[i]
            self.l1[i] += -4 * self.r2_coeff_i[i] - self.r1_coeff_i[i]
            self.c0[i] += 4 * self.r2_coeff_i[i] + 2 * self.r1_coeff_i[i]
            self.r1[i] += np.nan
            self.r2[i] += np.nan


# ---------------------------------------------------------------------------
# 2-D subclass replicating the old (pre-fix) shear-ghost treatment
# ---------------------------------------------------------------------------

class F2D_OldFreeEndBC(F2D):
    """F2D with the pre-2026 ghost-node treatment for zero_moment_zero_shear.

    Overrides _apply_bc_flexure: calls super() (which applies the corrected
    moment ghost) then undoes the moment-ghost modification at the first
    interior row/column and replaces it with the original staggered shear ghost.
    """

    def _apply_bc_flexure(self):
        super()._apply_bc_flexure()

        if self.bc_west == "zero_moment_zero_shear":
            # New code at j=1 added: cj_1i0 += 2*coeff, cj0i0 -= coeff
            # Old code at j=1:       cj_1i0 += 2*coeff, cj1i0 += -2*coeff, cj2i0 += coeff
            # (cj_2i0[:, 1] is already inf from both old and new)
            j = 1
            c = self.cj_2i0_coeff_ij[:, j]
            self.cj0i0[:, j] += c        # undo  cj0i0 -= c  from new code
            self.cj1i0[:, j] += -2 * c   # add old shear-ghost term
            self.cj2i0[:, j] += c        # add old shear-ghost term

        if self.bc_east == "zero_moment_zero_shear":
            # New code at j=-2 added: cj0i0 -= coeff, cj1i0 += 2*coeff
            # Old code at j=-2:       cj_2i0 += coeff, cj_1i0 += -2*coeff, cj1i0 += 2*coeff
            # (cj2i0[:, -2] is already inf from both old and new)
            j = -2
            c = self.cj2i0_coeff_ij[:, j]
            self.cj0i0[:, j] += c        # undo  cj0i0 -= c  from new code
            self.cj_2i0[:, j] += c       # add old shear-ghost term
            self.cj_1i0[:, j] += -2 * c  # add old shear-ghost term

        if self.bc_north == "zero_moment_zero_shear":
            # New code at i=1 added: cj0i_1 += 2*coeff, cj0i0 -= coeff
            # Old code at i=1:       cj0i_1 += 2*coeff, cj0i1 += -2*coeff, cj0i2 += coeff
            i = 1
            c = self.cj0i_2_coeff_ij[i, :]
            self.cj0i0[i, :] += c        # undo  cj0i0 -= c  from new code
            self.cj0i1[i, :] += -2 * c   # add old shear-ghost term
            self.cj0i2[i, :] += c        # add old shear-ghost term

        if self.bc_south == "zero_moment_zero_shear":
            # New code at i=-2 added: cj0i0 -= coeff, cj0i1 += 2*coeff
            # Old code at i=-2:       cj0i_2 += coeff, cj0i_1 += -2*coeff, cj0i1 += 2*coeff
            i = -2
            c = self.cj0i2_coeff_ij[i, :]
            self.cj0i0[i, :] += c        # undo  cj0i0 -= c  from new code
            self.cj0i_2[i, :] += c       # add old shear-ghost term
            self.cj0i_1[i, :] += -2 * c  # add old shear-ghost term


# ---------------------------------------------------------------------------
# Run helpers
# ---------------------------------------------------------------------------

def run_mms_1d(solver_class, nx, te, E, nu, rho_m, rho_fill, g, W0, L):
    """Solve the 1-D MMS problem; return (w_num, w_exact, L-inf rel error, dx, xi)."""
    dx = L / (nx - 1)
    xi = np.arange(nx) / (nx - 1)

    D      = E * te**3 / (12.0 * (1.0 - nu**2))
    drho_g = (rho_m - rho_fill) * g

    w_ex = -W0 * _f(xi)
    qs   = D * W0 / L**4 * _f4(xi) + drho_g * W0 * _f(xi)

    s = solver_class()
    s.dx       = dx
    s.te       = te
    s.E        = E
    s.nu       = nu
    s.rho_m    = rho_m
    s.rho_fill = rho_fill
    s.g        = g
    s.qs       = qs
    s.bc_west  = "zero_moment_zero_shear"
    s.bc_east  = "zero_moment_zero_shear"
    s.method   = "fd"
    s.quiet    = True
    s.verbose  = False
    s.debug    = False
    s.initialize()
    s.run()
    s.finalize()

    err = np.max(np.abs(s.w - w_ex)) / np.max(np.abs(w_ex))
    return s.w, w_ex, err, dx, xi


def run_mms_2d(solver_class, ny, nx, te, E, nu, rho_m, rho_fill, g, W0, L):
    """Solve the 2-D MMS problem; return (w_num, w_exact, L-inf rel error, dx)."""
    dx = L / (nx - 1)
    dy = L / (ny - 1)
    # xi: x-direction, shape (nx,); eta: y-direction, shape (ny,)
    xi  = np.arange(nx) / (nx - 1)
    eta = np.arange(ny) / (ny - 1)

    D      = E * te**3 / (12.0 * (1.0 - nu**2))
    drho_g = (rho_m - rho_fill) * g

    # w_ex[i, j] = -W0 * f(eta[i]) * f(xi[j]),  shape (ny, nx)
    w_ex = -W0 * _f(eta)[:, np.newaxis] * _f(xi)[np.newaxis, :]

    # qs[i, j] uses eta for rows and xi for columns, shape (ny, nx)
    fxi   = _f(xi)[np.newaxis, :]    # (1, nx)
    feta  = _f(eta)[:, np.newaxis]   # (ny, 1)
    f4xi  = _f4(xi)[np.newaxis, :]
    f4eta = _f4(eta)[:, np.newaxis]
    f2xi  = _f2(xi)[np.newaxis, :]
    f2eta = _f2(eta)[:, np.newaxis]
    qs = (D * W0 / L**4 * (f4xi * feta + 2 * f2xi * f2eta + fxi * f4eta)
          + drho_g * W0 * fxi * feta)

    s = solver_class()
    s.quiet    = True
    s.method   = "fd"
    s.solver   = "direct"
    s.g        = g
    s.E        = E
    s.nu       = nu
    s.rho_m    = rho_m
    s.rho_fill = rho_fill
    s.te       = te
    s.qs       = qs
    s.dx       = dx
    s.dy       = dy
    s.bc_west  = "zero_moment_zero_shear"
    s.bc_east  = "zero_moment_zero_shear"
    s.bc_north = "zero_moment_zero_shear"
    s.bc_south = "zero_moment_zero_shear"
    s.initialize()
    s.run()
    s.finalize()

    err = np.max(np.abs(s.w - w_ex)) / np.max(np.abs(w_ex))
    return s.w, w_ex, err, dx


# ---------------------------------------------------------------------------
# Physical parameters
# ---------------------------------------------------------------------------

te       = 30.0e3    # elastic thickness, m
E        = 65.0e9    # Young's modulus, Pa
nu       = 0.25
rho_m    = 3300.0    # mantle density, kg m-3
rho_fill = 0.0       # infill (air)
g        = 9.81      # m s-2
L        = 600.0e3   # plate length, m
# W0 so that max |w_exact| ≈ 100 m:
# max of f(t)=t^4(1-t)^4 at t=0.5 is (1/2)^8 = 1/256  → W0 = 100*256 = 25600
W0       = 25600.0

D      = E * te**3 / (12.0 * (1.0 - nu**2))
drho_g = (rho_m - rho_fill) * g
alpha  = (4.0 * D / drho_g) ** 0.25


# ===========================================================================
# 1-D analysis
# ===========================================================================

print("=" * 66)
print("  MMS error analysis: zero_moment_zero_shear fix — 1-D")
print("=" * 66)
print(f"  L = {L/1e3:.0f} km,  te = {te/1e3:.0f} km,  "
      f"alpha = {alpha/1e3:.0f} km")
print(f"  W0 = {W0:.0f} m  ->  |w_exact|_max = {W0/256:.1f} m")
print()

# 1. Single-resolution comparison
nx_ref = 201

w_new_1d, w_ex_1d, err_new_1d, dx_ref_1d, xi_ref = run_mms_1d(
    F1D, nx_ref, te, E, nu, rho_m, rho_fill, g, W0, L
)
w_old_1d, _,       err_old_1d, _,         _       = run_mms_1d(
    F1D_OldFreeEndBC, nx_ref, te, E, nu, rho_m, rho_fill, g, W0, L
)

print(f"  nx = {nx_ref},  dx = {dx_ref_1d/1e3:.2f} km")
print(f"  {'Implementation':<20}  {'L-inf rel. error':>18}")
print(f"  {'-'*20}  {'-'*18}")
print(f"  {'Old (pre-fix)':<20}  {err_old_1d:>18.4e}")
print(f"  {'New (corrected)':<20}  {err_new_1d:>18.4e}")
if err_new_1d > 0:
    print(f"  {'Improvement':<20}  {err_old_1d/err_new_1d:>17.1f}x")
print()

# 2. Convergence study
nx_vals_1d = [26, 51, 101, 201, 401, 801]

errs_new_1d = []
errs_old_1d = []
dx_vals_1d  = []

print(f"  {'nx':>4}  {'dx [km]':>8}  {'old error':>12}  {'new error':>12}")
print(f"  {'-'*4}  {'-'*8}  {'-'*12}  {'-'*12}")

for nx in nx_vals_1d:
    _, _, e_new, dx, _ = run_mms_1d(F1D,               nx, te, E, nu, rho_m, rho_fill, g, W0, L)
    _, _, e_old, _,  _ = run_mms_1d(F1D_OldFreeEndBC,  nx, te, E, nu, rho_m, rho_fill, g, W0, L)
    errs_new_1d.append(e_new)
    errs_old_1d.append(e_old)
    dx_vals_1d.append(dx)
    print(f"  {nx:>4}  {dx/1e3:>8.2f}  {e_old:>12.3e}  {e_new:>12.3e}")

dx_arr_1d       = np.array(dx_vals_1d)
errs_new_arr_1d = np.array(errs_new_1d)
errs_old_arr_1d = np.array(errs_old_1d)

n_fit = len(nx_vals_1d) // 2
slope_new_1d = np.polyfit(np.log(dx_arr_1d[-n_fit:]), np.log(errs_new_arr_1d[-n_fit:]), 1)[0]
slope_old_1d = np.polyfit(np.log(dx_arr_1d[-n_fit:]), np.log(errs_old_arr_1d[-n_fit:]), 1)[0]

print()
print(f"  Convergence order (finest {n_fit} points):")
print(f"    Old: O(dx^{slope_old_1d:.2f})")
print(f"    New: O(dx^{slope_new_1d:.2f})")
print()


# ===========================================================================
# 2-D analysis
# ===========================================================================

print("=" * 66)
print("  MMS error analysis: zero_moment_zero_shear fix — 2-D")
print("=" * 66)

n_vals_2d = [26, 51, 101, 201]

errs_new_2d = []
errs_old_2d = []
dx_vals_2d  = []

print(f"  {'n':>4}  {'dx [km]':>8}  {'old error':>12}  {'new error':>12}")
print(f"  {'-'*4}  {'-'*8}  {'-'*12}  {'-'*12}")

for n in n_vals_2d:
    _, _, e_new, dx = run_mms_2d(F2D,               n, n, te, E, nu, rho_m, rho_fill, g, W0, L)
    _, _, e_old, _  = run_mms_2d(F2D_OldFreeEndBC,  n, n, te, E, nu, rho_m, rho_fill, g, W0, L)
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

# 2-D reference run for figure
n_ref = 101
w_new_2d, w_ex_2d, err_new_2d_ref, dx_ref_2d = run_mms_2d(
    F2D,              n_ref, n_ref, te, E, nu, rho_m, rho_fill, g, W0, L
)
w_old_2d, _,       err_old_2d_ref, _         = run_mms_2d(
    F2D_OldFreeEndBC, n_ref, n_ref, te, E, nu, rho_m, rho_fill, g, W0, L
)

xi_ref_2d = np.arange(n_ref) / (n_ref - 1)
x_km_2d   = xi_ref_2d * L / 1e3


# ===========================================================================
# Figures
# ===========================================================================

x_km = xi_ref * L / 1e3

# 1-D figure
fig1, axes = plt.subplots(1, 3, figsize=(14, 4))

ax = axes[0]
ax.plot(x_km, w_ex_1d,  "k-",  lw=2,   label="Exact (MMS)")
ax.plot(x_km, w_new_1d, "b--", lw=1.5, label="New (corrected)")
ax.plot(x_km, w_old_1d, "r:",  lw=1.5, label="Old (pre-fix)")
ax.set_xlabel("x [km]")
ax.set_ylabel("w [m]")
ax.set_title("Deflection profiles")
ax.legend(fontsize=8)

ax = axes[1]
ax.plot(x_km, w_new_1d - w_ex_1d, "b-",  lw=1.5, label=f"New  (err={err_new_1d:.2e})")
ax.plot(x_km, w_old_1d - w_ex_1d, "r--", lw=1.5, label=f"Old  (err={err_old_1d:.2e})")
ax.axhline(0, color="k", lw=0.5)
ax.set_xlabel("x [km]")
ax.set_ylabel("w_num − w_exact [m]")
ax.set_title(f"Residuals at dx = {dx_ref_1d/1e3:.1f} km")
ax.legend(fontsize=8)

ax = axes[2]
ax.loglog(dx_arr_1d/1e3, errs_new_arr_1d, "bs-",  label=f"New  (O(dx^{slope_new_1d:.2f}))")
ax.loglog(dx_arr_1d/1e3, errs_old_arr_1d, "r^--", label=f"Old  (O(dx^{slope_old_1d:.2f}))")
ref_x = np.array([dx_arr_1d.min(), dx_arr_1d.max()]) / 1e3
ax.loglog(ref_x, 5e-5 * (ref_x / ref_x[0])**2, "k:",  lw=0.8, label="O(dx²)")
ax.set_xlabel("dx [km]")
ax.set_ylabel("L-inf relative error")
ax.set_title("Convergence with grid refinement")
ax.legend(fontsize=8)

fig1.suptitle(
    "MMS error analysis: zero_moment_zero_shear  (1-D)\n"
    f"te={te/1e3:.0f} km,  L={L/1e3:.0f} km,  "
    f"alpha={alpha/1e3:.0f} km,  W0={W0:.0f} m",
    fontsize=10,
)
fig1.tight_layout()
fig1.savefig("benchmarks/results/free_end_bc_error_1d.png", dpi=150, bbox_inches="tight")
print("  Figure saved to benchmarks/results/free_end_bc_error_1d.png")


# 2-D figure
fig2, axes2 = plt.subplots(1, 3, figsize=(14, 4))
mid = n_ref // 2

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
ax.loglog(ref_x2, 5e-5 * (ref_x2 / ref_x2[0])**2, "k:",  lw=0.8, label="O(dx²)")
ax.set_xlabel("dx [km]")
ax.set_ylabel("L-inf relative error")
ax.set_title("Convergence with grid refinement")
ax.legend(fontsize=8)

fig2.suptitle(
    "MMS error analysis: zero_moment_zero_shear  (2-D)\n"
    f"te={te/1e3:.0f} km,  L={L/1e3:.0f} km,  "
    f"alpha={alpha/1e3:.0f} km,  W0={W0:.0f} m",
    fontsize=10,
)
fig2.tight_layout()
fig2.savefig("benchmarks/results/free_end_bc_error_2d.png", dpi=150, bbox_inches="tight")
print("  Figure saved to benchmarks/results/free_end_bc_error_2d.png")
