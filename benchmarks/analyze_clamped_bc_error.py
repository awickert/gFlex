#!/usr/bin/env python3
"""
Error analysis: zero_displacement_zero_slope ghost-node fix (2026)

Quantifies the error of the original (pre-fix) and corrected ghost-node
implementations using a Method of Manufactured Solutions (MMS) exact solution
for a 1-D clamped-clamped plate on an elastic foundation.

Manufactured solution
---------------------
  xi = x / L,   L = (N-1) * dx

  w_exact(xi) = -W0 * xi**2 * (1 - xi)**2

satisfies both BCs at each end exactly:
  w(0) = 0,  w(L) = 0          (zero displacement)
  w'(0) = 0, w'(L) = 0         (zero slope)

The 4th derivative is constant:
  d4w/dx4 = -24 * W0 / L**4

so the load that produces this solution under the full flexure equation
  D * d4w/dx4 + drho*g * w = -qs
is:
  qs(xi) = 24*D*W0/L**4 + drho*g * W0 * xi**2 * (1 - xi)**2

Both terms are spatially varying (the second through the foundation restoring
force), making this a nontrivial test of the full governing equation.

Error metric
------------
  e_rel = max|w_num - w_exact| / max|w_exact|   (L-inf relative error)

The old implementation used implicit ghost = 0 and left the boundary row
coupled to interior nodes; this script shows how the relative error and its
convergence rate differ from the corrected version.
"""

import numpy as np
import matplotlib.pyplot as plt

from gflex.f1d import F1D


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
fig.savefig("benchmarks/results/clamped_bc_error.png", dpi=150, bbox_inches="tight")
print("  Figure saved to benchmarks/results/clamped_bc_error.png")
