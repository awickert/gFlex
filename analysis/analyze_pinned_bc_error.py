#!/usr/bin/env python3
"""
MMS convergence verification: zero_displacement_zero_moment BC (2026)

Verifies O(dx²) convergence of the pinned (simply-supported) boundary
condition for 1-D and 2-D elastic flexure on an elastic foundation.

Background
----------
zero_displacement_zero_moment enforces w = 0 (Dirichlet) and M = d²w/dx² = 0
(zero bending moment) at each boundary — the simply-supported or pinned
condition.  The implementation uses two distinct operations on each side:

  1. Boundary node: Dirichlet decoupling — all off-diagonal stencil terms
     are zeroed, leaving c0·w = q.  For q = 0 this gives w = 0 exactly.
  2. First interior node: odd-reflection ghost, w[ghost] = -w[interior],
     folds the zero-moment condition into the diagonal coefficient.

The 2-D implementation used Dirichlet from the outset (commit 4039013).
The 1-D code was reformulated in 678f5d4 to make the Dirichlet step
explicit; the old approach (subtracting l_coeff from r at the boundary
node) happened to give the same matrix because the mirror D-ghost makes
the stencil symmetric at the boundary node.

Manufactured solutions
----------------------
1-D  (xi = x/L, xi in [0,1]):

  w_exact(xi) = -W0 * sin(pi*xi)

  satisfies w = 0 and d²w/dx² = 0 at xi = 0 and xi = 1.
  Manufactured load:

    qs(xi) = W0 * (D*(pi/L)^4 + drho*g) * sin(pi*xi)

2-D  (xi = x/L, eta = y/L, square domain):

  w_exact(xi, eta) = -W0 * sin(pi*xi) * sin(pi*eta)

  satisfies w = 0 and M = 0 on all four edges.
  Manufactured load (from D*nabla^4 w + drho*g*w = -qs):

    qs(xi, eta) = W0 * (4*D*(pi/L)^4 + drho*g) * sin(pi*xi) * sin(pi*eta)

Error metric
------------
  e_rel = max|w_num - w_exact| / max|w_exact|   (L-inf relative error)
"""

import numpy as np
import matplotlib.pyplot as plt

from gflex.f1d import F1D
from gflex.f2d import F2D


# ---------------------------------------------------------------------------
# Physical parameters
# ---------------------------------------------------------------------------

E        = 65.0e9    # Young's modulus, Pa
nu       = 0.25
rho_m    = 3300.0    # mantle density, kg m-3
rho_fill = 0.0       # infill (air)
g        = 9.81      # m s-2
te       = 30.0e3    # elastic thickness, m
L        = 600.0e3   # plate length, m
W0       = 1600.0    # MMS amplitude; max|w_exact| = W0 at xi = 0.5

D      = E * te**3 / (12.0 * (1.0 - nu**2))
drho_g = (rho_m - rho_fill) * g
alpha  = (4.0 * D / drho_g) ** 0.25


# ---------------------------------------------------------------------------
# Run helpers
# ---------------------------------------------------------------------------

def run_1d(nx):
    """Return (w_num, w_exact, L-inf relative error, dx, xi) for 1-D MMS."""
    dx = L / (nx - 1)
    xi = np.arange(nx) / (nx - 1)

    w_ex = -W0 * np.sin(np.pi * xi)
    qs   = W0 * (D * (np.pi / L)**4 + drho_g) * np.sin(np.pi * xi)

    s = F1D()
    s.dx       = dx
    s.te       = te
    s.E        = E
    s.nu       = nu
    s.rho_m    = rho_m
    s.rho_fill = rho_fill
    s.g        = g
    s.qs       = qs
    s.bc_west  = "zero_displacement_zero_moment"
    s.bc_east  = "zero_displacement_zero_moment"
    s.method   = "fd"
    s.quiet    = True
    s.verbose  = False
    s.debug    = False
    s.initialize()
    s.run()
    s.finalize()

    err = np.max(np.abs(s.w - w_ex)) / np.max(np.abs(w_ex))
    return s.w, w_ex, err, dx, xi


def run_2d(ny, nx):
    """Return (w_num, w_exact, L-inf relative error, dx) for 2-D MMS."""
    dx = L / (nx - 1)
    dy = L / (ny - 1)
    xi  = np.arange(nx) / (nx - 1)
    eta = np.arange(ny) / (ny - 1)

    sin_xi  = np.sin(np.pi * xi)[np.newaxis, :]   # (1, nx)
    sin_eta = np.sin(np.pi * eta)[:, np.newaxis]   # (ny, 1)

    w_ex = -W0 * sin_eta * sin_xi
    qs   = W0 * (4.0 * D * (np.pi / L)**4 + drho_g) * sin_xi * sin_eta

    s = F2D()
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
    s.bc_west  = "zero_displacement_zero_moment"
    s.bc_east  = "zero_displacement_zero_moment"
    s.bc_north = "zero_displacement_zero_moment"
    s.bc_south = "zero_displacement_zero_moment"
    s.initialize()
    s.run()
    s.finalize()

    err = np.max(np.abs(s.w - w_ex)) / np.max(np.abs(w_ex))
    return s.w, w_ex, err, dx


# ===========================================================================
# 1-D analysis
# ===========================================================================

print("=" * 66)
print("  MMS convergence: zero_displacement_zero_moment — 1-D")
print("=" * 66)
print(f"  L = {L/1e3:.0f} km,  te = {te/1e3:.0f} km,  "
      f"alpha = {alpha/1e3:.0f} km")
print(f"  W0 = {W0:.0f} m  ->  |w_exact|_max = {W0:.1f} m")
print()

nx_vals_1d = [26, 51, 101, 201, 401]

errs_1d = []
dx_vals_1d = []

print(f"  {'nx':>4}  {'dx [km]':>8}  {'error':>12}")
print(f"  {'-'*4}  {'-'*8}  {'-'*12}")

for nx in nx_vals_1d:
    _, _, e, dx, _ = run_1d(nx)
    errs_1d.append(e)
    dx_vals_1d.append(dx)
    print(f"  {nx:>4}  {dx/1e3:>8.2f}  {e:>12.3e}")

dx_arr_1d   = np.array(dx_vals_1d)
errs_arr_1d = np.array(errs_1d)

n_fit = len(nx_vals_1d) // 2
slope_1d = np.polyfit(np.log(dx_arr_1d[-n_fit:]), np.log(errs_arr_1d[-n_fit:]), 1)[0]

print()
print(f"  Convergence order (finest {n_fit} points): O(dx^{slope_1d:.2f})")
print()

# Reference run for figure
nx_ref = 201
w_num_1d, w_ex_1d, err_ref_1d, dx_ref_1d, xi_ref = run_1d(nx_ref)
x_km = xi_ref * L / 1e3


# ===========================================================================
# 2-D analysis
# ===========================================================================

print("=" * 66)
print("  MMS convergence: zero_displacement_zero_moment — 2-D")
print("=" * 66)
print(f"  (2-D implementation used Dirichlet decoupling from the start.)")
print()

n_vals_2d = [26, 51, 101, 201]

errs_2d  = []
dx_vals_2d = []

print(f"  {'n':>4}  {'dx [km]':>8}  {'error':>12}")
print(f"  {'-'*4}  {'-'*8}  {'-'*12}")

for n in n_vals_2d:
    _, _, e, dx = run_2d(n, n)
    errs_2d.append(e)
    dx_vals_2d.append(dx)
    print(f"  {n:>4}  {dx/1e3:>8.2f}  {e:>12.3e}")

dx_arr_2d   = np.array(dx_vals_2d)
errs_arr_2d = np.array(errs_2d)

n_fit_2d = len(n_vals_2d) // 2
slope_2d = np.polyfit(np.log(dx_arr_2d[-n_fit_2d:]), np.log(errs_arr_2d[-n_fit_2d:]), 1)[0]

print()
print(f"  Convergence order (finest {n_fit_2d} points): O(dx^{slope_2d:.2f})")
print()

# Reference run for figure
n_ref = 101
w_num_2d, w_ex_2d, err_ref_2d, dx_ref_2d = run_2d(n_ref, n_ref)
xi_ref_2d = np.arange(n_ref) / (n_ref - 1)
x_km_2d   = xi_ref_2d * L / 1e3


# ===========================================================================
# Figures
# ===========================================================================

# 1-D figure
fig1, axes1 = plt.subplots(1, 3, figsize=(14, 4))

ax = axes1[0]
ax.plot(x_km, w_ex_1d,  "k-",  lw=2,   label="Exact (MMS)")
ax.plot(x_km, w_num_1d, "b--", lw=1.5, label=f"Numerical (err={err_ref_1d:.2e})")
ax.set_xlabel("x [km]")
ax.set_ylabel("w [m]")
ax.set_title("Deflection profile")
ax.legend(fontsize=8)

ax = axes1[1]
ax.plot(x_km, w_num_1d - w_ex_1d, "b-", lw=1.5)
ax.axhline(0, color="k", lw=0.5)
ax.set_xlabel("x [km]")
ax.set_ylabel("w_num − w_exact [m]")
ax.set_title(f"Residual at dx = {dx_ref_1d/1e3:.1f} km")

ax = axes1[2]
ax.loglog(dx_arr_1d/1e3, errs_arr_1d, "bs-", label=f"O(dx^{slope_1d:.2f})")
ref_x = np.array([dx_arr_1d.min(), dx_arr_1d.max()]) / 1e3
ax.loglog(ref_x, 3e-4 * (ref_x / ref_x[0])**2, "k:", lw=0.8, label="O(dx²)")
ax.set_xlabel("dx [km]")
ax.set_ylabel("L-inf relative error")
ax.set_title("Convergence with grid refinement")
ax.legend(fontsize=8)

fig1.suptitle(
    "MMS convergence: zero_displacement_zero_moment  (1-D)\n"
    f"te={te/1e3:.0f} km,  L={L/1e3:.0f} km,  "
    f"alpha={alpha/1e3:.0f} km,  W0={W0:.0f} m",
    fontsize=10,
)
fig1.tight_layout()
fig1.savefig("analysis/results/pinned_bc_error_1d.png", dpi=150, bbox_inches="tight")
print("  Figure saved to analysis/results/pinned_bc_error_1d.png")


# 2-D figure
fig2, axes2 = plt.subplots(1, 3, figsize=(14, 4))
mid = n_ref // 2

ax = axes2[0]
ax.plot(x_km_2d, w_ex_2d[mid, :],  "k-",  lw=2,   label="Exact (MMS)")
ax.plot(x_km_2d, w_num_2d[mid, :], "b--", lw=1.5, label=f"Numerical (err={err_ref_2d:.2e})")
ax.set_xlabel("x [km]")
ax.set_ylabel("w [m]  (centre row)")
ax.set_title("Deflection profile")
ax.legend(fontsize=8)

ax = axes2[1]
ax.plot(x_km_2d, (w_num_2d - w_ex_2d)[mid, :], "b-", lw=1.5)
ax.axhline(0, color="k", lw=0.5)
ax.set_xlabel("x [km]")
ax.set_ylabel("w_num − w_exact [m]  (centre row)")
ax.set_title(f"Residual at dx = {dx_ref_2d/1e3:.1f} km")

ax = axes2[2]
ax.loglog(dx_arr_2d/1e3, errs_arr_2d, "bs-", label=f"O(dx^{slope_2d:.2f})")
ref_x2 = np.array([dx_arr_2d.min(), dx_arr_2d.max()]) / 1e3
ax.loglog(ref_x2, 3e-4 * (ref_x2 / ref_x2[0])**2, "k:", lw=0.8, label="O(dx²)")
ax.set_xlabel("dx [km]")
ax.set_ylabel("L-inf relative error")
ax.set_title("Convergence with grid refinement")
ax.legend(fontsize=8)

fig2.suptitle(
    "MMS convergence: zero_displacement_zero_moment  (2-D)\n"
    f"te={te/1e3:.0f} km,  L={L/1e3:.0f} km,  "
    f"alpha={alpha/1e3:.0f} km,  W0={W0:.0f} m",
    fontsize=10,
)
fig2.tight_layout()
fig2.savefig("analysis/results/pinned_bc_error_2d.png", dpi=150, bbox_inches="tight")
print("  Figure saved to analysis/results/pinned_bc_error_2d.png")
