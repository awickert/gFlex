#!/usr/bin/env python
"""Broken-plate flexure: prescribed shear force at one end (1-D FD).

This script demonstrates the dict-style boundary condition introduced in
gFlex v2.0, which lets you prescribe non-zero values for displacement,
slope, bending moment, or shear force at any edge.

Physical scenario
-----------------
A semi-infinite oceanic plate is pulled down at the trench by a vertical
force V0 (slab pull).  The plate bends downward at the trench and a
forebulge (outer rise) develops at approximately π·α/2 from the trench.
This is the classical "broken-plate" scenario of Turcotte and Schubert (2002).

Boundary conditions
-------------------
  bc_west = {"moment": 0.0, "shear": V0}   -- zero moment, prescribed shear (trench end)
  bc_east = "zero_moment_zero_shear"        -- free far end

The dict BC uses the same finite-difference stencil as the equivalent string
BC ("zero_moment_zero_shear").  The prescribed V0 enters as an additive
correction to the right-hand side.  Setting V0 = 0 reproduces the string BC
exactly.

Analytical solution
-------------------
For a plate on a Winkler foundation with no distributed load and a shear
force V0 applied at x = 0:

    w(x) = [V0 / (2 D λ³)] · e^{-λx} · cos(λx)

    λ = (Δρ g / 4D)^{1/4},  α = 1/λ  (flexural parameter)

Key outputs (V0 = -1e8 N/m, Te = 30 km):
    w(0)           ≈ -93 mm   (downward — the "trench")
    forebulge peak ≈ +6 mm    at x ≈ 156 km  (~2.4 α from the trench)
    λ error        < 0.02 %   (L-inf vs. analytical)
"""

import warnings

import numpy as np

from gflex import F1D

# ---------------------------------------------------------------------------
# Physical parameters
# ---------------------------------------------------------------------------
E        = 65e9     # Pa — Young's modulus
nu       = 0.25
rho_m    = 3300.0   # kg m⁻³ — mantle density
rho_fill = 0.0      # kg m⁻³ — air infill
g        = 9.8      # m s⁻²
Te       = 30e3     # m  — elastic thickness

# Derived quantities
D     = E * Te**3 / (12.0 * (1.0 - nu**2))           # flexural rigidity [N·m]
lam   = ((rho_m - rho_fill) * g / (4.0 * D)) ** 0.25  # 1/α [m⁻¹]
alpha = 1.0 / lam                                      # flexural parameter [m]

print(f"Flexural parameter α = {alpha/1e3:.1f} km")
print(f"Flexural rigidity  D = {D:.3e} N·m")

# ---------------------------------------------------------------------------
# Grid: 10 α long, 401 nodes  (Δx ≈ α/40, so O(Δx²) error ≈ 0.06 %)
# ---------------------------------------------------------------------------
nx = 401
dx = 10.0 * alpha / (nx - 1)
x  = np.arange(nx) * dx
print(f"Grid: {nx} cells, Δx = {dx/1e3:.2f} km, domain = {x[-1]/1e3:.0f} km")

# ---------------------------------------------------------------------------
# Prescribed boundary load
# ---------------------------------------------------------------------------
V0 = -1e8   # N/m  — slab-pull shear force (negative = downward at west end)

# ---------------------------------------------------------------------------
# Run gFlex
# ---------------------------------------------------------------------------
flex = F1D()
flex.quiet    = True
flex.method   = "fd"
flex.g        = g
flex.E        = E
flex.nu       = nu
flex.rho_m    = rho_m
flex.rho_fill = rho_fill
flex.te       = Te
flex.qs       = np.zeros(nx)          # no distributed load
flex.dx       = dx
flex.bc_west  = {"moment": 0.0, "shear": V0}   # trench: zero moment, prescribed shear
flex.bc_east  = "zero_moment_zero_shear"         # free far end

with warnings.catch_warnings():
    warnings.simplefilter("ignore")   # suppress proximity warning for demonstration
    flex.initialize()
    flex.run()

w = flex.w.copy()
flex.finalize()

# ---------------------------------------------------------------------------
# Report results
# ---------------------------------------------------------------------------
i_bulge = np.argmax(w)
print(f"\nDeflection at loaded end (x=0): {w[0]*1e3:.1f} mm  (downward = trench)")
print(f"Forebulge peak:                  {w[i_bulge]*1e3:.1f} mm  at x = {x[i_bulge]/1e3:.0f} km")

# Compare with analytical solution
w_analytical = (V0 / (2.0 * D * lam**3)) * np.exp(-lam * x) * np.cos(lam * x)
l_inf_err = np.max(np.abs(w - w_analytical)) / np.max(np.abs(w_analytical))
print(f"L-inf error vs. analytical:      {l_inf_err:.4%}")

# ---------------------------------------------------------------------------
# Plot  (requires matplotlib)
# ---------------------------------------------------------------------------
try:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(x / 1e3, w * 1e3, "C0",     lw=2,   label="gFlex (FD)")
    ax.plot(x / 1e3, w_analytical * 1e3, "C1--", lw=1.5, label="Analytical")
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("Distance from trench [km]")
    ax.set_ylabel("Deflection [mm]")
    ax.set_title(
        f"Broken-plate flexure: V₀ = {V0:.0e} N/m, Tₑ = {Te/1e3:.0f} km"
    )
    ax.legend()
    plt.tight_layout()
    plt.show()

except ImportError:
    print("(matplotlib not available — skipping plot)")
