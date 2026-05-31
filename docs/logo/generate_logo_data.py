#!/usr/bin/env python
"""Generate gFlex deflection data for logo rendering.

Variable Te: smooth sigmoid transition from Te=15 km (west/soft) to
Te=35 km (east/stiff), so the deflection is asymmetric — deeper and
sharper on the soft side, broader and shallower with a more prominent
forebulge on the stiff side.

Domain is 150×150 cells at 5 km spacing (750 km), giving ≥7α coverage
on the stiff side and avoiding boundary-condition artefacts.
"""

import numpy as np
from gflex import F2D

# ── Grid ──────────────────────────────────────────────────────────────────────
nrows = ncols = 150
dx = dy = 5000.0   # 5 km → 750 km domain

# ── Physical parameters ───────────────────────────────────────────────────────
E        = 65e9
nu       = 0.25
rho_m    = 3300.0
rho_fill = 0.0
g        = 9.8

cx = cy = 0.5 * nrows * dx

xi = (np.arange(ncols) + 0.5) * dx
yi = (np.arange(nrows) + 0.5) * dy
XX, YY = np.meshgrid(xi, yi)

# ── Variable Te: sigmoid from 15 km (west) → 50 km (east) ────────────────────
Te_min, Te_max = 15e3, 35e3
x_norm = (XX - cx) / (8 * dx)                    # 8-cell (~40 km) transition
Te_grid = Te_min + (Te_max - Te_min) * 0.5 * (1.0 + np.tanh(x_norm))

D_grid = E * Te_grid**3 / (12.0 * (1.0 - nu**2))
drho   = rho_m - rho_fill

# Reference alpha using mean Te for domain-sizing guidance
Te_mean = Te_grid.mean()
D_mean  = E * Te_mean**3 / (12.0 * (1.0 - nu**2))
alpha   = (D_mean / (drho * g))**0.25
print(f"Mean Te = {Te_mean/1e3:.1f} km  |  alpha_2D (mean) = {alpha/1e3:.1f} km")
print(f"Domain = {nrows*dx/1e3:.0f} km = {nrows*dx/alpha:.1f} alpha")

# ── Load: circular plateau, radius ≈ 0.4 alpha ────────────────────────────────
load_radius = 0.4 * alpha
R = np.sqrt((XX - cx)**2 + (YY - cy)**2)
qs = np.where(R <= load_radius, 3e7, 0.0)   # 30 MPa flat-topped load

# ── Run gFlex ─────────────────────────────────────────────────────────────────
flex = F2D()
flex.quiet = True
flex.method = "FD"
flex.g = g;  flex.E = E;  flex.nu = nu
flex.rho_m = rho_m;  flex.rho_fill = rho_fill
flex.te = Te_grid
flex.qs = qs
flex.dx = dx;  flex.dy = dy
flex.bc_west = flex.bc_east = flex.bc_south = flex.bc_north = "zero_moment_zero_shear"
flex.initialize()
flex.run()
flex.finalize()

w = flex.w
print(f"Max subsidence : {w.min():+.2f} m")
print(f"Max uplift     : {w.max():+.2f} m")
print(f"Uplift/subsid  : {abs(w.max()/w.min())*100:.2f}%")

np.savez("/tmp/gflex_logo_data.npz",
         w=w, x=xi, y=yi, qs=qs, Te=Te_grid,
         cx=cx, cy=cy, load_radius=load_radius, alpha=alpha, dx=dx, dy=dy)
print("Saved → /tmp/gflex_logo_data.npz")

# ── Write pure-Python mesh file for Blender ───────────────────────────────────
PLATE_BU = 10.0
m_per_bu = (ncols * dx) / PLATE_BU
Z_EX     = 1000.0
z_scale  = Z_EX / m_per_bu
half     = PLATE_BU / 2.0

w_min_bu = float(w.min()) * z_scale
w_max_bu = float(w.max()) * z_scale

# Clip colormap so forebulge is clearly visible:
# show ±8× the forebulge amplitude (forebulge sits at 87.5% of the warm scale)
clip_bu = max(abs(w_max_bu) * 8.0, abs(w_min_bu) * 0.015)
print(f"Colormap clip  : ±{clip_bu:.4f} BU  (forebulge at {w_max_bu:.4f} BU)")

verts = []
for j in range(nrows):
    for i in range(ncols):
        bx = (xi[i] / (ncols * dx)) * PLATE_BU - half
        by = (yi[j] / (nrows * dy)) * PLATE_BU - half
        bz = float(w[j, i]) * z_scale
        verts.append((bx, by, bz))

faces = []
for j in range(nrows - 1):
    for i in range(ncols - 1):
        v0 = j * ncols + i
        faces.append((v0, v0+1, (j+1)*ncols+i+1, (j+1)*ncols+i))

load_bx     = (cx / (ncols * dx)) * PLATE_BU - half
load_by     = (cy / (nrows * dy)) * PLATE_BU - half
load_br     = (load_radius / (ncols * dx)) * PLATE_BU * 2.8   # wider visual
base_z      = float(w[nrows//2, ncols//2]) * z_scale
load_height = abs(w_min_bu) * 0.25   # short plateau

# Per-vertex Te in Blender units (needed for floor-plane and slab variants)
te_verts = [float(Te_grid[j, i]) for j in range(nrows) for i in range(ncols)]

with open("/tmp/gflex_logo_mesh.py", "w") as fh:
    fh.write(f"w_min_bu    = {w_min_bu}\n")
    fh.write(f"w_max_bu    = {w_max_bu}\n")
    fh.write(f"clip_bu     = {clip_bu}\n")
    fh.write(f"load_bx     = {load_bx}\n")
    fh.write(f"load_by     = {load_by}\n")
    fh.write(f"load_br     = {load_br}\n")
    fh.write(f"base_z      = {base_z}\n")
    fh.write(f"load_height = {load_height}\n")
    fh.write(f"nrows       = {nrows}\n")
    fh.write(f"ncols       = {ncols}\n")
    fh.write(f"te_min      = {float(Te_grid.min())}\n")
    fh.write(f"te_max      = {float(Te_grid.max())}\n")
    fh.write(f"vertices = {verts!r}\n")
    fh.write(f"faces    = {faces!r}\n")
    fh.write(f"te_verts = {te_verts!r}\n")

print(f"Saved → /tmp/gflex_logo_mesh.py  ({len(verts)} verts, {len(faces)} quads)")
