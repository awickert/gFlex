"""
Nested-domain flexure: Greenland ice sheet + hypothetical Irminger Sea seamount.

Demonstrates inhomogeneous (prescribed-value) boundary conditions to couple a
coarse regional model to a fine local model.

Coarse model (10 km):
  Ice-sheet loading over Greenland (BedMachine v6) with spatially variable Te
  (Steffen et al. 2026).  Clamped boundary conditions (zero_displacement_zero_slope)
  with domain padding.

Fine model (2 km, 300 × 300 km):
  Sub-domain centred at ~64°N 40°W on the SE Greenland coast, straddling the
  ice-sheet margin.  Boundary conditions are interpolated from the coarse model
  (displacement + slope).  The ice load within the fine domain is resampled
  from the coarse BedMachine grid to 2 km to drive the local plate response.

  A hypothetical Gaussian seamount (h_peak = 2 km, σ = 20 km) is placed at the
  domain centre to illustrate how coastal volcanism would perturb the plate on
  top of the existing ice-sheet deflection.

Usage
-----
Download BedMachineGreenland-v6.nc from NSIDC (Earthdata login required):
    https://nsidc.org/data/idbmg4

Then run from the repository root::

    python examples/greenland_volcano_nested.py

The Steffen et al. Te file is downloaded automatically from Zenodo on first run.

References
----------
Morlighem M. et al. (2017), BedMachine v3: Complete bed topography and ocean
    bathymetry mapping of Greenland, GRL 44, 11051–11061,
    doi:10.1002/2017GL074954.
Morlighem M. et al. (2025), IceBridge BedMachine Greenland, Version 6,
    NSIDC IDBMG4, doi:10.5067/6B6B225B8V2D.
Steffen R., Audet P., and Lund B. (2018), Weakened lithosphere beneath
    Greenland inferred from effective elastic thickness: A hot spot effect?,
    GRL 45(10), 4733–4742, doi:10.1029/2017GL076885.
Steffen R., Audet P., Strykowski G., and Lund B. (2026), Greenland Moho and
    effective elastic thickness models based on satellite gravity data,
    Zenodo, doi:10.5281/zenodo.18403685.
"""

import os
import urllib.request
import numpy as np
import xarray as xr
from scipy.interpolate import RegularGridInterpolator
from pyproj import Transformer
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from cmcrameri import cm as cmc

from gflex import F2D, pad_domain

# ---- Paths ----
BEDMACHINE = os.path.expanduser("~/Downloads/BedMachineGreenland-v6.nc")
_script_dir = os.path.dirname(os.path.abspath(__file__))
TE_FILE = os.path.join(_script_dir, "data", "Greenland_Te.dat")
TE_URL  = "https://zenodo.org/records/18403685/files/Greenland_Te.dat"

# ---- Physical parameters ----
G          = 9.81
RHO_ICE    = 917.0
RHO_MANTLE = 3300.0
RHO_FILL   = 0.0
RHO_BASALT = 2900.0   # submarine seamount (basalt)
RHO_SEA    = 1025.0   # seawater
E          = 65e9
NU         = 0.25

# ---- Fine sub-domain parameters (SE Greenland coast, ~64°N 40°W) ----
CX      = 250_000.0    # EPSG:3413 x-centre, m
CY      = -2_850_000.0 # EPSG:3413 y-centre, m
FINE_DX = 2_000.0      # fine grid spacing, m
N_FINE  = 150          # cells per side → 300 km domain

# ---- Seamount parameters ----
H_PEAK  = 2_000.0    # peak height above seafloor, m
SIGMA   = 20_000.0   # Gaussian σ, m

# ---- Step 1: Download Te if needed ----
os.makedirs(os.path.dirname(TE_FILE), exist_ok=True)
if not os.path.exists(TE_FILE):
    print("Downloading Greenland Te from Zenodo (doi:10.5281/zenodo.18403685)...")
    urllib.request.urlretrieve(TE_URL, TE_FILE)

# ---- Step 2: Load BedMachine at 10 km ----
TARGET_DX = 10_000.0
print("Loading BedMachine Greenland v6...")
with xr.open_dataset(BEDMACHINE) as ds:
    x_native = ds.x.values.astype(float)
    y_native = ds.y.values.astype(float)
    stride   = max(1, round(TARGET_DX / abs(x_native[1] - x_native[0])))
    x = x_native[::stride]
    y = y_native[::stride]
    sl = slice(None, None, stride)
    thickness = ds["thickness"].isel(y=sl, x=sl).values.astype(float)
    mask      = ds["mask"].isel(y=sl, x=sl).values

dx_c = abs(x[1] - x[0])   # coarse grid spacing (both axes)
dy_c = abs(y[1] - y[0])
ny_c, nx_c = len(y), len(x)
print(f"  Coarse grid: {ny_c}×{nx_c}, dx = dy = {dx_c:.0f} m")

ice = np.where(mask == 2, thickness, 0.0)   # grounded ice only
qs_coarse = RHO_ICE * G * ice

# ---- Step 3: Interpolate Steffen Te onto coarse grid ----
print("Interpolating elastic thickness onto coarse grid...")
te_data   = np.loadtxt(TE_FILE)
lon_u     = np.unique(te_data[:, 0])
lat_u     = np.unique(te_data[:, 1])
dlon = lon_u[1] - lon_u[0]
dlat = lat_u[1] - lat_u[0]
te_ll = np.full((len(lat_u), len(lon_u)), np.nan)
j_idx = np.round((te_data[:, 0] - lon_u[0]) / dlon).astype(int)
i_idx = np.round((te_data[:, 1] - lat_u[0]) / dlat).astype(int)
te_ll[i_idx, j_idx] = te_data[:, 2]
te_fill_km = float(np.nanmedian(te_ll))
te_ll = np.where(np.isnan(te_ll), te_fill_km, te_ll)

te_ll_interp = RegularGridInterpolator(
    (lat_u, lon_u), te_ll, method="linear",
    bounds_error=False, fill_value=te_fill_km)

transformer = Transformer.from_crs("EPSG:3413", "EPSG:4326", always_xy=True)
XX_c, YY_c = np.meshgrid(x, y)
lon_c, lat_c = transformer.transform(XX_c.ravel(), YY_c.ravel())
te_coarse = te_ll_interp(
    np.column_stack([lat_c, lon_c])).reshape(ny_c, nx_c) * 1000.0  # km → m

# ---- Step 4: Run coarse Greenland model ----
pad_kw_c = dict(dx=dx_c, dy=dy_c, E=E, nu=NU, rho_m=RHO_MANTLE,
                rho_fill=RHO_FILL, g=G)
te_pad, qs_pad, p = pad_domain(te_coarse, qs_coarse, **pad_kw_c)
print(f"Running coarse gFlex (padded grid "
      f"{te_pad.shape[0]}×{te_pad.shape[1]}, p={p})...")

flex_c = F2D()
flex_c.quiet    = True
flex_c.method   = "fd"
flex_c.g        = G
flex_c.E        = E
flex_c.nu       = NU
flex_c.rho_m    = RHO_MANTLE
flex_c.rho_fill = RHO_FILL
flex_c.te       = te_pad
flex_c.qs       = qs_pad
flex_c.dx       = dx_c
flex_c.dy       = dy_c
flex_c.bc_west = flex_c.bc_east = flex_c.bc_north = flex_c.bc_south = \
    "zero_displacement_zero_slope"
flex_c.initialize()
flex_c.run()
w_coarse = flex_c.w[p:-p, p:-p].copy()
flex_c.finalize()
print(f"  Coarse: max depression {-w_coarse.min():.0f} m, "
      f"forebulge {w_coarse.max():.0f} m")

# ---- Step 5: Define fine sub-domain grid ----
# Row 0 = northernmost, row N_FINE-1 = southernmost  (consistent with coarse).
# x increases eastward, y decreases southward (EPSG:3413 convention).
fine_x = CX + (np.arange(N_FINE) - N_FINE // 2) * FINE_DX     # increasing east
fine_y = CY + (N_FINE // 2 - np.arange(N_FINE)) * FINE_DX     # decreasing south

print(f"\nFine sub-domain (~64°N 40°W, SE Greenland coast):")
print(f"  x: {fine_x[0]/1e3:.0f} – {fine_x[-1]/1e3:.0f} km, "
      f"y: {fine_y[-1]/1e3:.0f} – {fine_y[0]/1e3:.0f} km")

# ---- Step 6: Interpolate Te onto fine grid ----
XX_f, YY_f = np.meshgrid(fine_x, fine_y)
lon_f, lat_f = transformer.transform(XX_f.ravel(), YY_f.ravel())
te_fine = te_ll_interp(
    np.column_stack([lat_f, lon_f])).reshape(N_FINE, N_FINE) * 1000.0  # km → m
print(f"  Fine Te: {te_fine.min()/1e3:.0f} – {te_fine.max()/1e3:.0f} km")

# ---- Step 6b: Resample coarse ice load to fine grid ----
# The fine sub-domain straddles the SE Greenland ice-sheet margin.  Including
# the ice load explicitly in the fine model drives the local plate response; the
# coarse BCs supply the far-field boundary constraint.
# (Note: y_inc is built in Step 8 below — build it here for reuse.)
y_inc = y[::-1]
_thick_rgi = RegularGridInterpolator(
    (y_inc, x), thickness[::-1, :],
    method="linear", bounds_error=False, fill_value=0.)
_mask_rgi  = RegularGridInterpolator(
    (y_inc, x), mask[::-1, :].astype(float),
    method="nearest", bounds_error=False, fill_value=0.)

_pts_f = np.column_stack([YY_f.ravel(), XX_f.ravel()])
thick_fine = _thick_rgi(_pts_f).reshape(N_FINE, N_FINE).clip(0)
mask_fine  = _mask_rgi(_pts_f).reshape(N_FINE, N_FINE).round().astype(int)
ice_fine   = np.where(mask_fine == 2, thick_fine, 0.)
qs_ice_fine = RHO_ICE * G * ice_fine
print(f"  Fine ice: {ice_fine.max():.0f} m max thickness, "
      f"{(mask_fine == 2).sum()} grounded-ice cells")

# ---- Step 7: Seamount load ----
R2 = (XX_f - CX)**2 + (YY_f - CY)**2
h_seamount = H_PEAK * np.exp(-R2 / (2.0 * SIGMA**2))
# Net load = excess weight of basalt over displaced seawater.
qs_seamount = (RHO_BASALT - RHO_SEA) * G * h_seamount

# ---- Step 8: Extract coarse-model displacement + slope BCs ----
#
# Convention (same as np.gradient / gFlex Dirichlet BC):
#   west / east:   slope = dw/dx  (positive = eastward, axis=1)
#   north / south: slope = dw/dy  (positive = southward, i.e. increasing
#                                  row-index direction, axis=0)
#
# Both dx_c and dy_c are positive magnitudes, so np.gradient with these
# spacings gives the gradient in the increasing-index direction (east / south).
# That matches gFlex's row-index convention for all four edges.
dw_dx_c = np.gradient(w_coarse, dx_c, axis=1)   # dw/dx, positive = east
dw_dy_c = np.gradient(w_coarse, dy_c, axis=0)   # dw/dy, positive = south

# RegularGridInterpolator requires monotonically increasing axes.
# Coarse y is decreasing → flip rows before building the interpolator.
# (y_inc already defined in Step 6b)
_interp_kw = dict(method="linear", bounds_error=False, fill_value=None)
w_rgi     = RegularGridInterpolator((y_inc, x), w_coarse[::-1, :],     **_interp_kw)
dwdx_rgi  = RegularGridInterpolator((y_inc, x), dw_dx_c[::-1, :],     **_interp_kw)
dwdy_rgi  = RegularGridInterpolator((y_inc, x), dw_dy_c[::-1, :],     **_interp_kw)


def _sample(pts_x, pts_y):
    """Sample coarse w, dw/dx, dw/dy at given (x, y) points."""
    pts = np.column_stack([pts_y, pts_x])   # (y, x) order for RGI
    return w_rgi(pts), dwdx_rgi(pts), dwdy_rgi(pts)


# West edge: x = fine_x[0], y varies north → south (row 0 → row N_FINE-1)
w_w, dwdx_w, _      = _sample(np.full(N_FINE, fine_x[0]),   fine_y)
# East edge: x = fine_x[-1]
w_e, dwdx_e, _      = _sample(np.full(N_FINE, fine_x[-1]),  fine_y)
# North edge: y = fine_y[0], x varies west → east (col 0 → col N_FINE-1)
w_n, _, dwdy_n      = _sample(fine_x, np.full(N_FINE, fine_y[0]))
# South edge: y = fine_y[-1]
w_s, _, dwdy_s      = _sample(fine_x, np.full(N_FINE, fine_y[-1]))

# ---- Step 9: Run fine F2D (helper) ----

def _run_fine(qs_load):
    flex_f = F2D()
    flex_f.quiet    = True
    flex_f.method   = "fd"
    flex_f.g        = G
    flex_f.E        = E
    flex_f.nu       = NU
    flex_f.rho_m    = RHO_MANTLE
    flex_f.rho_fill = RHO_FILL
    flex_f.te       = te_fine
    flex_f.qs       = qs_load
    flex_f.dx       = FINE_DX
    flex_f.dy       = FINE_DX
    flex_f.bc_west  = {"displacement": w_w, "slope": dwdx_w}
    flex_f.bc_east  = {"displacement": w_e, "slope": dwdx_e}
    flex_f.bc_north = {"displacement": w_n, "slope": dwdy_n}
    flex_f.bc_south = {"displacement": w_s, "slope": dwdy_s}
    flex_f.initialize()
    flex_f.run()
    w_out = flex_f.w.copy()
    flex_f.finalize()
    return w_out


print("Running fine F2D — background (ice-sheet load, no seamount)...")
w_bg = _run_fine(qs_ice_fine)

print("Running fine F2D — total (ice-sheet + seamount)...")
w_total = _run_fine(qs_ice_fine + qs_seamount)

w_volcanic = w_total - w_bg
print(f"  Background range:   {w_bg.min():.1f} – {w_bg.max():.1f} m")
print(f"  Total range:        {w_total.min():.1f} – {w_total.max():.1f} m")
print(f"  Seamount signal:    {-w_volcanic.min():.1f} m depression, "
      f"{w_volcanic.max():.1f} m uplift")

# ---- Step 10: Plot — Figure 1: regional context with fine sub-domain box ----
ext_c = [x[0]/1e3,      x[-1]/1e3,  y[-1]/1e3,      y[0]/1e3]
ext_f = [fine_x[0]/1e3, fine_x[-1]/1e3, fine_y[-1]/1e3, fine_y[0]/1e3]
kw_c  = dict(extent=ext_c, origin="upper", aspect="equal")
kw_f  = dict(extent=ext_f, origin="upper", aspect="equal")

# Size figure to the map's natural aspect ratio so the domain fills the axes.
map_w_km = (x[-1] - x[0]) / 1e3
map_h_km = (y[0] - y[-1]) / 1e3
fig_h = 9.0
fig_w = fig_h * map_w_km / map_h_km + 2.0   # extra width for colorbar

fig1, ax_ctx = plt.subplots(1, 1, figsize=(fig_w, fig_h))

norm_c = mcolors.TwoSlopeNorm(vcenter=0, vmin=w_coarse.min(), vmax=80)
im_ctx = ax_ctx.imshow(w_coarse, cmap=cmc.vik, norm=norm_c, **kw_c)
cb_ctx = plt.colorbar(im_ctx, ax=ax_ctx,
                      label="Deflection (m,  − = down)", shrink=0.85)
cb_ctx.set_ticks([-800, -600, -400, -200, 0, 20, 40, 60, 80])

rect = mpatches.Rectangle(
    (fine_x[0]/1e3, fine_y[-1]/1e3),
    (fine_x[-1] - fine_x[0]) / 1e3, (fine_y[0] - fine_y[-1]) / 1e3,
    linewidth=2, edgecolor="k", facecolor="none")
ax_ctx.add_patch(rect)

ax_ctx.set_xlabel("x (km, EPSG:3413)")
ax_ctx.set_ylabel("y (km, EPSG:3413)")
ax_ctx.set_title("Greenland lithospheric deflection — regional model\n"
                 "(box marks the 300 × 300 km fine sub-domain)", fontsize=12)
plt.tight_layout()

out1 = os.path.join(_script_dir, "greenland_volcano_context.png")
plt.savefig(out1, dpi=150, bbox_inches="tight")
print(f"\nSaved {out1}")
plt.close()

# ---- Figure 2: three fine-domain solution panels ----
# Panels 1–2 (background, total) share a colour scale so the seamount imprint
# is visible by direct comparison; panel 3 (seamount signal) has its own scale.
fine_abs = max(abs(w_bg.min()), abs(w_bg.max()),
               abs(w_total.min()), abs(w_total.max()), 1.0)
norm_fine = mcolors.TwoSlopeNorm(vcenter=0, vmin=-fine_abs, vmax=fine_abs)

fig2, axes2 = plt.subplots(1, 3, figsize=(15, 5))
ax_bg, ax_tot, ax_sig = axes2

im_bg = ax_bg.imshow(w_bg, cmap=cmc.vik, norm=norm_fine, **kw_f)
plt.colorbar(im_bg, ax=ax_bg, label="Deflection (m,  − = down)", shrink=0.85)
ax_bg.set_title("Background\n(ice-sheet only, no seamount)")

im_tot = ax_tot.imshow(w_total, cmap=cmc.vik, norm=norm_fine, **kw_f)
plt.colorbar(im_tot, ax=ax_tot, label="Deflection (m,  − = down)", shrink=0.85)
ax_tot.set_title("Total\n(ice-sheet + seamount)")

v_abs = max(abs(w_volcanic.min()), abs(w_volcanic.max()), 0.1)
im_sig = ax_sig.imshow(w_volcanic, cmap=cmc.vik, vmin=-v_abs, vmax=v_abs, **kw_f)
plt.colorbar(im_sig, ax=ax_sig, label="Δw seamount (m,  − = down)", shrink=0.85)
ax_sig.set_title("Seamount signal\n(total − background)")

for ax in axes2:
    ax.set_xlabel("x (km, EPSG:3413)")
    ax.set_ylabel("y (km, EPSG:3413)")

fig2.suptitle(
    "SE Greenland fine sub-domain — nested-domain gFlex (300 × 300 km at 2 km)",
    fontsize=12)
plt.tight_layout()

out2 = os.path.join(_script_dir, "greenland_volcano_nested.png")
plt.savefig(out2, dpi=150, bbox_inches="tight")
print(f"Saved {out2}")
plt.close()
