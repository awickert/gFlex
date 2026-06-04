"""
Greenland lithospheric flexure under ice sheet loading.

Ice thickness:      BedMachine Greenland v6 (Morlighem et al. 2017/2022,
                    NSIDC IDBMG4), downsampled from 150 m to ~10 km.
Elastic thickness:  Steffen et al. (Zenodo doi:10.5281/zenodo.18403685),
                    reprojected from geographic lon/lat to EPSG:3413.

Usage
-----
Download BedMachineGreenland-v6.nc from NSIDC (Earthdata login required):
    https://nsidc.org/data/idbmg4

Then run from the repository root::

    python examples/greenland_flexure.py

The Steffen et al. Te file is downloaded automatically from Zenodo on first run.

References
----------
Morlighem M. et al. (2017), BedMachine v3: Complete bed topography and ocean
    bathymetry mapping of Greenland, GRL 44, 11051–11061,
    doi:10.1002/2017GL074954.
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
from cmcrameri import cm as cmc

from gflex import F2D, pad_domain

# ---- Paths (adjust BEDMACHINE to your local copy) ----
BEDMACHINE = os.path.expanduser("~/Downloads/BedMachineGreenland-v6.nc")
_script_dir = os.path.dirname(os.path.abspath(__file__))
TE_FILE = os.path.join(_script_dir, "data", "Greenland_Te.dat")
TE_URL = "https://zenodo.org/records/18403685/files/Greenland_Te.dat"

# ---- Physical parameters ----
G = 9.81           # m s⁻²
RHO_ICE = 917.0    # kg m⁻³  (BedMachine nominal)
RHO_MANTLE = 3300.0
RHO_FILL = 0.0     # air: conservative (no isostatic infill above plate)
E = 65e9           # Pa  (Young's modulus)
NU = 0.25          # Poisson's ratio

# ---- Step 1: Download Te data if needed ----
os.makedirs(os.path.dirname(TE_FILE), exist_ok=True)
if not os.path.exists(TE_FILE):
    print("Downloading Greenland Te from Zenodo (doi:10.5281/zenodo.18403685)...")
    urllib.request.urlretrieve(TE_URL, TE_FILE)

# ---- Step 2: Load BedMachine, subsample to ~10 km ----
TARGET_DX = 10_000.0  # m

print("Loading BedMachine Greenland v6...")
with xr.open_dataset(BEDMACHINE) as ds:
    x_native = ds.x.values.astype(float)
    y_native = ds.y.values.astype(float)
    stride = max(1, round(TARGET_DX / abs(x_native[1] - x_native[0])))

    x = x_native[::stride]
    y = y_native[::stride]
    sl = slice(None, None, stride)
    thickness = ds["thickness"].isel(y=sl, x=sl).values.astype(float)
    mask = ds["mask"].isel(y=sl, x=sl).values

dx = abs(x[1] - x[0])
dy = abs(y[1] - y[0])
ny, nx = len(y), len(x)
print(f"  Grid: {ny}×{nx} cells, dx = dy = {dx:.0f} m")

# Grounded ice only (mask==2); floating ice is hydrostatically supported
ice = np.where(mask == 2, thickness, 0.0)
qs = RHO_ICE * G * ice   # Pa

# ---- Step 3: Interpolate Steffen Te onto the 10 km grid ----
print("Loading and interpolating elastic thickness...")
te_data = np.loadtxt(TE_FILE)   # columns: lon, lat, Te (km), uncertainty (km)
lon_pts, lat_pts, te_km = te_data[:, 0], te_data[:, 1], te_data[:, 2]

lon_u = np.unique(lon_pts)
lat_u = np.unique(lat_pts)
dlon = lon_u[1] - lon_u[0]
dlat = lat_u[1] - lat_u[0]

te_ll = np.full((len(lat_u), len(lon_u)), np.nan)
j_idx = np.round((lon_pts - lon_u[0]) / dlon).astype(int)
i_idx = np.round((lat_pts - lat_u[0]) / dlat).astype(int)
te_ll[i_idx, j_idx] = te_km

# Fill gaps (ocean / outside coverage) with the domain median
te_fill_km = float(np.nanmedian(te_ll))
te_ll = np.where(np.isnan(te_ll), te_fill_km, te_ll)

te_interp = RegularGridInterpolator(
    (lat_u, lon_u), te_ll,
    method="linear", bounds_error=False, fill_value=te_fill_km,
)

# Convert EPSG:3413 grid → geographic for Te lookup
XX, YY = np.meshgrid(x, y)
transformer = Transformer.from_crs("EPSG:3413", "EPSG:4326", always_xy=True)
lon_grid, lat_grid = transformer.transform(XX.ravel(), YY.ravel())
te_proj = te_interp(
    np.column_stack([lat_grid, lon_grid])
).reshape(ny, nx) * 1000.0   # km → m

print(f"  Te range: {te_proj.min()/1e3:.1f}–{te_proj.max()/1e3:.1f} km")

# ---- Step 4: Pad domain and run gFlex ----
pad_kw = dict(dx=dx, dy=dy, E=E, nu=NU, rho_m=RHO_MANTLE, rho_fill=RHO_FILL, g=G)
te_pad, qs_pad, p = pad_domain(te_proj, qs, **pad_kw)
print(f"Padding: {p} cells on each side → padded grid {te_pad.shape[0]}×{te_pad.shape[1]}")

print("Running gFlex F2D (fd solver)...")
flex = F2D()
flex.quiet = True
flex.method = "fd"
flex.g = G
flex.E = E
flex.nu = NU
flex.rho_m = RHO_MANTLE
flex.rho_fill = RHO_FILL
flex.te = te_pad
flex.qs = qs_pad
flex.dx = dx
flex.dy = dy
flex.bc_west = flex.bc_east = flex.bc_north = flex.bc_south = "zero_displacement_zero_slope"
flex.initialize()
flex.run()

w = flex.w[p:-p, p:-p].copy()   # trim padding; m, negative = downward
flex.finalize()

print(f"  Maximum depression:  {-w.min():.0f} m")
print(f"  Maximum forebulge:   {w.max():.0f} m")

# ---- Step 5: Plot ----
ext = [x[0] / 1e3, x[-1] / 1e3, y[-1] / 1e3, y[0] / 1e3]  # km
kw = dict(extent=ext, origin="upper", aspect="equal")

fig, axes = plt.subplots(1, 3, figsize=(17, 6))

# Te — Roma (diverging, dark-blue to dark-red; low Te = blue, high Te = red)
im0 = axes[0].imshow(te_proj / 1e3, cmap=cmc.roma, vmin=0, vmax=90, **kw)
plt.colorbar(im0, ax=axes[0], label="$T_e$ (km)", shrink=0.7)
axes[0].set_title("Effective elastic thickness\n(Steffen et al.)")

# Ice thickness — Devon, grounded ice only (non-ice masked → transparent)
ice_masked = np.ma.array(ice, mask=(mask != 2))
im1 = axes[1].imshow(ice_masked, cmap=cmc.devon, vmin=0, **kw)
plt.colorbar(im1, ax=axes[1], label="Ice thickness (m)", shrink=0.7)
axes[1].set_title("Ice thickness\n(BedMachine Greenland v6)")

# Deflection — Vik with TwoSlopeNorm: asymmetric limits, zero at colour centre
# Blue = subsidence (w < 0 downward), red = uplift (w > 0)
norm_w = mcolors.TwoSlopeNorm(vcenter=0, vmin=w.min(), vmax=80)
im2 = axes[2].imshow(w, cmap=cmc.vik, norm=norm_w, **kw)
cb2 = plt.colorbar(im2, ax=axes[2], label="Deflection (m,  − = down)", shrink=0.7)
cb2.set_ticks([-800, -600, -400, -200, 0, 20, 40, 60, 80])
axes[2].set_title("Lithospheric deflection")

for ax in axes:
    ax.set_xlabel("x (km, EPSG:3413)")
    ax.set_ylabel("y (km, EPSG:3413)")

fig.suptitle("Greenland glacial isostasy — gFlex", fontsize=13)
plt.tight_layout()

out = os.path.join(_script_dir, "greenland_flexure.png")
plt.savefig(out, dpi=150, bbox_inches="tight")
print(f"Saved {out}")
plt.close()

# ---- Step 6: Uniform-Te comparison ----
# Mean Te under grounded ice only
te_mean = float(np.mean(te_proj[mask == 2]))
print(f"\nMean Te under grounded ice: {te_mean/1e3:.1f} km")
te_uniform = np.full_like(te_proj, te_mean)

te_unif_pad, qs_unif_pad, p2 = pad_domain(te_uniform, qs, **pad_kw)

flex2 = F2D()
flex2.quiet = True
flex2.method = "fd"
flex2.g = G
flex2.E = E
flex2.nu = NU
flex2.rho_m = RHO_MANTLE
flex2.rho_fill = RHO_FILL
flex2.te = te_unif_pad
flex2.qs = qs_unif_pad
flex2.dx = dx
flex2.dy = dy
flex2.bc_west = flex2.bc_east = flex2.bc_north = flex2.bc_south = "zero_displacement_zero_slope"
flex2.initialize()
print("Solving uniform-Te case...")
flex2.run()

w_unif = flex2.w[p2:-p2, p2:-p2].copy()
flex2.finalize()

dw = w - w_unif   # positive where variable Te gives less subsidence than uniform

print(f"  Uniform-Te depression:  {-w_unif.min():.0f} m")
print(f"  Difference max/min:     {dw.max():.0f} / {dw.min():.0f} m")

fig2, axes2 = plt.subplots(1, 3, figsize=(17, 6))

norm_u = mcolors.TwoSlopeNorm(vcenter=0, vmin=w_unif.min(), vmax=80)
im_u = axes2[0].imshow(w_unif, cmap=cmc.vik, norm=norm_u, **kw)
cb_u = plt.colorbar(im_u, ax=axes2[0], label="Deflection (m,  − = down)", shrink=0.7)
cb_u.set_ticks([-800, -600, -400, -200, 0, 20, 40, 60, 80])
axes2[0].set_title(f"Uniform $T_e$ = {te_mean/1e3:.0f} km")

norm_v = mcolors.TwoSlopeNorm(vcenter=0, vmin=w.min(), vmax=80)
im_v = axes2[1].imshow(w, cmap=cmc.vik, norm=norm_v, **kw)
cb_v = plt.colorbar(im_v, ax=axes2[1], label="Deflection (m,  − = down)", shrink=0.7)
cb_v.set_ticks([-800, -600, -400, -200, 0, 20, 40, 60, 80])
axes2[1].set_title("Variable $T_e$ (Steffen et al.)")

dw_abs = max(abs(dw.min()), dw.max())
im_d = axes2[2].imshow(dw, cmap=cmc.vik, vmin=-dw_abs, vmax=dw_abs, **kw)
plt.colorbar(im_d, ax=axes2[2],
             label="Δw = variable − uniform (m,  − = more subsidence)", shrink=0.7)
axes2[2].set_title("Difference (variable − uniform $T_e$)")

for ax in axes2:
    ax.set_xlabel("x (km, EPSG:3413)")
    ax.set_ylabel("y (km, EPSG:3413)")

fig2.suptitle("Effect of spatially variable $T_e$ on Greenland flexure", fontsize=13)
plt.tight_layout()

out2 = os.path.join(_script_dir, "greenland_flexure_comparison.png")
plt.savefig(out2, dpi=150, bbox_inches="tight")
print(f"Saved {out2}")
plt.close()
