"""
Pre-process the Greenland teaching example: build gFlex input grids.

This script is run *once* to turn the raw source data (a large ice-sheet
netCDF and a scattered-point elastic-thickness file) into two ready-to-use
gFlex input grids, saved into ``examples/data/``:

    greenland_load.npy   surface load stress q0 = rho_ice * g * H_ice  [Pa]
    greenland_te.npy     effective elastic thickness T_e               [m]

Once these exist, the teaching example runs with nothing but gFlex and NumPy
(see ``greenland.yaml`` + ``run_greenland.py`` and ``greenland_script.py``);
learners do not need the 3 GB netCDF, xarray, pyproj, or scipy.

Both grids are on the same EPSG:3413 raster at ~10 km spacing.  The grid
spacing is written to ``greenland_grid.txt`` and hard-coded into the YAML.

Source data
-----------
Ice thickness:      BedMachine Greenland v6 (Morlighem et al. 2017/2025,
                    NSIDC IDBMG4), downsampled from 150 m to ~10 km.
Elastic thickness:  Steffen et al. (Zenodo doi:10.5281/zenodo.18403685),
                    reprojected from geographic lon/lat to EPSG:3413.

Usage
-----
Download BedMachineGreenland-v6.nc from NSIDC (Earthdata login required):
    https://nsidc.org/data/idbmg4
adjust ``BEDMACHINE`` below to point at your copy, then run from the repo root::

    python examples/preprocess_greenland.py

The Steffen et al. Te point file is downloaded automatically from Zenodo.
"""

import os
import urllib.request

import numpy as np
import xarray as xr
from pyproj import Transformer
from scipy.interpolate import RegularGridInterpolator

from gflex import pad_domain

# ---- Paths (adjust BEDMACHINE to your local copy) ----
BEDMACHINE = os.path.expanduser("~/Downloads/BedMachineGreenland-v6.nc")
_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
TE_POINTS = os.path.join(_dir, "Greenland_Te.dat")
TE_URL = "https://zenodo.org/records/18403685/files/Greenland_Te.dat"

# ---- Physical constants used to convert ice thickness -> load stress ----
G = 9.81           # m s^-2
RHO_ICE = 917.0    # kg m^-3  (BedMachine nominal)

# ---- Plate parameters (must match greenland.yaml; used only for padding) ----
E = 6.5e10         # Pa, Young's modulus
NU = 0.25          # Poisson's ratio
RHO_MANTLE = 3300.0
RHO_FILL = 0.0

TARGET_DX = 10_000.0  # m, subsampling target

os.makedirs(_dir, exist_ok=True)

# ---- Step 1: fetch the Steffen Te point file if needed ----
if not os.path.exists(TE_POINTS):
    print("Downloading Greenland Te from Zenodo (doi:10.5281/zenodo.18403685)...")
    urllib.request.urlretrieve(TE_URL, TE_POINTS)

# ---- Step 2: load BedMachine, subsample to ~10 km, build load stress ----
print("Loading BedMachine Greenland v6...")
with xr.open_dataset(BEDMACHINE) as ds:
    x_native = ds.x.values.astype(float)
    stride = max(1, round(TARGET_DX / abs(x_native[1] - x_native[0])))
    sl = slice(None, None, stride)
    x = ds.x.values[sl].astype(float)
    y = ds.y.values[sl].astype(float)
    thickness = ds["thickness"].isel(y=sl, x=sl).values.astype(float)
    mask = ds["mask"].isel(y=sl, x=sl).values

dx = abs(x[1] - x[0])
dy = abs(y[1] - y[0])
ny, nx = len(y), len(x)
print(f"  Grid: {ny} x {nx} cells, dx = {dx:.0f} m, dy = {dy:.0f} m")

# Grounded ice only (mask == 2); floating ice is hydrostatically supported.
ice = np.where(mask == 2, thickness, 0.0)
load = RHO_ICE * G * ice   # Pa

# ---- Step 3: interpolate Steffen Te onto the same grid ----
print("Loading and interpolating elastic thickness...")
te_data = np.loadtxt(TE_POINTS)   # columns: lon, lat, Te (km), uncertainty (km)
lon_pts, lat_pts, te_km = te_data[:, 0], te_data[:, 1], te_data[:, 2]

lon_u = np.unique(lon_pts)
lat_u = np.unique(lat_pts)
te_ll = np.full((len(lat_u), len(lon_u)), np.nan)
j = np.round((lon_pts - lon_u[0]) / (lon_u[1] - lon_u[0])).astype(int)
i = np.round((lat_pts - lat_u[0]) / (lat_u[1] - lat_u[0])).astype(int)
te_ll[i, j] = te_km

# Fill gaps (ocean / outside coverage) with the domain median.
te_fill_km = float(np.nanmedian(te_ll))
te_ll = np.where(np.isnan(te_ll), te_fill_km, te_ll)

te_interp = RegularGridInterpolator(
    (lat_u, lon_u), te_ll, method="linear",
    bounds_error=False, fill_value=te_fill_km,
)

XX, YY = np.meshgrid(x, y)
lon_grid, lat_grid = Transformer.from_crs(
    "EPSG:3413", "EPSG:4326", always_xy=True
).transform(XX.ravel(), YY.ravel())
te = te_interp(np.column_stack([lat_grid, lon_grid])).reshape(ny, nx) * 1000.0  # m

print(f"  Te range: {te.min() / 1e3:.1f}-{te.max() / 1e3:.1f} km")

# ---- Step 4: pad the domain, then save ----
# Greenland ice reaches near all four map edges, so the flexural forebulge
# would be suppressed by a boundary placed right at the coast.  pad_domain()
# grows the grid by ~one flexural wavelength (zero load, tapered Te) so the
# clamped boundaries in greenland.yaml sit far from the load.  We bake the
# padding into the saved grids, keeping the run scripts as short as possible.
te_pad, load_pad, p = pad_domain(
    te, load, dx=dx, dy=dy, E=E, nu=NU, rho_m=RHO_MANTLE, rho_fill=RHO_FILL, g=G
)
print(f"  Padded {p} cells per side: {ny} x {nx} -> "
      f"{load_pad.shape[0]} x {load_pad.shape[1]} cells")

np.save(os.path.join(_dir, "greenland_load.npy"), load_pad)
np.save(os.path.join(_dir, "greenland_te.npy"), te_pad)
with open(os.path.join(_dir, "greenland_grid.txt"), "w") as fh:
    fh.write("# gFlex Greenland example grid, EPSG:3413 (padded)\n")
    fh.write(f"grid_spacing_x {dx:.1f}\n")
    fh.write(f"grid_spacing_y {dy:.1f}\n")
    fh.write(f"shape {load_pad.shape[0]} {load_pad.shape[1]}\n")
    fh.write(f"pad_cells_per_side {p}\n")

print(f"Saved greenland_load.npy and greenland_te.npy to {_dir}")
print(f"Use grid_spacing_x = {dx:.1f} m, grid_spacing_y = {dy:.1f} m in the YAML.")
