# gFlex examples

Worked examples of lithospheric flexure with gFlex, built around the modern
Greenland Ice Sheet.

## Start here: the Greenland teaching example

Lithospheric deflection under the Greenland Ice Sheet, with a spatially
variable effective elastic thickness. The two input grids are pre-built, so
you can run gFlex with nothing but gFlex and NumPy.

Two equivalent ways to run the same model:

**1. From a configuration file** — `greenland.yaml` defines everything
(parameters, input grids, boundary conditions, plot). Run it with either:

```bash
python -m gflex examples/greenland.yaml     # gFlex's command-line runner
python examples/run_greenland.py            # a three-line Python launcher
```

**2. From a Python script** — `greenland_script.py` builds the same model by
setting attributes on an `F2D` object directly. Run it with:

```bash
python examples/greenland_script.py
```

Both produce gFlex's standard load + deflection plot.

### The input grids

The load and elastic-thickness grids in `data/` are the only pieces that need
geospatial preprocessing, so that work is done once, up front:

| file | quantity | units |
| --- | --- | --- |
| `data/greenland_load.npy` | surface load stress, `rho_ice * g * H_ice` | Pa |
| `data/greenland_te.npy` | effective elastic thickness `T_e` | m |

`preprocess_greenland.py` builds them from BedMachine Greenland v6 (ice
thickness) and the Steffen et al. elastic-thickness model, downsampled to
~10 km on an EPSG:3413 grid and padded by ~one flexural wavelength so the
boundaries sit far from the load. You only need to rerun it to change the
resolution or source data; it requires the 3 GB BedMachine netCDF plus
`xarray`, `pyproj`, and `scipy`.

## Going further

- `greenland_flexure.py` — full research-grade version of the same problem:
  downloads/reads the source data itself, and compares variable vs. uniform
  `T_e`.
- `greenland_volcano_nested.py` — nested coarse/fine domains coupled through
  prescribed-value (inhomogeneous) boundary conditions, with a hypothetical
  seamount load.
