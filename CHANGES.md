# Release Notes

## [Unreleased]

### New features

- **FFT spectral solver** for 1-D and 2-D problems.  Exact for periodic
  boundaries; other boundary conditions are handled with 4α zero-padding
  (NoOutsideLoads approximation).  Requires scalar (uniform) elastic
  thickness.  In-plane stresses (`sigma_xx`, `sigma_yy`, `sigma_xy`) are
  fully supported.
- **YAML configuration file format** alongside the legacy INI format.  Pass
  a `.yaml` / `.yml` file to the CLI or to the `F1D` / `F2D` constructor.
- **Domain-padding utilities** to reduce spurious boundary effects when
  using spatially variable elastic thickness:
  - 2-D: `pad_domain`, `smooth_pad_Te`, `recommended_pad_width`
  - 1-D: `pad_domain_1d`, `smooth_pad_Te_1d`, `recommended_pad_width_1d`
- **FD boundary-condition warnings** — `F1D` and `F2D` now issue
  `UserWarning` messages for `'0Moment0Shear'` (free broken end — verify
  a rifted margin is intended), `'0Slope0Shear'` (no clear geological
  analog), and when the nearest loaded cell is within one flexural
  wavelength of a `'0Displacement0Slope'` boundary (forebulge suppression).

### Improvements

- Iterative FD solver (`Solver = iterative`) upgraded from Jacobi to ILU
  preconditioning; falls back to a direct solve if LGMRES does not converge.
  The `ConvergenceTolerance` parameter is now passed as `rtol` (relative
  residual) to SciPy's LGMRES; default remains `1e-3`.
- In-plane membrane stresses (`sigma_xx`, `sigma_yy`, `sigma_xy`) now
  supported by all FD and FFT solvers in both 1-D and 2-D.

### Bug fixes

- Fixed 2-D FD constant-Te stencil: corrected dx/dy swap, σ swap, and
  missing `/dx²` factors that produced incorrect results for asymmetric grids.
- Fixed 2-D FFT `sigma_xy` assembly: double-wrapped corner entries and wrong
  coefficient in the Periodic right-roll buffer.
- Fixed `0Slope0Shear` description in configuration reference and README.

### Documentation

- New Sphinx / Read the Docs site with a theory page (governing equations,
  solution methods, in-plane stresses), an accuracy page, a configuration
  reference, API reference, and changelog.
- Theory page covers the full governing PDE with in-plane stress terms for
  both 1-D and 2-D, the variable-D FD expansion, and the FFT transfer
  function.

### Tests

- Comprehensive new test suites for 1-D SAS/SAS_NG, 2-D SAS/SAS_NG,
  1-D FFT, 2-D FFT (including sigma_xy), 1-D FD boundary conditions with
  analytical cross-validation, 2-D FD, FD boundary-condition warnings, and
  all domain-padding utilities.

## [1.3.0](https://github.com/awickert/gFlex/releases/tag/v1.3.0) - 2026-05-27

- New: `BmiGflex` — CSDMS Basic Model Interface (BMI) v2 implementation;
  `bmipy` is an optional dependency and gFlex works without it.
- New: `F1D`, `F2D`, and `BmiGflex` are now accessible directly from the
  `gflex` package namespace (e.g. `import gflex; gflex.F2D()`); previously
  this raised `AttributeError`.
- Fix: updated `scipy.sparse.linalg` import for compatibility with modern
  scipy.
- Fix: latent typo `self.T_e` → `self.Te` in the scalar-Te branch of
  `get_coeff_values` (introduced 2021-06-26); would produce incorrect
  results when tectonic stress terms (sigma_xx, sigma_yy, sigma_xy) are
  non-zero.
- Removed outdated root-level `gflex_bmi.py`.

## [1.2.0](https://github.com/awickert/gFlex/releases/tag/v1.2.0) - 2024-01-08

- Python 3 modernisation: black formatting, isort, updated imports,
  argparse CLI, `pyproject.toml` replacing `setup.py`.
- GitHub Actions CI replacing Travis CI.
- Eric W. H. Hutton added as author.

## [1.1.1](https://github.com/awickert/gFlex/releases/tag/v1.1.1) - 2021-06-26

- Updated PyPI support: twine upload, README.md
- Throw meaningful error if a nonuniform self.Te grid is used with the
  analytical solution

## [1.1.0](https://github.com/awickert/gFlex/releases/tag/v1.1.0) - 2018-05-28

- Support for both Python 2 and Python 3
- Main code
- Examples
- README.md updates
- PATH updates in
- Code tests included
- Code testing on commit by Travis
- Updated on PyPI

## [1.0.1](https://github.com/awickert/gFlex/releases/tag/v1.0.1) - 2018-05-25

- Final Python 2 (only) release
- Additional documentation added

## [1.0.0](https://github.com/awickert/gFlex/releases/tag/v1.0.0) - 2017-05-10

This is the update is what is available now from pypi.

- Minor updates to the main gFlex codes
- Addition of a missing "12" in the flexural wavelength calculator in the utilities.


## [1.0](https://github.com/awickert/gFlex/releases/tag/v1.0) - 2016-01-28

- First full release of gFlex in association with the now-accepted GMD paper,
"Open-source modular solutions for flexural isostasy: gFlex v1.0", by A. D. Wickert.


## [0.9](https://github.com/awickert/gFlex/releases/tag/v0.9) - 2015-03-05

- Release submitted to GMD and that appears on PyPI; this is the same as v1.0a
  on umn-earth-surface.

## [0.8.1](https://github.com/awickert/gFlex/releases/tag/v0.8.1) - 2015-03-04

- Fixed error in PyPI integration from v0.8 and updated README.md to include
  PyPI integration.


## [0.8](https://github.com/awickert/gFlex/releases/tag/v0.8) - 2015-03-04

- First fully-functional and error-checked release, and therefore the first
  true non-"pre-release".


## [0.7](https://github.com/awickert/gFlex/releases/tag/0.7) - 2015-01-21

This release may be short-lived, but is the turning point at which work
on gFlex will go from fixing deficiencies to adding new capabilities.
As such, it is receiving version 0.7 (the last version of its parent project,
"Flexure", was 0.6 back in 2012). It now seems to be functional and stable,
and thus this tagged release will be the basis for this stage of the model
to be pulled to @csdms-contrib and @umn-earth-surface.
