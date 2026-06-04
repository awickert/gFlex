# Release Notes

## [Unreleased]

### Breaking changes

- **`finalize()` now clears all model state** — ``w``, ``qs``, and the
  coefficient matrix are deleted.  Read ``w`` and call ``output()``
  **before** ``finalize()``.  The CLI (``gflex.py``) has been fixed to
  follow this order.

- **INI configuration file format removed** — only YAML (`.yaml` / `.yml`)
  configuration files are supported.  Passing a `.ini` file now raises
  `ValueError`.  Convert any remaining INI configs to YAML format.

- **Iterative FD solver removed** — `F1D` now always uses a direct sparse LU
  factorization.  Setting `Solver = 'iterative'` has no effect (a message is
  printed and the direct solver is used instead).  The `ConvergenceTolerance`
  configuration key is ignored.
- **G2009 plate solution removed** — the finite-difference solver now always
  uses the vWC1994 stencil.  Setting `PlateSolutionType = 'G2009'` previously
  selected a different stencil; it is now a no-op attribute (ignored silently).
- **`PlateSolutionType` removed from the public interface** — the attribute is
  no longer documented or read by any solver.  Existing scripts that set it
  will not raise an error, but the value has no effect.
- **BC string case normalised; v1.x PascalCase strings deprecated** — all
  boundary-condition strings are now lowercase (`"zero_displacement_zero_slope"`,
  `"zero_moment_zero_shear"`, `"mirror"`, `"periodic"`).  The old v1.x
  PascalCase names (`"0Displacement0Slope"`, `"0Moment0Shear"`, `"Mirror"`,
  etc.) are still accepted but trigger a `DeprecationWarning` and will be
  removed in a future release.  See the `boundary_conditions` page for the
  full mapping.

### New features

- **Inhomogeneous (prescribed-value) boundary conditions** in 1-D and 2-D FD
  solver.  Any edge can carry a combination of prescribed displacement, slope,
  moment, and/or shear by passing a dict instead of a string:

  ```python
  flex.bc_west = {"displacement": w_arr, "slope": dwdx_arr}
  flex.bc_east = {"moment": M0, "shear": V0}
  ```

  The primary use case is **nested-domain (coarse-to-fine) modelling**: extract
  `w` and `np.gradient`-based slopes at the boundary of a fine sub-domain from
  a coarse regional solve, impose them as BCs, and the fine solver inherits the
  full far-field loading through those four boundary arrays.  Also supports
  broken-plate problems with applied edge loads.  Both 1-D and 2-D; verified
  by a round-trip unit test (`TestNestedModelGradientRoundTrip`).

- **`"clamped"` and `"free"` BC aliases** — concise alternatives to the full
  canonical names: `"clamped"` normalises to `"zero_displacement_zero_slope"`;
  `"free"` normalises to `"zero_moment_zero_shear"`.  Both produce
  bit-identical results to their canonical names.

- **`"zero_slope_zero_shear"` normalises to `"mirror"`** — the two BCs are
  mathematically identical (same even-reflection ghost equations); the former
  name is accepted silently and collapses to `"mirror"` internally.

### Performance

- **LU factorization cache** — set ``flex.cache_factorization = True`` to
  cache the sparse-LU factorization of the FD coefficient matrix across
  ``run()`` calls.  The cache is reused automatically when the matrix is
  unchanged (Te, grid, BCs, and physical parameters all fixed) and
  invalidated by a hash comparison when any of those inputs change.  Set
  ``flex.cache_factorization = "no_check"`` to skip the hash and reuse the
  factorization unconditionally — maximum performance for coupling workflows
  where the user guarantees matrix stability.  Default ``False`` preserves
  existing behaviour.
- **Multithreaded FFT** — ``scipy.fft`` transforms now use all available CPU
  threads (``workers=-1``).  Benefit is marginal at typical grid sizes;
  meaningful at very large grids (≳ 2000 × 2000 cells).

### Improvements

- `twoSurfplots()` now shows three panels (load, elastic thickness, deflection)
  when `Te` is a 2-D array, using an asymmetric diverging colormap for the
  deflection.
- `export_for_blender()` writes load-cylinder geometry (radius from loaded
  area, height proportional to load magnitude) and per-vertex colour weights
  to the mesh file, enabling the companion Blender script to render a
  physically scaled load cylinder.
- Blender scene script: cylinder base deforms with the plate surface; camera
  frames the full scene vertically.

### Bug fixes

- Fixed `SyntaxWarning` for invalid escape sequences (`\sigma`) in the F2D
  class docstring (Python 3.12+).
- Fixed a longstanding ghost-node error in the 1-D
  `zero_displacement_zero_slope` FD boundary condition.  The original
  implementation (present since the first version) dropped ghost-node stencil
  terms rather than eliminating them via even reflection, and left the
  boundary row coupled to interior nodes rather than decoupling it as a
  Dirichlet constraint.  The practical effect was that the boundary deflection
  was not strictly zero and the zero-slope condition was not enforced; results
  were accurate only when the boundary was far enough from any load that
  deflection was naturally negligible there (the recommended use case, and the
  intent of the existing proximity warning).  The corrected implementation
  decouples the boundary row (enforcing w = 0 exactly) and folds the
  even-reflected ghost into the adjacent interior node (enforcing dw/dx = 0).
- Fixed the same two ghost-node bugs in the **2-D `zero_displacement_zero_slope`**
  boundary condition (boundary rows/columns not decoupled as Dirichlet
  constraints; even-reflection ghost absent at the first interior row/column on
  all four edges).  Convergence recovers from first order (O(Δx^0.92)) to
  second order (O(Δx^1.99)).
- Fixed a ghost-node inconsistency in the **`zero_moment_zero_shear`** FD
  boundary condition (1-D and 2-D).  The first interior node's stencil used a
  shear ghost evaluated at a staggered location one cell inward, inconsistent
  with the moment ghost used at the boundary node itself.  The corrected
  implementation uses the moment condition consistently at both rows, recovering
  second-order convergence (O(Δx^2.00) in 1-D, O(Δx^2.01) in 2-D) from first
  order.

### Documentation

- New `boundary_conditions` page: comprehensive reference covering all BC
  types with a summary table (structural-mechanics and geophysical names), SVG
  diagrams for each condition, physical guidance on when each is appropriate,
  and a deprecation table mapping v1.x PascalCase strings to their v2.0
  lowercase equivalents.
- New `greenland_example` page: a realistic Greenland ice-sheet flexure model
  (BedMachine v6 ice thickness; Steffen et al. spatially variable elastic
  thickness) and a hypothetical nested-domain seamount scenario that illustrates
  how inhomogeneous BCs couple a coarse regional model to a fine local one.
- Accuracy page extended with MMS convergence tables and figures for the 2-D
  `zero_displacement_zero_slope` ghost-node fix and the 1-D/2-D
  `zero_moment_zero_shear` ghost-node fix.

### Tests

- 325+ tests passing across 1-D and 2-D FD, FFT, SAS/SAS_NG solvers, all BC
  types, inhomogeneous BCs, domain padding, warnings, and BC aliases.
- `TestNestedModelGradientRoundTrip`: verifies that `np.gradient`-extracted
  slopes used as inhomogeneous Dirichlet BCs reproduce the full-domain interior
  to within 2 % of peak deflection, confirming the slope sign convention for
  all four edges in the nested-domain workflow.
- Parametrised alias tests confirm `"clamped"` and `"free"` produce
  bit-identical results to their canonical names in both 1-D and 2-D.

## [1.4.0](https://github.com/awickert/gFlex/releases/tag/v1.4.0) - 2026-05-29

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
