# Release Notes

## [Unreleased]

### Breaking changes

- **BMI load variable renamed** — `BmiGflex` input variable
  ``lithosphere__load_pressure`` has been renamed to
  ``load__normal_component_of_stress``.  The old name implied the load was a
  property of the lithosphere, but it is applied *to* the lithosphere from
  above and is material-agnostic (ice, sediment, water, topography, etc.
  all express as a surface-normal stress in Pa before being passed to gFlex).
  The new name follows CSDMS quantity vocabulary and matches the field name
  used by the Landlab gFlex component.  Update all ``set_value`` /
  ``get_value`` calls accordingly:
  ``set_value("load__normal_component_of_stress", q_array)``.

- **Matrix-determining array inputs are read-only** — when `T_e` or the
  in-plane stresses (`sigma_xx`, `sigma_yy`, `sigma_xy`) are set as NumPy
  arrays, they are now stored as read-only copies.  An in-place edit
  (``flex.T_e[i] = ...``) raises `ValueError`; change a value by reassigning
  the whole array (``flex.T_e = new_array``).  This closes a silent-error
  path: the cached coefficient matrix / LU factorisation is keyed to these
  inputs, and an in-place edit previously bypassed cache invalidation, so the
  next ``run()`` returned deflections computed from the *old* values.
  Reassignment invalidates the cache correctly.

### New features

- **BMI `lithosphere__elastic_thickness` input variable** — `BmiGflex` now
  exposes elastic thickness as a second input variable (CSDMS Standard Name
  ``lithosphere__elastic_thickness`` [m]).  Te is usually constant but can be
  updated between coupling steps via ``set_value()``; the new value is pushed
  to the solver immediately, invalidating the cached coefficient matrix so the
  next ``update()`` recomputes deflection with the revised rigidity field.

- **BMI scalar physical constants** — `BmiGflex` now exposes the five
  physical constants as single-element (grid 1) input variables, supporting
  introspection and ensemble initialisation without separate config files:

  | CSDMS Standard Name | gFlex attr | Units |
  |---|---|---|
  | `lithosphere__young_modulus` | `E` | Pa |
  | `lithosphere__poisson_ratio` | `nu` | 1 |
  | `mantle__mass-per-volume_density` | `rho_m` | kg m⁻³ |
  | `infill_material__mass-per-volume_density` | `rho_fill` | kg m⁻³ |
  | `planet_surface__gravitational_acceleration` | `g` | m s⁻² |

  Values are read from the config file at ``initialize()`` and exposed on
  grid 1 (shape ``(1,)``).  Updating any constant via ``set_value()``
  propagates to the solver immediately and invalidates the cached LU
  factorisation; the next ``update()`` rebuilds the stiffness matrix
  automatically.  `rho_fill` is the constant most likely to change at
  runtime (e.g. when a subaerial basin becomes subaqueous).

- **FD one-sided periodic warning** — `F1D` and `F2D` now issue a
  `UserWarning` when `method='fd'` and exactly one side of an opposite
  boundary pair (west/east or north/south) is `'periodic'`.  Periodic BCs
  tile opposite edges together, so a one-sided periodic is non-physical in
  most cases; the warning names both sides and suggests the fix.  The solver
  proceeds rather than raising an error, in case an unforeseen use case
  requires it.
- **FFT per-axis boundary conditions** — `F2D` now handles the west/east
  and north/south boundary-condition pairs independently for
  `method='fft'`.  Setting both sides of a pair to `'periodic'` makes
  that axis genuinely periodic; leaving them unset (or `'no_outside_loads'`)
  zero-pads that axis.  Mixing periodic and non-periodic axes is supported
  (e.g. `bc_west = bc_east = 'periodic'` with north/south unset gives an
  x-periodic, y-padded domain).  Setting only one side of a pair to
  `'periodic'` raises a `UserWarning` and falls back to zero-padding for
  that axis.  `F1D` behaviour is unchanged (single axis, both sides must
  match).
- **`fft_pad_n_alpha`** — new instance attribute on `F1D` and `F2D`
  (default ``4``) that controls the FFT zero-padding width: the load is
  zero-padded by `fft_pad_n_alpha × α` cells on each side, placing
  periodic images at least `2 × fft_pad_n_alpha × α` apart.  `F1D` uses
  α₁D = (4D/Δρg)^0.25; `F2D` uses α₂D = (D/Δρg)^0.25.
- **`'infinite'` alias** — added ``'infinite'`` as a short alias for the
  ``'no_outside_loads'`` boundary condition, consistent with the existing
  aliases (``'clamped'``, ``'pinned'``, ``'free'``, ``'mirror'``).

- **`'pinned'` BC alias** — added `'pinned'` as a short alias for
  `'zero_displacement_zero_moment'`, consistent with the existing
  `'clamped'`, `'free'`, and `'mirror'` aliases.
  (The `zero_displacement_zero_moment` BC itself was added in 2.0.0b1.)

- **FD `'no_outside_loads'` boundary condition** — `F1D` and `F2D` now
  accept `'no_outside_loads'` as a valid `bc_*` string when
  `method='fd'`.  On each side where it is set, the solver automatically
  pads the domain by one flexural wavelength with zero loads, applies a
  clamped (`zero_displacement_zero_slope`) outer boundary, solves, and
  crops `w` back to the original domain.  The padding is invisible to the
  caller: `w.shape == qs.shape` on return.  Sides with explicit BCs are
  unaffected, so mixing `'no_outside_loads'` with `'mirror'`,
  `'zero_displacement_zero_slope'`, etc. on opposite sides is supported.
  Variable-`Te` arrays are tapered smoothly into the padded region.

- **Domain-padding API unified** — `pad_domain(Te, qs, dx, ...)` now
  dispatches to 1-D or 2-D based on the shape of `qs`.  The former
  `pad_domain_1d` public function is now the private `_pad_domain_1d`
  (called internally by `pad_domain`).  `smooth_pad_Te`,
  `recommended_pad_width`, `smooth_pad_Te_1d`, and
  `recommended_pad_width_1d` remain public for callers that need finer
  control.

### Verbosity and logging

- **Default verbosity changed to quiet** — `F1D` and `F2D` now default to
  ``quiet=True, verbose=False, debug=False``.  In library and coupling
  contexts (CSDMS BMI, Landlab, GRASS GIS) the calling framework controls
  logging; gFlex should be silent unless asked.  To restore the previous
  behaviour set ``flex.quiet = False; flex.verbose = True`` after
  construction.
- **Verbose output restructured** — the ``quiet=False`` / ``verbose=True``
  output now follows the BMI lifecycle: an ``initialized`` block shows
  grid dimensions, T_e, and method; a ``run`` block shows BCs and
  per-phase timing (6 decimal places); a ``finalized`` block shows total
  runtime.  gFlex output is bracketed by blank lines to separate it from
  surrounding caller messages.
- **Timing attributes** — after each ``run()``, three instance attributes
  report wall-clock times measured with ``time.perf_counter()``:
  ``time_to_solve`` (total solve time; all solvers),
  ``coeff_creation_time`` (seconds to construct the sparse coefficient
  matrix; FD only), and ``linear_solve_time`` (seconds for the backsolve;
  FD only).  Non-FD solvers (SAS, FFT) do not set the latter two.  The
  FD split is useful in coupling loops with ``cache_factorization``
  enabled, where ``coeff_creation_time ≈ 0`` confirms the matrix was
  reused.
- **Python ``logging`` module** — all ``print()`` calls in the solver path
  have been replaced with ``logging`` calls on the ``gflex`` logger
  (``_logger.info``, ``_logger.warning``, ``_logger.debug``).  A
  ``NullHandler`` is registered on the ``gflex`` logger so the library
  emits no output by default.  Standalone scripts can enable console output
  by setting ``quiet=False`` (``logging.INFO``) or ``debug=True``
  (``logging.DEBUG``); coupling frameworks configure their own handlers on
  the ``"gflex"`` logger as needed.

### Input validation

- `initialize()`, `flexural_wavelengths()`, `recommended_pad_width_1d()`, and
  `recommended_pad_width()` now raise a clear `ValueError` when
  ``rho_fill >= rho_m``.  A non-positive density contrast means the foundation
  has no restoring force and the flexural parameter α is undefined; previously
  the code produced a ``ZeroDivisionError`` or silent ``nan``.

### Bug fixes

- Fixed ``bc_check()`` silently skipping BC validation when ``quiet=True``:
  the guard ``if not self.quiet: raise ValueError(...)`` meant a
  misconfigured analytical solver (e.g. SAS with FD boundary conditions set)
  would proceed without error in quiet mode.  Analytical solvers always use
  the infinite-plate assumption and ignore FD BCs regardless; the fix issues
  a ``UserWarning`` unconditionally so the user is informed without the solve
  being blocked.
- Fixed FFT crash when no boundary conditions are set: `bc_check()` was
  populating legacy uppercase attributes (`BC_E`, `BC_W`, `BC_N`, `BC_S`)
  instead of the current lowercase properties (`bc_east`, `bc_west`, etc.).
  `_solve_fft()` reads the lowercase properties; with BCs unset the
  property getter raised `AttributeError`.  Callers that did not explicitly
  assign a BC (e.g. GRASS GIS `r.flexure`) were affected.  Fix: default
  unset BCs to `""` via the property setter, mirroring the SAS path.
- Fixed 2-D FFT zero-padding width: was using the 1-D flexural parameter
  α₁D = (4D/Δρg)^0.25 instead of the 2-D parameter α₂D = (D/Δρg)^0.25
  (α₁D ≈ √2 × α₂D), causing ~41 % over-padding on every non-periodic 2-D
  FFT run.  Both `F1D` and `F2D` now derive α from `flexural_wavelengths()`
  rather than hardcoding the formula.
- Fixed stale `self.nx` (F1D) and `self.nx`/`self.ny` (F2D) after a
  `'no_outside_loads'` FD run.  The auto-padding path calls
  `gridded_x()` / `build_diagonals()` on the enlarged domain, which set
  these attributes to the padded size.  They were not restored when `w`
  and `qs` were cropped back to the original domain, leaving them
  inconsistent with `w.shape` until the next run.  Now explicitly
  restored after the crop.
- Fixed `RuntimeError` when a coupling loop re-runs the FD solver after a
  property setter (e.g. assigning `T_e`) invalidates the cached coefficient
  matrix between runs.  `bc_check()` seeds the internal `_bc_*_norm`
  attributes with the raw user-supplied BC string, then resolves short
  aliases (``'mirror'`` → ``'zero_slope_zero_shear'``, etc.) only inside
  the ``coeff_matrix is None`` guard.  If the matrix was valid at
  `bc_check()` time but later invalidated by a setter, the guard was
  skipped and `_apply_bc_rigidity()` received an unrecognised alias string.
  Fix: alias resolution now runs unconditionally before the guard.
- Fixed variable-Te 2-D MMS manufactured load: the cross-derivative term
  used the simplified biharmonic coefficient 2 rather than the van Wees &
  Cloetingh (1994) value 2(1−ν).  For bilinear D the missing
  2ν·∂²D/∂x∂y·∂²w/∂x∂y term capped apparent convergence at
  O(Δx^1.91) near n = 200–400 instead of the correct O(Δx^2.00).  The
  fix is in test code only (`_mms_2d_variable_te` in
  `tests/test_bc_mms.py`); the solver itself is unaffected.
- Fixed the sign of the right-hand-side correction for a prescribed nonzero
  slope on the *east* boundary in `F1D` dict-style BCs
  (`{'displacement': w, 'slope': theta}`).  The ghost-node constant was
  moved to the RHS without the required sign flip, producing a large error
  localized at the east boundary whenever a nonzero east slope was
  prescribed (the west boundary and the 2-D code were already correct).  A
  quadratic manufactured solution now reproduces the deflection to solver
  precision (`tests/test_1D_inhomogeneous_BC.py`, Case E).
- Fixed an `AttributeError` (`'float' object has no attribute 'shape'`) when
  running with `debug=True` and a scalar elastic thickness: the solver debug
  path logged `self.T_e.shape`, and logging evaluates its arguments
  regardless of level.  `F1D` and `F2D` now log `np.shape(self.T_e)`, which
  returns `()` for a scalar; the `F2D` `contextlib.suppress` band-aid was
  removed so the shape is actually logged.
- Fixed the 2-D inhomogeneous moment/shear boundary correction wrapping
  elastic thickness across the domain corners.  The cross-derivative ghost
  terms in `_apply_bc_rhs_inhomogeneous_2d` used `np.roll`, whose wrap-around
  made each corner equation read the diagonally-opposite corner's Δ₁.  With
  uniform Te this is harmless (Δ₁ constant), but with variable Te it injected
  the far corner's rigidity into the near corner, degrading convergence and
  making the deflection at one corner depend non-physically on the opposite
  corner's elastic thickness.  The edge endpoints are now clamped rather than
  wrapped (`_shift_clamp`).
- Fixed `F2D` leaving `self.T_e` as the internally padded ``(ny+2, nx+2)``
  array after a run with a 2-D Te grid.  The solvers now keep their padded,
  mutable working thickness in a private copy (`_te`) and never modify the
  user-facing `T_e`, so it retains the supplied ``(ny, nx)`` shape across
  solves.
- Fixed a load applied **on** a fixed-displacement boundary node not being
  reacted by the support.  For clamped (`zero_displacement_zero_slope`),
  pinned (`zero_displacement_zero_moment`), and dict ``{"displacement": ...}``
  boundaries, the boundary row is decoupled to enforce the prescribed value,
  but a load on that node leaked into its right-hand side, giving
  ``w = w0 − qs/c0`` instead of ``w0`` and corrupting the nearby interior.
  Such loads (common when a whole-grid load array is nonzero on a clamped
  edge, e.g. in a coupling loop) are now reacted by the support so the node
  holds its prescribed displacement exactly, in `F1D` and `F2D`.  A
  `UserWarning` is issued once when a load sits on a fixed-displacement
  boundary.  Free and mirror edges, which do not fix displacement, are
  unaffected.
- Fixed `F1D.fd_solve` raising `AttributeError` when a caller supplies a
  prebuilt `coeff_matrix` (coupled mode), which skips matrix assembly and the
  creation of `_bc_rhs_correction`.  The 1-D path now guards for the missing
  correction, matching `F2D`.
- Fixed the `'sandbox'` BC easter egg raising `AttributeError`: it called
  `finalize()` (which deletes `w`) before reading `solver.w`.  It now reads
  the deflection before finalising, so it renders as intended.

### Documentation

- New `theory_and_numerics` page: governing PDE, flexural rigidity, the
  full variable-D van Wees & Cloetingh (1994) FD expansion, and
  descriptions of all solution methods (SAS, FD, FFT).  Replaces the
  orphaned `theory.rst` and `numerical_methods.rst` pages from the beta.
- `boundary_conditions` page expanded: `versionadded` markers for
  `zero_displacement_zero_moment` and `no_outside_loads`; proximity
  warning documentation; prescribed-BC usage examples.
- `api.rst` expanded: BMI variables split into grid-0 (spatial) and
  grid-1 (scalar constants) tables; LU cache and timing attributes
  documented; coupling guide added covering load conversion, `rho_fill`
  semantics, and topography-update patterns; Landlab component section
  added with field-name table, BMI comparison, and coupling example.
- `index.rst`: interfaces table reorganised into three groups
  (Python/CLI, Modelling frameworks, GIS); Status column removed;
  GRASS repo links added.
- Accuracy page: variable-Te 2-D MMS section updated with corrected
  governing equation (2(1−ν) cross-derivative coefficient) and revised
  convergence results (O(Δx^2.00)).

### Tests

- 509 tests passing (up from 335 in 2.0.0b1), covering all solvers,
  all BC types (including pinned and prescribed-value), MMS convergence
  for mirror (1-D and 2-D) and variable-Te (2-D), domain padding, LU
  cache, warnings, and the full BMI variable set including the five
  scalar physical constants.

----

## [2.0.0b1] - 2026-06-04

### Breaking changes

- **`finalize()` now clears all model state** — ``w``, ``qs``, and the
  coefficient matrix are deleted.  Read ``w`` and call ``output()``
  **before** ``finalize()``.  The CLI (``gflex.py``) has been fixed to
  follow this order.

- **INI configuration file format removed** — only YAML (`.yaml` / `.yml`)
  configuration files are supported.  Passing a `.ini` file now raises
  `ValueError`.  Convert any remaining INI configs to YAML format.

- **Iterative FD solver removed** — `F1D` now always uses a direct sparse LU
  factorization.  The `ConvergenceTolerance` configuration key is ignored.
- **G2009 plate solution removed** — the finite-difference solver now always
  uses the vWC1994 stencil.  Setting `PlateSolutionType = 'G2009'` previously
  selected a different stencil; it is now a no-op attribute (ignored silently).
- **`PlateSolutionType` removed from the public interface** — the attribute is
  no longer documented or read by any solver.  Existing scripts that set it
  will not raise an error, but the value has no effect.
- **`te` attribute renamed to `T_e`** — the elastic-thickness attribute on
  `F1D` and `F2D` is now `T_e` (PEP 8 subscript notation; avoids ambiguity
  with a time variable `t` and a generic elastic parameter `e`).  `flex.te`
  now raises `AttributeError`.  Update all occurrences: `flex.te = …` → `flex.T_e = …`.

- **`"mirror"` is now an alias for `"zero_slope_zero_shear"`** — the canonical
  BC name is `"zero_slope_zero_shear"` (consistent with all other long-form names);
  `"mirror"` remains a permanent short alias (like `"clamped"` and `"free"`) and
  does not trigger any warning.

- **`solver` attribute now raises `ValueError` for unsupported values** —
  only `"direct"` is supported in this release.  Passing any other string
  (e.g. `"iterative"`) previously silently fell back to the direct solver;
  it now raises `ValueError` immediately.  An iterative solver may be added
  in a future version, at which point additional values will become valid.

- **BC string case normalised; v1.x PascalCase strings removed** — all
  boundary-condition strings are now lowercase (`"zero_displacement_zero_slope"`,
  `"zero_moment_zero_shear"`, `"zero_slope_zero_shear"`, `"periodic"`).  The old v1.x
  PascalCase names (`"0Displacement0Slope"`, `"0Moment0Shear"`, `"Mirror"`,
  etc.) now raise `ValueError`.  See the `boundary_conditions` page for the
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

- **`gflex.VALID_BC_STRINGS_1D` and `gflex.VALID_BC_STRINGS_2D`** — public
  `frozenset` constants listing every accepted BC string (canonical names and
  aliases).  Wrappers can validate user input against these sets instead of
  maintaining a parallel copy that drifts with new releases.

- **FD `'zero_displacement_zero_moment'` (simply-supported / pinned) boundary
  condition** — `F1D` and `F2D` now accept `'zero_displacement_zero_moment'`
  as a valid `bc_*` string when `method='fd'`.  The plate is held at zero
  deflection but is free to rotate; no bending moment is transmitted across
  the boundary.  Corresponds to a simply-supported (pinned) end in structural
  mechanics.  Implemented via odd-reflection ghost nodes — contrast with the
  even-reflection `mirror` condition.  See the `boundary_conditions` page for
  diagrams.

- **`"clamped"`, `"free"`, and `"mirror"` BC aliases** — concise alternatives
  to the full canonical names: `"clamped"` normalises to
  `"zero_displacement_zero_slope"`; `"free"` normalises to
  `"zero_moment_zero_shear"`; `"mirror"` normalises to
  `"zero_slope_zero_shear"`.  All three produce bit-identical results to their
  canonical names.  (The `"pinned"` alias was added in the subsequent
  post-beta commit; see [Unreleased].)

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

- 335 tests passing across 1-D and 2-D FD, FFT, SAS/SAS_NG solvers, all BC
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
