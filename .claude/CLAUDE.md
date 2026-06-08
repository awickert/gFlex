# gFlex project conventions for Claude

## Directory layout

- `analysis/` — accuracy analysis, MMS convergence studies, old-vs-new BC
  comparators. This is where scripts like `analyze_clamped_bc_error.py`,
  `analyze_free_end_bc_error.py`, and `compare_bc_versions.py` live.
- `benchmarks/` — speed benchmarks only (`bench_solvers.py`). Do NOT put
  accuracy or error-analysis scripts here.
- `docs/_static/` — figures committed for use in Sphinx documentation.
- `analysis/results/` — local output from analysis scripts (gitignored).
- `benchmarks/results/` — local output from speed benchmarks (gitignored).

## Public API is frozen

`F1D`, `F2D`, and the IRF interface are used by many downstream projects
(Landlab, GRASS GIS add-ons, user scripts). Never change attribute names,
method signatures, or BC string values in a way that breaks existing code.
Deprecation warnings are acceptable; silent breakage is not.

## Commit discipline

- Never commit without explicit per-commit approval from Andy. Show the diff
  summary and proposed message and wait for a yes.
- Commits must be granular: one logical change per commit.
- Commit messages describe *what changed and why*, not just what.
- When a message lists multiple renames, use one line per change, not commas.
- Never amend a published commit; always create a new one.

## Releases and pushing

- Never push to a shared remote, create a git tag, or submit to any external
  registry (PyPI, Zenodo, etc.) without explicit per-action authorization in
  the current message. "These are the last steps" is not authorization.
- Andy has explicitly forbidden pushing in this project.

## Documentation style

- Describe what *is*, not what might be. Never forward-reference unreleased
  versions or planned features.
- Preserve Andy's voice and phrasing. Add and amend; do not rewrite.
- `accuracy.rst` version attributions must match the actual release in which
  a fix landed — double-check before writing version numbers.

## GitHub issues

Never close a GitHub issue without explicit per-issue instruction from Andy.

## BC string values (case)

Canonical BC strings (all lowercase): `"zero_displacement_zero_slope"`,
`"zero_displacement_zero_moment"`, `"zero_moment_zero_shear"`,
`"zero_slope_zero_shear"`, `"periodic"`, `"no_outside_loads"`. Short aliases
(permanent, no warning): `"clamped"` → `zero_displacement_zero_slope`,
`"free"` → `zero_moment_zero_shear`, `"mirror"` → `zero_slope_zero_shear`,
`"infinite"` → `no_outside_loads`. Uppercase v1.x variants trigger a
DeprecationWarning and are normalized internally.

## Analysis scripts

The `analysis/` scripts run from the repo root. Results (figures, text) go
to `analysis/results/` which is gitignored. Figures used in documentation
are copied to `docs/_static/` and committed there separately.
