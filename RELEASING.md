# Release process

## Pre-release

1. **Confirm CSDMS BMI variable names** with Eric Hutton (currently pending for
   `mantle__mass-per-volume_density`, `infill_material__mass-per-volume_density`,
   `planet_surface__gravitational_acceleration`).  If names change, update
   `gflex/bmi.py`, `tests/test_bmi.py`, `docs/api.rst`, and `CHANGES.md`.

2. **Verify CI is green** on `master` before starting.

3. **Bump version** in `gflex/_version.py`:
   ```
   __version__ = "2.0.0"
   ```
   (`docs/conf.py` and `pyproject.toml` both pull from this file dynamically.)

4. **Check `CHANGES.md`** — the top entry must be `## [X.Y.Z] — YYYY-MM-DD`
   with no `[Unreleased]` header above it.

5. **Run the test suite** — all tests must pass:
   ```
   python -m pytest
   ```

6. **Build the docs** — zero warnings required:
   ```
   python -m sphinx docs docs/_build/html -q
   ```

7. **Commit the version bump** (one commit, e.g. `release: bump version to 2.0.0`).

## Release

8. **Push** to GitHub (`origin/master`).

9. **Wait for CI** (`test.yml` and `build-check.yml`) to go green.

10. **Create a GitHub Release** with tag `v2.0.0` (target: `master`).
    - Title: `gFlex 2.0.0`
    - Body: paste the `[2.0.0]` section from `CHANGES.md`
    - Publishing the release triggers `publish.yml` automatically.

11. **Approve the PyPI deployment** — the `pypi` environment requires manual
    approval in the GitHub Actions UI before the package is uploaded.

12. **Verify** on PyPI: `pip install gflex==2.0.0`

## Post-release

13. **Add `## [Unreleased]`** at the top of `CHANGES.md` for future patch releases.

14. **Zenodo** — a new version record is created automatically via the GitHub
    integration.  Once it appears, update the version-specific DOI in
    `docs/references.rst` (currently uses the concept DOI as placeholder).

15. **GRASS GIS** — proceed with the `r.flexure` / `v.flexure` update
    (issue #55; branches `update/r.flexure` and `update/v.flexure` in
    `~/models/grass-addons-gflex`).

16. **Landlab** — PR `landlab/landlab#2420` is unblocked once `gflex>=2.0.0`
    is on PyPI.  Close issue #73 once that PR merges.
