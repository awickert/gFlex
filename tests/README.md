# gFlex test suite

Tests live in this directory and are run with **pytest**:

```
pytest tests/
```

The test suite is split into four files by solver type.  All tests that
require Landlab are automatically skipped when Landlab is not installed.

---

## Files

| File | What it tests |
|------|---------------|
| `test_1D_FD.py` | 1-D finite-difference solver (`F1D`) |
| `test_1D_SAS.py` | 1-D spatial-domain analytical solution (`F1D`, SAS method) |
| `test_2D_FD.py` | 2-D finite-difference solver (`F2D`) |
| `test_landlab.py` | Landlab component integration (skipped without Landlab) |

---

## Landlab integration tests (`test_landlab.py`)

Two independent analytical benchmarks are provided.

### 1 — Uniform load / isostatic equilibrium

For a spatially uniform load *q* (Pa) with all-periodic boundary
conditions the biharmonic term D∇⁴w vanishes for constant *w*, so the
governing equation reduces to

    (ρ_m − ρ_fill) · g · w = −q

giving the exact analytical result

    w = −q / ((ρ_m − ρ_fill) · g)

This is independent of elastic thickness and domain size.

### 2 — Point load / Kelvin-function solution

The deflection due to a concentrated vertical load *P* (N) at the origin
of an infinite 2-D elastic plate is (Turcotte & Schubert, 2002):

    w(r) = P α² / (2π D) · kei(r / α)

where

    D = E Te³ / (12(1 − ν²))        flexural rigidity [N·m]
    α = (D / (Δρ g))^(1/4)          2-D flexural parameter [m]
    kei(x) = Im[ e^{−iπ/4} K₀(x e^{iπ/4}) ]   Kelvin function

Because kei(x) ≤ 0 for x ≥ 0, a positive (downward) load produces a
negative deflection w, consistent with gFlex's sign convention.

**Note on the flexural parameter:**
The 2-D parameter uses D (not 4D as in the 1-D case).  The 1-D
parameter is α₁D = (4D / (Δρ g))^(1/4) = √2 · α.

The test uses a 100 × 100 grid at dx = 5 km (Te = 10 km, α ≈ 21 km,
domain ≈ 24 α wide) with `zero_moment_zero_shear` boundary conditions.
Comparison points are sampled in the interior at radii 1.5–3.5 α from
the load centre and at least 3 α from every boundary, keeping the
comparison well below the kei zero crossing (forebulge onset at
r/α ≈ 3.91) where the relative-error metric is undefined.  Tolerance
is 5 %.

---

## Grid-convergence results (2026-05-27)

A convergence study with E = 65 GPa, Te = 10 km, ν = 0.25,
ρ_m = 3300 kg m⁻³, ρ_fill = 0, g = 9.8 m s⁻² confirmed that the 2-D
FD scheme is **second-order accurate** (O(dx²)).

Grid spacings tested: dx = 10 000, 5 000, 2 500, 1 250 m
(dx/α ≈ 0.49, 0.24, 0.12, 0.06).

Convergence orders measured between successive refinements at exact
integer-multiple radii (all divisible by 1250 m to avoid sampling
jitter):

| r/α  | coarse→medium | medium→fine | fine→finest |
|------|---------------|-------------|-------------|
| 0.97 | 2.20          | 2.07        | 2.01        |
| 1.46 | 2.68          | 2.18        | 2.04        |
| 1.95 | 1.07          | 1.88        | 1.97        |
| 2.92 | −0.07         | 1.77        | 1.95        |

The −0.07 at r = 2.92 α, dx = 10 km is a near-cancellation artefact
(the error changes sign between those two resolutions); it is not a real
anomaly.

At the finest resolution (dx = 1250 m, dx/α ≈ 0.06) relative errors
are **0.02–0.1 %**.  The 5 % test tolerance is conservative for the
dx = 5 km grid used in `test_landlab.py`.
