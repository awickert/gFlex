"""Benchmarks for the gFlex Landlab component against analytical solutions.

Skipped automatically when Landlab is not installed.

Detailed interface tests live in the Landlab repository at
tests/components/gflex/.

Two independent analytical benchmarks are provided:

1. Uniform load / isostatic equilibrium
   For a spatially uniform load q [Pa] with all-periodic boundary conditions
   the biharmonic term vanishes (∇⁴w = 0 for constant w) and the exact
   solution is simply isostatic equilibrium:

       w = -q / ((rho_m - rho_fill) * g)

2. Point load / Kelvin-function solution
   The deflection due to a concentrated vertical load P [N] at the origin
   of an infinite 2-D elastic plate is given by the Kelvin (Thomson) function:

       w(r) = P * alpha² / (2π * D) * kei(r / alpha)

   where
       D     = E * Te³ / (12 * (1 - nu²))   [N·m]  flexural rigidity
       alpha = (D / (drho * g))^(1/4)        [m]    2-D flexural parameter
       kei   = Im[ e^{-iπ/4} * K₀(r e^{iπ/4}) ]    Kelvin function

   Note: the 2-D flexural parameter uses D, not 4*D as in the 1-D case.

   Because kei(x) < 0 for x ≥ 0, the deflection w is negative (downward)
   for a positive (downward) load, consistent with the Landlab component's
   sign convention for lithosphere_surface__elevation_increment.

   The comparison is made at several points in the interior of the grid,
   far enough from both the load and the domain boundaries that (a) the
   FD discretisation error is small and (b) boundary reflections are
   negligible.  The upper radius is kept below the forebulge onset
   (kei zero crossing at r/alpha ≈ 3.91) where relative error is undefined.
"""

import numpy as np
import pytest
from scipy.special import kei

pytest.importorskip("landlab")


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _make_grid(nrows, ncols, spacing):
    from landlab import RasterModelGrid
    mg = RasterModelGrid((nrows, ncols), xy_spacing=spacing)
    mg.add_zeros("surface_load__stress", at="node")
    return mg



# ---------------------------------------------------------------------------
# Benchmark 1 — uniform load / isostatic equilibrium
# ---------------------------------------------------------------------------

def test_uniform_load_isostatic_deflection():
    """Uniform load with all-periodic BCs yields exact isostatic deflection.

    For uniform w, ∇⁴w = 0, so the governing equation reduces to:
        (rho_m - rho_fill) * g * w = -q
    This is independent of Te and the domain size.
    """
    from landlab.components.gflex.flexure import gFlex

    rho_m, rho_fill, g, q = 3300.0, 0.0, 9.8, 1e4   # Pa

    mg = _make_grid(20, 20, 25000.0)
    mg.at_node["surface_load__stress"][:] = q

    gf = gFlex(
        mg,
        rho_mantle=rho_m, rho_fill=rho_fill, g=g,
        BC_W="Periodic", BC_E="Periodic",
        BC_N="Periodic", BC_S="Periodic",
        quiet=True,
    )
    gf.run_one_step()

    w = mg.at_node["lithosphere_surface__elevation_increment"][mg.core_nodes]
    expected = -q / ((rho_m - rho_fill) * g)   # ≈ −0.309 m
    np.testing.assert_allclose(w, expected, rtol=1e-3)


# ---------------------------------------------------------------------------
# Benchmark 2 — point load / Kelvin-function analytical solution
# ---------------------------------------------------------------------------

def test_point_load_kelvin_function():
    """FD deflection matches the Kelvin-function infinite-plate solution.

    Set-up
    ------
    * 100 × 100 grid, dx = dy = 5 km  →  500 km domain
    * Te = 10 km  →  alpha ≈ 21 km  →  domain ≈ 24 alpha wide
    * Central cell loaded with q = 1e6 Pa
    * 0Moment0Shear BCs on all edges (free edges, minimal reflection)

    Comparison points are chosen at radii between 1.5 and 3.5 alpha from
    the load centre (below the forebulge onset at ~3.9 alpha) and at least
    3 alpha from every boundary, so that both near-field FD discretisation
    error and boundary reflections are small.  A tolerance of 5 % is used.
    """
    from landlab.components.gflex.flexure import gFlex

    # Physical parameters
    E   = 65e9      # Pa
    Te  = 10000.0   # m
    nu  = 0.25
    rho_m    = 3300.0
    rho_fill = 0.0
    g        = 9.8

    dx = dy = 5000.0        # m
    nrows = ncols = 100

    D     = E * Te**3 / (12.0 * (1.0 - nu**2))          # flexural rigidity [N·m]
    drho  = rho_m - rho_fill
    alpha = (D / (drho * g)) ** 0.25                    # 2-D flexural parameter [m]

    # Build grid and apply a point load at the central cell
    mg = _make_grid(nrows, ncols, dx)
    ci = nrows // 2   # centre row index
    cj = ncols // 2   # centre col index
    q_load = 1e6      # Pa

    mg.at_node["surface_load__stress"][mg.grid_coords_to_node_id(ci, cj)] = q_load

    gf = gFlex(
        mg,
        Youngs_modulus=E, Poissons_ratio=nu,
        rho_mantle=rho_m, rho_fill=rho_fill, g=g,
        elastic_thickness=Te,
        BC_W="0Moment0Shear", BC_E="0Moment0Shear",
        BC_N="0Moment0Shear", BC_S="0Moment0Shear",
        quiet=True,
    )
    gf.run_one_step()

    w_grid = mg.at_node["lithosphere_surface__elevation_increment"].reshape(
        mg.shape
    )

    # Total applied force
    P = q_load * dx * dy   # N

    # Analytical solution:  w(r) = P * alpha² / (2π D) * kei(r/alpha)
    # kei(x) < 0 for x ≥ 0, so w < 0 for a downward (positive) load.
    def w_analytical(r):
        return P * alpha**2 / (2.0 * np.pi * D) * kei(r / alpha)

    # Sample points on the grid row passing through the load centre (east).
    # Stay below 3.5 alpha to avoid the forebulge (kei zero at r/alpha ≈ 3.91).
    # Restrict to cells at least 3 alpha from the eastern boundary.
    min_r  = 1.5 * alpha
    max_r  = 3.5 * alpha
    min_from_boundary = 3.0 * alpha

    errors = []
    for dj in range(1, ncols - cj):
        r = dj * dx
        col_idx = cj + dj
        dist_from_east = (ncols - 1 - col_idx) * dx
        if r < min_r or r > max_r:
            continue
        if dist_from_east < min_from_boundary:
            continue

        w_fd  = w_grid[ci, col_idx]
        w_ana = w_analytical(r)
        errors.append(abs(w_fd - w_ana) / abs(w_ana))

    assert len(errors) > 0, "No comparison points satisfied the distance criteria"
    assert max(errors) < 0.05, (
        f"Largest relative error {max(errors):.1%} exceeds 5 % tolerance; "
        "FD solution does not match the Kelvin-function analytical result."
    )
