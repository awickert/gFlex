"""Tests for 2-D FD on non-square grids (nx ≠ ny, dx ≠ dy).

Three properties are checked:

1. MMS convergence with dx ≠ dy — the manufactured-solution test from
   test_2D_FD.py is extended to a 2:1 wide grid (Nx = 2 Ny) with the
   domain kept square (Lx = Ly), so dx = Ly / (2N) ≠ dy = Ly / N.
   Second-order convergence in max(dx, dy) must still hold.

2. Bilateral symmetry — a square load centred on a non-square Periodic
   domain must produce a deflection that is symmetric in both x and y.
   This catches axis-swap bugs in the stencil assembly.

3. FD–SAS interior agreement on wide and tall domains — confirms that
   the correct infinite-plate deflection is reproduced regardless of
   aspect ratio.  The domain is large enough that boundary effects at the
   comparison region are < 5 % of peak deflection.

Grid convention: qs.shape = (ny, nx); y is axis-0, x is axis-1.

Physical parameters (shared):
  E = 65 GPa, Te = 30 km, nu = 0.25, rho_m = 3300 kg/m³, rho_fill = 0,
  g = 9.8 m/s².
  Flexural parameter  alpha = (D / (drho·g))^0.25 ≈ 47 km.
  At dx = dy = 4 km:  alpha / dx ≈ 11.75 cells.
"""

import numpy as np
import pytest

from gflex.f2d import F2D


# ---------------------------------------------------------------------------
# Shared physical parameters
# ---------------------------------------------------------------------------

_E        = 65e9
_nu       = 0.25
_Te       = 30e3
_rho_m    = 3300.0
_rho_fill = 0.0
_g        = 9.8

_D    = _E * _Te**3 / (12.0 * (1.0 - _nu**2))
_drho = _rho_m - _rho_fill
_alpha = (_D / (_drho * _g)) ** 0.25   # ≈ 47 km


def _run(ny, nx, qs, dx, dy, method="FD", bc="0Moment0Shear"):
    flex = F2D()
    flex.Quiet = True
    flex.Method = method
    flex.g = _g
    flex.E = _E
    flex.nu = _nu
    flex.rho_m = _rho_m
    flex.rho_fill = _rho_fill
    flex.Te = _Te
    flex.qs = qs.copy()
    flex.dx = dx
    flex.dy = dy
    flex.BC_W = flex.BC_E = flex.BC_N = flex.BC_S = bc
    flex.initialize()
    flex.run()
    return flex


# ---------------------------------------------------------------------------
# 1. MMS convergence with dx ≠ dy
# ---------------------------------------------------------------------------

def test_2d_fd_nonsquare_convergence_anisotropic():
    """FD achieves O(dx²) convergence on a 2:1 non-square grid with dx ≠ dy.

    Manufactured solution:  w_exact = cos(β_x · x) · cos(β_y · y)
    with β_x = 2π / Lx, β_y = 2π / Ly.

    The biharmonic identity
        ∇⁴[cos(β_x x) cos(β_y y)] = (β_x² + β_y²)² cos(β_x x) cos(β_y y)
    gives the manufactured load
        q = (D · (β_x² + β_y²)² + k) · cos(β_x x) · cos(β_y y).

    Periodic BCs on all sides eliminate boundary truncation error so the
    measured rate reflects only the interior stencil accuracy.

    Grid: Nx = 2·N cells, Ny = N cells, square domain Lx = Ly = L.
    This gives dx = L / (2N), dy = L / N, so dx ≠ dy at every refinement.
    The coarser spacing dy sets the asymptotic convergence rate.  Three
    refinements (N = 20, 40, 80) must each show rate > 1.8.
    """
    k = _drho * _g

    # Use the flexural parameter to fix domain size so both D and k contribute
    L = 2.0 * _alpha
    beta_x = 2.0 * np.pi / L
    beta_y = 2.0 * np.pi / L   # same wavelength; the key test is dx ≠ dy

    q_factor = _D * (beta_x**2 + beta_y**2)**2 + k

    def q_mms(X, Y):
        return q_factor * np.cos(beta_x * X) * np.cos(beta_y * Y)

    def w_exact(X, Y):
        # gFlex sign convention: positive load → negative w
        return -np.cos(beta_x * X) * np.cos(beta_y * Y)

    Ns = [20, 40, 80]
    errors = []

    for N in Ns:
        Nx, Ny = 2 * N, N       # 2:1 aspect ratio
        dx = L / Nx
        dy = L / Ny             # dy = 2 · dx

        x = (np.arange(Nx) + 0.5) * dx
        y = (np.arange(Ny) + 0.5) * dy
        X, Y = np.meshgrid(x, y)   # shape (Ny, Nx)

        qs = q_mms(X, Y)
        flex = _run(Ny, Nx, qs, dx, dy, bc="Periodic")

        err = np.max(np.abs(flex.w - w_exact(X, Y)))
        errors.append(err)

    for i in range(len(Ns) - 1):
        rate = np.log(errors[i] / errors[i + 1]) / np.log(Ns[i + 1] / Ns[i])
        assert rate > 1.8, (
            f"Expected O(dy²) convergence (rate ≥ 1.8) on 2:1 non-square grid; "
            f"got {rate:.2f} between N={Ns[i]} and N={Ns[i+1]} "
            f"(errors: {errors[i]:.3g} → {errors[i+1]:.3g} m)"
        )


# ---------------------------------------------------------------------------
# 2. Bilateral symmetry on non-square Periodic domain
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("nx,ny", [
    (120, 60),   # 2:1 wide
    (60, 120),   # 1:2 tall
    (160, 40),   # 4:1 wide
    (40, 160),   # 1:4 tall
])
def test_2d_fd_nonsquare_symmetry(nx, ny):
    """Symmetric load on a non-square Periodic domain → symmetric deflection.

    A square 6×6-cell load is placed at the domain centre.  Periodic BCs
    extend the solution periodically, so x-reflection about the load centre
    (w[i, j] = w[i, -j]) and y-reflection (w[i, j] = w[-i, j]) must hold
    to within floating-point rounding.

    This catches axis-swap bugs in the stencil (e.g. dx used where dy is
    needed or nx/ny transposed in array construction).
    """
    dx = dy = 4000.0
    qs = np.zeros((ny, nx))

    cy, cx = ny // 2, nx // 2
    qs[cy - 3 : cy + 3, cx - 3 : cx + 3] = 1e6

    flex = _run(ny, nx, qs, dx, dy, bc="Periodic")
    w = flex.w

    # x-symmetry: w[i, j] == w[i, nx-1-j]  (load is centred on nx/2)
    np.testing.assert_allclose(w, w[:, ::-1], rtol=1e-8, atol=1e-4,
                               err_msg=f"x-symmetry failed for {nx}×{ny} grid")

    # y-symmetry: w[i, j] == w[ny-1-i, j]
    np.testing.assert_allclose(w, w[::-1, :], rtol=1e-8, atol=1e-4,
                               err_msg=f"y-symmetry failed for {nx}×{ny} grid")


# ---------------------------------------------------------------------------
# 3. FD–SAS agreement near the load on wide and tall domains
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("nx,ny", [
    (200, 120),   # 5:3 wide;  load center 240 km ≈ 5.1 α from every edge
    (120, 200),   # 3:5 tall;  same geometry transposed
    (300, 120),   # 5:2 wide;  short axis still 240 km = 5.1 α
    (120, 300),   # 2:5 tall;  same
])
def test_2d_fd_nonsquare_agrees_with_sas(nx, ny):
    """FD on non-square grids reproduces the SAS solution near the load.

    The comparison is restricted to a ±10-cell window centred on the load,
    where deflections are ≥ 40 % of peak and BC corrections are negligible.

    Geometry: the load centre is nx//2, ny//2.  For all test cases the
    short axis has ny//2 (or nx//2) ≥ 60 cells from the domain edge, giving
    a distance of 240 km ≈ 5.1 α.  At 5.1 α the BC correction decays to
    exp(−5.1) ≈ 0.6 % of peak; within the ±10-cell window it stays below
    1 % everywhere.  The 5 % tolerance has ample headroom.
    """
    dx = dy = 4000.0
    qs = np.zeros((ny, nx))
    cy, cx = ny // 2, nx // 2
    half_load = 3
    qs[cy - half_load : cy + half_load, cx - half_load : cx + half_load] = 1e6

    flex_fd  = _run(ny, nx, qs, dx, dy, method="FD",  bc="0Moment0Shear")
    flex_sas = _run(ny, nx, qs, dx, dy, method="SAS", bc="0Moment0Shear")

    half = 10
    w_fd  = flex_fd.w [cy - half : cy + half, cx - half : cx + half]
    w_sas = flex_sas.w[cy - half : cy + half, cx - half : cx + half]

    np.testing.assert_allclose(
        w_fd, w_sas, rtol=0.05,
        err_msg=f"FD–SAS mismatch on {nx}×{ny} grid"
    )
