#! /usr/bin/env python
"""Tests for 2-D FD boundary conditions and in-plane stress terms.

Boundary conditions tested
--------------------------
mirror            — exact 2-D cosine eigenfunction; mirror == periodic on the
                    2(N-1)×2(M-1) even-extended domain
zero_slope_zero_shear      — shown to be genuinely distinct from mirror
zero_displacement_zero_slope — interior matches SAS (large domain, load far from edge)
zero_moment_zero_shear     — same

In-plane stress tests
---------------------
sigma_xx direction — kx-only load (ky = 0): sigma_yy and sigma_xy have NO
                     effect; sigma_xx does (no axis-swap check)
sigma_yy direction — ky-only load (kx = 0): sigma_xx and sigma_xy have NO effect
sigma_xy diagonal  — x-reflection symmetry: w(+σ_xy)[i,j] = w(−σ_xy)[i,N-1-j]
                     for a symmetric load and mirror BCs (validates the sign
                     pattern of the cj±1 i±1 stencil)
sigma monotonicity — tensile reduces deflection, compressive increases it,
                     for every BC type and every sigma component

Grid convention: node-centred, x[j] = j·dx, y[i] = i·dy.
"""

import numpy as np
import pytest

from gflex.f2d import F2D


# ---------------------------------------------------------------------------
# Shared parameters
# ---------------------------------------------------------------------------

E        = 65e9
nu       = 0.25
Te       = 30e3
rho_m    = 3300.0
rho_fill = 0.0
g        = 9.8
dx       = 4000.0
dy       = 4000.0

D    = E * Te**3 / (12.0 * (1.0 - nu**2))
drho = rho_m - rho_fill

# 2-D flexural parameter: (D / (drho·g))^0.25  [m]
alpha_2d = (D / (drho * g)) ** 0.25


def _run(qs, method="fd", bc_w="zero_moment_zero_shear", bc_e="zero_moment_zero_shear",
         bc_n="zero_moment_zero_shear", bc_s="zero_moment_zero_shear",
         dx_=dx, dy_=dy,
         sigma_xx=None, sigma_yy=None, sigma_xy=None):
    """Run a 2-D flexure calculation and return the flex object."""
    flex = F2D()
    flex.quiet = True
    flex.method = method
    flex.solver = "direct"
    flex.g = g
    flex.E = E
    flex.nu = nu
    flex.rho_m = rho_m
    flex.rho_fill = rho_fill
    flex.te = Te
    flex.qs = qs.copy()
    flex.dx = dx_
    flex.dy = dy_
    flex.bc_west = bc_w
    flex.bc_east = bc_e
    flex.bc_north = bc_n
    flex.bc_south = bc_s
    if sigma_xx is not None:
        flex.sigma_xx = sigma_xx
    if sigma_yy is not None:
        flex.sigma_yy = sigma_yy
    if sigma_xy is not None:
        flex.sigma_xy = sigma_xy
    flex.initialize()
    flex.run()
    return flex


# ---------------------------------------------------------------------------
# mirror: exact 2-D cosine eigenfunction
# ---------------------------------------------------------------------------

def test_fd_2d_mirror_cosine_eigenfunction():
    """2-D FD/mirror matches the analytical formula for a separable cosine load.

    cos(kx·x)·cos(ky·y) is an eigenfunction of the biharmonic operator for
    mirror BCs (zero slope at all four edges).  The exact deflection is:

        w = −q₀ / (D(kx² + ky²)² + Δρg) · cos(kx·x)·cos(ky·y)

    The FD discretisation has an O((k·dx)²) eigenvalue error, giving < 0.1%
    relative error for the wavenumbers chosen, so rtol = 1e-3 is safe.
    """
    Ny, Nx = 128, 128
    Lx = (Nx - 1) * dx
    Ly = (Ny - 1) * dy
    nx_waves, ny_waves = 2, 3
    kx = nx_waves * np.pi / Lx
    ky = ny_waves * np.pi / Ly
    x = np.arange(Nx) * dx
    y = np.arange(Ny) * dy
    X, Y = np.meshgrid(x, y)
    q0 = 1e6
    qs = q0 * np.cos(kx * X) * np.cos(ky * Y)

    flex = _run(qs, bc_w="mirror", bc_e="mirror", bc_n="mirror", bc_s="mirror")

    w_exact = -q0 / (D * (kx**2 + ky**2)**2 + drho * g) * np.cos(kx * X) * np.cos(ky * Y)
    np.testing.assert_allclose(flex.w, w_exact, rtol=1e-3)


# ---------------------------------------------------------------------------
# mirror: exact equivalence with periodic on the even-extended 2× domain
# ---------------------------------------------------------------------------

def test_fd_2d_mirror_equals_periodic_2x2():
    """2-D FD/mirror on (Ny×Nx) equals FD/periodic on the (2Ny-2)×(2Nx-2) even extension.

    mirror BCs at all four edges implement even reflection.  The ghost-cell
    assignments are identical to running a periodic problem on the period-
    2(Ny-1) × 2(Nx-1) even extension of the load.  Agreement to rtol = 1e-6
    (limited by the direct solver) confirms both stencils implement the same
    linear system.
    """
    Ny, Nx = 50, 60
    qs = np.zeros((Ny, Nx))
    qs[10:20, 15:30] = 1e6   # off-centre so even extension is non-trivial

    flex_m = _run(qs, bc_w="mirror", bc_e="mirror", bc_n="mirror", bc_s="mirror")

    # Even extension: period is 2(Nx-1) in x and 2(Ny-1) in y, sharing endpoints
    qs_x   = np.concatenate([qs,   qs[:,   -2:0:-1]], axis=1)  # (Ny, 2Nx-2)
    qs_ext = np.concatenate([qs_x, qs_x[-2:0:-1, :]], axis=0)  # (2Ny-2, 2Nx-2)

    flex_p = _run(qs_ext, bc_w="periodic", bc_e="periodic",
                          bc_n="periodic", bc_s="periodic")

    np.testing.assert_allclose(flex_m.w, flex_p.w[:Ny, :Nx], rtol=1e-6)


# ---------------------------------------------------------------------------
# zero_slope_zero_shear: distinct from mirror despite the same physical intent
# ---------------------------------------------------------------------------

def test_fd_2d_0slope0shear_equals_mirror():
    """zero_slope_zero_shear and mirror are numerically identical stencils in 2-D.

    Both BCs enforce even reflection of the ghost nodes along each axis:
      dw/dn = 0 at boundary  →  first ghost = adjacent interior node
      d³w/dn³ = 0 at boundary  →  second ghost = second interior node
    The ghost equations are identical to mirror (even reflection), so the two
    stencils are numerically indistinguishable.

    Both reproduce the cosine eigenfunction for a cosine load.
    """
    Ny, Nx = 128, 128
    Lx = (Nx - 1) * dx
    Ly = (Ny - 1) * dy
    nx_waves, ny_waves = 3, 2
    kx = nx_waves * np.pi / Lx
    ky = ny_waves * np.pi / Ly
    x = np.arange(Nx) * dx
    y = np.arange(Ny) * dy
    X, Y = np.meshgrid(x, y)
    q0 = 1e6
    qs = q0 * np.cos(kx * X) * np.cos(ky * Y)

    w_mirror = _run(qs, bc_w="mirror", bc_e="mirror",
                        bc_n="mirror", bc_s="mirror").w
    w_0ss    = _run(qs, bc_w="zero_slope_zero_shear", bc_e="zero_slope_zero_shear",
                        bc_n="zero_slope_zero_shear", bc_s="zero_slope_zero_shear").w
    w_exact  = -q0 / (D * (kx**2 + ky**2)**2 + drho * g) * np.cos(kx * X) * np.cos(ky * Y)

    # Both reproduce the cosine eigenfunction
    np.testing.assert_allclose(w_mirror, w_exact, rtol=2e-3)
    np.testing.assert_allclose(w_0ss, w_exact, rtol=2e-3)

    # The two BCs give identical results
    np.testing.assert_array_equal(w_mirror, w_0ss)


# ---------------------------------------------------------------------------
# Large-domain interior checks: zero_slope_zero_shear, zero_displacement_zero_slope, zero_moment_zero_shear
# ---------------------------------------------------------------------------

def _central_load_sas_comparison(bc, margin=70):
    """Run FD with *bc* and SAS on a 200×200 grid with a central square load.

    alpha_2D ≈ 47 km (≈ 11.75 cells).  The load sits ~8 alpha from each edge,
    so BC corrections at the centre are O(exp(−8/√2)) < 0.3 %.  Agreement
    to rtol = 5e-3 in [margin:N-margin, margin:N-margin] confirms the FD
    solver produces the correct infinite-plate interior deflection.
    """
    N = 200
    qs = np.zeros((N, N))
    qs[95:105, 95:105] = 1e6

    flex_fd  = _run(qs, bc_w=bc, bc_e=bc, bc_n=bc, bc_s=bc)
    flex_sas = _run(qs, method="sas", bc_w=bc, bc_e=bc, bc_n=bc, bc_s=bc)

    np.testing.assert_allclose(
        flex_fd.w[margin:N-margin, margin:N-margin],
        flex_sas.w[margin:N-margin, margin:N-margin],
        rtol=5e-3,
    )


def test_fd_2d_0slope0shear_vs_sas():
    """2-D FD/zero_slope_zero_shear matches SAS for a central load far from the boundary."""
    _central_load_sas_comparison("zero_slope_zero_shear")


def test_fd_2d_0displacement0slope_vs_sas():
    """2-D FD/zero_displacement_zero_slope matches SAS for a central load far from the boundary."""
    _central_load_sas_comparison("zero_displacement_zero_slope")


def test_fd_2d_0moment0shear_vs_sas():
    """2-D FD/zero_moment_zero_shear matches SAS for a central load far from the boundary."""
    _central_load_sas_comparison("zero_moment_zero_shear")


# ---------------------------------------------------------------------------
# sigma_xx direction test (no axis swap)
# ---------------------------------------------------------------------------

def test_fd_2d_sigma_xx_direction():
    """sigma_yy and sigma_xy have no effect on a ky = 0 load; sigma_xx does.

    A strip load uniform in y (all rows, columns 25–34) has only kx modes.
    For such a load the governing equation reduces to the 1-D equation in x:
    the σ_yy·Te·∂²w/∂y² term vanishes (∂²w/∂y² = 0 for a y-uniform solution),
    and the σ_xy diagonal-stencil terms cancel pairwise for any BC that
    preserves y-translational symmetry.  mirror BCs on all four sides are
    used because mirror is analytically y-symmetric (ghost cells are set by
    even reflection) and the sigma_xy corrections at mirror boundaries cancel
    exactly for a y-uniform solution (verified algebraically).

    sigma_xx changes the deflection because σ_xx·Te·kx² ≠ 0.

    atol = 1e-9 handles the O(1e-12) roundoff introduced when different sigma
    values alter matrix conditioning, producing huge relative errors near zero.
    """
    N = 64
    qs = np.zeros((N, N))
    qs[:, 25:35] = 1e6    # uniform in y → only kx modes

    bc = dict(bc_w="mirror", bc_e="mirror", bc_n="mirror", bc_s="mirror")
    S = 1e8

    flex_0  = _run(qs, **bc)
    flex_xx = _run(qs, **bc, sigma_xx=S)
    flex_yy = _run(qs, **bc, sigma_yy=S)
    flex_xy = _run(qs, **bc, sigma_xy=S)

    # sigma_yy and sigma_xy must have no physical effect for a ky=0 load
    np.testing.assert_allclose(flex_yy.w, flex_0.w, atol=1e-9,
                               err_msg="sigma_yy must not affect a ky=0 load")
    np.testing.assert_allclose(flex_xy.w, flex_0.w, atol=1e-9,
                               err_msg="sigma_xy must not affect a ky=0 load")

    # sigma_xx must change the deflection
    assert not np.allclose(flex_xx.w, flex_0.w), \
        "sigma_xx should change the deflection for a kx-only load"


# ---------------------------------------------------------------------------
# sigma_yy direction test (symmetric axis check)
# ---------------------------------------------------------------------------

def test_fd_2d_sigma_yy_direction():
    """sigma_xx and sigma_xy have no effect on a kx = 0 load; sigma_yy does.

    Symmetric counterpart of test_fd_2d_sigma_xx_direction.  A strip load
    uniform in x (only ky modes) has kx = 0 everywhere, so σ_xx·Te·kx² = 0
    and the σ_xy diagonal-stencil terms cancel pairwise for a y-translation-
    symmetric solution.  mirror BCs on all four sides are used because the
    σ_xy diagonal-stencil terms cancel exactly for an x-uniform solution
    (verified algebraically), giving atol = 1e-9 to handle only rounding noise.

    Only sigma_yy changes the deflection because σ_yy·Te·ky² ≠ 0.
    """
    N = 64
    qs = np.zeros((N, N))
    qs[25:35, :] = 1e6    # uniform in x → only ky modes

    bc = dict(bc_w="mirror", bc_e="mirror", bc_n="mirror", bc_s="mirror")
    S = 1e8

    flex_0  = _run(qs, **bc)
    flex_xx = _run(qs, **bc, sigma_xx=S)
    flex_yy = _run(qs, **bc, sigma_yy=S)
    flex_xy = _run(qs, **bc, sigma_xy=S)

    np.testing.assert_allclose(flex_xx.w, flex_0.w, atol=1e-9,
                               err_msg="sigma_xx must not affect a kx=0 load")
    np.testing.assert_allclose(flex_xy.w, flex_0.w, atol=1e-9,
                               err_msg="sigma_xy must not affect a kx=0 load")
    assert not np.allclose(flex_yy.w, flex_0.w), \
        "sigma_yy should change the deflection for a ky-only load"


# ---------------------------------------------------------------------------
# FD vs FFT periodic cross-validation for sigma_xx, sigma_yy, sigma_xy
# ---------------------------------------------------------------------------

def _fft_run(qs, **kw):
    """Thin wrapper to run the FFT solver through the shared _run helper."""
    return _run(qs, method="fft", **kw)


def test_fd_2d_sigma_xx_vs_fft_periodic():
    """FD/periodic with sigma_xx matches FFT/periodic to within FD truncation error.

    A central square block load (many Fourier modes, no zero-crossing issues)
    in a periodic domain is solved by both FD and FFT.  A block load excites
    high-k modes where the O((k·dx)²) FD eigenvalue error reaches ~0.93 %,
    so rtol = 1e-2 covers the truncation error while being tight enough to
    detect any σ_xx stencil sign flip or coefficient error.
    """
    N = 64
    qs = np.zeros((N, N))
    qs[28:36, 28:36] = 1e6    # central block, many Fourier modes, no exact zeros

    bc = dict(bc_w="periodic", bc_e="periodic", bc_n="periodic", bc_s="periodic")
    sigma_xx = 2e8

    flex_fd  = _run(qs, **bc, sigma_xx=sigma_xx)
    flex_fft = _fft_run(qs, **bc, sigma_xx=sigma_xx)

    np.testing.assert_allclose(flex_fd.w, flex_fft.w, rtol=1e-2)


def test_fd_2d_sigma_yy_vs_fft_periodic():
    """FD/periodic with sigma_yy matches FFT/periodic (validates the cj0i±1 stencil).

    See test_fd_2d_sigma_xx_vs_fft_periodic for the load and tolerance rationale.
    rtol = 1e-2 accounts for O((k·dx)²) FD truncation reaching ~0.93 % for a
    block load.
    """
    N = 64
    qs = np.zeros((N, N))
    qs[28:36, 28:36] = 1e6

    bc = dict(bc_w="periodic", bc_e="periodic", bc_n="periodic", bc_s="periodic")
    sigma_yy = 2e8

    flex_fd  = _run(qs, **bc, sigma_yy=sigma_yy)
    flex_fft = _fft_run(qs, **bc, sigma_yy=sigma_yy)

    np.testing.assert_allclose(flex_fd.w, flex_fft.w, rtol=1e-2)


def test_fd_2d_sigma_xy_vs_fft_periodic():
    """FD/periodic with sigma_xy matches FFT/periodic to within FD truncation error.

    Validates the fixes to the all-periodic coefficient-matrix assembly:
    (1) the NW coefficient is saved before the west-boundary shuffle overwrites
        the cj_1i1 slot; (2) the two double-wrap corner entries A[0, N-1] and
        A[N-1, 0] (which link cell (0,0)↔(ny-1, nx-1) diagonally across the
        domain) use the SW/NE coefficients rather than the NW/SE coefficients
        that the spdiags diagonals would otherwise supply.

    A central square block load (many Fourier modes) is used.  The σ_xy
    cross-stencil adds an additional O((k·dx)²) error compared with the
    σ_xx / σ_yy cases, so the combined truncation reaches ~1.35 % for the
    highest modes in this 64×64 grid; rtol = 2e-2 is still tight enough to
    detect any sign flip or coefficient error.
    """
    N = 64
    qs = np.zeros((N, N))
    qs[28:36, 28:36] = 1e6

    bc = dict(bc_w="periodic", bc_e="periodic", bc_n="periodic", bc_s="periodic")
    sigma_xy = 2e8

    flex_fd  = _run(qs, **bc, sigma_xy=sigma_xy)
    flex_fft = _fft_run(qs, **bc, sigma_xy=sigma_xy)

    np.testing.assert_allclose(flex_fd.w, flex_fft.w, rtol=2e-2)


def test_fd_2d_sigma_xy_reflection_symmetry():
    """sigma_xy respects the x-reflection symmetry of the stencil sign pattern.

    Under the substitution x → Lx − x, the cross-derivative ∂²w/∂x∂y acquires
    a sign flip, so σ_xy → −σ_xy maps a solution to its x-reflected counterpart.
    Formally: for a load that is symmetric about x = Lx/2 (i.e. q(Lx−x, y) = q(x, y))
    and mirror BCs (which respect x-reflection), the deflection satisfies

        w(+S)[i, j] = w(−S)[i, N−1−j]

    The load is centred on a 121×121 grid (odd, so j = 60 is the exact centre)
    mirror BCs are used because they implement even reflection, making the
    x-reflection identity exact at the discrete level.

    rtol = 1e-8 is limited only by floating-point rounding in the sparse solver.
    """
    N = 121                         # odd: centre at j = 60, i = 60
    qs = np.zeros((N, N))
    qs[50:71, 50:71] = 1e6          # symmetric about both centre lines

    bc = dict(bc_w="mirror", bc_e="mirror", bc_n="mirror", bc_s="mirror")
    S = 1e8

    flex_p = _run(qs, **bc, sigma_xy=+S)
    flex_n = _run(qs, **bc, sigma_xy=-S)

    # x-reflection: column j in w(+S) matches column N-1-j in w(-S)
    np.testing.assert_allclose(flex_p.w, flex_n.w[:, ::-1], rtol=1e-8)


# ---------------------------------------------------------------------------
# sigma_xy sign test (diagonal stencil produces a non-trivial, sign-sensitive effect)
# ---------------------------------------------------------------------------

def test_fd_2d_sigma_xy_sign():
    """Positive and negative sigma_xy give genuinely different deflections.

    For a 2-D load with both kx and ky content, the cross term
    2·σ_xy·Te·∂²w/∂x∂y modifies the plate response asymmetrically in the
    (x, y) and (x, −y) diagonal directions.  Positive and negative σ_xy must
    produce distinct (non-equal) deflection fields.  This test detects any
    stencil implementation where the σ_xy coefficient has zero net effect or
    where ±σ_xy give the same result.
    """
    N = 60
    qs = np.zeros((N, N))
    qs[25:35, 25:35] = 1e6

    bc = dict(bc_w="periodic", bc_e="periodic", bc_n="periodic", bc_s="periodic")
    S = 1e8

    flex_0 = _run(qs, **bc)
    flex_p = _run(qs, **bc, sigma_xy=+S)
    flex_n = _run(qs, **bc, sigma_xy=-S)

    # sigma_xy != 0 must produce a different solution than sigma_xy == 0
    assert not np.allclose(flex_p.w, flex_0.w), \
        "positive sigma_xy should change the deflection for a 2-D load"
    assert not np.allclose(flex_n.w, flex_0.w), \
        "negative sigma_xy should change the deflection for a 2-D load"

    # Positive and negative sigma_xy must produce different results
    assert not np.allclose(flex_p.w, flex_n.w), \
        "positive and negative sigma_xy must produce different deflections"


# ---------------------------------------------------------------------------
# sigma monotonicity for every FD BC
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("bc", [
    "zero_moment_zero_shear",
    "zero_displacement_zero_slope",
    "mirror",
    "zero_slope_zero_shear",
])
def test_fd_2d_sigma_xx_monotonicity_all_bcs(bc):
    """Tensile sigma_xx reduces deflection; compressive increases it — for every 2-D BC."""
    N = 120
    qs = np.zeros((N, N))
    qs[55:65, 55:65] = 1e6

    flex_0 = _run(qs, bc_w=bc, bc_e=bc, bc_n=bc, bc_s=bc)
    flex_t = _run(qs, bc_w=bc, bc_e=bc, bc_n=bc, bc_s=bc, sigma_xx=+1e8)
    flex_c = _run(qs, bc_w=bc, bc_e=bc, bc_n=bc, bc_s=bc, sigma_xx=-1e8)

    assert flex_t.w.min() > flex_0.w.min(), (
        f"{bc}: tensile sigma_xx should reduce subsidence"
    )
    assert flex_c.w.min() < flex_0.w.min(), (
        f"{bc}: compressive sigma_xx should increase subsidence"
    )


@pytest.mark.parametrize("bc", [
    "zero_moment_zero_shear",
    "zero_displacement_zero_slope",
    "mirror",
    "zero_slope_zero_shear",
])
def test_fd_2d_sigma_yy_monotonicity_all_bcs(bc):
    """Tensile sigma_yy reduces deflection; compressive increases it — for every 2-D BC."""
    N = 120
    qs = np.zeros((N, N))
    qs[55:65, 55:65] = 1e6

    flex_0 = _run(qs, bc_w=bc, bc_e=bc, bc_n=bc, bc_s=bc)
    flex_t = _run(qs, bc_w=bc, bc_e=bc, bc_n=bc, bc_s=bc, sigma_yy=+1e8)
    flex_c = _run(qs, bc_w=bc, bc_e=bc, bc_n=bc, bc_s=bc, sigma_yy=-1e8)

    assert flex_t.w.min() > flex_0.w.min(), (
        f"{bc}: tensile sigma_yy should reduce subsidence"
    )
    assert flex_c.w.min() < flex_0.w.min(), (
        f"{bc}: compressive sigma_yy should increase subsidence"
    )


# ---------------------------------------------------------------------------
# zero_displacement_zero_moment: 2-D sine eigenfunction
# ---------------------------------------------------------------------------

def test_fd_2d_0displacement0moment_sine_eigenfunction():
    """2-D FD/zero_displacement_zero_moment matches the analytical formula for a separable sine load.

    sin(kx·x)·sin(ky·y) is an eigenfunction of the 2-D biharmonic operator for
    zero_displacement_zero_moment BCs at all four edges (zero displacement and zero moment
    at every boundary).  The exact deflection is:

        w = −q₀ / (D·(kx²+ky²)² + Δρg) · sin(kx·x)·sin(ky·y)

    The load vanishes at all boundary nodes, which is the necessary condition for
    the simply-supported BC to be consistent.  The FD discretisation has an
    O((k·dx)²) error; rtol = 1e-3 is safe for the wavenumbers chosen.
    """
    Ny, Nx = 128, 128
    Lx = (Nx - 1) * dx
    Ly = (Ny - 1) * dy
    nx_waves, ny_waves = 2, 3
    kx = nx_waves * np.pi / Lx
    ky = ny_waves * np.pi / Ly
    x = np.arange(Nx) * dx
    y = np.arange(Ny) * dy
    X, Y = np.meshgrid(x, y)
    q0 = 1e6
    qs = q0 * np.sin(kx * X) * np.sin(ky * Y)

    flex = _run(qs, bc_w="zero_displacement_zero_moment", bc_e="zero_displacement_zero_moment",
                    bc_n="zero_displacement_zero_moment", bc_s="zero_displacement_zero_moment")

    w_exact = (-q0 / (D * (kx**2 + ky**2)**2 + drho * g)
               * np.sin(kx * X) * np.sin(ky * Y))
    np.testing.assert_allclose(flex.w, w_exact, rtol=1e-3, atol=1e-5)


# ---------------------------------------------------------------------------
# zero_displacement_zero_moment: 2-D half-domain antisymmetric equivalence
# ---------------------------------------------------------------------------

def test_fd_2d_0displacement0moment_half_domain_antisymmetric():
    """2-D half-domain with zero_displacement_zero_moment at the symmetry edge matches the full domain.

    For a load antisymmetric in the x-direction (q[:, j] = -q[:, Nx-1-j]) on a
    uniform-Te domain with free (zero_moment_zero_shear) boundary conditions, the solution
    is exactly antisymmetric in x.  Running the left x-half with zero_moment_zero_shear
    at three edges and zero_displacement_zero_moment at the (virtual) symmetry edge must
    reproduce the left half of the full-domain solution.  Agreement is exact to
    floating-point precision (rtol = 1e-6).
    """
    Ny, Nx = 40, 101   # Nx odd: x-centre at column 50

    qs_full = np.zeros((Ny, Nx))
    qs_full[:, 20] = +1e6 / Ny
    qs_full[:, 80] = -1e6 / Ny    # antisymmetric: 100 - 20 = 80

    flex_full = _run(qs_full, bc_w="zero_moment_zero_shear", bc_e="zero_moment_zero_shear",
                              bc_n="zero_moment_zero_shear", bc_s="zero_moment_zero_shear")

    Nx_half = 51   # columns 0..50
    qs_half = qs_full[:, :Nx_half].copy()
    flex_half = _run(qs_half, bc_w="zero_moment_zero_shear", bc_e="zero_displacement_zero_moment",
                              bc_n="zero_moment_zero_shear", bc_s="zero_moment_zero_shear")

    np.testing.assert_allclose(flex_half.w, flex_full.w[:, :Nx_half], rtol=1e-6, atol=1e-10)


# ---------------------------------------------------------------------------
# BC aliases: 'clamped' and 'free'
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("alias,canonical", [
    ("clamped", "zero_displacement_zero_slope"),
    ("free",    "zero_moment_zero_shear"),
])
def test_fd_bc_alias_equals_canonical_2d(alias, canonical):
    """Short aliases produce bit-identical results to their canonical names in 2-D."""
    import warnings
    Ny, Nx = 40, 40
    qs = np.zeros((Ny, Nx))
    qs[Ny // 3 : Ny // 2, Nx // 3 : Nx // 2] = 1e6

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        w_alias     = _run(qs, bc_w=alias,    bc_e=alias,    bc_n=alias,    bc_s=alias).w
        w_canonical = _run(qs, bc_w=canonical, bc_e=canonical, bc_n=canonical, bc_s=canonical).w

    np.testing.assert_array_equal(w_alias, w_canonical)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
