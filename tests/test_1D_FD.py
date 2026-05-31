#! /usr/bin/env python
"""Tests for the 1-D finite-difference and FFT solvers, including end-load (sigma_xx) support."""

import numpy as np
import pytest

from gflex.f1d import F1D


def _run_flex_1d(Te, qs, dx, method="fd", bc_w="zero_moment_zero_shear",
                 bc_e="zero_moment_zero_shear", sigma_xx=None, solver="direct"):
    """Helper: run a 1-D flexure calculation and return the flex object."""
    flex = F1D()
    flex.quiet = True
    flex.method = method
    flex.solver = solver
    flex.g = 9.8
    flex.E = 65e9
    flex.nu = 0.25
    flex.rho_m = 3300.0
    flex.rho_fill = 0.0
    flex.te = Te
    flex.qs = qs.copy()
    flex.dx = dx
    flex.bc_west = bc_w
    flex.bc_east = bc_e
    if sigma_xx is not None:
        flex.sigma_xx = sigma_xx
    flex.initialize()
    flex.run()
    flex.finalize()
    return flex


# ---------------------------------------------------------------------------
# Basic smoke test
# ---------------------------------------------------------------------------

def test_1d_fd_runs():
    """1-D FD solver produces a non-trivial, finite deflection."""
    N = 200
    dx = 4000.0
    qs = np.zeros(N)
    qs[90:110] = 1e6
    flex = _run_flex_1d(Te=30e3, qs=qs, dx=dx)
    assert flex.w.ndim == 1
    assert flex.w.size == N
    assert not np.any(np.isnan(flex.w))
    assert flex.w.min() < 0


# ---------------------------------------------------------------------------
# End-load (sigma_xx) correctness tests
# ---------------------------------------------------------------------------

def test_1d_sigma_xx_zero_matches_default():
    """sigma_xx=0 gives identical results to not setting sigma_xx."""
    N = 150
    dx = 4000.0
    qs = np.zeros(N)
    qs[60:90] = 1e6
    flex_default = _run_flex_1d(Te=30e3, qs=qs, dx=dx,
                                bc_w="periodic", bc_e="periodic")
    flex_zero = _run_flex_1d(Te=30e3, qs=qs, dx=dx,
                             bc_w="periodic", bc_e="periodic",
                             sigma_xx=0.0)
    np.testing.assert_array_equal(flex_default.w, flex_zero.w)


def test_1d_end_load_fft_vs_fd():
    """FFT and FD solutions agree for sigma_xx with periodic BCs.

    The FFT method applies the exact spectral transfer function
    W(k) = -Q(k) / (D k⁴ + σ_xx T_e k² + Δρg), making it the analytical
    reference for the FD end-load implementation.  Agreement to within FD
    truncation error (~rtol 1e-4 for this grid) confirms sign, scaling, and
    placement of the σ_xx term in the FD stencil.
    """
    dx = 4000.0
    N = 256
    L = N * dx
    Te = 30e3
    sigma_xx = 2e8  # tensile [Pa]

    n_waves = 3
    k = 2.0 * np.pi * n_waves / L
    x = (np.arange(N) + 0.5) * dx
    qs = np.sin(k * x) * 1e6

    flex_fft = _run_flex_1d(Te=Te, qs=qs, dx=dx, method="fft",
                             bc_w="periodic", bc_e="periodic",
                             sigma_xx=sigma_xx)
    flex_fd = _run_flex_1d(Te=Te, qs=qs, dx=dx, method="fd",
                            bc_w="periodic", bc_e="periodic",
                            sigma_xx=sigma_xx)

    # rtol=5e-4: FD discretisation error in k⁴ is O((kdx)²) ≈ 0.05 % for
    # this grid; FFT is the exact spectral answer.
    np.testing.assert_allclose(flex_fd.w, flex_fft.w, rtol=5e-4)


def test_1d_end_load_monotonicity():
    """Tensile sigma_xx reduces deflection; compressive increases it."""
    N = 200
    dx = 4000.0
    qs = np.zeros(N)
    qs[90:110] = 1e6

    flex_0 = _run_flex_1d(Te=30e3, qs=qs, dx=dx)
    flex_t = _run_flex_1d(Te=30e3, qs=qs, dx=dx, sigma_xx=+1e8)   # tensile
    flex_c = _run_flex_1d(Te=30e3, qs=qs, dx=dx, sigma_xx=-1e8)   # compressive

    assert flex_t.w.min() > flex_0.w.min(), "tensile end load should reduce subsidence"
    assert flex_c.w.min() < flex_0.w.min(), "compressive end load should increase subsidence"


# ---------------------------------------------------------------------------
# Variable elastic thickness
# ---------------------------------------------------------------------------

def test_1d_fd_uniform_te_array_equals_scalar():
    """Uniform Te array gives the same result as scalar Te.

    BC_Rigidity broadcasts a scalar D to a uniform array before building the
    stencil, so a uniform Te array should produce bit-for-bit identical output.
    """
    N  = 100
    dx = 4000.0
    qs = np.zeros(N)
    qs[45:55] = 1e6

    flex_scalar = _run_flex_1d(Te=30e3, qs=qs, dx=dx)
    flex_array  = _run_flex_1d(Te=np.full(N, 30e3), qs=qs, dx=dx)

    np.testing.assert_allclose(flex_array.w, flex_scalar.w, rtol=1e-12)


def test_1d_fd_variable_te_monotonicity():
    """Higher Te everywhere reduces deflection magnitude.

    Comparing two uniform-array Te runs exercises the variable-Te code path
    (scalar D is broadcast to an array in BC_Rigidity) and confirms the
    stencil uses D correctly.
    """
    N  = 200
    dx = 4000.0
    qs = np.zeros(N)
    qs[90:110] = 1e6

    flex_lo = _run_flex_1d(Te=np.full(N, 20e3), qs=qs, dx=dx)
    flex_hi = _run_flex_1d(Te=np.full(N, 40e3), qs=qs, dx=dx)

    assert flex_hi.w.min() > flex_lo.w.min(), (
        "higher Te should produce less subsidence"
    )


def test_1d_fd_variable_te_asymmetric_deflection():
    """Step change in Te produces an asymmetric deflection profile.

    A load centred on a domain where the left half has low Te and the right
    half has high Te must deflect more deeply on the soft (left) side than
    on the stiff (right) side, and the profile must be asymmetric.
    """
    N  = 200
    dx = 4000.0
    qs = np.zeros(N)
    qs[95:105] = 1e6   # load at centre

    Te_arr = np.full(N, 20e3)
    Te_arr[100:] = 40e3   # right half stiffer

    flex = _run_flex_1d(Te=Te_arr, qs=qs, dx=dx)

    # Left half (soft) should deflect more than the mirror of the right half
    assert flex.w[:100].min() < flex.w[100:][::-1].min(), (
        "deflection should be deeper on the low-Te (left) side"
    )
    # Profile must differ from any symmetric (uniform Te) solution
    assert not np.allclose(flex.w[:100], flex.w[100:][::-1]), (
        "step-Te deflection must be asymmetric"
    )


# ---------------------------------------------------------------------------
# Convergence order (Method of Manufactured Solutions)
# ---------------------------------------------------------------------------

def test_1d_fd_convergence_order():
    """1-D FD solver achieves second-order (O(dx²)) spatial convergence.

    Uses the Method of Manufactured Solutions (MMS) with:

        w_exact(x) = −cos(2π x / L)

    on a periodic domain of width L = 2α (α = (4D/Δρg)^0.25).  The
    manufactured load is:

        q_mms(x) = (D kx⁴ + Δρg) · cos(kx x),   kx = 2π/L

    Three grids (N = 50, 100, 200 cells) must show a convergence rate > 1.8
    between each successive pair, confirming second-order behaviour with
    modest tolerance for pre-asymptotic effects.
    """
    E_mms  = 65e9
    nu_mms = 0.25
    Te_mms = 30e3
    rho_m_mms  = 3300.0
    rho_f_mms  = 0.0
    g_mms  = 9.8

    D_mms  = E_mms * Te_mms**3 / (12.0 * (1.0 - nu_mms**2))
    k_mms  = (rho_m_mms - rho_f_mms) * g_mms
    alpha_mms = (4.0 * D_mms / k_mms) ** 0.25
    L      = 2.0 * alpha_mms
    kx     = 2.0 * np.pi / L
    q_fac  = D_mms * kx**4 + k_mms

    Ns     = [50, 100, 200]
    dxs    = []
    errors = []

    for N in Ns:
        dxi    = L / N
        x      = (np.arange(N) + 0.5) * dxi
        qs_mms = q_fac * np.cos(kx * x)
        flex   = _run_flex_1d(Te=Te_mms, qs=qs_mms, dx=dxi,
                               bc_w="periodic", bc_e="periodic")
        err    = np.max(np.abs(flex.w - (-np.cos(kx * x))))
        dxs.append(dxi)
        errors.append(err)

    for i in range(len(Ns) - 1):
        rate = np.log(errors[i] / errors[i + 1]) / np.log(dxs[i] / dxs[i + 1])
        assert rate > 1.8, (
            f"Expected O(dx²) convergence (rate ≥ 1.8), "
            f"got {rate:.2f} between N={Ns[i]} and N={Ns[i+1]} "
            f"(errors: {errors[i]:.3g} → {errors[i+1]:.3g} m)"
        )



if __name__ == "__main__":
    pytest.main([__file__, "-v"])
