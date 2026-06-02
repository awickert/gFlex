#! /usr/bin/env python
"""Tests for 1-D domain-padding utilities.

Unit tests
----------
recommended_pad_width_1d  — positive integer, scales with Te and n_wavelengths
smooth_pad_Te_1d          — shape, inner domain preserved, outer edge, taper
                            monotonicity, error conditions
pad_domain_1d             — shapes, inner domain, zero-padded load ring

Integration test
----------------
Padded FD agrees better with the SAS (infinite-plate) reference than unpadded
FD when the load is within one flexural wavelength of a zero_displacement_zero_slope
boundary.
"""

import warnings

import numpy as np
import pytest

from gflex.f1d import F1D, pad_domain_1d, recommended_pad_width_1d, smooth_pad_Te_1d

# ---------------------------------------------------------------------------
# Shared physical parameters
# Te = 35 km → alpha_1D ≈ 74.4 km → flexural wavelength ≈ 467 km
# ---------------------------------------------------------------------------

E = 65e9
nu = 0.25
rho_m = 3300.0
rho_fill = 0.0
g = 9.8
Te = 35e3
dx = 5000.0

# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _run_1d(Te_in, qs, bc="zero_displacement_zero_slope", method="fd"):
    flex = F1D()
    flex.quiet = True
    flex.method = method
    flex.solver = "direct"
    flex.g = g
    flex.E = E
    flex.nu = nu
    flex.rho_m = rho_m
    flex.rho_fill = rho_fill
    flex.te = Te_in.copy() if isinstance(Te_in, np.ndarray) else Te_in
    flex.qs = qs.copy()
    flex.dx = dx
    flex.bc_west = bc
    flex.bc_east = bc
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        flex.initialize()
        flex.run()
        w = flex.w
        flex.finalize()
    return w


# ---------------------------------------------------------------------------
# recommended_pad_width_1d
# ---------------------------------------------------------------------------

def test_recommended_pad_width_1d_type_and_sign():
    p = recommended_pad_width_1d(Te=35e3, dx=5000.0)
    assert isinstance(p, int) and p > 0


def test_recommended_pad_width_1d_scales_with_Te():
    """Stiffer plate → longer wavelength → wider padding."""
    p_thin = recommended_pad_width_1d(Te=20e3, dx=5000.0)
    p_thick = recommended_pad_width_1d(Te=50e3, dx=5000.0)
    assert p_thick > p_thin


def test_recommended_pad_width_1d_scales_with_n_wavelengths():
    p1 = recommended_pad_width_1d(Te=35e3, dx=5000.0, n_wavelengths=1.0)
    p2 = recommended_pad_width_1d(Te=35e3, dx=5000.0, n_wavelengths=2.0)
    assert p2 > p1


def test_recommended_pad_width_1d_decreases_with_dx():
    """Finer grid → same wavelength requires more cells."""
    p_coarse = recommended_pad_width_1d(Te=35e3, dx=10000.0)
    p_fine = recommended_pad_width_1d(Te=35e3, dx=5000.0)
    assert p_fine > p_coarse


# ---------------------------------------------------------------------------
# smooth_pad_Te_1d
# ---------------------------------------------------------------------------

def test_smooth_pad_Te_1d_shape():
    Te_arr = Te * np.ones(20)
    p = 5
    Te_pad = smooth_pad_Te_1d(Te_arr, pad_width=p)
    assert Te_pad.shape == (20 + 2 * p,)


def test_smooth_pad_Te_1d_inner_preserved():
    Te_arr = np.linspace(20e3, 50e3, 20)
    p = 6
    Te_pad = smooth_pad_Te_1d(Te_arr, pad_width=p)
    np.testing.assert_array_equal(Te_pad[p:-p], Te_arr)


def test_smooth_pad_Te_1d_outer_edge_equals_Te_out():
    Te_arr = Te * np.ones(20)
    Te_out = 20e3
    p = 5
    Te_pad = smooth_pad_Te_1d(Te_arr, pad_width=p, Te_out=Te_out)
    assert Te_pad[0] == Te_out
    assert Te_pad[-1] == Te_out


def test_smooth_pad_Te_1d_default_Te_out_is_mean():
    Te_arr = np.array([20e3, 30e3, 40e3, 50e3])
    Te_pad = smooth_pad_Te_1d(Te_arr, pad_width=3)
    assert Te_pad[0] == pytest.approx(Te_arr.mean())


def test_smooth_pad_Te_1d_left_taper_monotonic():
    """Left padding increases monotonically from outer edge (Te_out) to inner domain."""
    Te_arr = Te * np.ones(20)
    Te_out = 20e3  # less than Te_arr, so taper increases inward
    p = 8
    Te_pad = smooth_pad_Te_1d(Te_arr, pad_width=p, Te_out=Te_out)
    assert np.all(np.diff(Te_pad[:p]) >= 0)


def test_smooth_pad_Te_1d_right_taper_monotonic():
    """Right padding decreases monotonically from inner domain to outer edge (Te_out)."""
    Te_arr = Te * np.ones(20)
    Te_out = 20e3
    p = 8
    Te_pad = smooth_pad_Te_1d(Te_arr, pad_width=p, Te_out=Te_out)
    assert np.all(np.diff(Te_pad[-p:]) <= 0)


def test_smooth_pad_Te_1d_error_wrong_ndim():
    with pytest.raises(ValueError, match="1-D"):
        smooth_pad_Te_1d(np.ones((5, 5)), pad_width=2)


def test_smooth_pad_Te_1d_error_pad_width_zero():
    with pytest.raises(ValueError, match="pad_width"):
        smooth_pad_Te_1d(np.ones(10), pad_width=0)


# ---------------------------------------------------------------------------
# pad_domain_1d
# ---------------------------------------------------------------------------

def test_pad_domain_1d_shapes():
    nx = 20
    Te_arr = Te * np.ones(nx)
    qs = np.ones(nx)
    Te_pad, qs_pad, p = pad_domain_1d(Te_arr, qs, dx=dx)
    assert isinstance(p, int) and p > 0
    assert Te_pad.shape == (nx + 2 * p,)
    assert qs_pad.shape == (nx + 2 * p,)


def test_pad_domain_1d_inner_preserved():
    nx = 20
    Te_arr = Te * np.ones(nx)
    qs = np.ones(nx)
    Te_pad, qs_pad, p = pad_domain_1d(Te_arr, qs, dx=dx)
    np.testing.assert_array_equal(Te_pad[p:-p], Te_arr)
    np.testing.assert_array_equal(qs_pad[p:-p], qs)


def test_pad_domain_1d_qs_padding_zeros():
    nx = 20
    Te_arr = Te * np.ones(nx)
    qs = np.ones(nx)
    Te_pad, qs_pad, p = pad_domain_1d(Te_arr, qs, dx=dx)
    assert np.all(qs_pad[:p] == 0.0)
    assert np.all(qs_pad[-p:] == 0.0)


# ---------------------------------------------------------------------------
# Integration: padded FD agrees better with SAS than unpadded FD
# ---------------------------------------------------------------------------

def test_pad_domain_1d_improves_fd_accuracy():
    """Padded FD solution agrees better with SAS (infinite-plate) than unpadded FD.

    Load at cell 2 → 12.5 km from the W boundary, well within one flexural
    wavelength (≈ 467 km at Te = 35 km).  The zero_displacement_zero_slope BC suppresses
    the flexural forebulge and distorts the deflection.  pad_domain_1d pushes
    the effective boundary ≈ one wavelength away, recovering the infinite-plate
    response.
    """
    N = 40
    qs = np.zeros(N)
    qs[2] = 1e6
    Te_arr = np.full(N, Te)

    # SAS: infinite-plate reference (BCs irrelevant for SAS)
    w_sas = _run_1d(Te, qs, bc="zero_moment_zero_shear", method="sas")

    # Unpadded FD: boundary close to load
    w_unpadded = _run_1d(Te, qs, bc="zero_displacement_zero_slope", method="fd")

    # Padded FD: boundary pushed away; trim to inner domain
    Te_pad, qs_pad, p = pad_domain_1d(Te_arr, qs, dx=dx)
    w_padded = _run_1d(Te_pad, qs_pad, bc="zero_displacement_zero_slope", method="fd")[p:-p]

    # Compare over the western half where the W boundary effect is strongest
    inner = slice(0, N // 2)
    err_unpadded = np.max(np.abs(w_unpadded[inner] - w_sas[inner]))
    err_padded = np.max(np.abs(w_padded[inner] - w_sas[inner]))

    assert err_padded < err_unpadded, (
        f"Padded FD should agree better with SAS in the inner domain: "
        f"unpadded max err = {err_unpadded:.4g} m, "
        f"padded max err = {err_padded:.4g} m"
    )
