#! /usr/bin/env python

import numpy as np
import pytest

from utilities.flexural_wavelength_calculator import flexural_wavelengths


# Standard geophysical parameters used throughout
_PARAMS = dict(rho_m=3300.0, rho_fill=0.0, E=65e9, nu=0.25, g=9.81)


def test_lambda1D_is_sqrt2_times_lambda2D():
    """alpha_1D = (4D/drho/g)^0.25, alpha_2D = (D/drho/g)^0.25, so lambda_1D/lambda_2D = 4^0.25 = sqrt(2)."""
    r = flexural_wavelengths(Te=30e3, **_PARAMS)
    np.testing.assert_allclose(r["lambda_1D"] / r["lambda_2D"], np.sqrt(2), rtol=1e-10)


def test_zero_crossing_is_0375_wavelength():
    r = flexural_wavelengths(Te=30e3, **_PARAMS)
    np.testing.assert_allclose(r["zero_crossing_1D"], 0.375 * r["lambda_1D"], rtol=1e-10)
    np.testing.assert_allclose(r["zero_crossing_2D"], 0.375 * r["lambda_2D"], rtol=1e-10)


def test_wavelength_increases_with_Te():
    """Stiffer plate → larger D → longer flexural wavelength."""
    r_thin = flexural_wavelengths(Te=20e3, **_PARAMS)
    r_thick = flexural_wavelengths(Te=50e3, **_PARAMS)
    assert r_thick["lambda_2D"] > r_thin["lambda_2D"]
    assert r_thick["lambda_1D"] > r_thin["lambda_1D"]


def test_wavelength_decreases_with_drho():
    """Larger density contrast → shorter flexural wavelength."""
    r_low = flexural_wavelengths(Te=30e3, rho_m=3300.0, rho_fill=1000.0, E=65e9, nu=0.25, g=9.81)
    r_high = flexural_wavelengths(Te=30e3, rho_m=3300.0, rho_fill=0.0,   E=65e9, nu=0.25, g=9.81)
    assert r_high["lambda_2D"] < r_low["lambda_2D"]


def test_known_value_2D_alpha():
    """Check alpha_2D against a hand-computed reference value."""
    Te = 30e3
    E, nu, rho_m, rho_fill, g = 65e9, 0.25, 3300.0, 0.0, 9.81
    D = E * Te**3 / (12 * (1 - nu**2))
    alpha_expected = (D / ((rho_m - rho_fill) * g)) ** 0.25
    r = flexural_wavelengths(Te=Te, rho_m=rho_m, rho_fill=rho_fill, E=E, nu=nu, g=g)
    np.testing.assert_allclose(r["alpha_2D"], alpha_expected, rtol=1e-12)
