#! /usr/bin/env python
"""Compute flexural wavelength and related length scales for a thin elastic plate."""

import argparse

import numpy as np


def flexural_wavelengths(Te, rho_m, rho_fill, E, nu, g):
    """
    Compute flexural parameters and wavelengths for a thin elastic plate.

    Parameters
    ----------
    Te : float
        Elastic thickness [m].
    rho_m : float
        Mantle density [kg m^-3].
    rho_fill : float
        Infill density [kg m^-3] (e.g. 0 for air, 1000 for water).
    E : float
        Young's modulus [Pa].
    nu : float
        Poisson's ratio.
    g : float
        Gravitational acceleration [m s^-2].

    Returns
    -------
    dict with keys:
        alpha_1D, lambda_1D, zero_crossing_1D  (all in metres)
        alpha_2D, lambda_2D, zero_crossing_2D
    """
    drho = rho_m - rho_fill
    D = (E * Te**3) / (12.0 * (1.0 - nu**2))
    alpha_1D = (4.0 * D / (drho * g)) ** 0.25
    alpha_2D = (D / (drho * g)) ** 0.25
    lambda_1D = 2.0 * np.pi * alpha_1D
    lambda_2D = 2.0 * np.pi * alpha_2D
    return {
        "alpha_1D": alpha_1D,
        "lambda_1D": lambda_1D,
        "zero_crossing_1D": 0.375 * lambda_1D,
        "alpha_2D": alpha_2D,
        "lambda_2D": lambda_2D,
        "zero_crossing_2D": 0.375 * lambda_2D,
    }


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--Te", type=float, required=True,
                        metavar="m", help="Elastic thickness [m]")
    parser.add_argument("--rho-m", type=float, required=True,
                        metavar="kg/m3", help="Mantle density [kg m⁻³]")
    parser.add_argument("--rho-fill", type=float, required=True,
                        metavar="kg/m3", help="Infill density [kg m⁻³] (e.g. 0 for air, 1000 for water)")
    parser.add_argument("--E", type=float, required=True,
                        metavar="Pa", help="Young's modulus [Pa]")
    parser.add_argument("--nu", type=float, required=True,
                        help="Poisson's ratio")
    parser.add_argument("--g", type=float, required=True,
                        metavar="m/s2", help="Gravitational acceleration [m s⁻²]")
    args = parser.parse_args()

    r = flexural_wavelengths(
        Te=args.Te, rho_m=args.rho_m, rho_fill=args.rho_fill,
        E=args.E, nu=args.nu, g=args.g,
    )

    print()
    print("1D:")
    print("  Flexural wavelength:             ", r["lambda_1D"] / 1000, "km")
    print("  Distance to first zero-crossing: ", r["zero_crossing_1D"] / 1000, "km")
    print("  Flexural parameter:              ", r["alpha_1D"] / 1000, "km")
    print()
    print("2D:")
    print("  Flexural wavelength:             ", r["lambda_2D"] / 1000, "km")
    print("  Distance to first zero-crossing: ", r["zero_crossing_2D"] / 1000, "km")
    print("  Flexural parameter:              ", r["alpha_2D"] / 1000, "km")
    print()


if __name__ == "__main__":
    main()
