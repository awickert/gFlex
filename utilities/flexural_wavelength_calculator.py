#! /usr/bin/env python
"""Compute flexural wavelength and related length scales for a thin elastic plate."""

import argparse

import numpy as np


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--Te", type=float, default=30000.0,
                        metavar="m", help="Elastic thickness [m]")
    parser.add_argument("--rho-m", type=float, default=3300.0,
                        metavar="kg/m3", help="Mantle density [kg m⁻³]")
    parser.add_argument("--rho-fill", type=float, default=1000.0,
                        metavar="kg/m3", help="Infill density [kg m⁻³] (e.g. 1000 for water, 0 for air)")
    parser.add_argument("--E", type=float, default=65e9,
                        metavar="Pa", help="Young's modulus [Pa]")
    parser.add_argument("--nu", type=float, default=0.25,
                        help="Poisson's ratio")
    parser.add_argument("--g", type=float, default=9.8,
                        metavar="m/s2", help="Gravitational acceleration [m s⁻²]")
    args = parser.parse_args()

    drho = args.rho_m - args.rho_fill
    D = (args.E * args.Te**3) / (12 * (1 - args.nu**2))

    alpha1D = (4 * D / (drho * args.g)) ** 0.25
    alpha2D = (D / (drho * args.g)) ** 0.25

    lambda1D = alpha1D * 2 * np.pi
    lambda2D = alpha2D * 2 * np.pi

    print()
    print("1D:")
    print("  Flexural wavelength:              ", lambda1D / 1000, "km")
    print("  Distance to first zero-crossing:  ", 0.375 * lambda1D / 1000, "km")
    print("  Flexural parameter:               ", alpha1D / 1000, "km")
    print()
    print("2D:")
    print("  Flexural wavelength:              ", lambda2D / 1000, "km")
    print("  Distance to first zero-crossing:  ", 0.375 * lambda2D / 1000, "km")
    print("  Flexural parameter:               ", alpha2D / 1000, "km")
    print()


if __name__ == "__main__":
    main()
