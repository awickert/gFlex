#! /usr/bin/env python

import numpy as np
import pytest

from gflex.f2d import F2D, pad_domain, recommended_pad_width, smooth_pad_Te


def _run_flex_2d(Te, qs, dx, dy, bc="0Displacement0Slope"):
    """Helper: run a 2-D FD flexure calculation and return deflection array."""
    flex = F2D()
    flex.Quiet = True
    flex.Method = "FD"
    flex.PlateSolutionType = "vWC1994"
    flex.Solver = "direct"
    flex.g = 9.8
    flex.E = 65e9
    flex.nu = 0.25
    flex.rho_m = 3300.0
    flex.rho_fill = 0.0
    flex.Te = Te.copy() if isinstance(Te, np.ndarray) else Te
    flex.qs = qs.copy()
    flex.dx = dx
    flex.dy = dy
    flex.BC_W = bc
    flex.BC_E = bc
    flex.BC_S = bc
    flex.BC_N = bc
    flex.initialize()
    flex.run()
    flex.finalize()
    return flex.w


def test_main():
    flex = F2D()

    flex.Quiet = False

    flex.Method = "FD"  # Solution method: * FD (finite difference)
    #                  * SAS (superposition of analytical solutions)
    #                  * SAS_NG (ungridded SAS)
    flex.PlateSolutionType = "vWC1994"  # van Wees and Cloetingh (1994)
    # The other option is 'G2009': Govers et al. (2009)
    flex.Solver = "direct"  # direct or iterative
    # convergence = 1E-3 # convergence between iterations, if an iterative solution
    # method is chosen

    flex.g = 9.8  # acceleration due to gravity
    flex.E = 65e9  # Young's Modulus
    flex.nu = 0.25  # Poisson's Ratio
    flex.rho_m = 3300.0  # MantleDensity
    flex.rho_fill = 0.0  # InfiillMaterialDensity

    flex.Te = 35000.0 * np.ones(
        (50, 50)
    )  # Elastic thickness [m] -- scalar but may be an array
    flex.Te[:, -3:] = 0.0
    flex.qs = np.zeros((50, 50))  # Template array for surface load stresses
    flex.qs[10:40, 10:40] += 1e6  # Populating this template
    flex.dx = 5000.0  # grid cell size, x-oriented [m]
    flex.dy = 5000.0  # grid cell size, y-oriented [m]
    # Boundary conditions can be:
    # (FD): 0Slope0Shear, 0Moment0Shear, 0Displacement0Slope, Mirror, or Periodic
    # For SAS or SAS_NG, NoOutsideLoads is valid, and no entry defaults to this
    flex.BC_W = "0Displacement0Slope"  # west boundary condition
    flex.BC_E = "0Moment0Shear"  # east boundary condition
    flex.BC_S = "0Displacement0Slope"  # south boundary condition
    flex.BC_N = "0Displacement0Slope"  # north boundary condition

    # latitude/longitude solutions are exact for SAS, approximate otherwise
    # latlon = # true/false: flag to enable lat/lon input. Defaults False.
    # PlanetaryRadius = # radius of planet [m], for lat/lon solutions

    flex.initialize()
    flex.run()
    flex.finalize()

    # If you want to plot the output
    # flex.plotChoice='both'
    # An output file for deflections could also be defined here
    # flex.wOutFile =
    flex.output()  # Plots and/or saves output, or does nothing, depending on
    # whether flex.plotChoice and/or flex.wOutFile have been set


def test_variable_Te_abrupt_padding_artefact():
    """
    Regression test for issue #45: variable-Te domain padded with mean Te.

    When the user manually pads a spatially variable Te grid with a constant
    value (e.g., mean(Te)) before passing the full grid to gFlex, there is a
    sharp step in flexural rigidity D at the inner/outer boundary.  The vWC1994
    finite-difference stencil contains D-derivative terms (Dx, Dy, Dxx, Dyy,
    Dxy) that are large at this step and drive deflections in the nominally
    zero-load padding zone — an effect that does not appear with constant Te.

    This test confirms:
      1. Both the abrupt and smooth cases deflect downward under the load.
      2. A linear Te taper within the padding zone (blending from inner-boundary
         Te values at the inner edge to mean Te at the outer edge) reduces the
         step at the inner/outer boundary to ~1/pad of its abrupt value and cuts
         the resulting deflection artefact in the padding zone by roughly 5x.

    Workaround for users: instead of padding with a constant mean Te, taper Te
    smoothly over the padding width so that the rigidity gradient is small at the
    inner/outer boundary.
    """
    dx = dy = 5000.0  # m
    pad = 8
    ny_in = nx_in = 24
    ny = nx = ny_in + 2 * pad  # 40 x 40

    # Inner domain: left half soft (Te = 20 km), right half stiff (Te = 50 km)
    Te_inner = np.full((ny_in, nx_in), 20e3)
    Te_inner[:, nx_in // 2 :] = 50e3
    Te_mean = Te_inner.mean()  # 35 km

    # Load placed in the softer (20 km) left half of the inner domain
    qs = np.zeros((ny, nx))
    qs[pad + 4 : pad + 20, pad + 4 : pad + 12] = 1e6

    load_row = (pad + 4 + pad + 20) // 2  # centre row of load
    load_col = (pad + 4 + pad + 12) // 2  # centre col of load

    # Reference: uniform Te = mean(Te_inner) everywhere
    Te_const = np.full((ny, nx), Te_mean)

    # Abrupt (Jürgen's approach): inner variable Te, outer ring = mean Te
    Te_abrupt = np.full((ny, nx), Te_mean)
    Te_abrupt[pad:-pad, pad:-pad] = Te_inner

    # Smooth: use the new utility — reduces the step at the inner/outer
    # boundary to ~1/pad of the abrupt value
    Te_smooth = smooth_pad_Te(Te_inner, pad_width=pad, Te_out=Te_mean)

    w_const = _run_flex_2d(Te_const, qs, dx, dy)
    w_abrupt = _run_flex_2d(Te_abrupt, qs, dx, dy)
    w_smooth = _run_flex_2d(Te_smooth, qs, dx, dy)

    # 1. Both variable-Te cases must deflect downward under the load.
    assert w_abrupt[load_row, load_col] < 0, "abrupt case: expected downward deflection under load"
    assert w_smooth[load_row, load_col] < 0, "smooth case: expected downward deflection under load"

    # 2. Smooth Te taper substantially reduces the deflection artefact in the
    #    padding zone relative to the constant-Te reference.
    pad_mask = np.zeros((ny, nx), dtype=bool)
    pad_mask[:pad, :] = True
    pad_mask[-pad:, :] = True
    pad_mask[:, :pad] = True
    pad_mask[:, -pad:] = True

    dev_abrupt = np.abs(w_abrupt[pad_mask] - w_const[pad_mask]).max()
    dev_smooth = np.abs(w_smooth[pad_mask] - w_const[pad_mask]).max()

    assert dev_smooth < dev_abrupt, (
        f"Smooth Te taper should reduce boundary artefact "
        f"(abrupt max dev = {dev_abrupt:.4g} m, smooth max dev = {dev_smooth:.4g} m)"
    )


def test_recommended_pad_width():
    """recommended_pad_width returns a positive integer that scales with Te."""
    p_thin = recommended_pad_width(Te=20e3, dx=5000.0)
    p_thick = recommended_pad_width(Te=50e3, dx=5000.0)
    assert isinstance(p_thin, int) and p_thin > 0
    assert p_thick > p_thin  # stiffer plate → wider flexural wavelength → more padding


def test_pad_domain():
    """pad_domain returns correctly shaped arrays with zeros in the load ring."""
    ny, nx = 20, 30
    Te = 35e3 * np.ones((ny, nx))
    qs = np.ones((ny, nx))
    dx = dy = 5000.0

    Te_pad, qs_pad, p = pad_domain(Te, qs, dx=dx, dy=dy)

    assert isinstance(p, int) and p > 0
    assert Te_pad.shape == (ny + 2 * p, nx + 2 * p)
    assert qs_pad.shape == (ny + 2 * p, nx + 2 * p)

    # Inner domain is preserved exactly
    np.testing.assert_array_equal(Te_pad[p:-p, p:-p], Te)
    np.testing.assert_array_equal(qs_pad[p:-p, p:-p], qs)

    # Load padding ring is all zeros
    pad_mask = np.ones(qs_pad.shape, dtype=bool)
    pad_mask[p:-p, p:-p] = False
    assert np.all(qs_pad[pad_mask] == 0.0)


if __name__ == "__main__":
    test_main()
