#! /usr/bin/env python
"""Tests for FD boundary-condition warnings in F1D and F2D.

Warning types covered:
  0Moment0Shear   — fires for any side carrying that BC
  0Slope0Shear    — fires for any side carrying that BC
  proximity       — fires when nearest loaded cell is within one flexural
                    wavelength of a 0Displacement0Slope boundary; absent for
                    Mirror / Periodic, and when the load is far enough away
"""

import warnings

import numpy as np
import pytest

from gflex.f1d import F1D
from gflex.f2d import F2D

# ---------------------------------------------------------------------------
# Shared physical parameters
# ---------------------------------------------------------------------------

E = 65e9
nu = 0.25
rho_m = 3300.0
rho_fill = 0.0
g = 9.8

# Te = 35 km → alpha_1D ≈ 74.4 km → flexural wavelength ≈ 467 km
Te = 35e3
dx = 5000.0

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _run_1d(qs, bc_w, bc_e):
    flex = F1D()
    flex.Quiet = True
    flex.Method = "FD"
    flex.Solver = "direct"
    flex.g = g
    flex.E = E
    flex.nu = nu
    flex.rho_m = rho_m
    flex.rho_fill = rho_fill
    flex.Te = Te
    flex.qs = qs.copy()
    flex.dx = dx
    flex.BC_W = bc_w
    flex.BC_E = bc_e
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        flex.initialize()
        flex.run()
        flex.finalize()
    return [str(x.message) for x in w if issubclass(x.category, UserWarning)]


def _run_2d(qs, bc_w, bc_e, bc_n, bc_s):
    flex = F2D()
    flex.Quiet = True
    flex.Method = "FD"
    flex.PlateSolutionType = "vWC1994"
    flex.Solver = "direct"
    flex.g = g
    flex.E = E
    flex.nu = nu
    flex.rho_m = rho_m
    flex.rho_fill = rho_fill
    flex.Te = Te
    flex.qs = qs.copy()
    flex.dx = dx
    flex.dy = dx
    flex.BC_W = bc_w
    flex.BC_E = bc_e
    flex.BC_N = bc_n
    flex.BC_S = bc_s
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        flex.initialize()
        flex.run()
        flex.finalize()
    return [str(x.message) for x in w if issubclass(x.category, UserWarning)]


# ---------------------------------------------------------------------------
# 1-D: BC-type warnings
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("side", ["W", "E"])
def test_1d_moment_shear_warning_fires(side):
    """0Moment0Shear triggers a warning on the named side."""
    qs = np.zeros(80)
    qs[40] = 1e6
    bcs = {"W": "Mirror", "E": "Mirror"}
    bcs[side] = "0Moment0Shear"
    msgs = _run_1d(qs, bcs["W"], bcs["E"])
    assert any(f"BC_{side} = '0Moment0Shear'" in m for m in msgs)


@pytest.mark.parametrize("side", ["W", "E"])
def test_1d_slope_shear_warning_fires(side):
    """0Slope0Shear triggers a warning on the named side."""
    qs = np.zeros(80)
    qs[40] = 1e6
    bcs = {"W": "Mirror", "E": "Mirror"}
    bcs[side] = "0Slope0Shear"
    msgs = _run_1d(qs, bcs["W"], bcs["E"])
    assert any(f"BC_{side} = '0Slope0Shear'" in m for m in msgs)


@pytest.mark.parametrize("bc", ["Mirror", "Periodic", "0Displacement0Slope"])
def test_1d_no_bc_type_warning(bc):
    """Mirror, Periodic, and 0Displacement0Slope do not trigger BC-type warnings."""
    qs = np.zeros(300)
    qs[150] = 1e6
    msgs = _run_1d(qs, bc, bc)
    assert not any("0Moment0Shear" in m or "0Slope0Shear" in m for m in msgs)


# ---------------------------------------------------------------------------
# 1-D: proximity warnings
# ---------------------------------------------------------------------------

def test_1d_proximity_warning_fires():
    """Load within one wavelength of a 0Displacement0Slope boundary warns."""
    # cell 2 → 12.5 km from W boundary; wavelength ≈ 467 km at Te=35 km
    qs = np.zeros(40)
    qs[2] = 1e6
    msgs = _run_1d(qs, "0Displacement0Slope", "0Displacement0Slope")
    assert any("flexural wavelength" in m for m in msgs)


def test_1d_proximity_warning_absent_when_far():
    """Load more than one wavelength from both boundaries: no proximity warning."""
    # N=300 cells, load at centre → 752.5 km from W, 747.5 km from E >> 467 km
    qs = np.zeros(300)
    qs[150] = 1e6
    msgs = _run_1d(qs, "0Displacement0Slope", "0Displacement0Slope")
    assert not any("flexural wavelength" in m for m in msgs)


@pytest.mark.parametrize("bc", ["Mirror", "Periodic"])
def test_1d_proximity_warning_not_for_other_bcs(bc):
    """Proximity warning only fires for 0Displacement0Slope, not Mirror or Periodic."""
    qs = np.zeros(40)
    qs[2] = 1e6
    msgs = _run_1d(qs, bc, bc)
    assert not any("flexural wavelength" in m for m in msgs)


def test_1d_proximity_no_warning_empty_load():
    """All-zero load array: no proximity warning."""
    qs = np.zeros(40)
    msgs = _run_1d(qs, "0Displacement0Slope", "0Displacement0Slope")
    assert not any("flexural wavelength" in m for m in msgs)


def test_1d_proximity_warning_message_content():
    """Proximity warning message contains wavelength fraction and pad_domain_1d hint."""
    qs = np.zeros(40)
    qs[2] = 1e6
    msgs = _run_1d(qs, "0Displacement0Slope", "0Displacement0Slope")
    prox = [m for m in msgs if "flexural wavelength" in m]
    assert prox, "expected at least one proximity warning"
    msg = prox[0]
    assert "flexural wavelengths" in msg
    assert "pad_domain_1d()" in msg


# ---------------------------------------------------------------------------
# 2-D: BC-type warnings
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("side", ["W", "E", "N", "S"])
def test_2d_moment_shear_warning_fires(side):
    """0Moment0Shear triggers a warning on the named side."""
    qs = np.zeros((80, 80))
    qs[40, 40] = 1e6
    bcs = {"W": "Mirror", "E": "Mirror", "N": "Mirror", "S": "Mirror"}
    bcs[side] = "0Moment0Shear"
    msgs = _run_2d(qs, bcs["W"], bcs["E"], bcs["N"], bcs["S"])
    assert any(f"BC_{side} = '0Moment0Shear'" in m for m in msgs)


@pytest.mark.parametrize("side", ["W", "E", "N", "S"])
def test_2d_slope_shear_warning_fires(side):
    """0Slope0Shear triggers a warning on the named side."""
    qs = np.zeros((80, 80))
    qs[40, 40] = 1e6
    bcs = {"W": "Mirror", "E": "Mirror", "N": "Mirror", "S": "Mirror"}
    bcs[side] = "0Slope0Shear"
    msgs = _run_2d(qs, bcs["W"], bcs["E"], bcs["N"], bcs["S"])
    assert any(f"BC_{side} = '0Slope0Shear'" in m for m in msgs)


@pytest.mark.parametrize("bc", ["Mirror", "Periodic", "0Displacement0Slope"])
def test_2d_no_bc_type_warning(bc):
    """Mirror, Periodic, and 0Displacement0Slope do not trigger BC-type warnings."""
    qs = np.zeros((300, 300))
    qs[150, 150] = 1e6
    msgs = _run_2d(qs, bc, bc, bc, bc)
    assert not any("0Moment0Shear" in m or "0Slope0Shear" in m for m in msgs)


# ---------------------------------------------------------------------------
# 2-D: proximity warnings
# ---------------------------------------------------------------------------

def test_2d_proximity_warning_fires():
    """Load within one wavelength of a 0Displacement0Slope boundary warns."""
    # row 2 → 12.5 km from N boundary; wavelength ≈ 467 km at Te=35 km
    qs = np.zeros((40, 40))
    qs[2, 20] = 1e6
    msgs = _run_2d(qs, "0Displacement0Slope", "0Displacement0Slope",
                   "0Displacement0Slope", "0Displacement0Slope")
    assert any("flexural wavelength" in m for m in msgs)


def test_2d_proximity_warning_absent_when_far():
    """Load more than one wavelength from all boundaries: no proximity warning."""
    # N=300, load at centre → 752.5 km from each boundary >> 467 km
    qs = np.zeros((300, 300))
    qs[150, 150] = 1e6
    msgs = _run_2d(qs, "0Displacement0Slope", "0Displacement0Slope",
                   "0Displacement0Slope", "0Displacement0Slope")
    assert not any("flexural wavelength" in m for m in msgs)


@pytest.mark.parametrize("bc", ["Mirror", "Periodic"])
def test_2d_proximity_warning_not_for_other_bcs(bc):
    """Proximity warning only fires for 0Displacement0Slope, not Mirror or Periodic."""
    qs = np.zeros((40, 40))
    qs[2, 20] = 1e6
    msgs = _run_2d(qs, bc, bc, bc, bc)
    assert not any("flexural wavelength" in m for m in msgs)


def test_2d_proximity_no_warning_empty_load():
    """All-zero load array: no proximity warning."""
    qs = np.zeros((40, 40))
    msgs = _run_2d(qs, "0Displacement0Slope", "0Displacement0Slope",
                   "0Displacement0Slope", "0Displacement0Slope")
    assert not any("flexural wavelength" in m for m in msgs)


def test_2d_proximity_warning_message_content():
    """Proximity warning message contains wavelength fraction and pad_domain() hint."""
    qs = np.zeros((40, 40))
    qs[2, 20] = 1e6
    msgs = _run_2d(qs, "0Displacement0Slope", "0Displacement0Slope",
                   "0Displacement0Slope", "0Displacement0Slope")
    prox = [m for m in msgs if "flexural wavelength" in m]
    assert prox, "expected at least one proximity warning"
    msg = prox[0]
    assert "flexural wavelengths" in msg
    assert "pad_domain()" in msg
