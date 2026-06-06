#! /usr/bin/env python
"""Tests for FD boundary-condition warnings in F1D and F2D.

Warning types covered:
  zero_moment_zero_shear   — UserWarning fires for any side carrying that BC
  proximity       — fires when nearest loaded cell is within one flexural
                    wavelength of a zero_displacement_zero_slope boundary; absent for
                    zero_slope_zero_shear / mirror / periodic, and when the load is far enough away
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
    flex.quiet = True
    flex.method = "fd"
    flex.solver = "direct"
    flex.g = g
    flex.E = E
    flex.nu = nu
    flex.rho_m = rho_m
    flex.rho_fill = rho_fill
    flex.T_e = Te
    flex.qs = qs.copy()
    flex.dx = dx
    flex.bc_west = bc_w
    flex.bc_east = bc_e
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        flex.initialize()
        flex.run()
        flex.finalize()
    return [str(x.message) for x in w if issubclass(x.category, UserWarning)]


def _run_2d(qs, bc_w, bc_e, bc_n, bc_s):
    flex = F2D()
    flex.quiet = True
    flex.method = "fd"
    flex.solver = "direct"
    flex.g = g
    flex.E = E
    flex.nu = nu
    flex.rho_m = rho_m
    flex.rho_fill = rho_fill
    flex.T_e = Te
    flex.qs = qs.copy()
    flex.dx = dx
    flex.dy = dx
    flex.bc_west = bc_w
    flex.bc_east = bc_e
    flex.bc_north = bc_n
    flex.bc_south = bc_s
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
    """zero_moment_zero_shear triggers a warning on the named side."""
    qs = np.zeros(80)
    qs[40] = 1e6
    bcs = {"W": "mirror", "E": "mirror"}
    bcs[side] = "zero_moment_zero_shear"
    msgs = _run_1d(qs, bcs["W"], bcs["E"])
    assert any(f"BC_{side} = 'zero_moment_zero_shear'" in m for m in msgs)


@pytest.mark.parametrize("bc", ["mirror", "zero_slope_zero_shear", "periodic", "zero_displacement_zero_slope"])
def test_1d_no_bc_type_warning(bc):
    """mirror, zero_slope_zero_shear, periodic, and zero_displacement_zero_slope do not trigger BC-type warnings."""
    qs = np.zeros(300)
    qs[150] = 1e6
    msgs = _run_1d(qs, bc, bc)
    assert not any("zero_moment_zero_shear" in m or "zero_slope_zero_shear" in m for m in msgs)


# ---------------------------------------------------------------------------
# 1-D: proximity warnings
# ---------------------------------------------------------------------------

def test_1d_proximity_warning_fires():
    """Load within one wavelength of a zero_displacement_zero_slope boundary warns."""
    # cell 2 → 12.5 km from W boundary; wavelength ≈ 467 km at Te=35 km
    qs = np.zeros(40)
    qs[2] = 1e6
    msgs = _run_1d(qs, "zero_displacement_zero_slope", "zero_displacement_zero_slope")
    assert any("flexural wavelength" in m for m in msgs)


def test_1d_proximity_warning_absent_when_far():
    """Load more than one wavelength from both boundaries: no proximity warning."""
    # N=300 cells, load at centre → 752.5 km from W, 747.5 km from E >> 467 km
    qs = np.zeros(300)
    qs[150] = 1e6
    msgs = _run_1d(qs, "zero_displacement_zero_slope", "zero_displacement_zero_slope")
    assert not any("flexural wavelength" in m for m in msgs)


@pytest.mark.parametrize("bc", ["mirror", "periodic"])
def test_1d_proximity_warning_not_for_other_bcs(bc):
    """Proximity warning only fires for zero_displacement_zero_slope, not mirror or periodic."""
    qs = np.zeros(40)
    qs[2] = 1e6
    msgs = _run_1d(qs, bc, bc)
    assert not any("flexural wavelength" in m for m in msgs)


def test_1d_proximity_no_warning_empty_load():
    """All-zero load array: no proximity warning."""
    qs = np.zeros(40)
    msgs = _run_1d(qs, "zero_displacement_zero_slope", "zero_displacement_zero_slope")
    assert not any("flexural wavelength" in m for m in msgs)


def test_1d_proximity_warning_message_content():
    """Proximity warning message contains wavelength fraction and pad_domain() hint."""
    qs = np.zeros(40)
    qs[2] = 1e6
    msgs = _run_1d(qs, "zero_displacement_zero_slope", "zero_displacement_zero_slope")
    prox = [m for m in msgs if "flexural wavelength" in m]
    assert prox, "expected at least one proximity warning"
    msg = prox[0]
    assert "flexural wavelengths" in msg
    assert "pad_domain()" in msg


# ---------------------------------------------------------------------------
# 2-D: BC-type warnings
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("side", ["W", "E", "N", "S"])
def test_2d_moment_shear_warning_fires(side):
    """zero_moment_zero_shear triggers a warning on the named side."""
    qs = np.zeros((80, 80))
    qs[40, 40] = 1e6
    bcs = {"W": "mirror", "E": "mirror", "N": "mirror", "S": "mirror"}
    bcs[side] = "zero_moment_zero_shear"
    msgs = _run_2d(qs, bcs["W"], bcs["E"], bcs["N"], bcs["S"])
    assert any(f"BC_{side} = 'zero_moment_zero_shear'" in m for m in msgs)


@pytest.mark.parametrize("bc", ["mirror", "zero_slope_zero_shear", "periodic", "zero_displacement_zero_slope"])
def test_2d_no_bc_type_warning(bc):
    """mirror, zero_slope_zero_shear, periodic, and zero_displacement_zero_slope do not trigger BC-type warnings."""
    qs = np.zeros((300, 300))
    qs[150, 150] = 1e6
    msgs = _run_2d(qs, bc, bc, bc, bc)
    assert not any("zero_moment_zero_shear" in m or "zero_slope_zero_shear" in m for m in msgs)


# ---------------------------------------------------------------------------
# 2-D: proximity warnings
# ---------------------------------------------------------------------------

def test_2d_proximity_warning_fires():
    """Load within one wavelength of a zero_displacement_zero_slope boundary warns."""
    # row 2 → 12.5 km from N boundary; wavelength ≈ 467 km at Te=35 km
    qs = np.zeros((40, 40))
    qs[2, 20] = 1e6
    msgs = _run_2d(qs, "zero_displacement_zero_slope", "zero_displacement_zero_slope",
                   "zero_displacement_zero_slope", "zero_displacement_zero_slope")
    assert any("flexural wavelength" in m for m in msgs)


def test_2d_proximity_warning_absent_when_far():
    """Load more than one wavelength from all boundaries: no proximity warning."""
    # N=300, load at centre → 752.5 km from each boundary >> 467 km
    qs = np.zeros((300, 300))
    qs[150, 150] = 1e6
    msgs = _run_2d(qs, "zero_displacement_zero_slope", "zero_displacement_zero_slope",
                   "zero_displacement_zero_slope", "zero_displacement_zero_slope")
    assert not any("flexural wavelength" in m for m in msgs)


@pytest.mark.parametrize("bc", ["mirror", "periodic"])
def test_2d_proximity_warning_not_for_other_bcs(bc):
    """Proximity warning only fires for zero_displacement_zero_slope, not mirror or periodic."""
    qs = np.zeros((40, 40))
    qs[2, 20] = 1e6
    msgs = _run_2d(qs, bc, bc, bc, bc)
    assert not any("flexural wavelength" in m for m in msgs)


def test_2d_proximity_no_warning_empty_load():
    """All-zero load array: no proximity warning."""
    qs = np.zeros((40, 40))
    msgs = _run_2d(qs, "zero_displacement_zero_slope", "zero_displacement_zero_slope",
                   "zero_displacement_zero_slope", "zero_displacement_zero_slope")
    assert not any("flexural wavelength" in m for m in msgs)


def test_2d_proximity_warning_message_content():
    """Proximity warning message contains wavelength fraction and pad_domain() hint."""
    qs = np.zeros((40, 40))
    qs[2, 20] = 1e6
    msgs = _run_2d(qs, "zero_displacement_zero_slope", "zero_displacement_zero_slope",
                   "zero_displacement_zero_slope", "zero_displacement_zero_slope")
    prox = [m for m in msgs if "flexural wavelength" in m]
    assert prox, "expected at least one proximity warning"
    msg = prox[0]
    assert "flexural wavelengths" in msg
    assert "pad_domain()" in msg


# ---------------------------------------------------------------------------
# Dict-style BC rejection on non-FD solvers (#66)
# ---------------------------------------------------------------------------

def _make_2d_base():
    """Return a minimally configured F2D with shared physical parameters."""
    flex = F2D()
    flex.E        = E
    flex.nu       = nu
    flex.rho_m    = rho_m
    flex.rho_fill = rho_fill
    flex.g        = g
    flex.T_e      = 35e3
    flex.dx = flex.dy = 10e3
    flex.qs = np.zeros((20, 20))
    return flex


@pytest.mark.parametrize("method", ["sas", "sas_ng", "fft"])
def test_dict_bc_rejected_on_non_fd_solver(method):
    """Dict-style BCs on non-FD solvers raise ValueError before any solve."""
    flex = _make_2d_base()
    flex.method   = method
    flex.bc_west  = {"moment": 1e12, "shear": 0.0}
    flex.bc_east  = flex.bc_north = flex.bc_south = ""
    flex.dimension = 2
    with pytest.raises(ValueError, match="method='fd'"):
        flex.bc_check()


@pytest.mark.parametrize("method", ["sas", "sas_ng", "fft"])
def test_string_bc_not_rejected_on_non_fd_solver(method):
    """String BCs on non-FD solvers do not trigger the dict guard."""
    flex = _make_2d_base()
    flex.method   = method
    flex.bc_west  = flex.bc_east = flex.bc_north = flex.bc_south = ""
    flex.dimension = 2
    # Should not raise — the guard must be silent for string BCs.
    # (bc_check may still exit for other reasons; we only care that
    #  ValueError is not raised by the dict guard.)
    try:
        flex.bc_check()
    except ValueError:
        pytest.fail("ValueError raised for string BCs on non-FD solver")


# ---------------------------------------------------------------------------
# FFT: partial-periodic BC warning
# ---------------------------------------------------------------------------

def _run_1d_fft(qs, bc_w, bc_e):
    flex = F1D()
    flex.quiet = True
    flex.method = "fft"
    flex.solver = "direct"
    flex.g = g
    flex.E = E
    flex.nu = nu
    flex.rho_m = rho_m
    flex.rho_fill = rho_fill
    flex.T_e = Te
    flex.qs = qs.copy()
    flex.dx = dx
    flex.bc_west = bc_w
    flex.bc_east = bc_e
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        flex.initialize()
        flex.run()
    return [str(x.message) for x in w if issubclass(x.category, UserWarning)]


def _run_2d_fft(qs, bc_w, bc_e, bc_n, bc_s):
    flex = F2D()
    flex.quiet = True
    flex.method = "fft"
    flex.solver = "direct"
    flex.g = g
    flex.E = E
    flex.nu = nu
    flex.rho_m = rho_m
    flex.rho_fill = rho_fill
    flex.T_e = Te
    flex.qs = qs.copy()
    flex.dx = dx
    flex.dy = dx
    flex.bc_west = bc_w
    flex.bc_east = bc_e
    flex.bc_north = bc_n
    flex.bc_south = bc_s
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        flex.initialize()
        flex.run()
    return [str(x.message) for x in w if issubclass(x.category, UserWarning)]


def test_fft_1d_partial_periodic_warns():
    """FFT 1-D: one periodic BC and one non-periodic raises UserWarning."""
    qs = np.zeros(100)
    qs[50] = 1e6
    msgs = _run_1d_fft(qs, bc_w="periodic", bc_e="")
    assert any("no_outside_loads" in m for m in msgs), msgs
    assert any("bc_east" in m for m in msgs), msgs


def test_fft_1d_all_periodic_no_warn():
    """FFT 1-D: both periodic → no partial-periodic warning."""
    qs = np.zeros(100)
    qs[50] = 1e6
    msgs = _run_1d_fft(qs, bc_w="periodic", bc_e="periodic")
    assert not any("no_outside_loads" in m for m in msgs), msgs


def test_fft_1d_none_periodic_no_warn():
    """FFT 1-D: no periodic BCs → no partial-periodic warning."""
    qs = np.zeros(100)
    qs[50] = 1e6
    msgs = _run_1d_fft(qs, bc_w="", bc_e="")
    assert not any("no_outside_loads" in m for m in msgs), msgs


def test_fft_2d_partial_periodic_warns():
    """FFT 2-D: one periodic BC and three non-periodic raises UserWarning."""
    qs = np.zeros((60, 60))
    qs[30, 30] = 1e6
    msgs = _run_2d_fft(qs, bc_w="periodic", bc_e="", bc_n="", bc_s="")
    assert any("no_outside_loads" in m for m in msgs), msgs
    assert any("bc_east" in m for m in msgs), msgs


def test_fft_2d_all_periodic_no_warn():
    """FFT 2-D: all four periodic → no partial-periodic warning."""
    qs = np.zeros((60, 60))
    qs[30, 30] = 1e6
    msgs = _run_2d_fft(qs, bc_w="periodic", bc_e="periodic",
                       bc_n="periodic", bc_s="periodic")
    assert not any("no_outside_loads" in m for m in msgs), msgs


def test_fft_2d_none_periodic_no_warn():
    """FFT 2-D: no periodic BCs → no partial-periodic warning."""
    qs = np.zeros((60, 60))
    qs[30, 30] = 1e6
    msgs = _run_2d_fft(qs, bc_w="", bc_e="", bc_n="", bc_s="")
    assert not any("no_outside_loads" in m for m in msgs), msgs


def test_fft_2d_partial_periodic_names_non_periodic_bcs():
    """FFT 2-D warning names the unpaired BC side."""
    qs = np.zeros((60, 60))
    qs[30, 30] = 1e6
    msgs = _run_2d_fft(qs, bc_w="periodic", bc_e="periodic", bc_n="periodic", bc_s="")
    partial = [m for m in msgs if "no_outside_loads" in m]
    assert partial, "expected a partial-periodic warning"
    assert "bc_south" in partial[0]
    assert "bc_west" not in partial[0]


def test_fft_2d_mixed_axes_no_warn():
    """FFT 2-D: paired x-periodic + unset y is valid — no warning."""
    qs = np.zeros((60, 60))
    qs[30, 30] = 1e6
    msgs = _run_2d_fft(qs, bc_w="periodic", bc_e="periodic", bc_n="", bc_s="")
    assert not any("no_outside_loads" in m for m in msgs), msgs


def test_fft_2d_mixed_axes_warns_for_unpaired_x():
    """FFT 2-D: bc_west='periodic' but bc_east unset → warning naming bc_east."""
    qs = np.zeros((60, 60))
    qs[30, 30] = 1e6
    msgs = _run_2d_fft(qs, bc_w="periodic", bc_e="", bc_n="periodic", bc_s="periodic")
    partial = [m for m in msgs if "no_outside_loads" in m]
    assert partial, "expected a partial-periodic warning"
    assert "bc_east" in partial[0]


def test_fft_2d_two_unpaired_warns_twice():
    """FFT 2-D: both pairs one-sided → two separate warnings."""
    qs = np.zeros((60, 60))
    qs[30, 30] = 1e6
    msgs = _run_2d_fft(qs, bc_w="periodic", bc_e="", bc_n="periodic", bc_s="")
    partial = [m for m in msgs if "no_outside_loads" in m]
    assert len(partial) == 2, f"expected 2 partial-periodic warnings, got {len(partial)}"


# ---------------------------------------------------------------------------
# FD: one-sided periodic BC warning
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("bc_w,bc_e", [
    ("periodic", "zero_displacement_zero_slope"),
    ("zero_displacement_zero_slope", "periodic"),
])
def test_fd_1d_one_sided_periodic_warns(bc_w, bc_e):
    """FD 1-D: exactly one of west/east 'periodic' raises UserWarning."""
    qs = np.zeros(80)
    qs[40] = 1e6
    msgs = _run_1d(qs, bc_w, bc_e)
    assert any("non-physical" in m for m in msgs), msgs


def test_fd_1d_both_periodic_no_warn():
    """FD 1-D: both periodic → no one-sided-periodic warning."""
    qs = np.zeros(80)
    qs[40] = 1e6
    msgs = _run_1d(qs, "periodic", "periodic")
    assert not any("non-physical" in m for m in msgs), msgs


def test_fd_1d_neither_periodic_no_warn():
    """FD 1-D: neither periodic → no one-sided-periodic warning."""
    qs = np.zeros(80)
    qs[40] = 1e6
    msgs = _run_1d(qs, "zero_displacement_zero_slope", "zero_displacement_zero_slope")
    assert not any("non-physical" in m for m in msgs), msgs


@pytest.mark.parametrize("pair", [
    ("periodic", "zero_displacement_zero_slope", "mirror", "mirror"),   # W only
    ("mirror", "periodic", "mirror", "mirror"),                          # E only
    ("mirror", "mirror", "periodic", "zero_displacement_zero_slope"),   # N only
    ("mirror", "mirror", "mirror", "periodic"),                          # S only
])
def test_fd_2d_one_sided_periodic_warns(pair):
    """FD 2-D: exactly one side of a pair 'periodic' raises UserWarning."""
    qs = np.zeros((80, 80))
    qs[40, 40] = 1e6
    bc_w, bc_e, bc_n, bc_s = pair
    msgs = _run_2d(qs, bc_w, bc_e, bc_n, bc_s)
    assert any("non-physical" in m for m in msgs), msgs


def test_fd_2d_all_periodic_no_warn():
    """FD 2-D: all four periodic → no one-sided-periodic warning."""
    qs = np.zeros((80, 80))
    qs[40, 40] = 1e6
    msgs = _run_2d(qs, "periodic", "periodic", "periodic", "periodic")
    assert not any("non-physical" in m for m in msgs), msgs


def test_fd_2d_one_sided_periodic_names_sides():
    """FD 2-D warning message names the periodic side and the missing side."""
    qs = np.zeros((80, 80))
    qs[40, 40] = 1e6
    msgs = _run_2d(qs, "periodic", "mirror", "mirror", "mirror")
    onesided = [m for m in msgs if "non-physical" in m]
    assert onesided, "expected a one-sided-periodic warning"
    assert "bc_west" in onesided[0]
    assert "bc_east" in onesided[0]
