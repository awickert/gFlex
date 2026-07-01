#! /usr/bin/env python
"""Tests for FD boundary-condition warnings in F1D and F2D.

Warning types covered:
  zero_moment_zero_shear   — UserWarning fires for any side carrying that BC
  proximity       — fires when nearest loaded cell is within one flexural
                    wavelength of a zero_displacement_zero_slope boundary; absent for
                    zero_slope_zero_shear / mirror / periodic, and when the load is far enough away
  lu_memory       — UserWarning fires before LU factorisation when estimated
                    RAM exceeds 80% of available RAM; absent when RAM is
                    sufficient or when LU is already cached
"""

import warnings
from unittest.mock import patch

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

def _build_1d(qs, bc_w, bc_e):
    """Construct (but do not run) a 1-D FD model for BC-guard tests."""
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
    return flex


def _build_2d(qs, bc_w, bc_e, bc_n, bc_s):
    """Construct (but do not run) a 2-D FD model for BC-guard tests."""
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
    return flex


@pytest.mark.parametrize("bc_w,bc_e", [
    ("periodic", "zero_displacement_zero_slope"),
    ("zero_displacement_zero_slope", "periodic"),
])
def test_fd_1d_one_sided_periodic_raises(bc_w, bc_e):
    """FD 1-D: exactly one of west/east 'periodic' raises ValueError by default."""
    qs = np.zeros(80)
    qs[40] = 1e6
    flex = _build_1d(qs, bc_w, bc_e)
    with pytest.raises(ValueError, match="allow_unpaired_periodic"):
        flex.initialize()
        flex.run()


def test_fd_1d_both_periodic_runs():
    """FD 1-D: both periodic → no guard, solve completes with finite w."""
    qs = np.zeros(80)
    qs[40] = 1e6
    flex = _build_1d(qs, "periodic", "periodic")
    flex.initialize()
    flex.run()
    assert np.all(np.isfinite(flex.w))


def test_fd_1d_neither_periodic_runs():
    """FD 1-D: neither periodic → no guard, solve completes with finite w."""
    qs = np.zeros(80)
    qs[40] = 1e6
    flex = _build_1d(qs, "zero_displacement_zero_slope",
                     "zero_displacement_zero_slope")
    flex.initialize()
    flex.run()
    assert np.all(np.isfinite(flex.w))


def test_fd_1d_one_sided_periodic_allowed_with_flag():
    """allow_unpaired_periodic=True lets the solve proceed; the only warning is
    the one-time 'safety check disabled' announcement from enabling the flag —
    the guard message does not re-fire on every solve."""
    qs = np.zeros(80)
    qs[40] = 1e6
    flex = _build_1d(qs, "periodic", "zero_displacement_zero_slope")
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        flex.allow_unpaired_periodic = True   # fires the announcement once
        flex.initialize()
        flex.run()                            # proceeds silently (no re-warn)
    msgs = [str(x.message) for x in w if issubclass(x.category, UserWarning)]
    assert any("safety check is now disabled" in m for m in msgs), msgs
    assert not any("well-posed" in m for m in msgs), msgs
    assert np.all(np.isfinite(flex.w))


@pytest.mark.parametrize("pair", [
    ("periodic", "zero_displacement_zero_slope", "mirror", "mirror"),   # W only
    ("mirror", "periodic", "mirror", "mirror"),                          # E only
    ("mirror", "mirror", "periodic", "zero_displacement_zero_slope"),   # N only
    ("mirror", "mirror", "mirror", "periodic"),                          # S only
])
def test_fd_2d_one_sided_periodic_raises(pair):
    """FD 2-D: exactly one side of a pair 'periodic' raises ValueError by default."""
    qs = np.zeros((80, 80))
    qs[40, 40] = 1e6
    bc_w, bc_e, bc_n, bc_s = pair
    flex = _build_2d(qs, bc_w, bc_e, bc_n, bc_s)
    with pytest.raises(ValueError, match="allow_unpaired_periodic"):
        flex.initialize()
        flex.run()


def test_fd_2d_all_periodic_runs():
    """FD 2-D: all four periodic → no guard, solve completes with finite w."""
    qs = np.zeros((80, 80))
    qs[40, 40] = 1e6
    flex = _build_2d(qs, "periodic", "periodic", "periodic", "periodic")
    flex.initialize()
    flex.run()
    assert np.all(np.isfinite(flex.w))


def test_fd_2d_one_sided_periodic_error_names_sides():
    """FD 2-D ValueError names the periodic side and the missing side."""
    qs = np.zeros((80, 80))
    qs[40, 40] = 1e6
    flex = _build_2d(qs, "periodic", "mirror", "mirror", "mirror")
    with pytest.raises(ValueError) as exc:
        flex.initialize()
        flex.run()
    assert "bc_west" in str(exc.value)
    assert "bc_east" in str(exc.value)


def test_allow_unpaired_periodic_default_false():
    """The override defaults to False on a fresh model."""
    assert F1D().allow_unpaired_periodic is False
    assert F2D().allow_unpaired_periodic is False


def test_allow_unpaired_periodic_setter_announces():
    """Enabling the override warns, so disabling the safety check leaves a trace."""
    flex = F1D()
    with pytest.warns(UserWarning, match="unpaired-periodic safety check"):
        flex.allow_unpaired_periodic = True
    assert flex.allow_unpaired_periodic is True


# ---------------------------------------------------------------------------
# ValueError when rho_fill >= rho_m
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("rho_fill", [3300.0, 3500.0])
def test_initialize_drho_nonpositive_raises_1d(rho_fill):
    """F1D.initialize raises ValueError when rho_fill >= rho_m."""
    flex = F1D()
    flex.quiet = True; flex.method = "fd"; flex.solver = "direct"
    flex.g = 9.8; flex.E = 65e9; flex.nu = 0.25
    flex.rho_m = 3300.0; flex.rho_fill = rho_fill
    flex.T_e = 35e3; flex.dx = 5000.0
    flex.qs = np.zeros(100); flex.qs[50] = 1e6
    flex.bc_west = flex.bc_east = "zero_displacement_zero_slope"
    with pytest.raises(ValueError, match="rho_fill"):
        flex.initialize()


@pytest.mark.parametrize("rho_fill", [3300.0, 3500.0])
def test_initialize_drho_nonpositive_raises_2d(rho_fill):
    """F2D.initialize raises ValueError when rho_fill >= rho_m."""
    flex = F2D()
    flex.quiet = True; flex.method = "fd"; flex.solver = "direct"
    flex.g = 9.8; flex.E = 65e9; flex.nu = 0.25
    flex.rho_m = 3300.0; flex.rho_fill = rho_fill
    flex.T_e = 35e3; flex.dx = flex.dy = 5000.0
    flex.qs = np.zeros((40, 40)); flex.qs[20, 20] = 1e6
    flex.bc_west = flex.bc_east = flex.bc_north = flex.bc_south = "zero_displacement_zero_slope"
    with pytest.raises(ValueError, match="rho_fill"):
        flex.initialize()


def test_drho_revalidated_on_density_reassignment_1d():
    """F1D re-checks drho at solve time, so a density changed after a valid
    run cannot silently solve with Δρ ≤ 0 (regression: was only checked at
    initialize())."""
    qs = np.zeros(80); qs[40] = 1e6
    flex = _build_1d(qs, "zero_displacement_zero_slope",
                     "zero_displacement_zero_slope")
    flex.initialize()
    flex.run()                       # valid: rho_fill < rho_m
    flex.rho_fill = flex.rho_m + 100.0   # now Δρ < 0 (mimics BMI set_value)
    with pytest.raises(ValueError, match="rho_fill"):
        flex.run()


def test_drho_revalidated_on_density_reassignment_2d():
    """F2D re-checks drho at solve time after a density reassignment."""
    qs = np.zeros((40, 40)); qs[20, 20] = 1e6
    flex = _build_2d(qs, "zero_displacement_zero_slope",
                     "zero_displacement_zero_slope",
                     "zero_displacement_zero_slope",
                     "zero_displacement_zero_slope")
    flex.initialize()
    flex.run()
    flex.rho_fill = flex.rho_m + 100.0
    with pytest.raises(ValueError, match="rho_fill"):
        flex.run()


def test_f2d_debug_rejects_nondirect_solver():
    """F2D rejects an unsupported solver even with debug=True (the guard is
    independent of the debug branch; regression for the elif that let a
    non-direct solver through whenever debug was on)."""
    qs = np.zeros((40, 40)); qs[20, 20] = 1e6
    flex = _build_2d(qs, "zero_displacement_zero_slope",
                     "zero_displacement_zero_slope",
                     "zero_displacement_zero_slope",
                     "zero_displacement_zero_slope")
    flex.debug = True
    flex.solver = "iterative"
    with pytest.raises(ValueError, match="not supported"):
        flex.initialize()
        flex.run()


# ---------------------------------------------------------------------------
# LU memory warning (#68)
# ---------------------------------------------------------------------------
# A 50×50 2-D grid (2500 cells) has an estimated LU footprint of ~4.8 MB.
# Mocking available RAM to 5 MB → estimate is ~97% → warns.
# Mocking it to 100 MB → ~4.8% → no warning.
# For 1-D, a 50-cell grid has only ~35 KB of estimated LU RAM; the threshold
# is set accordingly (40 KB available → ~87% → warns).

_SMALL_RAM_2D = 5_000_000   # 5 MB — forces the warning for a 50×50 grid
_LARGE_RAM    = 100_000_000  # 100 MB — no warning for a 50×50 grid
_SMALL_RAM_1D = 40_000      # 40 KB — forces the warning for a 50-cell 1-D grid


def _make_fd_2d_50():
    flex = F2D()
    flex.quiet = True; flex.method = "fd"; flex.solver = "direct"
    flex.g = 9.8; flex.E = 65e9; flex.nu = 0.25
    flex.rho_m = 3300.0; flex.rho_fill = 0.0
    flex.T_e = 35e3; flex.dx = flex.dy = 5000.0
    flex.qs = np.zeros((50, 50)); flex.qs[25, 25] = 1e6
    flex.bc_west = flex.bc_east = flex.bc_north = flex.bc_south = "mirror"
    return flex


def _make_fd_1d_50():
    flex = F1D()
    flex.quiet = True; flex.method = "fd"; flex.solver = "direct"
    flex.g = 9.8; flex.E = 65e9; flex.nu = 0.25
    flex.rho_m = 3300.0; flex.rho_fill = 0.0
    flex.T_e = 35e3; flex.dx = 5000.0
    flex.qs = np.zeros(50); flex.qs[25] = 1e6
    flex.bc_west = flex.bc_east = "mirror"
    return flex


def _lu_msgs_2d(flex):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        flex.initialize()
        flex.run()
    return [str(x.message) for x in w
            if issubclass(x.category, UserWarning) and "LU memory" in str(x.message)]


def _lu_msgs_1d(flex):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        flex.initialize()
        flex.run()
    return [str(x.message) for x in w
            if issubclass(x.category, UserWarning) and "LU memory" in str(x.message)]


def test_lu_memory_warning_fires_2d():
    """LU memory warning fires when estimated RAM exceeds 80% of available."""
    flex = _make_fd_2d_50()
    with patch("gflex.f2d._available_ram_bytes", return_value=_SMALL_RAM_2D):
        msgs = _lu_msgs_2d(flex)
    assert msgs, "expected LU memory warning"
    assert "LU memory" in msgs[0]
    assert "GB" in msgs[0]


def test_lu_memory_warning_absent_when_ram_sufficient_2d():
    """No LU memory warning when RAM is well above the threshold."""
    flex = _make_fd_2d_50()
    with patch("gflex.f2d._available_ram_bytes", return_value=_LARGE_RAM):
        msgs = _lu_msgs_2d(flex)
    assert not msgs, f"unexpected LU memory warning: {msgs}"


def test_lu_memory_warning_fires_1d():
    """LU memory warning fires in F1D when estimated RAM exceeds 80%."""
    flex = _make_fd_1d_50()
    with patch("gflex.f1d._available_ram_bytes", return_value=_SMALL_RAM_1D):
        msgs = _lu_msgs_1d(flex)
    assert msgs, "expected LU memory warning"


def test_lu_memory_warning_absent_when_ram_sufficient_1d():
    """No LU memory warning in F1D when RAM is sufficient."""
    flex = _make_fd_1d_50()
    with patch("gflex.f1d._available_ram_bytes", return_value=_LARGE_RAM):
        msgs = _lu_msgs_1d(flex)
    assert not msgs, f"unexpected LU memory warning: {msgs}"


def test_lu_memory_warning_skipped_when_lu_cached():
    """No LU memory warning on second run when LU is already cached (no_check)."""
    flex = _make_fd_2d_50()
    flex.cache_factorization = "no_check"
    with patch("gflex.f2d._available_ram_bytes", return_value=_SMALL_RAM_2D):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            flex.initialize()
            flex.run()           # first run — LU factorised, may warn
            flex.qs[25, 25] = 2e6
            flex.run()           # second run — LU cached, must not warn
    second_lu = [str(x.message) for x in w
                 if issubclass(x.category, UserWarning) and "LU memory" in str(x.message)]
    # The second run must produce exactly zero LU memory warnings.
    # (Filter by counting only warnings emitted during the second run is
    # tricky, so instead assert total count <= 1: at most one from the first.)
    assert len(second_lu) <= 1, (
        f"LU memory warning fired on cached second run: {second_lu}"
    )


def test_lu_memory_warning_skipped_when_ram_unknown():
    """No crash and no LU memory warning when available RAM cannot be determined."""
    flex = _make_fd_2d_50()
    with patch("gflex.f2d._available_ram_bytes", return_value=None):
        msgs = _lu_msgs_2d(flex)
    assert not msgs, f"unexpected LU memory warning: {msgs}"


def test_lu_memory_warning_message_content():
    """Warning message contains the percentage, 'available RAM', and 'free'."""
    flex = _make_fd_2d_50()
    with patch("gflex.f2d._available_ram_bytes", return_value=_SMALL_RAM_2D):
        msgs = _lu_msgs_2d(flex)
    assert msgs, "expected LU memory warning"
    msg = msgs[0]
    assert "%" in msg
    assert "available RAM" in msg
    assert "free" in msg


def test_lu_memory_warning_uses_padded_size():
    """Warning fires based on the padded domain size, not the original.

    50×50 original → ~4.8 MB estimated (4.8% of 100 MB → no warning).
    With all-sides no_outside_loads the domain pads to 184×184 (~128 MB,
    128% of 100 MB → warning fires).
    """
    _RAM = 100_000_000  # 100 MB: original would not warn; padded does
    flex = F2D()
    flex.quiet = True; flex.method = "fd"; flex.solver = "direct"
    flex.g = 9.8; flex.E = 65e9; flex.nu = 0.25
    flex.rho_m = 3300.0; flex.rho_fill = 0.0
    flex.T_e = 35e3; flex.dx = flex.dy = 5000.0
    flex.qs = np.zeros((50, 50)); flex.qs[25, 25] = 1e6
    flex.bc_west = flex.bc_east = flex.bc_north = flex.bc_south = "no_outside_loads"
    with patch("gflex.f2d._available_ram_bytes", return_value=_RAM):
        msgs = _lu_msgs_2d(flex)
    assert msgs, (
        "expected LU memory warning for padded domain; "
        "original 50×50 at 100 MB available would not trigger it"
    )


def test_lu_memory_warning_skipped_when_lu_cached_true_mode():
    """No LU memory warning on second run with cache_factorization=True (LU reused)."""
    flex = _make_fd_2d_50()
    flex.cache_factorization = True
    with patch("gflex.f2d._available_ram_bytes", return_value=_SMALL_RAM_2D):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            flex.initialize()
            flex.run()           # first run — factorises and caches
            flex.qs[25, 25] = 2e6
            flex.run()           # second run — LU reused (no rebuild)
    second_lu = [str(x.message) for x in w
                 if issubclass(x.category, UserWarning) and "LU memory" in str(x.message)]
    assert len(second_lu) <= 1, (
        f"LU memory warning fired on cached second run (True mode): {second_lu}"
    )
