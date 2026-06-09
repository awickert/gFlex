#! /usr/bin/env python

import logging
import time
import warnings

import numpy as np
import pytest

from gflex import F1D, F2D


# ── Helpers ────────────────────────────────────────────────────────────────────

def _make_fd_1d(**kwargs):
    flex = F1D()
    flex.quiet = True
    flex.method = "fd"
    flex.bc_west = flex.bc_east = "zero_moment_zero_shear"
    flex.E = 65e9; flex.nu = 0.25
    flex.rho_m = 3300.0; flex.rho_fill = 0.0; flex.g = 9.8
    flex.T_e = 30e3; flex.dx = 10e3
    n = 51
    flex.qs = np.zeros(n); flex.qs[n // 2] = 1e6
    for k, v in kwargs.items():
        setattr(flex, k, v)
    return flex


def _run_minimal_1d(quiet=True, verbose=False, debug=False):
    flex = F1D()
    flex.quiet = quiet
    flex.verbose = verbose
    flex.debug = debug
    flex.method = "sas"
    flex.bc_west = ""
    flex.bc_east = ""
    flex.E = 65e9
    flex.nu = 0.25
    flex.rho_m = 3300.0
    flex.rho_fill = 0.0
    flex.g = 9.8
    flex.T_e = 30e3
    n = 51
    flex.dx = 10e3
    flex.qs = np.zeros(n)
    flex.qs[n // 2] = 1e6
    flex.initialize()
    flex.run()
    flex.finalize()


# ── 1. Default verbosity flags ─────────────────────────────────────────────────

def test_default_flags():
    flex = F1D()
    assert flex.quiet is True
    assert flex.verbose is False
    assert flex.debug is False


# ── 2. No console output by default ───────────────────────────────────────────

def test_default_produces_no_stdout(capfd):
    _run_minimal_1d(quiet=True)
    out, err = capfd.readouterr()
    assert out == ""
    assert err == ""


# ── 3. verbose=True produces log records ──────────────────────────────────────

def test_verbose_produces_log_records(caplog):
    with caplog.at_level(logging.INFO, logger="gflex"):
        _run_minimal_1d(quiet=False, verbose=True)
    assert len(caplog.records) > 0


# ── 4. Flag interactions ───────────────────────────────────────────────────────

def test_quiet_overrides_verbose(caplog):
    """quiet=True must suppress INFO records even when verbose=True."""
    with caplog.at_level(logging.INFO, logger="gflex"):
        _run_minimal_1d(quiet=True, verbose=True)
    info_records = [r for r in caplog.records if r.levelno == logging.INFO]
    assert info_records == []


def test_debug_flag_enables_debug_level():
    flex = F1D()
    flex.debug = True
    assert logging.getLogger("gflex").level == logging.DEBUG


def test_quiet_flag_sets_warning_level():
    flex = F1D()
    flex.quiet = True
    assert logging.getLogger("gflex").level == logging.WARNING


# ── 5. Timing attributes ───────────────────────────────────────────────────────

def test_fd_sets_linear_solve_time():
    """FD run sets linear_solve_time; it is positive and <= time_to_solve."""
    flex = _make_fd_1d()
    flex.initialize()
    flex.run()
    assert hasattr(flex, "linear_solve_time")
    assert flex.linear_solve_time >= 0
    assert flex.linear_solve_time <= flex.time_to_solve


def test_sas_does_not_set_linear_solve_time():
    """SAS run does not set linear_solve_time (no matrix factorisation)."""
    flex = F1D()
    flex.quiet = True
    flex.method = "sas"
    flex.bc_west = flex.bc_east = ""
    flex.E = 65e9; flex.nu = 0.25
    flex.rho_m = 3300.0; flex.rho_fill = 0.0; flex.g = 9.8
    flex.T_e = 30e3; flex.dx = 10e3
    n = 51
    flex.qs = np.zeros(n); flex.qs[n // 2] = 1e6
    flex.initialize()
    flex.run()
    assert not hasattr(flex, "linear_solve_time")


def test_2d_fd_sets_linear_solve_time():
    """2D FD run sets linear_solve_time >= 0 and <= time_to_solve."""
    flex = F2D()
    flex.quiet = True
    flex.method = "fd"
    flex.bc_north = flex.bc_south = "zero_moment_zero_shear"
    flex.bc_west  = flex.bc_east  = "zero_moment_zero_shear"
    flex.E = 65e9; flex.nu = 0.25
    flex.rho_m = 3300.0; flex.rho_fill = 0.0; flex.g = 9.8
    flex.T_e = 30e3; flex.dx = flex.dy = 10e3
    n = 20
    flex.qs = np.zeros((n, n)); flex.qs[n // 2, n // 2] = 1e6
    flex.initialize()
    flex.run()
    assert hasattr(flex, "linear_solve_time")
    assert flex.linear_solve_time >= 0
    assert flex.linear_solve_time <= flex.time_to_solve


def test_fft_does_not_set_linear_solve_time():
    """FFT run does not set linear_solve_time (no matrix factorisation)."""
    flex = F1D()
    flex.quiet = True
    flex.method = "fft"
    flex.bc_west = flex.bc_east = ""
    flex.E = 65e9; flex.nu = 0.25
    flex.rho_m = 3300.0; flex.rho_fill = 0.0; flex.g = 9.8
    flex.T_e = 30e3; flex.dx = 10e3
    n = 51
    flex.qs = np.zeros(n); flex.qs[n // 2] = 1e6
    flex.initialize()
    flex.run()
    assert not hasattr(flex, "linear_solve_time")


def test_total_start_time_precedes_solver():
    """_total_start_time is set before the solver starts."""
    flex = _make_fd_1d()
    flex.initialize()
    assert hasattr(flex, "_total_start_time")
    flex.run()
    assert flex._total_start_time <= time.time() - flex.time_to_solve


# ── 6. bc_check UserWarning for FD BCs passed to analytical solver ─────────────

def test_bc_check_warns_fd_bcs_on_sas():
    """bc_check issues UserWarning when FD BCs are passed to SAS solver."""
    flex = F1D()
    flex.quiet = True
    flex.method = "sas"
    flex.bc_west = flex.bc_east = "zero_moment_zero_shear"
    flex.E = 65e9; flex.nu = 0.25
    flex.rho_m = 3300.0; flex.rho_fill = 0.0; flex.g = 9.8
    flex.T_e = 30e3; flex.dx = 10e3
    n = 51
    flex.qs = np.zeros(n); flex.qs[n // 2] = 1e6
    flex.initialize()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        flex.run()
    messages = [str(w.message) for w in caught if issubclass(w.category, UserWarning)]
    assert any("no_outside_loads" in m for m in messages)
