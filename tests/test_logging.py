#! /usr/bin/env python

import logging

import numpy as np
import pytest

from gflex import F1D


# ── Helpers ────────────────────────────────────────────────────────────────────

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
