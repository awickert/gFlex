#!/usr/bin/env python3
"""Smoke tests for framework-style usage that relies on gFlex defaults.

These reproduce the minimal setup pattern used by downstream callers such as
GRASS GIS r.flexure: only required physical parameters are assigned; optional
attributes (boundary conditions, verbosity flags) are left completely unset.
Several bugs found during GRASS integration only surfaced under this pattern —
notably the FFT AttributeError when bc_west / bc_east were never assigned.
"""

import warnings

import numpy as np

from gflex import F1D, F2D


def _base_1d(n=51):
    """Minimum required F1D setup — no BC attributes assigned."""
    flex = F1D()
    flex.E        = 65e9
    flex.nu       = 0.25
    flex.rho_m    = 3300.0
    flex.rho_fill = 0.0
    flex.g        = 9.8
    flex.T_e      = 30e3
    flex.dx       = 10e3
    flex.qs       = np.zeros(n)
    flex.qs[n // 2] = 1e6
    return flex


def _base_2d(n=20):
    """Minimum required F2D setup — no BC attributes assigned."""
    flex = F2D()
    flex.E        = 65e9
    flex.nu       = 0.25
    flex.rho_m    = 3300.0
    flex.rho_fill = 0.0
    flex.g        = 9.8
    flex.T_e      = 30e3
    flex.dx       = flex.dy = 10e3
    flex.qs       = np.zeros((n, n))
    flex.qs[n // 2, n // 2] = 1e6
    return flex


# ---------------------------------------------------------------------------
# FFT — no BCs ever assigned (the GRASS r.flexure pattern before the fix)
# ---------------------------------------------------------------------------

class TestFftNoBCsAssigned:
    """FFT must run when bc_west / bc_east / bc_north / bc_south are never set.

    The original bug: bc_check() set BC_E / BC_W (legacy uppercase attributes)
    instead of bc_east / bc_west; _solve_fft() then raised AttributeError.
    """

    def test_1d_fft_no_bcs_runs(self):
        flex = _base_1d()
        flex.method = "fft"
        flex.initialize()
        flex.run()
        assert np.any(flex.w != 0.0)
        flex.finalize()

    def test_2d_fft_no_bcs_runs(self):
        flex = _base_2d()
        flex.method = "fft"
        flex.initialize()
        flex.run()
        assert np.any(flex.w != 0.0)
        flex.finalize()

    def test_1d_fft_no_bcs_output_shape(self):
        flex = _base_1d()
        flex.method = "fft"
        flex.initialize()
        flex.run()
        assert flex.w.shape == flex.qs.shape
        flex.finalize()

    def test_2d_fft_no_bcs_output_shape(self):
        flex = _base_2d()
        flex.method = "fft"
        flex.initialize()
        flex.run()
        assert flex.w.shape == flex.qs.shape
        flex.finalize()


# ---------------------------------------------------------------------------
# SAS — no BCs ever assigned
# ---------------------------------------------------------------------------

class TestSasNoBCsAssigned:
    """SAS must run (and produce finite deflection) when BCs are never set."""

    def test_1d_sas_no_bcs_runs(self):
        flex = _base_1d()
        flex.method = "sas"
        flex.initialize()
        flex.run()
        assert np.any(flex.w != 0.0)
        assert np.all(np.isfinite(flex.w))
        flex.finalize()

    def test_2d_sas_no_bcs_runs(self):
        flex = _base_2d()
        flex.method = "sas"
        flex.initialize()
        flex.run()
        assert np.any(flex.w != 0.0)
        assert np.all(np.isfinite(flex.w))
        flex.finalize()


# ---------------------------------------------------------------------------
# quiet=True default — bc_check must still issue UserWarnings
# ---------------------------------------------------------------------------

class TestQuietDefaultDoesNotSuppressWarnings:
    """quiet=True (default) must not silently swallow bc_check UserWarnings.

    The original bug: bc_check guarded the warning emission with
    ``if not self.quiet``, so passing FD BCs to SAS in quiet mode (the
    default since v2.0.0) produced no diagnostic at all.
    """

    def test_sas_with_fd_bcs_warns_when_quiet(self):
        flex = _base_1d()
        flex.method  = "sas"
        flex.bc_west = flex.bc_east = "zero_moment_zero_shear"
        assert flex.quiet is True          # default — must not suppress the warn
        flex.initialize()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            flex.run()
        messages = [str(w.message) for w in caught if issubclass(w.category, UserWarning)]
        assert any("no_outside_loads" in m for m in messages), (
            "bc_check must warn about FD BCs on SAS even in the default quiet=True mode"
        )
