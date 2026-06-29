"""The 'sandbox' BC easter egg must run without crashing.

It previously called ``finalize()`` (which deletes ``w``) before reading
``solver.w``, raising ``AttributeError`` every time it was triggered.
"""

import contextlib
import io
import warnings

import numpy as np
import pytest

from gflex.f1d import F1D


def test_sandbox_bc_runs_and_produces_output():
    f = F1D()
    f.quiet = True
    f.verbose = False
    f.debug = False
    f.method = "fd"
    f.g, f.E, f.nu, f.rho_m, f.rho_fill = 9.8, 65e9, 0.25, 3300.0, 0.0
    f.dx = 30.0e3
    f.T_e = 30.0e3
    f.qs = np.zeros(40)
    f.bc_west = "sandbox"
    f.bc_east = "zero_displacement_zero_slope"
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        f.initialize()
        # The easter egg prints its scene, then exits cleanly (by design).
        # The bug was an AttributeError before it could print anything.
        with pytest.raises(SystemExit) as exc:
            f.run()
    assert exc.value.code == 0
    assert buf.getvalue().strip(), "sandbox easter egg produced no output"
