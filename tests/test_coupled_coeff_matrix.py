"""Supplying a prebuilt coefficient matrix (coupled mode) must not crash.

When a caller assigns ``flex.coeff_matrix`` directly, matrix assembly — and
with it the creation of ``_bc_rhs_correction`` — is skipped.  ``fd_solve``
must tolerate the missing correction (as the F2D path already did) rather
than raising ``AttributeError``.
"""

import warnings

import numpy as np

from gflex.f1d import F1D


def _make():
    f = F1D()
    f.quiet = True
    f.verbose = False
    f.debug = False
    f.method = "fd"
    f.g, f.E, f.nu, f.rho_m, f.rho_fill = 9.8, 65e9, 0.25, 3300.0, 0.0
    f.dx = 30.0e3
    f.T_e = 30.0e3
    qs = np.zeros(40)
    qs[20] = 1.0e8
    f.qs = qs
    f.bc_west = f.bc_east = "zero_displacement_zero_slope"
    return f


def test_1d_external_coeff_matrix_does_not_crash():
    ref = _make()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ref.initialize()
        ref.run()
    w_ref = ref.w.copy()
    matrix = ref.coeff_matrix

    # Fresh solver: supply the prebuilt matrix so assembly is skipped.
    flex = _make()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        flex.initialize()
    flex.coeff_matrix = matrix
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        flex.run()
    np.testing.assert_allclose(flex.w, w_ref, rtol=1e-12, atol=0)
