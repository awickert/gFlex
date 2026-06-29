"""v1.x PascalCase boundary-condition strings are rejected in v2.0 (clean break).

They are not normalised or deprecation-warned: they raise ``ValueError``.
Use the lowercase canonical names or the short aliases instead.
"""

import warnings

import numpy as np
import pytest

from gflex.f1d import F1D
from gflex.f2d import F2D

_V1X = ["0Displacement0Slope", "0Moment0Shear", "Mirror", "Periodic", "NoOutsideLoads"]


def _run_1d(bc_west):
    f = F1D()
    f.quiet = True
    f.verbose = False
    f.debug = False
    f.method = "fd"
    f.g, f.E, f.nu, f.rho_m, f.rho_fill = 9.8, 65e9, 0.25, 3300.0, 0.0
    f.dx = 30.0e3
    f.T_e = 30.0e3
    f.qs = np.zeros(20)
    f.bc_west = bc_west
    f.bc_east = "zero_displacement_zero_slope"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        f.initialize()
        f.run()


@pytest.mark.parametrize("bc", _V1X)
def test_1d_v1x_pascalcase_bc_raises(bc):
    with pytest.raises(ValueError):
        _run_1d(bc)


def test_2d_v1x_pascalcase_bc_raises():
    f = F2D()
    f.quiet = True
    f.verbose = False
    f.debug = False
    f.method = "fd"
    f.solver = "direct"
    f.g, f.E, f.nu, f.rho_m, f.rho_fill = 9.8, 65e9, 0.25, 3300.0, 0.0
    f.dx = f.dy = 30.0e3
    f.T_e = 30.0e3
    f.qs = np.zeros((20, 20))
    f.bc_west = "0Moment0Shear"
    f.bc_east = f.bc_north = f.bc_south = "zero_displacement_zero_slope"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with pytest.raises(ValueError):
            f.initialize()
            f.run()


def test_lowercase_canonical_and_aliases_accepted():
    # Sanity: the supported forms still work.
    _run_1d("mirror")          # short alias
    _run_1d("zero_moment_zero_shear")   # canonical
