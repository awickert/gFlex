#!/usr/bin/env python3
"""
MMS runner for BC comparison.  Imported by compare_bc_versions.py.

Runs the zero_displacement_zero_slope MMS convergence study using
whatever gflex package is first on sys.path, and writes a JSON
result to stdout.

Manufactured solution
---------------------
  xi = x / L,  L = (N-1)*dx

  w_exact(xi) = -W0 * xi**2 * (1 - xi)**2

  qs(xi)      = 24*D*W0/L**4 + drho_g * W0 * xi**2 * (1 - xi)**2

Both BCs satisfied exactly at xi=0 and xi=1:
  w=0, dw/dx=0  (zero displacement and slope at each end)
"""

import json
import sys
import warnings

import numpy as np

import gflex as _gflex_mod
from gflex.f1d import F1D

# Physical parameters — identical to benchmarks/analyze_clamped_bc_error.py
TE      = 30.0e3    # elastic thickness, m
E       = 65.0e9    # Young's modulus, Pa
NU      = 0.25
RHO_M   = 3300.0    # mantle density, kg m-3
RHO_F   = 0.0       # infill (air)
G       = 9.81      # m s-2
L       = 600.0e3   # plate length, m
W0      = 1600.0    # MMS amplitude, m  (max deflection = W0/16 = 100 m)

NX_VALS = [26, 51, 101, 201, 401, 801]


def _run_once(nx):
    dx    = L / (nx - 1)
    xi    = np.arange(nx) / (nx - 1)
    D     = E * TE**3 / (12.0 * (1.0 - NU**2))
    drho_g = (RHO_M - RHO_F) * G
    w_ex  = -W0 * xi**2 * (1 - xi)**2
    qs    = 24.0 * D * W0 / L**4 + drho_g * W0 * xi**2 * (1 - xi)**2

    s = F1D()
    s.dx      = dx
    s.te      = TE
    s.E       = E
    s.nu      = NU
    s.rho_m   = RHO_M
    s.rho_fill = RHO_F
    s.g       = G
    s.qs      = qs
    s.bc_west = "zero_displacement_zero_slope"
    s.bc_east = "zero_displacement_zero_slope"
    s.method  = "fd"
    s.quiet   = True
    s.verbose = False
    s.debug   = False
    s.initialize()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s.run()
    s.finalize()

    err = np.max(np.abs(s.w - w_ex)) / np.max(np.abs(w_ex))
    return float(err), float(s.w[0]), float(s.w[-1])


errors = []
w0s    = []
for nx in NX_VALS:
    err, w0, wN = _run_once(nx)
    errors.append(err)
    w0s.append(w0)

dx_vals = [L / (nx - 1) for nx in NX_VALS]

# Convergence slope over finest half
n_fit = len(NX_VALS) // 2
slope = float(np.polyfit(
    np.log(dx_vals[-n_fit:]), np.log(errors[-n_fit:]), 1
)[0])

result = {
    "gflex_version": _gflex_mod.__version__,
    "nx_vals":  NX_VALS,
    "dx_km":    [dx / 1e3 for dx in dx_vals],
    "errors":   errors,
    "w0_at_finest": w0s[-1],
    "convergence_slope": slope,
}

json.dump(result, sys.stdout, indent=2)
sys.stdout.write("\n")
