"""
Greenland ice-sheet flexure -- the same run as greenland.yaml, but built
entirely in Python instead of from a configuration file.

The two input grids are pre-built by examples/preprocess_greenland.py; this
script just loads them, sets the plate parameters, runs gFlex, and makes the
standard load + deflection plot.

    python examples/greenland_script.py
"""

import os

import numpy as np

import gflex

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

flex = gflex.F2D()
flex.method = "fd"

# Plate parameters (SI units)
flex.g = 9.81            # gravitational acceleration
flex.E = 6.5e10          # Young's modulus
flex.nu = 0.25           # Poisson's ratio
flex.rho_m = 3300.0      # mantle density
flex.rho_fill = 0.0      # infill density (air)

# Pre-built input grids (padded; see preprocess_greenland.py)
flex.qs = np.load(os.path.join(DATA, "greenland_load.npy"))   # load stress [Pa]
flex.T_e = np.load(os.path.join(DATA, "greenland_te.npy"))    # elastic thickness [m]
flex.dx = flex.dy = 10050.0                                   # grid spacing [m]

# Clamped boundaries, held far from the load by the padding
flex.bc_west = flex.bc_east = "zero_displacement_zero_slope"
flex.bc_north = flex.bc_south = "zero_displacement_zero_slope"

flex.initialize()
flex.run()

flex.plot_choice = "both"    # standard gFlex load + deflection plot
flex.output()
flex.finalize()
