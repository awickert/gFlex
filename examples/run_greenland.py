"""
Run the Greenland ice-sheet flexure example from its YAML configuration.

Everything -- parameters, input grids, boundary conditions, and the plot --
is defined in ``greenland.yaml``.  This launcher just hands that file to
gFlex and lets it do the work, including the standard load + deflection plot.

    python examples/run_greenland.py

Equivalent command-line form:

    python -m gflex examples/greenland.yaml
"""

import os

import gflex

CONFIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "greenland.yaml")

flex = gflex.F2D(CONFIG)
flex.initialize(CONFIG)
flex.run()
flex.output()    # standard gFlex plot (plot: both) + optional file output
flex.finalize()
