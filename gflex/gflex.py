#! /usr/bin/env python
"""
Multiple methods to solve elastic plate flexure, designed for applications
to Earth's lithosphere.


To generate an input file, please see the examples in the "input"
directory of this install.

To run in a Python script or shell, follow this general pattern:")

```
import gflex
flex = gflex.F1D()
flex.method = ...
# ...more variable setting...
# see the 'input' directory for examples
```
"""

import argparse
import sys

from gflex.base import WhichModel
from gflex.f1d import F1D
from gflex.f2d import F2D

from ._version import __version__

LICENSE = """
This file is part of gFlex.
gFlex computes lithospheric flexural isostasy with heterogeneous rigidity
Copyright (C) 2010-2018 Andrew D. Wickert

gFlex is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

gFlex is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with gFlex.  If not, see <http://www.gnu.org/licenses/>.
"""


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("filename", nargs=1, help="gflex configuration file.")
    parser.add_argument("--version", action="version", version=f"gflex {__version__}")
    parser.add_argument(
        "--verbose", action="store_true", help="print debugging information."
    )
    parser.add_argument("--silent", action="store_true", help="print minimal output.")

    args = parser.parse_args()

    if not args.silent:
        print(LICENSE, file=sys.stderr)

    obj = WhichModel(args.filename)

    obj.Debug = args.verbose
    obj.Quiet = args.silent

    ########################################
    ## SET MODEL TYPE AND DIMENSIONS HERE ##
    ########################################

    if obj.dimension == 1:
        obj = F1D(args.filename)
    elif obj.dimension == 2:
        obj = F2D(args.filename)

    obj.initialize(args.filename)

    ############################################
    ##       SET MODEL PARAMETERS HERE        ##
    ## (if not defined in configuration file) ##
    ############################################
    # obj.set_value('method','FD') # for example

    obj.run()
    obj.finalize()

    obj.output()  # Not part of IRF or BMI: Does standalone plotting and file output

    #####################
    ## GET VALUES HERE ##
    ##   (if desired)  ##
    #####################
    # wout = obj.get_value('Deflection') # for example


if __name__ == "__main__":
    main()
