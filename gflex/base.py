"""
This file is part of gFlex.
gFlex computes lithospheric flexural isostasy with heterogeneous rigidity
Copyright (C) 2010-2026 Andrew D. Wickert

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

import contextlib
import os

import warnings
from enum import Enum

import numpy as np
from matplotlib import pyplot as plt

from ._version import __version__

# Scientific colour maps (optional but strongly recommended).
# Install with:  pip install cmcrameri
try:
    from cmcrameri import cm as _cmc
    _cmap_deflection = _cmc.vik        # blue → white → red; diverging, zero-centred
    _cmap_load = _cmc.lajolla_r       # pale cream (zero) → deep warm brown; heavier = darker
    _cmap_te   = _cmc.roma            # warm red (thin) → cool blue-green (thick)
    _HAVE_CRAMERI = True
except ImportError:
    _cmap_deflection = plt.cm.RdBu
    _cmap_load = plt.cm.YlOrBr
    _cmap_te   = plt.cm.RdBu_r        # warm red (thin) → cool blue (thick)
    _HAVE_CRAMERI = False
    warnings.warn(
        "cmcrameri is not installed — using matplotlib fallback colormaps for gFlex "
        "plots.\nFor perceptually uniform scientific colormaps (vik for deflection, "
        "lajolla for load), install it with:\n    pip install cmcrameri",
        UserWarning,
        stacklevel=2,
    )


class _RigidityBC(Enum):
    """Internal enum for flexural-rigidity boundary condition type."""
    PERIODIC = "periodic"
    MIRROR = "mirror_symmetry"
    ZERO_CURVATURE = "zero_curvature"


# ---------------------------------------------------------------------------
# Inhomogeneous (prescribed-value) BC support  (#52)
# ---------------------------------------------------------------------------

_VALID_BC_KEYS = frozenset({"displacement", "slope", "moment", "shear"})

# Maps the frozenset of two prescribed keys → canonical BC type string.
# Each entry corresponds to one existing homogeneous BC type whose stencil
# structure is re-used when nonzero values are prescribed.
_BC_DICT_TO_STRING = {
    frozenset({"displacement", "slope"}):  "zero_displacement_zero_slope",
    frozenset({"displacement", "moment"}): "zero_displacement_zero_moment",
    frozenset({"moment", "shear"}):        "zero_moment_zero_shear",
    frozenset({"slope", "shear"}):         "zero_slope_zero_shear",
}


def _resolve_bc(bc, edge_label, dimension):
    """Return the canonical BC type string for *bc* (str or dict).

    String presets are returned unchanged.  Dict BCs are validated and
    mapped to their equivalent canonical string so that downstream stencil
    and rigidity code never sees a raw dict.  Raises ``ValueError`` or
    ``TypeError`` with a descriptive message on any error.
    """
    if isinstance(bc, str):
        return bc

    if not isinstance(bc, dict):
        raise TypeError(
            f"Boundary condition at {edge_label!r} must be a string or dict, "
            f"got {type(bc).__name__!r}."
        )
    keys = frozenset(bc)
    unknown = keys - _VALID_BC_KEYS
    if unknown:
        raise ValueError(
            f"Boundary condition at {edge_label!r}: unknown key(s) "
            f"{sorted(unknown)}.  Valid keys: {sorted(_VALID_BC_KEYS)}."
        )
    if len(keys) != 2:
        raise ValueError(
            f"Boundary condition at {edge_label!r}: dict must have exactly 2 keys, "
            f"got {len(keys)} ({sorted(keys)})."
        )
    _ILL_POSED = {
        frozenset({"displacement", "shear"}): (
            "'displacement' and 'shear' are work-conjugate (both relate to "
            "transverse displacement): prescribing both on the same edge is "
            "ill-posed.  Pair 'displacement' with 'slope' or 'moment', or "
            "pair 'shear' with 'slope' or 'moment'."
        ),
        frozenset({"slope", "moment"}): (
            "'slope' and 'moment' are work-conjugate (both relate to "
            "rotation): prescribing both on the same edge is ill-posed.  "
            "Pair 'slope' with 'displacement' or 'shear', or pair 'moment' "
            "with 'displacement' or 'shear'."
        ),
    }
    if keys not in _BC_DICT_TO_STRING:
        if keys in _ILL_POSED:
            raise ValueError(
                f"Boundary condition at {edge_label!r}: {sorted(keys)!r} is "
                f"ill-posed.  {_ILL_POSED[keys]}"
            )
        raise ValueError(
            f"Boundary condition at {edge_label!r}: {sorted(keys)!r} is not a valid "
            f"pair.  Valid pairs: "
            f"{[sorted(p) for p in _BC_DICT_TO_STRING]}."
        )
    for k, v in bc.items():
        if isinstance(v, np.ndarray):
            if v.ndim != 1:
                raise ValueError(
                    f"Boundary condition at {edge_label!r}: value for {k!r} must be "
                    f"a scalar or 1-D array, got ndim={v.ndim}."
                )
        else:
            try:
                float(v)
            except (TypeError, ValueError):
                raise TypeError(
                    f"Boundary condition at {edge_label!r}: value for {k!r} must be "
                    f"a scalar or 1-D array, got {type(v).__name__!r}."
                )
    return _BC_DICT_TO_STRING[keys]


class Utility:

    """
    Generic utility functions
    """

    def configGet(
        self, vartype, category, name, optional=False, specialReturnMessage=None
    ):
        """
        Reads a configuration value from self.config (a nested dict loaded from YAML).

        vartype can be 'float', 'str' or 'string' (str and string are the same),
        'int' or 'integer' (also the same), or 'bool' or 'boolean'.

        optional=True: if the key is absent, return None instead of raising ValueError.

        specialReturnMessage: extra text appended to the error message on failure.
        """

        try:
            raw = self.config[category][name]
            if vartype == "float":
                var = float(raw)
            elif vartype == "string" or vartype == "str":
                var = "" if raw is None else str(raw)
                if var == "" and not optional:
                    # "" is acceptable for boundary conditions
                    if name[:18] != "boundary_condition":
                        if not self.quiet:
                            print(
                                "An empty input string here is not an acceptable"
                                " option."
                            )
                            print(name, "is not optional.")
                            print("Program crash likely to occur.")
            elif vartype == "integer" or vartype == "int":
                var = int(raw)
            elif vartype == "boolean" or vartype == "bool":
                if isinstance(raw, bool):
                    var = raw
                elif isinstance(raw, int):
                    var = bool(raw)
                else:
                    var = str(raw).lower() in ("true", "yes", "1", "on")
            else:
                raise ValueError(
                    "Please enter 'float', 'string' (or 'str'), 'integer' (or 'int'),"
                    " or 'boolean' (or 'bool') for vartype"
                )
            return var
        except Exception:
            if optional:
                # Carry on if the variable is optional
                var = None
                if self.verbose or self.debug:
                    if not self.grass:
                        print("")
                        print('No value entered for optional parameter "' + name + '"')
                        print('in category "' + category + '" in configuration file.')
                        print(
                            "No action related to this optional parameter will be taken."
                        )
                        print("")
            else:
                msg = (
                    f'Problem loading {vartype} "{name}" in category '
                    f'"{category}" from configuration file.'
                )
                if specialReturnMessage:
                    msg += f"\n{specialReturnMessage}"
                raise ValueError(msg)

    def _load_config(self, filename):
        """Read a YAML configuration file and return a nested dict.

        Only YAML (``.yaml`` / ``.yml``) configuration files are supported.
        The file must use the standard section names (``mode``, ``parameter``,
        ``input``, ``output``, ``numerical``, ``numerical2D``, ``verbosity``).
        """
        ext = os.path.splitext(filename)[1].lower()
        if ext not in (".yaml", ".yml"):
            raise ValueError(
                f"Configuration file '{filename}' does not have a .yaml or .yml "
                "extension. Only YAML configuration files are supported."
            )
        try:
            import yaml
        except ImportError:
            raise ImportError(
                "PyYAML is required to read YAML configuration files. "
                "Install it with: pip install pyyaml"
            )
        with open(filename) as fh:
            data = yaml.safe_load(fh)
        return {k: v for k, v in data.items() if isinstance(v, dict)}

    def readyCoeff(self):
        from scipy import sparse

        if sparse.issparse(self.coeff_matrix):
            pass  # Good type
        else:
            try:
                self.coeff_matrix = sparse.dia_matrix(self.coeff_matrix)
            except Exception:
                raise ValueError(
                    "Failed to make a sparse array or load a sparse matrix from the input."
                )

    def greatCircleDistance(self, lat1, long1, lat2, long2, radius):
        """
        Returns the great circle distance between two points.
        Useful when using the SAS_NG solution in lat/lon coordinates
        Modified from http://www.johndcook.com/blog/python_longitude_latitude/
        It should be able to take numpy arrays.
        """

        # Convert latitude and longitude to
        # spherical coordinates in radians.
        degrees_to_radians = np.pi / 180.0

        # theta = colatitude = 90 - latitude
        theta1rad = (90.0 - lat1) * degrees_to_radians
        theta2rad = (90.0 - lat2) * degrees_to_radians

        # lambda = longitude
        lambda1rad = long1 * degrees_to_radians
        lambda2rad = long2 * degrees_to_radians

        # Compute spherical distance from spherical coordinates.

        # For two locations in spherical coordinates
        # (1, theta, phi) and (1, theta, phi)
        # cosine( arc length ) =
        #    sin(theta) * sin(theta') * cos(theta-theta') + cos(phi) * cos(phi')
        # distance = radius * arc length

        cos_arc_length = np.sin(theta1rad) * np.sin(theta2rad) * np.cos(
            lambda1rad - lambda2rad
        ) + np.cos(theta1rad) * np.cos(theta2rad)
        arc = np.arccos(cos_arc_length)

        great_circle_distance = radius * arc

        return great_circle_distance

    def define_points_grid(self):
        """
        This is experimental code that could be used in the spatial_domain_no_grid
        section to build a grid of points on which to generate the solution.
        However, the current development plan (as of 27 Jan 2015) is to have the
        end user supply the list of points where they want a solution (and/or for
        it to be provided in a more automated way by GRASS GIS). But because this
        (untested) code may still be useful, it will remain as its own function
        here.
        It used to be in f2d.py.
        """
        # Grid making step
        # In this case, an output at different (x,y), e.g., on a grid, is desired
        # First, see if there is a need for a grid, and then make it
        # latlon arrays must have a pre-set grid
        if not self.latlon:
            # Warn that any existing grid will be overwritten
            try:
                self.dx
            except AttributeError:
                try:
                    self.dy
                except AttributeError:
                    pass
                else:
                    if not self.quiet:
                        print("dx and dy being overwritten -- supply a full grid")
            else:
                if not self.quiet:
                    print("dx and dy being overwritten -- supply a full grid")
            # Boundaries
            n = np.max(self.y) + self.alpha
            s = np.min(self.y) - self.alpha
            w = np.min(self.x) + self.alpha
            e = np.max(self.x) - self.alpha
            # Grid spacing
            dxprelim = self.alpha / 50.0  # x or y
            nx = np.ceil((e - w) / dxprelim)
            ny = np.ceil((n - s) / dxprelim)
            dx = (e - w) / nx
            dy = (n - s) / ny
            self.dx = self.dy = (dx + dy) / 2.0  # Average of these to create a
            # square grid for more compatibility
            self.xw = np.linspace(w, e, nx)
            self.yw = np.linspace(s, n, ny)
        else:
            raise ValueError(
                "Lat/lon xw and yw must be pre-set: grid will not be square "
                "and may run into issues with poles, so to ensure the proper "
                "output points are chosen, the end user should do this."
            )

    def loadFile(self, var, close_on_fail=True):
        """
        A special function to replace a variable name that is a string file path
        with the loaded file.
        var is a string on input
        output is a numpy array or a None-type object (success vs. failure)
        """
        out = None
        try:
            # First see if it is a full path or directly links from the current
            # working directory
            out = np.load(var)
        except Exception:
            try:
                out = np.loadtxt(var)
            except Exception:
                # Then see if it is relative to the location of the configuration file
                try:
                    out = np.load(self.inpath + var)
                except Exception:
                    try:
                        out = np.loadtxt(self.inpath + var)
                    # If failure
                    except Exception:
                        pass
                    else:
                        format_name = "ASCII"
                else:
                    format_name = "numpy binary"
            else:
                format_name = "ASCII"
        else:
            format_name = "numpy binary"

        if out is None and close_on_fail:
            raise FileNotFoundError(
                f"Cannot find {var} file. "
                f"Looked relative to model Python files and relative to "
                f"configuration file path {self.inpath!r}."
            )
        elif out is not None and self.verbose:
            print(f"Loading {var} from {format_name}")

        return out


class Plotting:
    """
    Mixin that provides quick-look plotting for F1D and F2D results.

    Methods are called by :meth:`Flexure.output`; they are not intended to
    be called directly by users.  Set ``self.plot_choice`` before calling
    :meth:`output` to trigger a plot.

    .. note::
        Plot, if desired.
        1D all here, 2D in functions — just because there is often more code
        in 2D plotting functions.  Also, yes, this portion of the code is NOT
        efficient or elegant in how it handles functions.  But it's just a
        simple way to visualize results easily!  And not too hard to improve
        with a bit of time.  Anyway, the main goal here is the visualization,
        not the beauty of the code :)
    """

    def plotting(self):
        """
        Dispatch to the appropriate plot routine based on ``self.plot_choice``.

        Valid values for ``plot_choice``:

        * ``'q'`` — surface load only
        * ``'w'`` — deflection only
        * ``'both'`` — load and deflection in separate subplots (1D and 2D)
        * ``'combo'`` — load and deflection overlaid (1D only)
        """
        # try:
        #  self.plot_choice
        # except:
        #  self.plot_choice = None
        if self.plot_choice:
            if self.verbose:
                print("Starting to plot " + self.plot_choice)
            if self.dimension == 1:
                if self.plot_choice == "q":
                    plt.figure(1)
                    if self.method == "sas_ng":
                        plt.plot(self.x / 1000.0, self.q / (self.rho_m * self.g), "ko-")
                        plt.ylabel(
                            "Load volume, mantle equivalent [m$^3$]",
                            fontsize=12,
                            fontweight="bold",
                        )
                    else:
                        plt.plot(self.x / 1000.0, self.qs / (self.rho_m * self.g), "k-")
                        plt.ylabel(
                            "Load thickness, mantle equivalent [km]",
                            fontsize=12,
                            fontweight="bold",
                        )
                    plt.xlabel(
                        "Distance along profile [km]", fontsize=12, fontweight="bold"
                    )
                    plt.tight_layout()
                    plt.show()
                elif self.plot_choice == "w":
                    plt.figure(1)
                    if self.method == "sas_ng":
                        plt.plot(self.xw / 1000.0, self.w, "k-")
                    else:
                        plt.plot(self.x / 1000.0, self.w, "k-")
                    plt.ylabel("Deflection [m]", fontsize=12, fontweight="bold")
                    plt.xlabel(
                        "Distance along profile [km]", fontsize=12, fontweight="bold"
                    )
                    plt.tight_layout()
                    plt.show()
                elif self.plot_choice == "both":
                    plt.figure(1, figsize=(6, 9))
                    ax = plt.subplot(212)
                    if self.method == "sas_ng":
                        ax.plot(self.xw / 1000.0, self.w, "k-")
                    else:
                        ax.plot(self.x / 1000.0, self.w, "k-")
                    ax.set_ylabel("Deflection [m]", fontsize=12, fontweight="bold")
                    ax.set_xlabel(
                        "Distance along profile [m]", fontsize=12, fontweight="bold"
                    )
                    plt.subplot(211)
                    plt.title("Loads and Lithospheric Deflections", fontsize=16)
                    if self.method == "sas_ng":
                        plt.plot(self.x / 1000.0, self.q / (self.rho_m * self.g), "ko-")
                        plt.ylabel(
                            "Load volume, mantle equivalent [m$^3$]",
                            fontsize=12,
                            fontweight="bold",
                        )
                        plt.xlim(ax.get_xlim())
                    else:
                        plt.plot(self.x / 1000.0, self.qs / (self.rho_m * self.g), "k-")
                        plt.ylabel(
                            "Load thickness, mantle equivalent [m]",
                            fontsize=12,
                            fontweight="bold",
                        )
                    plt.xlabel(
                        "Distance along profile [km]", fontsize=12, fontweight="bold"
                    )
                    plt.tight_layout()
                    plt.show()
                elif self.plot_choice == "combo":
                    fig = plt.figure(1, figsize=(10, 6))
                    titletext = "Loads and Lithospheric Deflections"
                    ax = fig.add_subplot(1, 1, 1)
                    # Plot undeflected load
                    if self.method == "sas_ng":
                        if not self.quiet:
                            print(
                                "Combo plot can't work with SAS_NG! Don't have mechanism"
                                " in place to calculate load width."
                            )
                            print(
                                "Big problem -- what is the area represented by the loads"
                                " at the extreme ends of the array?"
                            )
                    else:
                        ax.plot(
                            self.x / 1000.0,
                            self.qs / (self.rho_m * self.g),
                            "g--",
                            linewidth=2,
                            label="Load thickness [m mantle equivalent]",
                        )
                    # Plot deflected load
                    if self.method == "sas_ng":
                        pass
                        # ax.plot(
                        #     self.x / 1000.0,
                        #     self.q / (self.rho_m * self.g) + self.w,
                        #     "go-",
                        #     linewidth=2,
                        #     label="Load volume [m^3] mantle equivalent]",
                        # )
                    else:
                        ax.plot(
                            self.x / 1000.0,
                            self.qs / (self.rho_m * self.g) + self.w,
                            "g-",
                            linewidth=2,
                            label="Deflection [m] + load thickness [m mantle equivalent]",
                        )
                    # Plot deflection
                    if self.method == "sas_ng":
                        ax.plot(
                            self.xw / 1000.0,
                            self.w,
                            "ko-",
                            linewidth=2,
                            label="Deflection [m]",
                        )
                    else:
                        ax.plot(
                            self.x / 1000.0,
                            self.w,
                            "k-",
                            linewidth=2,
                            label="Deflection [m]",
                        )
                    # Set y min to equal to the absolute value maximum of y max and y min
                    # (and therefore show isostasy better)
                    yabsmax = max(abs(np.array(plt.ylim())))
                    # Y axis label
                    plt.ylim((-yabsmax, yabsmax))
                    # Plot title selector -- be infomrative
                    try:
                        self.te
                        if self.method == "fd":
                            if type(self.te) is np.ndarray:
                                if (self.te != (self.te).mean()).any():
                                    plt.title(titletext, fontsize=16)
                                else:
                                    plt.title(
                                        titletext
                                        + ", $T_e$ = "
                                        + str((self.te / 1000).mean())
                                        + " km",
                                        fontsize=16,
                                    )
                            else:
                                plt.title(
                                    titletext
                                    + ", $T_e$ = "
                                    + str(self.te / 1000)
                                    + " km",
                                    fontsize=16,
                                )
                        else:
                            plt.title(
                                titletext + ", $T_e$ = " + str(self.te / 1000) + " km",
                                fontsize=16,
                            )
                    except AttributeError:
                        plt.title(titletext, fontsize=16)
                    # x and y labels
                    plt.ylabel("Loads and flexural response [m]", fontsize=16)
                    plt.xlabel("Distance along profile [km]", fontsize=16)
                    # legend -- based on lables
                    plt.legend(loc=0, numpoints=1, fancybox=True)
                    plt.tight_layout()
                    plt.show()
                else:
                    if not self.quiet:
                        print(
                            'Incorrect plot_choice input, "'
                            + self.plot_choice
                            + '" provided.'
                        )
                        print(
                            "Possible input strings are: q, w, both, and (for 1D) combo"
                        )
                        print("Unable to produce plot.")
            elif self.dimension == 2:
                if self.plot_choice == "q":
                    fig = plt.figure(1, figsize=(8, 6))
                    if self.method != "sas_ng":
                        self.surfplot(
                            self.qs / (self.rho_m * self.g),
                            "Load thickness, mantle equivalent [m]",
                            cmap=_cmap_load,
                        )
                        plt.show()
                    else:
                        self.xyzinterp(
                            self.x,
                            self.y,
                            self.q,
                            "Load volume, mantle equivalent [m$^3$]",
                            cmap=_cmap_load,
                        )
                    plt.tight_layout()
                    plt.show()
                elif self.plot_choice == "w":
                    fig = plt.figure(1, figsize=(8, 6))
                    w_abs = float(np.abs(self.w).max())
                    if self.method != "sas_ng":
                        self.surfplot(self.w, "Deflection [m]",
                                      cmap=_cmap_deflection, vmin=-w_abs, vmax=w_abs)
                        plt.show()
                    else:
                        self.xyzinterp(self.xw, self.yw, self.w, "Deflection [m]",
                                       cmap=_cmap_deflection, vmin=-w_abs, vmax=w_abs)
                    plt.tight_layout()
                    plt.show()
                elif self.plot_choice == "both":
                    plt.figure(1, figsize=(6, 9))
                    if self.method != "sas_ng":
                        self.twoSurfplots()
                        plt.show()
                    else:
                        w_abs = float(np.abs(self.w).max())
                        plt.subplot(211)
                        self.xyzinterp(
                            self.x,
                            self.y,
                            self.q,
                            "Load volume, mantle equivalent [m$^3$]",
                            cmap=_cmap_load,
                        )
                        plt.subplot(212)
                        self.xyzinterp(self.xw, self.yw, self.w, "Deflection [m]",
                                       cmap=_cmap_deflection, vmin=-w_abs, vmax=w_abs)
                        plt.tight_layout()
                        plt.show()
                else:
                    if not self.quiet:
                        print(
                            'Incorrect plot_choice input, "'
                            + self.plot_choice
                            + '" provided.'
                        )
                        print(
                            "Possible input strings are: q, w, both, and (for 1D) combo"
                        )
                        print("Unable to produce plot.")

    def surfplot(self, z, titletext, cmap=None, vmin=None, vmax=None):
        """
        Plot if you want to - for troubleshooting - 1 figure
        """
        if cmap is None:
            cmap = _cmap_load
        if self.latlon:
            plt.imshow(
                z, extent=(0, self.dx * z.shape[1], self.dy * z.shape[0], 0),
                cmap=cmap, vmin=vmin, vmax=vmax,
            )  # ,interpolation='nearest'
            plt.xlabel("longitude [deg E]", fontsize=12, fontweight="bold")
            plt.ylabel("latitude [deg N]", fontsize=12, fontweight="bold")
        else:
            plt.imshow(
                z,
                extent=(
                    0,
                    self.dx / 1000.0 * z.shape[1],
                    self.dy / 1000.0 * z.shape[0],
                    0,
                ),
                cmap=cmap, vmin=vmin, vmax=vmax,
            )  # ,interpolation='nearest'
            plt.xlabel("x [km]", fontsize=12, fontweight="bold")
            plt.ylabel("y [km]", fontsize=12, fontweight="bold")
        plt.colorbar()

        plt.title(titletext, fontsize=16)

    def twoSurfplots(self):
        """
        Plot multiple subplot figure for 2D array
        """
        # Could more elegantly just call surfplot twice
        # And also could include xyzinterp as an option inside surfplot.
        # Noted here in case anyone wants to take that on in the future...

        from matplotlib.colors import TwoSlopeNorm

        w_min = float(self.w.min())
        w_max = float(self.w.max())
        w_abs = max(abs(w_min), abs(w_max))

        # When forebulge << subsidence the forebulge is invisible with a
        # symmetric clim.  TwoSlopeNorm keeps white at zero while giving the
        # positive arm enough range to show the forebulge.
        if w_max > 0 and w_max < 0.2 * abs(w_min):
            _defl_norm = TwoSlopeNorm(
                vmin=w_min, vcenter=0.0,
                vmax=max(w_max, 0.1 * abs(w_min)),
            )
            _defl_vmin = _defl_vmax = None
        else:
            _defl_norm = None
            _defl_vmin, _defl_vmax = -w_abs, w_abs

        _has_te_grid = isinstance(self.te, np.ndarray) and self.te.ndim == 2

        xlabel = "longitude [deg E]" if self.latlon else "x [km]"
        ylabel = "latitude [deg N]" if self.latlon else "y [km]"

        def _ext(arr):
            if self.latlon:
                return (0, self.dx * arr.shape[1], self.dy * arr.shape[0], 0)
            return (0, self.dx / 1000. * arr.shape[1],
                       self.dy / 1000. * arr.shape[0], 0)

        if _has_te_grid:
            # Three-panel horizontal layout: load | Te | deflection
            plt.gcf().set_size_inches(14, 4.5)
            ax_load = plt.subplot(1, 3, 1)
            ax_te   = plt.subplot(1, 3, 2)
            ax_defl = plt.subplot(1, 3, 3)
        else:
            ax_load = plt.subplot(2, 1, 1)
            ax_te   = None
            ax_defl = plt.subplot(2, 1, 2)

        # Load panel
        ax_load.set_title("Load thickness, mantle equivalent [m]", fontsize=16)
        im_load = ax_load.imshow(
            self.qs / (self.rho_m * self.g), extent=_ext(self.qs), cmap=_cmap_load,
        )
        ax_load.set_xlabel(xlabel, fontsize=12, fontweight="bold")
        ax_load.set_ylabel(ylabel, fontsize=12, fontweight="bold")
        plt.colorbar(im_load, ax=ax_load)

        # Elastic thickness panel (only when Te is a 2-D array)
        if _has_te_grid:
            ax_te.set_title(r"Elastic thickness $T_e$ [km]", fontsize=16)
            im_te = ax_te.imshow(
                self.te / 1e3, extent=_ext(self.te), cmap=_cmap_te,
            )
            ax_te.set_xlabel(xlabel, fontsize=12, fontweight="bold")
            ax_te.set_ylabel(ylabel, fontsize=12, fontweight="bold")
            plt.colorbar(im_te, ax=ax_te)

        # Deflection panel
        ax_defl.set_title("Deflection [m]", fontsize=16)
        im_defl = ax_defl.imshow(
            self.w, extent=_ext(self.w), cmap=_cmap_deflection,
            norm=_defl_norm, vmin=_defl_vmin, vmax=_defl_vmax,
        )
        ax_defl.set_xlabel(xlabel, fontsize=12, fontweight="bold")
        ax_defl.set_ylabel(ylabel, fontsize=12, fontweight="bold")
        cb = plt.colorbar(im_defl, ax=ax_defl)
        if _defl_norm is not None:
            step = 10 ** int(np.log10(abs(w_min)))
            neg_ticks = list(range(0, int(w_min) - 1, -step))
            vmax_norm = _defl_norm.vmax
            pos_ticks = list(range(10, int(vmax_norm), 10))
            cb.set_ticks(sorted(neg_ticks + pos_ticks))

        plt.tight_layout()

    def xyzinterp(self, x, y, z, titletext, cmap=None, vmin=None, vmax=None):
        """
        Interpolates and plots ungridded model outputs from SAS_NG solution
        """
        # Help from http://wiki.scipy.org/Cookbook/Matplotlib/Gridding_irregularly_spaced_data

        if self.verbose:
            print("Starting to interpolate grid for plotting -- can be a slow process!")

        from scipy.interpolate import griddata

        # define grid.
        xmin = np.min(self.xw)
        # xmean = np.mean(self.xw)  # not used right now
        xmax = np.max(self.xw)
        ymin = np.min(self.yw)
        # ymean = np.mean(self.yw)  # not used right now
        ymax = np.max(self.yw)
        # x_range = xmax - xmin
        # y_range = ymax - ymin

        # x and y grids
        # 100 cells on each side -- just for plotting, not so important
        # to optimize with how many points are plotted
        # xi = np.linspace(xmin-.05*x_range, xmax+.05*x_range, 200)
        # yi = np.linspace(ymin-.05*y_range, ymax+.05*y_range, 200)
        xi = np.linspace(xmin, xmax, 200)
        yi = np.linspace(ymin, ymax, 200)
        # grid the z-axis
        zi = griddata((x, y), z, (xi[None, :], yi[:, None]), method="cubic")
        # turn nan into 0 -- this will just be outside computation area for q
        zi[np.isnan(zi)] = 0
        # contour the gridded outputs, plotting dots at the randomly spaced data points.
        # CS = plt.contour(xi,yi,zi,15,linewidths=0.5,colors='k') -- don't need lines
        if cmap is None:
            cmap = _cmap_load
        if self.latlon:
            plt.contourf(xi, yi, zi, 100, cmap=cmap, vmin=vmin, vmax=vmax)
        else:
            plt.contourf(xi / 1000.0, yi / 1000.0, zi, 100, cmap=cmap, vmin=vmin, vmax=vmax)
        plt.colorbar()  # draw colorbar
        # plot model points.
        # Computed at
        if self.latlon:
            plt.plot(
                x, y, "o", markerfacecolor=".6", markeredgecolor=".6", markersize=1
            )
            plt.plot(
                self.x,
                self.y,
                "o",
                markerfacecolor=".2",
                markeredgecolor=".2",
                markersize=1,
            )
        else:
            plt.plot(
                x / 1000.0,
                y / 1000.0,
                "o",
                markerfacecolor=".6",
                markeredgecolor=".6",
                markersize=1,
            )
            # Load sources (overlay computed at)
            plt.plot(
                self.x / 1000.0,
                self.y / 1000.0,
                "o",
                markerfacecolor=".2",
                markeredgecolor=".2",
                markersize=1,
            )
        if self.latlon:
            plt.xlabel("longitude [deg E]", fontsize=12, fontweight="bold")
            plt.ylabel("latitude [deg N]", fontsize=12, fontweight="bold")
        else:
            plt.xlabel("x [km]", fontsize=12, fontweight="bold")
            plt.ylabel("y [km]", fontsize=12, fontweight="bold")
        # Limits -- to not get messed up by points (view wants to be wider so whole
        # point visible)
        if self.latlon:
            plt.xlim((xi[0], xi[-1]))
            plt.ylim((yi[0], yi[-1]))
        else:
            plt.xlim((xi[0] / 1000.0, xi[-1] / 1000.0))
            plt.ylim((yi[0] / 1000.0, yi[-1] / 1000.0))
        # Title
        plt.title(titletext, fontsize=16)


class WhichModel(Utility):
    """
    Reads a configuration file to determine which model class (F1D or F2D)
    to instantiate.  Used internally by the file-driven workflow; not needed
    when setting attributes programmatically.
    """

    def __init__(self, filename=None):
        """
        Read ``dimension`` from the ``[mode]`` section of the configuration
        file so the calling code knows whether to create an :class:`F1D` or
        :class:`F2D` instance.
        """
        self.filename = filename
        if self.filename:
            try:
                # only let this function imoprt things once
                self.whichModel_AlreadyRun
            except AttributeError:
                # Open parser and get what kind of model
                try:
                    self.config = self._load_config(filename)
                    # Need to change this and all slashes to be Windows compatible
                    self.inpath = os.path.dirname(os.path.realpath(filename)) + "/"
                    # Need to have these guys inside "try" to make sure it is set up OK
                    # (at least for them)
                    self.dimension = self.configGet("integer", "mode", "dimension")
                    self.whichModel_AlreadyRun = True
                except Exception:
                    raise FileNotFoundError(
                        "Cannot locate specified configuration file."
                    )


class Flexure(Utility, Plotting):
    """
    Solves flexural isostasy both analytically (for constant flexural rigidity)
    and numerically (for either variable or constant flexural rigidity).

    Analytical solutions are by superposition of analytical solutions
    in the spatial domain (i.e. a sum of Green's functions)

    Numerical solutions are finite difference by a direct sparse matrix solver.
    """

    def __init__(self, filename=None):
        """Set verbosity defaults and record the optional config filename.

        All parameter loading and solver setup is deferred to :meth:`initialize`
        so that attributes can be set programmatically between construction and
        the initialize–run–finalize call sequence.
        """
        # 17 Nov 2014: Splitting out initialize from __init__ to allow space
        # to use getters and setters to define values

        # Use standard routine to pull out values
        # If no filename provided, will not initialize configuration file.
        self.filename = filename

        # DEFAULT VERBOSITY
        # Set default "quiet" to False, unless set by setter or overwritten by
        # the configuration file.
        self.quiet = False
        # And also set default verbosity
        self.verbose = True
        self.debug = False

        # x and y to None for checks
        self.x = None
        self.y = None

        # Set GRASS GIS usage flag: if GRASS is used, don't display error
        # messages related to unset options. This sets it to False if it
        # hasn't already been set (and it can be set after this too)
        # (Though since this is __init__, would have to go through WhichModel
        # for some reason to define self.grass before this
        try:
            self.grass
        except AttributeError:
            self.grass = False

        # Default solver; may be overridden programmatically or via config file
        self.solver = "direct"

        # Backing store for the method property (see below)
        self._method = None

        # Default values for lat/lon usage -- defaulting not to use it
        try:
            self.latlon
        except AttributeError:
            self.latlon = False
        try:
            self.planetary_radius
        except AttributeError:
            self.planetary_radius = None

    @property
    def method(self):
        """Solution method: ``'fd'``, ``'fft'``, ``'sas'``, or ``'sas_ng'``."""
        return self._method

    @method.setter
    def method(self, value):
        if isinstance(value, str):
            lo = value.lower()
            if value != lo:
                import warnings
                warnings.warn(
                    f"method='{value}' is deprecated and will likely be removed in v2.0; use '{lo}' instead.",
                    DeprecationWarning,
                    stacklevel=2,
                )
            self._method = lo
        else:
            self._method = value

    def initialize(self, filename=None):
        """
        Load parameters and prepare internal state.

        If a configuration file path is available, reads all physical and
        numerical parameters from it (YAML format).  Otherwise expects the
        caller to have set attributes directly (programmatic use).  Called
        automatically by :meth:`F1D.initialize` and :meth:`F2D.initialize`.

        Parameters
        ----------
        filename : str, optional
            Path to a gFlex YAML configuration file (``.yaml`` / ``.yml``).
        """
        # Values from configuration file

        # If a filename is provided here, overwrite any prior value
        if filename:
            if self.filename:
                pass  # Don't overwrite if filename is None-type
                # "debug" not yet defined.
                # if self.debug:
                #  print("Overwriting filename from '__init__' step with that from\n"+\
                #        "initialize step."
            else:
                # Update to new filename
                self.filename = filename

        if self.filename:
            try:
                self.config = self._load_config(self.filename)
                self.inpath = os.path.dirname(os.path.realpath(self.filename)) + "/"
                # Need to have these guys inside "try" to make sure it is set up OK
                # (at least for them)
                self.dimension = self.configGet("integer", "mode", "dimension")
                self.whichModel_AlreadyRun = True
            except Exception:
                raise FileNotFoundError(
                    "No configuration file at specified path, or configuration file"
                    " configured incorrectly"
                )

            # Set verbosity for model run
            # Default is "verbose" with no debug or quiet
            # Verbose
            self.verbose = self.configGet(
                "bool", "verbosity", "verbose", optional=True
            ) or self.verbose
            # Debug means that whole arrays, etc., can be printed
            self.debug = self.configGet(
                "bool", "verbosity", "debug", optional=True
            ) or self.debug
            # Quiet suppresses all output
            self.quiet = self.configGet(
                "bool", "verbosity", "quiet", optional=True
            ) or self.quiet
        # Quiet overrides all others
        if self.quiet:
            self.debug = False
            self.verbose = False

        # Introduce model
        # After configuration file can define "quiet", and getter/setter should be done
        # by this point if we are going that way.
        if not self.quiet:
            print("")  # Blank line at start of run
            print("")
            print("****************************" + "*" * len(__version__))
            print("*** Initializing gFlex v" + __version__ + " ***")
            print("****************************" + "*" * len(__version__))
            print("")
            print("Open-source licensed under GNU GPL v3")
            print("")

        if self.filename:
            # Set clocks to None so if they are called by the getter before the
            # calculation is performed, there won't be an error
            self.coeff_creation_time = None
            self.time_to_solve = None

            self.method = self.configGet("string", "mode", "method")
            # Boundary conditions
            # This used to be nested inside an "if self.method == 'FD'", but it seems
            # better to define these to ensure there aren't mistaken impressions
            # about what they do for the SAS case
            # Not optional: flexural solutions can be very sensitive to b.c.'s
            self.bc_east = self.configGet(
                "string", "numerical", "boundary_condition_east", optional=False
            )
            self.bc_west = self.configGet(
                "string", "numerical", "boundary_condition_west", optional=False
            )
            if self.dimension == 2:
                self.bc_north = self.configGet(
                    "string", "numerical2D", "boundary_condition_north", optional=False
                )
                self.bc_south = self.configGet(
                    "string", "numerical2D", "boundary_condition_south", optional=False
                )

            # Parameters
            self.g = self.configGet("float", "parameter", "gravitational_acceleration")
            self.rho_m = self.configGet("float", "parameter", "mantle_density")
            self.rho_fill = self.configGet(
                "float", "parameter", "infill_material_density"
            )

            # Grid spacing
            if self.method != "sas_ng":
                # No meaning for ungridded superimposed analytical solutions
                # From configuration file
                self.dx = self.configGet("float", "numerical", "grid_spacing_x")
                if self.dimension == 2:
                    self.dy = self.configGet("float", "numerical2D", "grid_spacing_y")

            # Loading grid
            # q0 is either a load array or an x,y,q array.
            # Therefore q_0, initial q, before figuring out what it really is
            # for grid, q0 could also be written as $q_\sigma$ or q/(dx*(dy))
            # it is a surface normal stress that is h_load * rho_load * g
            # it later is combined with dx and (if 2D) dy for FD cases
            # for point loads, need mass: q0 should be written as [x, (y), force])
            self.q0 = self.configGet("string", "input", "loads")

        # Parameters -- rho_m and rho_fill defined, so this outside
        # of if-statement (to work with getters/setters as well)
        self.drho = self.rho_m - self.rho_fill
        if self.filename:
            self.E = self.configGet("float", "parameter", "youngs_modulus")
            self.nu = self.configGet("float", "parameter", "poissons_ratio")

        # Stop program if there is no q0 defined or if it is None-type
        try:
            self.q0
        except AttributeError:
            try:
                self.q
            except AttributeError:
                try:
                    self.qs
                except AttributeError:
                    raise RuntimeError(
                        "Must define q0, q, or qs by this stage in the initialization"
                        " step from either configuration file (string) or direct array"
                        " import."
                    )
        else:
            # Stop program if q0 is None-type
            if self.q0 is None:  # if is None type, just be patient
                raise RuntimeError(
                    "Must define non-None-type q0 by this stage in the initialization"
                    " step from either configuration file (string) or direct array"
                    " import."
                )

        # Ignore this if no q0 set
        try:
            self.q0
        except AttributeError:
            self.q0 = None
        if self.q0 == "":
            self.q0 = None
        if isinstance(self.q0, str):
            self.q0 = self.loadFile(self.q0)  # Won't do this if q0 is None

        # Check consistency of dimensions
        if self.q0 is not None:
            if self.method != "sas_ng":
                if self.q0.ndim != self.dimension:
                    raise ValueError(
                        f"Number of dimensions in loads array ({self.q0.ndim}D) is "
                        f"inconsistent with the solution dimension ({self.dimension}D)."
                    )

        # Plotting selection
        self.plot_choice = self.configGet("string", "output", "plot", optional=True)

        # Ensure that Te is of floating-point type to avoid integer math
        # and floor division
        try:
            self.te = self.te.astype(float)  # array
        except AttributeError:
            # Integer scalar Te does not seem to be a problem, but taking this step
            # anyway for consistency
            try:
                self.te = float(self.te)  # integer
            except (AttributeError, ValueError, TypeError):
                # If not already defined, then an input file is being used, and this
                # code should bring the grid in as floating point type... just later.
                pass

        # Check for end loads; otherwise set as 0
        # Do this for 2D; in the 1D case, xy and yy will just not be used
        try:
            self.sigma_xx
        except AttributeError:
            self.sigma_xx = 0
        else:
            if self.method not in ("fd", "fft"):
                warnings.warn(
                    "End loads have been set but will not be implemented because the"
                    " solution method is not finite difference or FFT",
                    category=RuntimeWarning,
                    stacklevel=2,
                )
        try:
            self.sigma_xy
        except AttributeError:
            self.sigma_xy = 0
        else:
            if self.method not in ("fd", "fft"):
                warnings.warn(
                    "End loads have been set but will not be implemented because the"
                    " solution method is not finite difference or FFT",
                    category=RuntimeWarning,
                    stacklevel=2,
                )
        try:
            self.sigma_yy
        except AttributeError:
            self.sigma_yy = 0
        else:
            if self.method not in ("fd", "fft"):
                warnings.warn(
                    "End loads have been set but will not be implemented because the"
                    " solution method is not finite difference or FFT",
                    category=RuntimeWarning,
                    stacklevel=2,
                )

    # Finalize
    def finalize(self):
        """
        Release cached solver state.

        Deletes the coefficient matrix so it is rebuilt fresh on the next
        call, which matters when running iteratively with changing inputs.
        Called automatically by :meth:`F1D.finalize` and
        :meth:`F2D.finalize`.
        """
        # Can include an option for this later, but for the moment, this will
        # clear the coefficient array so it doens't cause problems for model runs
        # searching for the proper rigidity
        with contextlib.suppress(AttributeError):
            del self.coeff_matrix
        if not self.quiet:
            print("")

    # SAVING TO FILE AND PLOTTING STEPS

    # Output: One of the functions run by isostasy.py; not part of IRF
    # (for standalone model use)
    def output(self):
        """
        Save deflection to file and/or plot, based on optional attributes.

        Does nothing if neither ``w_out_file`` nor ``plot_choice`` has been
        set.  Set ``w_out_file`` to a path ending in ``'.npy'`` for a binary
        NumPy array, or any other extension for an ASCII grid.  Set
        ``plot_choice`` to ``'q'``, ``'w'``, ``'both'``, or (1D only)
        ``'combo'`` to display plots.
        """
        if self.verbose:
            print("Output step")
        self.outputDeflections()
        self.plotting()

    # Save output deflections to file, if desired
    def outputDeflections(self):
        """
        Outputs a grid of deflections if an output directory is defined in the
        configuration file

        If the filename given in the configuration file ends in ".npy", then a binary
        numpy grid will be exported.

        Otherwise, an ASCII grid will be exported.
        """
        try:
            # If w_out_file exists, has already been set by a setter
            self.w_out_file
        except AttributeError:
            # Otherwise, set from config; None means no output file
            self.w_out_file = self.configGet(
                "string", "output", "deflection_out", optional=True
            )
        else:
            if self.verbose:
                print("Output filename provided.")
        if self.w_out_file:
            if self.w_out_file[-4:] == ".npy":
                from numpy import save

                save(self.w_out_file, self.w)
            else:
                from numpy import savetxt

                # Shouldn't need more than mm precision, at very most
                savetxt(self.w_out_file, self.w, fmt="%.3f")
                if self.verbose:
                    print("Saving deflections --> " + self.w_out_file)

    def bc_check(self):
        """
        Validate boundary conditions and prepare BC state for the chosen solver.

        For ``FD``: checks that every edge BC is one of the five accepted
        strings (``'zero_displacement_zero_slope'``, ``'zero_moment_zero_shear'``,
        ``'zero_slope_zero_shear'``, ``'mirror'``, ``'periodic'``) and exits with an
        informative message if not.  A pre-built ``coeff_matrix`` bypasses
        the check (coupled-model use case).

        For ``FFT``: ensures the BC attributes exist; the FFT solver
        interprets ``'periodic'`` on all sides as an exact periodic transform
        and treats any other value as a request for zero-padded
        ``no_outside_loads`` behaviour.

        For ``SAS`` / ``SAS_NG``: sets missing BC attributes to an empty
        string, which the analytical solvers treat as ``'no_outside_loads'``.
        """
        # Reject dict-style BCs on non-FD solvers before any method-specific work.
        if self.method != "fd":
            dict_edges = [
                name
                for name, attr in (
                    ("bc_west",  getattr(self, "bc_west",  None)),
                    ("bc_east",  getattr(self, "bc_east",  None)),
                    ("bc_north", getattr(self, "bc_north", None)),
                    ("bc_south", getattr(self, "bc_south", None)),
                )
                if isinstance(attr, dict)
            ]
            if dict_edges:
                raise ValueError(
                    f"Inhomogeneous (dict-style) boundary conditions require "
                    f"method='fd'. The {self.method!r} solver assumes periodic or "
                    f"zero-padded boundaries and cannot apply prescribed moments, "
                    f"shears, displacements, or slopes. "
                    f"Affected edges: {dict_edges}."
                )

        # Check that boundary conditions are acceptable with code implementation
        # Acceptable b.c.'s
        if self.method == "fft":
            # Ensure BC attributes exist; FFT handles them internally
            # 'periodic' → exact transform; anything else → zero-padded (no_outside_loads)
            for attr in ("BC_E", "BC_W"):
                if not hasattr(self, attr):
                    setattr(self, attr, "")
            if self.dimension == 2:
                for attr in ("BC_N", "BC_S"):
                    if not hasattr(self, attr):
                        setattr(self, attr, "")
            else:
                self.bc_south = None
                self.bc_north = None
        elif self.method == "fd":
            # Check if a coefficient array has been defined
            # It would only be by a getter or setter;
            # no way to do I/O with this with present configuration files
            # Define as None for use later.
            try:
                self.coeff_matrix
            except AttributeError:
                self.coeff_matrix = None

            # Seed _bc_*_norm to the raw user-set values so that coupled-mode
            # callers (coeff_matrix already provided) get working defaults even
            # if the validation block below is skipped.
            self._bc_west_norm = self.bc_west
            self._bc_east_norm = self.bc_east
            # _bc_*_values stores the user-supplied dict for dict BCs (used by
            # the inhomogeneous RHS assembly in #51); None for string presets.
            self._bc_west_values = self.bc_west if isinstance(self.bc_west, dict) else None
            self._bc_east_values = self.bc_east if isinstance(self.bc_east, dict) else None
            if self.dimension == 2:
                self._bc_north_norm = self.bc_north
                self._bc_south_norm = self.bc_south
                self._bc_north_values = self.bc_north if isinstance(self.bc_north, dict) else None
                self._bc_south_values = self.bc_south if isinstance(self.bc_south, dict) else None

            # No need to create a coeff_matrix if one already exists
            if self.coeff_matrix is None:
                # Acceptable string boundary conditions
                self.bc1D = np.array(
                    [
                        "zero_displacement_zero_slope",
                        "zero_displacement_zero_moment",
                        "periodic",
                        "mirror",
                        "zero_moment_zero_shear",
                        "zero_slope_zero_shear",
                        "sandbox",
                    ]
                )
                self.bc2D = np.array(
                    [
                        "zero_displacement_zero_slope",
                        "zero_displacement_zero_moment",
                        "periodic",
                        "mirror",
                        "zero_moment_zero_shear",
                        "zero_slope_zero_shear",
                    ]
                )
                # Boundary conditions should be defined by this point -- whether via
                # the configuration file or the getters and setters.
                # Validate each edge BC and compute its canonical type string
                # (_bc_*_norm).  Dict BCs are mapped to their equivalent string;
                # string presets are validated against bc1D / bc2D and kept as-is.
                edge_labels = ["east", "west"]
                bcs = [self.bc_east, self.bc_west]
                if self.dimension == 2:
                    edge_labels += ["north", "south"]
                    bcs += [self.bc_north, self.bc_south]
                for bc, edge in zip(bcs, edge_labels):
                    if isinstance(bc, dict):
                        norm = _resolve_bc(bc, edge, self.dimension)
                    elif self.dimension == 1:
                        if bc not in self.bc1D:
                            raise ValueError(
                                f"{bc!r} is not an acceptable 1D finite difference "
                                "boundary condition. Acceptable boundary conditions are: "
                                f"{', '.join(repr(b) for b in self.bc1D)}."
                            )
                        norm = bc
                    elif self.dimension == 2:
                        if bc not in self.bc2D:
                            raise ValueError(
                                f"{bc!r} is not an acceptable 2D finite difference "
                                "boundary condition. Acceptable boundary conditions are: "
                                f"{', '.join(repr(b) for b in self.bc2D)}."
                            )
                        norm = bc
                    else:
                        raise ValueError(
                            "For a flexural solution, grid must be 1D or 2D."
                        )
                    if norm == "zero_slope_zero_shear":
                        norm = "mirror"
                    setattr(self, f"_bc_{edge}_norm", norm)
                # Validate array BC value lengths against edge dimensions.
                # Only when dict BCs are present (config-file paths use string BCs
                # and qs may not be set yet at this point).
                if self.dimension == 2:
                    _any_dict = any(
                        getattr(self, f"_bc_{e}_values", None) is not None
                        for e in ("west", "east", "north", "south")
                    )
                    if _any_dict:
                        ny, nx = self.qs.shape
                        _edge_lengths = {"west": ny, "east": ny, "north": nx, "south": nx}
                        for _edge, _length in _edge_lengths.items():
                            _bv = getattr(self, f"_bc_{_edge}_values", None)
                            if _bv is not None:
                                for _k, _v in _bv.items():
                                    if isinstance(_v, np.ndarray) and len(_v) != _length:
                                        raise ValueError(
                                            f"Boundary condition at {_edge!r}: value for "
                                            f"{_k!r} has length {len(_v)}, expected "
                                            f"{_length} ({_length} nodes on that edge)."
                                        )
        else:
            # Analytical solution boundary conditions
            # If they aren't set, it is because no input file has been used
            # Just set them to an empty string (like input file would do)
            try:
                self.bc_east
            except AttributeError:
                self.bc_east = ""
            try:
                self.bc_west
            except AttributeError:
                self.bc_west = ""
            if self.dimension == 2:
                try:
                    self.bc_south
                except AttributeError:
                    self.bc_south = ""
                try:
                    self.bc_north
                except AttributeError:
                    self.bc_north = ""
            else:
                # Simplifies flow control a few lines down to define these as None-type
                self.bc_south = None
                self.bc_north = None
            if (
                self.bc_east == "no_outside_loads"
                or self.bc_east == ""
                and self.bc_west == "no_outside_loads"
                or self.bc_west == ""
            ) and (
                self.dimension != 2
                or (
                    self.bc_east == "no_outside_loads"
                    or self.bc_east == ""
                    and self.bc_west == "no_outside_loads"
                    or self.bc_west == ""
                )
            ):
                if (
                    self.bc_east == ""
                    or self.bc_west == ""
                    or self.bc_south == ""
                    or self.bc_north == ""
                ):
                    if self.verbose:
                        print(
                            "Assuming no_outside_loads boundary condition, as this is"
                            " implicit in the superposition-based analytical solution"
                        )
            else:
                if not self.quiet:
                    raise ValueError(
                        "Boundary conditions improperly defined for an analytical "
                        "solution. Analytical solvers require 'no_outside_loads' (or "
                        "blank) on all boundaries. FD boundary conditions such as "
                        "'zero_moment_zero_shear' cannot be applied to analytical "
                        "solutions."
                    )

    def coeffArraySizeCheck(self):
        """
        Make sure that q0 and coefficient array are the right size compared to
        each other (for finite difference if loading a pre-build coefficient
        array). Otherwise, exit.
        """
        if np.prod(self.coeff_matrix.shape) != np.long(
            np.prod(np.array(self.qs.shape, dtype=np.int64) + 2) ** 2
        ):
            raise ValueError(
                "Inconsistent size of q0 array and coefficient matrix."
            )

    def te_array_size_check(self):
        """
        Checks that Te and q0 array sizes are compatible
        For finite difference solution.
        """
        # Only if they are both defined and are arrays
        # Both being arrays is a possible bug in this check routine that I have
        # intentionally introduced
        if isinstance(self.te, np.ndarray) and isinstance(self.qs, np.ndarray):
            # Doesn't touch non-arrays or 1D arrays
            if type(self.te) is np.ndarray:
                if (np.array(self.te.shape) != np.array(self.qs.shape)).any():
                    raise ValueError(
                        f"q0 and Te arrays have incompatible shapes: "
                        f"qs={self.qs.shape}, te={self.te.shape}."
                    )
            else:
                if self.debug:
                    print("Te and qs array sizes pass consistency check")

    ### need to determine its interface, it is best to have a uniform interface
    ### no matter it is 1D or 2D; but if it can't be that way, we can set up a
    ### variable-length arguments, which is the way how Python overloads functions.

    def _solve_fd(self):
        """
        Set-up for the finite difference solution method
        """
        if self.verbose:
            print("Finite Difference Solution Technique")
        # Used to check for coeff_matrix here, but now doing so in self.bc_check()
        # called by f1d and f2d at the start
        #
        # Define a stress-based qs = q0
        # But only if the latter has not already been defined
        # (e.g., by the getters and setters)
        try:
            self.qs
        except AttributeError:
            self.qs = self.q0.copy()
            # Remove self.q0 to avoid issues with multiply-defined inputs
            # q0 is the parsable input to either a qs grid or contains (x,(y),q)
            del self.q0
        # Give it x and y dimensions for help with plotting tools
        # (not implemented internally, but a help with external methods)
        self.x = np.arange(self.dx / 2.0, self.dx * self.qs.shape[0], self.dx)
        if self.dimension == 2:
            self.y = np.arange(self.dy / 2.0, self.dy * self.qs.shape[1], self.dy)
        # Config file may override the default solver
        if self.filename:
            _solver = self.configGet("string", "numerical", "solver", optional=True)
            if _solver is not None:
                self.solver = _solver
        # Check consistency of size if coeff array was loaded
        if self.filename:
            # Try to import Te grid or scalar for the finite difference solution
            # Try float first (scalar Te); fall back to string (file path)
            self.te = self.configGet(
                "float", "input", "elastic_thickness", optional=True
            )
            if self.te is None:
                Tepath = self.configGet(
                    "string", "input", "elastic_thickness", optional=False
                )
                self.te = Tepath
            else:
                Tepath = None
            if self.te is None:
                if self.coeff_matrix is not None:
                    pass
                else:
                    # Have to bring this out here in case it was discovered in the
                    # try statement that there is no value given
                    raise ValueError(
                        "No input elastic thickness or coefficient matrix supplied."
                    )
        # or if getter/setter
        if isinstance(self.te, str):
            # Try to import Te grid or scalar for the finite difference solution
            Tepath = self.te
        else:
            Tepath = None  # in case no self.filename present (like for GRASS GIS)
        # If there is a Tepath, import Te
        # Assume that even if a coeff_matrix is defined
        # That the user wants Te if they gave the path
        if Tepath:
            self.te = self.loadFile(self.te, close_on_fail=False)
            if self.te is None:
                print("Requested Te file is provided but cannot be located.")
                print("No scalar elastic thickness is provided in configuration file")
                print("(Typo in path to input Te grid?)")
                if self.coeff_matrix is not None:
                    print("But a coefficient matrix has been found.")
                    print("Calculations will be carried forward using it.")
                else:
                    raise FileNotFoundError(
                        "Requested Te file cannot be located and no scalar elastic "
                        "thickness or coefficient matrix is available."
                    )

            # Check that Te is the proper size if it was loaded
            # Will be array if it was loaded
            if self.te.any():
                self.te_array_size_check()

    def _solve_fft(self):
        """
        Set-up for the FFT spectral solution method.

        Resolves the load array from ``q0`` if needed, assigns the ``x``
        (and, in two dimensions, ``y``) coordinate arrays for plotting
        utilities, and loads a scalar elastic thickness from the
        configuration file when one is present.  Mirrors the set-up
        performed by :meth:`FD` and :meth:`SAS`; the scalar-:math:`T_e`
        requirement is enforced in :class:`F1D` and :class:`F2D`.
        """
        if self.verbose:
            print("FFT Spectral Solution Technique")
        # Define qs from q0 if not already set by a getter/setter
        try:
            self.qs
        except AttributeError:
            self.qs = self.q0.copy()
            # Remove self.q0 to avoid issues with multiply-defined inputs
            del self.q0
        # Give it x (and y) dimensions for help with plotting tools
        self.x = np.arange(self.dx / 2.0, self.dx * self.qs.shape[0], self.dx)
        if self.dimension == 2:
            self.y = np.arange(self.dy / 2.0, self.dy * self.qs.shape[1], self.dy)
        # Config-file parameter loading
        # FFT requires scalar (uniform) Te; that check is performed in F1D/F2D
        if self.filename:
            self.te = self.configGet("float", "input", "elastic_thickness")
            # qs may still be coming from q0 when driven by a config file
            try:
                self.qs
            except AttributeError:
                self.qs = self.q0.copy()
                del self.q0

    # SAS and SAS_NG are the exact same here; leaving separate just for symmetry
    # with other functions

    def _solve_sas(self):
        """
        Set-up for the rectangularly-gridded superposition of analytical solutions
        method for solving flexure
        """
        if self.x is None:
            self.x = np.arange(self.dx / 2.0, self.dx * self.qs.shape[0], self.dx)
        if self.filename:
            # Define the (scalar) elastic thickness
            self.te = self.configGet("float", "input", "elastic_thickness")
            # Define a stress-based qs = q0
            self.qs = self.q0.copy()
            # Remove self.q0 to avoid issues with multiply-defined inputs
            # q0 is the parsable input to either a qs grid or contains (x,(y),q)
            del self.q0
        if self.dimension == 2:
            if self.y is None:
                self.y = np.arange(self.dy / 2.0, self.dy * self.qs.shape[0], self.dy)
            # Define a stress-based qs = q0
            # But only if the latter has not already been defined
            # (e.g., by the getters and setters)
            try:
                self.qs
            except AttributeError:
                self.qs = self.q0.copy()
                # Remove self.q0 to avoid issues with multiply-defined inputs
                # q0 is the parsable input to either a qs grid or contains (x,(y),q)
                del self.q0

    def _solve_sas_ng(self):
        """
        Set-up for the ungridded superposition of analytical solutions
        method for solving flexure
        """
        if self.filename:
            # Define the (scalar) elastic thickness
            self.te = self.configGet("float", "input", "elastic_thickness")
            # See if it wants to be run in lat/lon
            # Could put under in 2D if-statement, but could imagine an eventual desire
            # to change this and have 1D lat/lon profiles as well.
            # So while the options will be under "numerical2D", this place here will
            # remain held for an eventual future.
            self.latlon = self.configGet(
                "string", "numerical2D", "latlon", optional=True
            )
            self.planetary_radius = self.configGet(
                "float", "numerical2D", "planetary_radius", optional=True
            )
        # Parse out input q0 into variables of imoprtance for solution
        if self.dimension == 1:
            try:
                # If these have already been set, e.g., by getters/setters, great!
                self.x
                self.q
            except AttributeError:
                # Using [x, y, w] configuration file
                if self.q0.shape[1] == 2:
                    self.x = self.q0[:, 0]
                    self.q = self.q0[:, 1]
                else:
                    raise ValueError(
                        "For 1D (ungridded) SAS_NG, need an [x, w] array. "
                        f"Got shape: {self.q0.shape}."
                    )
        else:
            try:
                # If these have already been set, e.g., by getters/setters, great!
                self.x
                self.y
                self.q
            except AttributeError:
                # Using [x, y, w] configuration file
                if self.q0.shape[1] == 3:
                    self.x = self.q0[:, 0]
                    self.y = self.q0[:, 1]
                    self.q = self.q0[:, 2]
                else:
                    raise ValueError(
                        "For 2D (ungridded) SAS_NG, need an [x, y, w] array. "
                        f"Got shape: {self.q0.shape}."
                    )
        # x, y are in absolute coordinates. Create a local grid reference to
        # these. This local grid, which starts at (0,0), is defined just so that
        # we have a way of running the model without defined real-world
        # coordinates
        self.x = self.x
        if self.dimension == 2:
            self.y = self.y
        # Remove self.q0 to avoid issues with multiply-defined inputs
        # q0 is the parsable input to either a qs grid or contains (x,(y),q)
        del self.q0

        # Check if a seperate output set of x,y points has been defined
        # otherwise, set those values to None
        # First, try to load the arrays
        try:
            self.xw
        except AttributeError:
            self.xw = self.configGet("string", "input", "xw", optional=True) or None
        # If strings, load arrays
        if isinstance(self.xw, str):
            self.xw = self.loadFile(self.xw)
        if self.dimension == 2:
            try:
                # already set by setter?
                self.yw
            except AttributeError:
                self.yw = self.configGet("string", "input", "yw", optional=True) or None
            # At this point, can check if we have both None or both defined
            if (self.xw is not None and self.yw is None) or (
                self.xw is None and self.yw is not None
            ):
                raise ValueError(
                    "SAS_NG output at specified points requires both xw and yw to be defined"
                )
            # All right, now just finish defining
            if isinstance(self.yw, str):
                self.yw = self.loadFile(self.yw)
            elif self.yw is None:
                self.yw = self.y.copy()
        if self.xw is None:
            self.xw = self.x.copy()
