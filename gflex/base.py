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
import logging
import os

import warnings
from enum import Enum

import numpy as np
from matplotlib import pyplot as plt

from ._version import __version__

_logger = logging.getLogger(__name__)

# Module-level console handler for the convenience quiet=False mode.
_console_handler: logging.StreamHandler | None = None


def _apply_log_level(level: int) -> None:
    """Set the gflex logger level and manage the convenience console handler.

    Called by the quiet/verbose/debug property setters.  In a coupling
    context the application configures its own handlers and log level for
    the 'gflex' logger; this function is a convenience for standalone use.
    """
    global _console_handler
    gflex_log = logging.getLogger("gflex")
    if level < logging.WARNING:
        if _console_handler is None:
            _console_handler = logging.StreamHandler()
            _console_handler.setFormatter(logging.Formatter("%(message)s"))
            gflex_log.addHandler(_console_handler)
        gflex_log.setLevel(level)
    else:
        if _console_handler is not None:
            gflex_log.removeHandler(_console_handler)
            _console_handler = None
        gflex_log.setLevel(logging.WARNING)


# Sentinel for "attribute not yet assigned" — used by property setters to
# distinguish a first-time assignment (no cache to invalidate) from a genuine
# value change.
_UNSET = object()


def _value_changed(current, new):
    """Return True if *new* differs meaningfully from *current*.

    Handles None, Python scalars, NumPy scalars, NumPy arrays, strings, and
    dicts (inhomogeneous BC specifications).  Returns False when *current* is
    ``_UNSET`` so that first-time assignments never trigger cache invalidation.
    """
    if current is _UNSET:
        return False
    if isinstance(current, dict) or isinstance(new, dict):
        return current != new
    try:
        return not np.array_equal(current, new)
    except Exception:
        return True


def _readonly_array_copy(value):
    """Store a matrix-determining input, locking array values read-only.

    If *value* is a NumPy array, return an owned, read-only ``float`` copy;
    otherwise return *value* unchanged.  The coefficient-matrix / LU cache is
    keyed to these inputs, so an in-place edit would silently desynchronise it.
    Storing arrays read-only makes such edits raise ``ValueError`` instead;
    the input is changed by reassignment, which invalidates the cache.
    """
    if isinstance(value, np.ndarray):
        value = np.array(value, dtype=float)
        value.flags.writeable = False
    return value


def _normalize_cache_factorization(value):
    """Validate ``cache_factorization``; map the deprecated ``'no_check'``.

    Valid values are ``False`` (no LU cache) and ``True`` (cache the LU
    factorisation and trust setter-based invalidation).  ``'no_check'`` is a
    deprecated alias for ``True`` — they are now equivalent, because
    matrix-determining inputs invalidate the cache on reassignment and array
    inputs are read-only, so the cached matrix cannot silently desynchronise.
    """
    if value == "no_check":
        warnings.warn(
            "cache_factorization='no_check' is deprecated; use True.  The "
            "cached factorisation is now always trusted (matrix-determining "
            "inputs invalidate the cache on change, and array inputs are "
            "read-only), so True and 'no_check' are equivalent.",
            DeprecationWarning,
            stacklevel=3,
        )
        return True
    if value not in (False, True):
        raise ValueError(
            "cache_factorization must be False or True ('no_check' is a "
            f"deprecated alias for True); got {value!r}"
        )
    return value


# Boundary conditions that fix the boundary-node displacement (Dirichlet in w).
# A load on such a node is reacted by the support, not a deflection.
_FIXED_DISPLACEMENT_BCS = (
    "zero_displacement_zero_slope",
    "zero_displacement_zero_moment",
)


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


# ── LU memory estimate helpers ────────────────────────────────────────────────
# Empirical power-law fit to SuperLU/COLAMD on a gFlex-style 13-diagonal
# sparse matrix: M_MiB ≈ C × N^exp.  Used to warn before large allocations.

_LU_RAM_C   = 2.41e-4  # fit coefficient (MiB per cell^_LU_RAM_EXP)
_LU_RAM_EXP = 1.26     # empirical log-log slope


def _available_ram_bytes():
    """Return available system RAM in bytes, or None if indeterminate."""
    try:
        import psutil
        return psutil.virtual_memory().available
    except ImportError:
        pass
    try:
        with open("/proc/meminfo") as fh:
            for line in fh:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) * 1024
    except OSError:
        pass
    return None


def _estimate_lu_ram_bytes(n_cells):
    """Estimated peak LU RAM in bytes for a *n_cells*-cell FD grid.

    Uses M ≈ 2.41×10⁻⁴ · N^1.26 MiB, fit to SuperLU/COLAMD on a gFlex
    13-diagonal sparse matrix.  Accurate to within ~20 % for square 2-D grids;
    conservative (overestimates) for 1-D problems.
    """
    return _LU_RAM_C * (n_cells ** _LU_RAM_EXP) * 1024 * 1024


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

VALID_BC_STRINGS_2D = frozenset({
    "zero_displacement_zero_slope", "clamped",
    "zero_displacement_zero_moment", "pinned",
    "zero_moment_zero_shear", "free",
    "zero_slope_zero_shear", "mirror",
    "periodic",
    "no_outside_loads", "infinite",
})
"""All boundary-condition strings accepted by :class:`~gflex.F2D`'s FD solver.

Canonical names (``zero_displacement_zero_slope``, ``zero_displacement_zero_moment``,
``zero_moment_zero_shear``, ``zero_slope_zero_shear``, ``periodic``,
``no_outside_loads``) plus their short aliases (``clamped``, ``pinned``, ``free``,
``mirror``, ``infinite``).  Wrappers can validate user input against this set instead
of maintaining a parallel copy.
"""

VALID_BC_STRINGS_1D = VALID_BC_STRINGS_2D | frozenset({"sandbox"})
"""All boundary-condition strings accepted by :class:`~gflex.F1D`'s FD solver.

A superset of :data:`VALID_BC_STRINGS_2D` that additionally includes
``"sandbox"`` (free end, 1-D only).
"""


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
                        _logger.warning(
                            "An empty input string here is not an acceptable"
                            " option."
                        )
                        _logger.warning("%s is not optional.", name)
                        _logger.warning("Program crash likely to occur.")
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
                if not self.grass:
                    _logger.info('No value entered for optional parameter "%s"', name)
                    _logger.info('in category "%s" in configuration file.', category)
                    _logger.info(
                        "No action related to this optional parameter will be taken."
                    )
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
                    _logger.info("dx and dy being overwritten -- supply a full grid")
            else:
                _logger.info("dx and dy being overwritten -- supply a full grid")
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
        elif out is not None:
            _logger.info("Loading %s from %s", var, format_name)

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
            _logger.info("Starting to plot %s", self.plot_choice)
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
                        _logger.warning(
                                "Combo plot can't work with SAS_NG! Don't have mechanism"
                                " in place to calculate load width."
                            )
                        _logger.warning(
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
                        self.T_e
                        if self.method == "fd":
                            if type(self.T_e) is np.ndarray:
                                if (self.T_e != (self.T_e).mean()).any():
                                    plt.title(titletext, fontsize=16)
                                else:
                                    plt.title(
                                        titletext
                                        + ", $T_e$ = "
                                        + str((self.T_e / 1000).mean())
                                        + " km",
                                        fontsize=16,
                                    )
                            else:
                                plt.title(
                                    titletext
                                    + ", $T_e$ = "
                                    + str(self.T_e / 1000)
                                    + " km",
                                    fontsize=16,
                                )
                        else:
                            plt.title(
                                titletext + ", $T_e$ = " + str(self.T_e / 1000) + " km",
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
                    _logger.warning(
                            'Incorrect plot_choice input, "%s" provided.',
                            self.plot_choice,
                        )
                    _logger.warning(
                            "Possible input strings are: q, w, both, and (for 1D) combo"
                        )
                    _logger.warning("Unable to produce plot.")
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
                    _logger.warning(
                            'Incorrect plot_choice input, "%s" provided.',
                            self.plot_choice,
                        )
                    _logger.warning(
                            "Possible input strings are: q, w, both, and (for 1D) combo"
                        )
                    _logger.warning("Unable to produce plot.")

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

        _has_te_grid = isinstance(self.T_e, np.ndarray) and self.T_e.ndim == 2

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
                self.T_e / 1e3, extent=_ext(self.T_e), cmap=_cmap_te,
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

        _logger.info("Starting to interpolate grid for plotting -- can be a slow process!")

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
        # Default to quiet for library / coupling use.  Set quiet=False,
        # verbose=True, or debug=True to enable console output.
        self.quiet = True
        self.verbose = False
        self.debug = False

        # FFT padding: number of flexural-parameter units to zero-pad on each
        # side for non-periodic runs.  F1D uses α₁D = (4D/Δρg)^0.25;
        # F2D uses α₂D = (D/Δρg)^0.25.  Periodic images are separated by
        # 2 × fft_pad_n_alpha × α, so the default of 4 gives 8α separation.
        # Only used when method='fft' and BCs are not all 'periodic'.
        self.fft_pad_n_alpha = 4

        # x and y to None for checks
        self.x = None
        self.y = None

        # Output defaults — overridden by config file or explicit setter
        self.plot_choice = None
        self.w_out_file = None

        # LU factorization cache.
        # False  → no cache; spsolve on every run() (default).
        # True   → cache the factorized callable; reuse when the coefficient
        #           matrix hash is unchanged.
        # "no_check" → reuse the cached callable unconditionally; user
        #           guarantees the matrix is stable across run() calls.
        self.cache_factorization = False
        self._lu = None

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

    def _invalidate_matrix_cache(self):
        """Clear the FD coefficient matrix and LU factorisation caches.

        Called automatically by property setters whenever a matrix-determining
        input (``T_e``, ``E``, ``nu``, boundary conditions, grid spacing, …) is
        reassigned to a different value.  The caches are rebuilt transparently
        on the next ``run()`` call.

        Invalidation happens only on (re)assignment, so matrix-determining
        inputs must be **reassigned**, not mutated in place.  To enforce this
        for the most error-prone case, ``T_e`` arrays are stored read-only:
        an in-place edit (``flex.T_e[5] = 40e3``) raises ``ValueError`` rather
        than silently desynchronising the cache.  Reassign the full array
        (``flex.T_e = new_array``) to change it.
        """
        d = self.__dict__
        if "coeff_matrix" in d:
            d["coeff_matrix"] = None
        if "_lu" in d:
            d["_lu"] = None

    def _edge_fixes_displacement(self, norm, values):
        """True if an edge's BC prescribes its boundary-node displacement.

        Covers the homogeneous clamped/pinned BCs (``zero_displacement_*``)
        and dict-style ``{"displacement": ...}`` BCs.
        """
        if values is not None and "displacement" in values:
            return True
        return norm in _FIXED_DISPLACEMENT_BCS

    def _cancel_load_on_fixed_nodes(self, rhs):
        """Cancel any load sitting on a fixed-displacement boundary node.

        The boundary row of such a node is decoupled (``c0·w = rhs``) with
        ``rhs = -qs + correction``.  A load there is reacted by the support, so
        it must not deflect the node; adding ``qs`` back makes the row solve
        ``c0·w = c0·w0`` (``w0 = 0`` for the homogeneous clamped/pinned BCs,
        the prescribed value for dict displacement), enforcing ``w = w0``
        exactly.  Warns once if a load is actually present on a fixed node.

        *rhs* is the flat right-hand side; the boolean mask of fixed nodes is
        provided per-dimension by ``_fixed_displacement_mask``.
        """
        mask = self._fixed_displacement_mask()
        if not mask.any():
            return rhs
        mflat = mask.reshape(-1, order="C")
        qflat = np.asarray(self.qs).reshape(-1, order="C")
        if np.any(qflat[mflat] != 0.0):
            warnings.warn(
                "Load applied on a fixed-displacement boundary node (a "
                "clamped, pinned, or prescribed-displacement edge).  Such a "
                "load is reacted by the support and does not deflect the node; "
                "its displacement is held at the prescribed value.  Move the "
                "load off the boundary, or change the boundary condition, if "
                "this is unintended.",
                UserWarning,
                stacklevel=4,
            )
        rhs[mflat] += qflat[mflat]
        return rhs

    # ── Verbosity properties ───────────────────────────────────────────────────

    def _sync_log_level(self) -> None:
        """Recompute and apply the effective log level from the three flags."""
        if getattr(self, "_debug", False):
            _apply_log_level(logging.DEBUG)
        elif getattr(self, "_verbose", False) and not getattr(self, "_quiet", True):
            _apply_log_level(logging.INFO)
        else:
            _apply_log_level(logging.WARNING)

    @property
    def quiet(self):
        return self._quiet

    @quiet.setter
    def quiet(self, value):
        self._quiet = bool(value)
        self._sync_log_level()

    @property
    def verbose(self):
        return self._verbose

    @verbose.setter
    def verbose(self, value):
        self._verbose = bool(value)
        self._sync_log_level()

    @property
    def debug(self):
        return self._debug

    @debug.setter
    def debug(self, value):
        self._debug = bool(value)
        self._sync_log_level()

    # ── Properties: matrix-determining inputs ─────────────────────────────────
    # Each setter invalidates the coefficient matrix and LU cache when the
    # value actually changes.  Repeated run() calls with identical physics
    # therefore reuse the cached factorisation at no extra cost.

    @property
    def T_e(self):
        """Elastic thickness [m] (scalar or array).

        The pristine, user-supplied grid.  The setter stores it in ``_T_e``;
        the solution methods work on ``_te``, a separate writeable copy that
        may be padded or otherwise modified during a solve.

        An array assigned here is stored as a **read-only** copy: the cache of
        the coefficient matrix / LU factorisation is keyed to this value, and
        an in-place edit (``flex.T_e[i] = ...``) would silently desynchronise
        it.  Such edits raise ``ValueError``; reassign the whole array
        (``flex.T_e = new_array``) to change it, which invalidates the cache.
        """
        return self._T_e

    @T_e.setter
    def T_e(self, value):
        if _value_changed(self.__dict__.get("_T_e", _UNSET), value):
            self._invalidate_matrix_cache()
        self._T_e = _readonly_array_copy(value)

    @property
    def E(self):
        """Young's modulus [Pa]."""
        return self._E

    @E.setter
    def E(self, value):
        if _value_changed(self.__dict__.get("_E", _UNSET), value):
            self._invalidate_matrix_cache()
        self._E = value

    @property
    def nu(self):
        """Poisson's ratio [-]."""
        return self._nu

    @nu.setter
    def nu(self, value):
        if _value_changed(self.__dict__.get("_nu", _UNSET), value):
            self._invalidate_matrix_cache()
        self._nu = value

    @property
    def g(self):
        """Gravitational acceleration [m s⁻²]."""
        return self._g

    @g.setter
    def g(self, value):
        if _value_changed(self.__dict__.get("_g", _UNSET), value):
            self._invalidate_matrix_cache()
        self._g = value

    @property
    def rho_m(self):
        """Mantle density [kg m⁻³]."""
        return self._rho_m

    @rho_m.setter
    def rho_m(self, value):
        changed = _value_changed(self.__dict__.get("_rho_m", _UNSET), value)
        self._rho_m = value
        if "_rho_fill" in self.__dict__:
            self.__dict__["drho"] = value - self._rho_fill
        if changed:
            self._invalidate_matrix_cache()

    @property
    def rho_fill(self):
        """Infill density [kg m⁻³]."""
        return self._rho_fill

    @rho_fill.setter
    def rho_fill(self, value):
        changed = _value_changed(self.__dict__.get("_rho_fill", _UNSET), value)
        self._rho_fill = value
        if "_rho_m" in self.__dict__:
            self.__dict__["drho"] = self._rho_m - value
        if changed:
            self._invalidate_matrix_cache()

    def _validate_drho(self):
        """Check that the mantle/infill density difference is positive.

        ``drho = rho_m - rho_fill`` sets the sign of the restoring force
        ``drho·g`` and hence the flexural parameter α; a non-positive value
        is unphysical and makes α undefined.  This runs both at
        initialisation *and* before every solve (see :meth:`run`), so a
        density reassigned after ``initialize()`` — e.g. through the BMI
        ``set_value`` interface — cannot silently solve with a non-positive
        restoring force.  Validating at solve time rather than in the density
        setters keeps a transient invalid state legal during a two-step
        reassignment (e.g. lowering ``rho_m`` below the old ``rho_fill``
        before lowering ``rho_fill``).
        """
        self.drho = self.rho_m - self.rho_fill
        if self.drho <= 0:
            raise ValueError(
                f"rho_fill ({self.rho_fill} kg m⁻³) must be strictly less than "
                f"rho_m ({self.rho_m} kg m⁻³). "
                f"drho = {self.drho} kg m⁻³; the restoring force drho·g is "
                "non-positive and the flexural parameter α is undefined."
            )

    def _warn_if_lu_ram_high(self):
        """Warn if the estimated LU-factorisation memory exceeds 80% of free RAM.

        Shared by the F1D and F2D FD solvers (``self.qs.size`` is the cell
        count in either dimension).  A no-op when free RAM cannot be
        determined.
        """
        avail = _available_ram_bytes()
        if avail is None:
            return
        est = _estimate_lu_ram_bytes(self.qs.size)
        pct = 100.0 * est / avail
        if pct > 80:
            warnings.warn(
                f"FD solver: estimated LU memory requirement "
                f"~{est / 1e9:.1f} GB out of {avail / 1e9:.1f} GB free.\n"
                f"{pct:.0f}% of available RAM may be needed.\n"
                "The solve may exhaust memory.\n"
                "Consider method='fft' (requires uniform Te) "
                "or reducing grid resolution.",
                UserWarning,
                stacklevel=4,
            )

    @property
    def allow_unpaired_periodic(self):
        """Permit a one-sided ('unpaired') periodic boundary condition.

        Periodic boundary conditions connect opposite edges, so setting
        ``'periodic'`` on only one side of a pair is not a well-posed
        boundary: the finite-difference matrix carries a spurious partial
        wrap-around and the deflection near that edge is not a valid
        solution.  :meth:`bc_check` therefore raises a ``ValueError`` for an
        unpaired periodic BC by default.  Set this to ``True`` to override the
        guard and solve anyway — a deliberate, use-at-your-own-risk escape
        hatch.  Defaults to ``False``.
        """
        return self.__dict__.get("_allow_unpaired_periodic", False)

    @allow_unpaired_periodic.setter
    def allow_unpaired_periodic(self, value):
        value = bool(value)
        if value and not self.__dict__.get("_allow_unpaired_periodic", False):
            warnings.warn(
                "allow_unpaired_periodic = True: the unpaired-periodic safety "
                "check is now disabled.",
                UserWarning, stacklevel=2,
            )
        self.__dict__["_allow_unpaired_periodic"] = value

    @property
    def dx(self):
        """Grid spacing in x [m]."""
        return self._dx

    @dx.setter
    def dx(self, value):
        if _value_changed(self.__dict__.get("_dx", _UNSET), value):
            self._invalidate_matrix_cache()
        self._dx = value

    @property
    def dy(self):
        """Grid spacing in y [m]."""
        return self._dy

    @dy.setter
    def dy(self, value):
        if _value_changed(self.__dict__.get("_dy", _UNSET), value):
            self._invalidate_matrix_cache()
        self._dy = value

    @property
    def bc_west(self):
        """West boundary condition."""
        return self._bc_west

    @bc_west.setter
    def bc_west(self, value):
        if _value_changed(self.__dict__.get("_bc_west", _UNSET), value):
            self._invalidate_matrix_cache()
        self._bc_west = value

    @property
    def bc_east(self):
        """East boundary condition."""
        return self._bc_east

    @bc_east.setter
    def bc_east(self, value):
        if _value_changed(self.__dict__.get("_bc_east", _UNSET), value):
            self._invalidate_matrix_cache()
        self._bc_east = value

    @property
    def bc_north(self):
        """North boundary condition (2-D only)."""
        return self._bc_north

    @bc_north.setter
    def bc_north(self, value):
        if _value_changed(self.__dict__.get("_bc_north", _UNSET), value):
            self._invalidate_matrix_cache()
        self._bc_north = value

    @property
    def bc_south(self):
        """South boundary condition (2-D only)."""
        return self._bc_south

    @bc_south.setter
    def bc_south(self, value):
        if _value_changed(self.__dict__.get("_bc_south", _UNSET), value):
            self._invalidate_matrix_cache()
        self._bc_south = value

    @property
    def sigma_xx(self):
        """In-plane normal stress in x [Pa]."""
        return self._sigma_xx

    @sigma_xx.setter
    def sigma_xx(self, value):
        if _value_changed(self.__dict__.get("_sigma_xx", _UNSET), value):
            self._invalidate_matrix_cache()
        self._sigma_xx = _readonly_array_copy(value)

    @property
    def sigma_yy(self):
        """In-plane normal stress in y [Pa]."""
        return self._sigma_yy

    @sigma_yy.setter
    def sigma_yy(self, value):
        if _value_changed(self.__dict__.get("_sigma_yy", _UNSET), value):
            self._invalidate_matrix_cache()
        self._sigma_yy = _readonly_array_copy(value)

    @property
    def sigma_xy(self):
        """In-plane shear stress [Pa]."""
        return self._sigma_xy

    @sigma_xy.setter
    def sigma_xy(self, value):
        if _value_changed(self.__dict__.get("_sigma_xy", _UNSET), value):
            self._invalidate_matrix_cache()
        self._sigma_xy = _readonly_array_copy(value)

    @property
    def planetary_radius(self):
        """Planetary radius for spherical-domain corrections [m], or None."""
        return self.__dict__.get("_planetary_radius")

    @planetary_radius.setter
    def planetary_radius(self, value):
        if _value_changed(self.__dict__.get("_planetary_radius", _UNSET), value):
            self._invalidate_matrix_cache()
        self._planetary_radius = value

    # ── End matrix-determining properties ─────────────────────────────────────

    @property
    def method(self):
        """Solution method: ``'fd'``, ``'fft'``, ``'sas'``, or ``'sas_ng'``."""
        return self._method

    @method.setter
    def method(self, value):
        if isinstance(value, str):
            self._method = value.lower()
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
        # _sync_log_level() in the property setters already enforces that quiet
        # takes priority; no separate override block needed here.

        # Introduce model
        # After configuration file can define "quiet", and getter/setter should be done
        # by this point if we are going that way.
        _logger.info("")
        _logger.info("****************************" + "*" * len(__version__))
        _logger.info("*** Initializing gFlex v" + __version__ + " ***")
        _logger.info("****************************" + "*" * len(__version__))
        _logger.info("Open-source licensed under GNU GPL v3")
        _logger.info("")

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
        self._validate_drho()
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
            self.T_e = self.T_e.astype(float)  # array
        except AttributeError:
            # Integer scalar Te does not seem to be a problem, but taking this step
            # anyway for consistency
            try:
                self.T_e = float(self.T_e)  # integer
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
        Release all model state.

        Deletes the deflection result (``w``), the load array (``qs``), and
        the cached coefficient matrix.  Read ``w`` and save any output
        **before** calling ``finalize``; accessing ``w`` afterwards will raise
        ``AttributeError``.  Called automatically by :meth:`F1D.finalize` and
        :meth:`F2D.finalize`.
        """
        for _attr in ("w", "qs", "coeff_matrix", "_lu"):
            with contextlib.suppress(AttributeError):
                delattr(self, _attr)

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
        _logger.info("Output step")
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
            _logger.info("Output filename provided.")
        if self.w_out_file:
            if self.w_out_file[-4:] == ".npy":
                from numpy import save

                save(self.w_out_file, self.w)
            else:
                from numpy import savetxt

                # Shouldn't need more than mm precision, at very most
                savetxt(self.w_out_file, self.w, fmt="%.3f")
                _logger.info("Saving deflections --> %s", self.w_out_file)

    def bc_check(self):
        """
        Validate boundary conditions and prepare BC state for the chosen solver.

        For ``FD``: checks that every edge BC is one of the accepted strings
        (``'zero_displacement_zero_slope'`` / ``'clamped'``,
        ``'zero_moment_zero_shear'`` / ``'free'``,
        ``'mirror'``, ``'periodic'``, etc.) and exits with an informative
        message if not.  A pre-built ``coeff_matrix`` bypasses
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
            # Ensure BC attributes exist; FFT handles them internally.
            # 'periodic' → exact transform; anything else → zero-padded (no_outside_loads).
            for edge in ("bc_east", "bc_west"):
                if not hasattr(self, edge):
                    setattr(self, edge, "")
            if self.dimension == 2:
                for edge in ("bc_north", "bc_south"):
                    if not hasattr(self, edge):
                        setattr(self, edge, "")
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

            # Resolve short string aliases unconditionally.  The coeff_matrix
            # guard below may be skipped on a re-run (e.g. when a property
            # setter invalidates the cached matrix after bc_check() has already
            # seeded _bc_*_norm with the raw user value), so _apply_bc_rigidity
            # must never see un-aliased names such as "mirror" or "free".
            _BC_ALIASES = {
                "clamped":  "zero_displacement_zero_slope",
                "pinned":   "zero_displacement_zero_moment",
                "free":     "zero_moment_zero_shear",
                "mirror":   "zero_slope_zero_shear",
                "infinite": "no_outside_loads",
            }
            _bc_norm_attrs = ["_bc_west_norm", "_bc_east_norm"]
            if self.dimension == 2:
                _bc_norm_attrs += ["_bc_north_norm", "_bc_south_norm"]
            for _attr in _bc_norm_attrs:
                _raw = getattr(self, _attr)
                if isinstance(_raw, str):
                    setattr(self, _attr, _BC_ALIASES.get(_raw, _raw))

            # No need to create a coeff_matrix if one already exists
            if self.coeff_matrix is None:
                # Acceptable string boundary conditions
                self.bc1D = np.array(
                    [
                        "zero_displacement_zero_slope",
                        "clamped",
                        "zero_displacement_zero_moment",
                        "pinned",
                        "periodic",
                        "zero_slope_zero_shear",
                        "mirror",
                        "zero_moment_zero_shear",
                        "free",
                        "sandbox",
                        "no_outside_loads",
                        "infinite",
                    ]
                )
                self.bc2D = np.array(
                    [
                        "zero_displacement_zero_slope",
                        "clamped",
                        "zero_displacement_zero_moment",
                        "pinned",
                        "periodic",
                        "zero_slope_zero_shear",
                        "mirror",
                        "zero_moment_zero_shear",
                        "free",
                        "no_outside_loads",
                        "infinite",
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
                    # Resolve any short alias to its canonical name using the
                    # same _BC_ALIASES map applied above — one source of truth.
                    norm = _BC_ALIASES.get(norm, norm)
                    setattr(self, f"_bc_{edge}_norm", norm)
                # Reject an unpaired ('one-sided') periodic boundary by
                # default.  Periodic BCs connect opposite edges, so setting
                # 'periodic' on only one side of a pair is not a well-posed
                # boundary: the assembled FD matrix carries a spurious partial
                # wrap-around and the deflection near that edge is not a valid
                # solution.  A user who nonetheless wants this can set
                # allow_unpaired_periodic = True to override the guard.
                for side_a, side_b in (("west", "east"), ("north", "south")):
                    if self.dimension == 1 and side_a == "north":
                        continue
                    norm_a = getattr(self, f"_bc_{side_a}_norm", None)
                    norm_b = getattr(self, f"_bc_{side_b}_norm", None)
                    if (norm_a == "periodic") != (norm_b == "periodic"):
                        periodic_side  = side_a if norm_a == "periodic" else side_b
                        missing_side   = side_b if norm_a == "periodic" else side_a
                        message = (
                            f"FD method: bc_{periodic_side} is 'periodic' but "
                            f"bc_{missing_side} is not. Periodic boundary "
                            f"conditions connect opposite edges, so a one-sided "
                            f"periodic is not well-posed. If you have good "
                            f"reason to want such a pair of boundary "
                            f"conditions, set allow_unpaired_periodic = True."
                        )
                        if not self.allow_unpaired_periodic:
                            raise ValueError(message)
                        # Flag is set: the user took responsibility when they
                        # enabled it (the setter announced it once), so proceed
                        # silently rather than re-warning on every solve.
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
                    _logger.info(
                            "Assuming no_outside_loads boundary condition, as this is"
                            " implicit in the superposition-based analytical solution"
                        )
            else:
                warnings.warn(
                    "Boundary conditions improperly defined for an analytical "
                    "solution. Analytical solvers require 'no_outside_loads' (or "
                    "blank) on all boundaries. FD boundary conditions such as "
                    "'zero_moment_zero_shear' are ignored by analytical solvers; "
                    "the infinite-plate (no_outside_loads) assumption will be used.",
                    UserWarning,
                    stacklevel=3,
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

    def T_e_array_size_check(self):
        """
        Checks that Te and q0 array sizes are compatible
        For finite difference solution.
        """
        # Only if they are both defined and are arrays
        # Both being arrays is a possible bug in this check routine that I have
        # intentionally introduced
        if isinstance(self.T_e, np.ndarray) and isinstance(self.qs, np.ndarray):
            # Doesn't touch non-arrays or 1D arrays
            if type(self.T_e) is np.ndarray:
                if (np.array(self.T_e.shape) != np.array(self.qs.shape)).any():
                    raise ValueError(
                        f"q0 and Te arrays have incompatible shapes: "
                        f"qs={self.qs.shape}, te={self.T_e.shape}."
                    )
            else:
                _logger.debug("Te and qs array sizes pass consistency check")

    ### need to determine its interface, it is best to have a uniform interface
    ### no matter it is 1D or 2D; but if it can't be that way, we can set up a
    ### variable-length arguments, which is the way how Python overloads functions.

    def _solve_fd(self):
        """
        Set-up for the finite difference solution method
        """
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
            self.T_e = self.configGet(
                "float", "input", "elastic_thickness", optional=True
            )
            if self.T_e is None:
                Tepath = self.configGet(
                    "string", "input", "elastic_thickness", optional=False
                )
                self.T_e = Tepath
            else:
                Tepath = None
            if self.T_e is None:
                if self.coeff_matrix is not None:
                    pass
                else:
                    # Have to bring this out here in case it was discovered in the
                    # try statement that there is no value given
                    raise ValueError(
                        "No input elastic thickness or coefficient matrix supplied."
                    )
        # or if getter/setter
        if isinstance(self.T_e, str):
            # Try to import Te grid or scalar for the finite difference solution
            Tepath = self.T_e
        else:
            Tepath = None  # in case no self.filename present (like for GRASS GIS)
        # If there is a Tepath, import Te
        # Assume that even if a coeff_matrix is defined
        # That the user wants Te if they gave the path
        if Tepath:
            self.T_e = self.loadFile(self.T_e, close_on_fail=False)
            if self.T_e is None:
                _logger.warning("Requested Te file is provided but cannot be located.")
                _logger.warning("No scalar elastic thickness is provided in configuration file")
                _logger.warning("(Typo in path to input Te grid?)")
                if self.coeff_matrix is not None:
                    _logger.warning("But a coefficient matrix has been found.")
                    _logger.warning("Calculations will be carried forward using it.")
                else:
                    raise FileNotFoundError(
                        "Requested Te file cannot be located and no scalar elastic "
                        "thickness or coefficient matrix is available."
                    )

            # Check that Te is the proper size if it was loaded
            # Will be array if it was loaded
            if self.T_e.any():
                self.T_e_array_size_check()

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
            self.T_e = self.configGet("float", "input", "elastic_thickness")
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
            self.T_e = self.configGet("float", "input", "elastic_thickness")
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
            self.T_e = self.configGet("float", "input", "elastic_thickness")
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
