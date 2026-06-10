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
import time
import warnings

_logger = logging.getLogger(__name__)

import numpy as np
import scipy
import scipy.fft
from scipy.signal import fftconvolve
from scipy.special import kei
from scipy.sparse.linalg import factorized, spsolve

from gflex.base import Flexure, _RigidityBC, _matrix_hash


def flexural_wavelengths(Te, E, nu, rho_m, rho_fill, g):
    """
    Compute flexural parameters and wavelengths for a thin elastic plate.

    Parameters
    ----------
    Te : float
        Elastic thickness [m].
    E : float
        Young's modulus [Pa].
    nu : float
        Poisson's ratio.
    rho_m : float
        Mantle density [kg m^-3].
    rho_fill : float
        Infill density [kg m^-3] (e.g. 0 for air, 1000 for water).
    g : float
        Gravitational acceleration [m s^-2].

    Returns
    -------
    dict
        Keys ``alpha_1D``, ``lambda_1D``, ``zero_crossing_1D``,
        ``alpha_2D``, ``lambda_2D``, ``zero_crossing_2D`` — all in metres.

    Examples
    --------
    >>> from gflex import flexural_wavelengths
    >>> r = flexural_wavelengths(Te=30e3, E=65e9, nu=0.25,
    ...                          rho_m=3300., rho_fill=0., g=9.8)
    >>> round(r["alpha_2D"] / 1e3, 1)  # km
    46.9
    >>> round(r["lambda_1D"] / r["lambda_2D"], 6)
    1.414214
    """
    drho = rho_m - rho_fill
    if drho <= 0:
        raise ValueError(
            f"rho_fill ({rho_fill} kg m⁻³) must be strictly less than "
            f"rho_m ({rho_m} kg m⁻³). "
            f"drho = {drho} kg m⁻³; the flexural parameter α is "
            "undefined for non-positive drho·g."
        )
    D = (E * float(Te) ** 3) / (12.0 * (1.0 - nu**2))
    alpha_1D = (4.0 * D / (drho * g)) ** 0.25
    alpha_2D = (D / (drho * g)) ** 0.25
    lambda_1D = 2.0 * np.pi * alpha_1D
    lambda_2D = 2.0 * np.pi * alpha_2D
    return {
        "alpha_1D": alpha_1D,
        "lambda_1D": lambda_1D,
        "zero_crossing_1D": 0.375 * lambda_1D,
        "alpha_2D": alpha_2D,
        "lambda_2D": lambda_2D,
        "zero_crossing_2D": 0.375 * lambda_2D,
    }


def recommended_pad_width(Te, dx, E=65e9, nu=0.25, rho_m=3300.0, rho_fill=0.0,
                            g=9.8, n_wavelengths=1.0):
    """
    Return the recommended padding width in grid cells for a variable-Te run.

    The padded domain boundary should be at least one flexural wavelength from
    the load so that the plate's response to the load is negligible at the
    boundary.  Half a wavelength is often sufficient in practice.

    The 2-D flexural wavelength is computed from the maximum Te value via
    :func:`flexural_wavelengths`, giving the most conservative (widest)
    padding estimate.

    Parameters
    ----------
    Te : scalar or 2-D array
        Elastic thickness [m].  The maximum value is used.
    dx : float
        Grid cell size [m].  Use the smaller of dx and dy if they differ.
    E : float, optional
        Young's modulus [Pa].  Default 65 GPa.
    nu : float, optional
        Poisson's ratio.  Default 0.25.
    rho_m : float, optional
        Mantle density [kg m^-3].  Default 3300.
    rho_fill : float, optional
        Infill density [kg m^-3].  Default 0 (air).
    g : float, optional
        Gravitational acceleration [m s^-2].  Default 9.8.
    n_wavelengths : float, optional
        Number of flexural wavelengths to use as the padding width.
        Default 1.0.  Use 0.5 for a less conservative estimate.

    Returns
    -------
    int
        Recommended padding width in grid cells.
    """
    r = flexural_wavelengths(
        Te=float(np.max(Te)), rho_m=rho_m, rho_fill=rho_fill, E=E, nu=nu, g=g,
    )
    return int(np.ceil(n_wavelengths * r["lambda_2D"] / dx))


def smooth_pad_Te(Te, pad_width, Te_out=None):
    """
    Pad a 2-D elastic thickness array with a smooth linear taper.

    When a spatially variable Te grid is padded with a constant value before
    being passed to :class:`F2D`, the abrupt step in flexural rigidity D at
    the inner/outer boundary drives spurious deflections via the D-derivative
    terms in the vWC1994 stencil (issue #45).

    This function eliminates that step by linearly blending the inner-domain
    edge values toward *Te_out* across the padding ring, reducing the
    rigidity gradient at the inner/outer boundary by a factor of ~pad_width
    compared with an abrupt step.

    The corresponding surface load array should be padded with zeros, e.g.::

        qs_padded = numpy.pad(qs, pad_width, mode='constant')

    Parameters
    ----------
    Te : (M, N) array
        Elastic thickness [m] for the inner domain.
    pad_width : int
        Width of the padding ring in grid cells.
    Te_out : float, optional
        Te value at the outer edge of the padding ring.
        Defaults to ``Te.mean()``.

    Returns
    -------
    Te_padded : (M + 2*pad_width, N + 2*pad_width) array
        Padded elastic thickness with a smooth linear taper.
    """
    Te = np.asarray(Te, dtype=float)
    if Te.ndim != 2:
        raise ValueError("Te must be a 2-D array")
    if pad_width < 1:
        raise ValueError("pad_width must be >= 1")

    if Te_out is None:
        Te_out = float(Te.mean())
    Te_out = float(Te_out)

    ny, nx = Te.shape
    p = pad_width
    ny_p, nx_p = ny + 2 * p, nx + 2 * p

    Te_padded = np.full((ny_p, nx_p), Te_out)
    Te_padded[p:-p, p:-p] = Te

    for k in range(p):
        # Linear weight: 0 at the outermost cell, (p-1)/p at the innermost.
        # At k = p-1 the step from padding to inner domain is reduced by 1/p.
        alpha = k / p
        # Top and bottom rows (inner columns only; corners handled below)
        Te_padded[k,          p:-p] = (1.0 - alpha) * Te_out + alpha * Te[0, :]
        Te_padded[ny_p-1-k,   p:-p] = (1.0 - alpha) * Te_out + alpha * Te[-1, :]
        # Left and right columns (inner rows only)
        Te_padded[p:-p,       k   ] = (1.0 - alpha) * Te_out + alpha * Te[:, 0]
        Te_padded[p:-p, nx_p-1-k  ] = (1.0 - alpha) * Te_out + alpha * Te[:, -1]

    # Corner cells: blend toward the nearest inner-domain corner value using
    # the L-infinity distance so the taper connects smoothly to both edges.
    for ky in range(p):
        for kx in range(p):
            alpha = max(ky, kx) / p
            Te_padded[ky,          kx         ] = (1-alpha)*Te_out + alpha*Te[0,  0 ]
            Te_padded[ky,          nx_p-1-kx  ] = (1-alpha)*Te_out + alpha*Te[0,  -1]
            Te_padded[ny_p-1-ky,   kx         ] = (1-alpha)*Te_out + alpha*Te[-1, 0 ]
            Te_padded[ny_p-1-ky,   nx_p-1-kx  ] = (1-alpha)*Te_out + alpha*Te[-1, -1]

    return Te_padded


def _auto_pad_Te_2d(Te, pn, ps, pw, pe, Te_out=None):
    """Asymmetric smooth Te taper for F2D auto-padding.

    Like :func:`smooth_pad_Te` but supports independent widths on each side.
    Each padded side is tapered linearly from the inner-domain edge value to
    *Te_out*.  Corner cells use an L-infinity blend of the two adjacent edge
    tapers.  Sides with a zero pad width are left at *Te_out*.
    """
    Te = np.asarray(Te, dtype=float)
    if Te_out is None:
        Te_out = float(Te.mean())
    ny_i, nx_i = Te.shape
    ny_o, nx_o = ny_i + pn + ps, nx_i + pw + pe
    out = np.full((ny_o, nx_o), Te_out)
    out[pn : pn + ny_i, pw : pw + nx_i] = Te
    # Edge tapers (inner column/row band only; corners overwritten below)
    for k in range(pn):
        out[k, pw : pw + nx_i] = (1 - k / pn) * Te_out + (k / pn) * Te[0, :]
    for k in range(ps):
        out[ny_o - 1 - k, pw : pw + nx_i] = (
            (1 - k / ps) * Te_out + (k / ps) * Te[-1, :]
        )
    for k in range(pw):
        out[pn : pn + ny_i, k] = (1 - k / pw) * Te_out + (k / pw) * Te[:, 0]
    for k in range(pe):
        out[pn : pn + ny_i, nx_o - 1 - k] = (
            (1 - k / pe) * Te_out + (k / pe) * Te[:, -1]
        )
    # Corner blending (L-infinity distance)
    for ky in range(pn):
        ay = ky / pn
        for kx in range(pw):
            out[ky, kx] = (1 - max(ay, kx / pw)) * Te_out + max(ay, kx / pw) * Te[0, 0]
        for kx in range(pe):
            out[ky, nx_o - 1 - kx] = (
                (1 - max(ay, kx / pe)) * Te_out + max(ay, kx / pe) * Te[0, -1]
            )
    for ky in range(ps):
        ay = ky / ps
        for kx in range(pw):
            out[ny_o - 1 - ky, kx] = (
                (1 - max(ay, kx / pw)) * Te_out + max(ay, kx / pw) * Te[-1, 0]
            )
        for kx in range(pe):
            out[ny_o - 1 - ky, nx_o - 1 - kx] = (
                (1 - max(ay, kx / pe)) * Te_out + max(ay, kx / pe) * Te[-1, -1]
            )
    return out


def _pad_domain_2d(Te, qs, dx, dy=None, n_wavelengths=1.0, Te_out=None,
                   E=65e9, nu=0.25, rho_m=3300.0, rho_fill=0.0, g=9.8):
    """Pad a 2-D domain. Called by :func:`pad_domain`; see that function for docs."""
    if dy is None:
        dy = dx
    p = recommended_pad_width(
        Te, min(dx, dy), E=E, nu=nu, rho_m=rho_m, rho_fill=rho_fill,
        g=g, n_wavelengths=n_wavelengths,
    )
    Te_arr = np.asarray(Te, dtype=float)
    if Te_arr.ndim == 0:
        Te_padded = float(Te_arr)
    else:
        Te_padded = smooth_pad_Te(Te_arr, p, Te_out=Te_out)
    qs_padded = np.pad(qs, p, mode="constant")
    return Te_padded, qs_padded, p


def pad_domain(Te, qs, dx, dy=None, n_wavelengths=1.0, Te_out=None,
               E=65e9, nu=0.25, rho_m=3300.0, rho_fill=0.0, g=9.8):
    """
    Pad a flexure domain for use with :class:`F1D` or :class:`F2D`.

    Dispatches to a 1-D or 2-D implementation based on the shape of *qs*.
    For array *Te*, the padding region is tapered from the inner-domain edge
    values toward *Te_out* to avoid an abrupt rigidity step.  For scalar *Te*,
    only the load array is zero-padded; *Te* is returned unchanged as a float.

    The returned pad width *p* can be used to trim the deflection output after
    the run::

        # 2-D
        w_inner = flex.w[p:-p, p:-p]
        # 1-D
        w_inner = flex.w[p:-p]

    Parameters
    ----------
    Te : scalar or array
        Elastic thickness [m].  A scalar is broadcast to the full padded grid
        internally by :class:`F1D` / :class:`F2D`; no tapering is applied.
        A 1-D array is expected when *qs* is 1-D; a 2-D array when *qs* is 2-D.
    qs : 1-D or 2-D array
        Surface load [Pa] for the inner domain.  Shape determines whether
        1-D or 2-D padding is applied.
    dx : float
        Grid cell size [m].  For 2-D grids, this is the x-direction spacing.
    dy : float, optional
        Grid cell size in the y-direction [m].  2-D only; ignored for 1-D.
        Defaults to *dx*.
    n_wavelengths : float, optional
        Padding width expressed as a number of flexural wavelengths.
        Default 1.0; use 0.5 for a less conservative (narrower) padding.
    Te_out : float, optional
        Te value at the outer edge of the padding region.
        Only used for array *Te*; defaults to ``Te.mean()``.
    E : float, optional
        Young's modulus [Pa].  Default 65 GPa.
    nu : float, optional
        Poisson's ratio.  Default 0.25.
    rho_m : float, optional
        Mantle density [kg m^-3].  Default 3300.
    rho_fill : float, optional
        Infill density [kg m^-3].  Default 0 (air).
    g : float, optional
        Gravitational acceleration [m s^-2].  Default 9.8.

    Returns
    -------
    Te_padded : float or array
        Elastic thickness for the padded domain.  Float when *Te* is scalar;
        array of shape ``(len(qs) + 2p,)`` (1-D) or ``(M+2p, N+2p)`` (2-D)
        when *Te* is an array.
    qs_padded : array
        Surface load zero-padded to the padded domain shape.
    p : int
        Pad width in grid cells (same on both ends / all four sides).

    Examples
    --------
    2-D scalar Te:

    >>> import numpy as np
    >>> from gflex import pad_domain
    >>> qs = np.zeros((5, 5))
    >>> Te_pad, qs_pad, p = pad_domain(10e3, qs, dx=10000.,
    ...     E=65e9, nu=0.25, rho_m=3300., rho_fill=0., g=9.8)
    >>> p
    13
    >>> qs_pad.shape
    (31, 31)
    >>> Te_pad   # scalar returned unchanged
    10000.0

    2-D array Te:

    >>> Te = 10e3 * np.ones((5, 5))
    >>> Te_pad, qs_pad, p = pad_domain(Te, qs, dx=10000.,
    ...     E=65e9, nu=0.25, rho_m=3300., rho_fill=0., g=9.8)
    >>> Te_pad.shape
    (31, 31)

    1-D:

    >>> qs_1d = np.zeros(5)
    >>> Te_pad, qs_pad, p = pad_domain(10e3, qs_1d, dx=10000.)
    >>> qs_pad.shape
    (43,)
    """
    qs = np.asarray(qs)
    if qs.ndim == 1:
        from gflex.f1d import _pad_domain_1d
        return _pad_domain_1d(
            Te, qs, dx, n_wavelengths=n_wavelengths, Te_out=Te_out,
            E=E, nu=nu, rho_m=rho_m, rho_fill=rho_fill, g=g,
        )
    return _pad_domain_2d(
        Te, qs, dx, dy=dy, n_wavelengths=n_wavelengths, Te_out=Te_out,
        E=E, nu=nu, rho_m=rho_m, rho_fill=rho_fill, g=g,
    )


class F2D(Flexure):
    """
    Two-dimensional lithospheric flexure solver.

    Computes the deflection *w(x, y)* of a thin elastic plate overlying an
    inviscid fluid (mantle) given a surface load stress *qs*.  Supports
    spatially variable elastic thickness *Te*.

    Set instance attributes, then call :meth:`initialize`, :meth:`run`, and
    :meth:`finalize` in sequence.  Read ``flex.w`` **before** calling
    :meth:`finalize`; finalize clears all model state including ``w``.

    Attributes
    ----------
    method : str
        Solution method.  ``'fd'`` (finite difference, supports variable
        *Te*), ``'fft'`` (spectral, requires scalar *Te*; 2-D only),
        ``'sas'`` (superposition of analytical solutions, constant *Te*
        only), or ``'sas_ng'`` (SAS on an ungridded point cloud).
    solver : str
        Linear solver: ``'direct'`` (sparse LU, default).
    g : float
        Gravitational acceleration [m s⁻²].
    E : float
        Young's modulus [Pa].
    nu : float
        Poisson's ratio.
    rho_m : float
        Mantle density [kg m⁻³].
    rho_fill : float
        Infill material density [kg m⁻³] (0 for air, ~1000 for water,
        ~2700 for rock).
    T_e : float or ndarray of shape (M, N)
        Elastic thickness [m].  A scalar is broadcast to the full grid.
    qs : ndarray of shape (M, N)
        Surface load stress [Pa].
    dx : float
        Grid spacing in the x (column) direction [m].
    dy : float
        Grid spacing in the y (row) direction [m].
    bc_west, bc_east, bc_north, bc_south : str
        Boundary conditions on the west, east, north, and south edges.
        FD options: ``'zero_displacement_zero_slope'`` (alias ``'clamped'``),
        ``'zero_displacement_zero_moment'`` (alias ``'pinned'``),
        ``'zero_moment_zero_shear'`` (alias ``'free'``),
        ``'zero_slope_zero_shear'`` (alias ``'mirror'``), ``'periodic'``,
        ``'no_outside_loads'`` (auto-pad by one flexural wavelength and apply
        ``'zero_displacement_zero_slope'`` at the new outer edges; ``self.w``
        is trimmed to the original domain).
        SAS option: ``'no_outside_loads'`` (the default when unset).
        FFT: each opposite pair (west/east, north/south) is treated
        independently — set both sides of a pair to ``'periodic'`` for
        exact periodicity along that axis, or leave them unset (or set to
        ``'no_outside_loads'``) to zero-pad that axis.  Mixing periodic
        and non-periodic axes is valid.  Setting only one side of a pair
        to ``'periodic'`` raises a ``UserWarning`` and falls back to
        zero-padding for that axis.
    sigma_xx : float, optional
        Normal in-plane stress in the x-direction :math:`\\sigma_{xx}` [Pa].
        Supported by ``fd`` and ``fft``.  Default ``0``.
    sigma_yy : float, optional
        Normal in-plane stress in the y-direction :math:`\\sigma_{yy}` [Pa].
        Supported by ``fd`` and ``fft``.  Default ``0``.
    sigma_xy : float, optional
        In-plane shear stress :math:`\\sigma_{xy}` [Pa].
        Supported by ``fd`` and ``fft``.  Default ``0``.
    fft_pad_n_alpha : int or float
        Number of 2-D flexural-parameter units (α₂D = (D/Δρg)^0.25) to
        zero-pad on each side for non-periodic FFT runs.  Periodic images
        of the load are separated by 2 × ``fft_pad_n_alpha`` × α₂D.
        Default ``4`` (8α₂D total separation).  Ignored when
        ``method != 'fft'`` or when both sides of the relevant axis are
        ``'periodic'``.
    quiet : bool
        Suppress timing output.  Default ``False``.
    verbose : bool
        Print progress messages.  Default ``True``.

    Examples
    --------
    Minimal finite-difference run::

        import numpy as np
        from gflex import F2D

        flex = F2D()
        flex.quiet = True
        flex.method = 'fd'
        flex.solver = 'direct'
        flex.g = 9.8
        flex.E = 65e9
        flex.nu = 0.25
        flex.rho_m = 3300.
        flex.rho_fill = 0.
        flex.T_e = 30e3 * np.ones((50, 50))
        flex.qs = np.zeros((50, 50))
        flex.qs[20:30, 20:30] = 1e6   # 10 × 10 cell load
        flex.dx = flex.dy = 5000.     # 5 km grid
        flex.bc_west = flex.bc_east = flex.bc_south = flex.bc_north = 'zero_moment_zero_shear'
        flex.initialize()
        flex.run()
        flex.finalize()
        deflection = flex.w           # (50, 50) array, negative downward
    """

    def initialize(self, filename=None):
        """
        Validate inputs and prepare the solver.

        Must be called once before :meth:`run`.  If a configuration-file
        path was passed to the constructor (or to this method), parameters
        are read from that file; otherwise they are taken from the instance
        attributes set by the caller.

        Parameters
        ----------
        filename : str, optional
            Path to a gFlex YAML configuration file.  Overrides
            any filename supplied to the constructor.
        """
        self.dimension = 2  # Set it here in case it wasn't set for selection before
        self._total_start_time = time.time()
        super().initialize()
        _logger.info("F2D initialized")
        if hasattr(self, "qs") and hasattr(self, "dx") and hasattr(self, "dy"):
            _logger.info(
                "  Grid: ny=%d, nx=%d, dy=%.1f m, dx=%.1f m",
                self.qs.shape[0], self.qs.shape[1], self.dy, self.dx,
            )
        if hasattr(self, "T_e"):
            if np.ndim(self.T_e) == 0:
                _logger.info("  T_e: %.1f m", float(self.T_e))
            else:
                _logger.info(
                    "  T_e: array, shape=%s, range=[%.1f, %.1f] m",
                    self.T_e.shape, np.min(self.T_e), np.max(self.T_e),
                )
        if hasattr(self, "method"):
            _logger.info("  Method: %s", self.method)
        _logger.info("")

    def run(self):
        """
        Execute the flexural solution.

        Selects and runs the method specified by ``self.method``.  The
        deflection array is stored in ``self.w`` on return.  Call
        :meth:`finalize` afterwards to restore any internally modified
        state.

        For repeated solves (e.g. a coupling loop), set
        ``cache_factorization = True`` before :meth:`initialize` to reuse the
        LU factorisation when only ``qs`` changes; use ``"no_check"`` for
        maximum throughput when the coefficient matrix is guaranteed stable.
        """
        self.bc_check()
        self.solver_start_time = time.perf_counter()

        if self.method == "fd":
            # Finite difference
            super()._solve_fd()
            self.method_func = self._solve_fd
        elif self.method == "fft":
            # Fast Fourier transform
            super()._solve_fft()
            self.method_func = self._solve_fft
        elif self.method == "sas":
            # Superposition of analytical solutions
            super()._solve_sas()
            self.method_func = self._solve_sas
        elif self.method == "sas_ng":
            # Superposition of analytical solutions,
            # nonuniform points (no grid)
            super()._solve_sas_ng()
            self.method_func = self._solve_sas_ng
        else:
            raise ValueError('method must be "fd", "fft", "sas", or "sas_ng"')

        _logger.info("F2D run")
        self.method_func()

        self.time_to_solve = time.perf_counter() - self.solver_start_time
        if hasattr(self, "linear_solve_time"):
            _logger.info("  Time for linear solve [s]: %.6f", self.linear_solve_time)
        else:
            _logger.info("  Time to solve [s]: %.6f", self.time_to_solve)

    def finalize(self):
        """
        Release all model state.

        Restores ``self.T_e`` to its pre-run value if gFlex padded it
        internally.  Then calls the base ``finalize``, which deletes
        ``w``, ``qs``, and the cached coefficient matrix.  Read ``w``
        before calling this method.
        """
        with contextlib.suppress(AttributeError):
            self.T_e = self.T_e_unpadded
        _logger.info("")
        _logger.info("F2D finalized")
        if hasattr(self, "_total_start_time"):
            _logger.info("  Total runtime [s]: %.6f", time.time() - self._total_start_time)
        _logger.info("")
        super().finalize()

    ########################################
    ## FUNCTIONS FOR EACH SOLUTION METHOD ##
    ########################################

    def _check_warnings_FD(self):
        """Issue UserWarnings for potentially problematic FD boundary conditions.

        Two categories of warning are raised:

        **BC-type warnings** — fired for boundary types whose physical meaning
        deserves verification: ``'zero_moment_zero_shear'`` (free broken end; valid for
        rifted/passive margins, subduction trenches with an edge load, and
        broken-plate flexure).

        **Proximity warnings** — fired for ``'zero_displacement_zero_slope'`` boundaries
        when the nearest loaded cell is within one flexural wavelength
        (:math:`\\lambda = 2\\pi\\alpha`, :math:`\\alpha = (4D/\\Delta\\rho g)^{1/4}`)
        of that boundary.  Within this distance the flexural forebulge is
        suppressed by the zero-displacement condition.  The warning message
        reports the distance as a fraction of the local flexural wavelength and
        points to the domain-padding utilities.
        """
        bc_sides = {"W": self.bc_west, "E": self.bc_east, "N": self.bc_north, "S": self.bc_south}
        for side, bc in bc_sides.items():
            if bc == "zero_moment_zero_shear":
                warnings.warn(
                    f"BC_{side} = 'zero_moment_zero_shear': assumes a free broken plate end "
                    "(zero moment and shear force). Valid for rifted/passive margins, "
                    "subduction trenches with an applied edge load, and broken-plate "
                    "flexure (Turcotte & Schubert). Verify this is physically "
                    "appropriate for your setup.",
                    UserWarning,
                    stacklevel=4,
                )
        # Load-proximity warning: only for zero_displacement_zero_slope boundaries
        loaded = np.argwhere(self.qs != 0)
        if loaded.size == 0:
            return
        ny, nx = self.qs.shape
        Te_arr = (
            self.T_e
            if isinstance(self.T_e, np.ndarray)
            else np.full((ny, nx), float(self.T_e))
        )
        rows, cols = loaded[:, 0], loaded[:, 1]
        Te_loaded = Te_arr[rows, cols]
        D_loaded = self.E * Te_loaded**3 / (12 * (1 - self.nu**2))
        alpha_loaded = (4 * D_loaded / (self.drho * self.g)) ** 0.25
        wavelength_loaded = 2 * np.pi * alpha_loaded
        dist_fns = {
            "W": lambda: (cols + 0.5) * self.dx,
            "E": lambda: (nx - cols - 0.5) * self.dx,
            "N": lambda: (rows + 0.5) * self.dy,
            "S": lambda: (ny - rows - 0.5) * self.dy,
        }
        for side, bc in bc_sides.items():
            if bc != "zero_displacement_zero_slope":
                continue
            distances = dist_fns[side]()
            frac = distances / wavelength_loaded
            worst = np.argmin(frac)
            if frac[worst] < 1.0:
                r, c = rows[worst], cols[worst]
                warnings.warn(
                    f"BC_{side} = 'zero_displacement_zero_slope': nearest loaded cell "
                    f"(row {r}, col {c}) is {distances[worst]/1e3:.1f} km from the boundary "
                    f"({frac[worst]:.2f} flexural wavelengths; "
                    f"wavelength ≈ {wavelength_loaded[worst]/1e3:.1f} km "
                    f"at Te = {Te_loaded[worst]/1e3:.1f} km). "
                    "The flexural forebulge peaks near one wavelength from the load "
                    "and will be suppressed by this boundary. "
                    "Use pad_domain() to extend the domain before solving, "
                    "then trim w to the original extent.",
                    UserWarning,
                    stacklevel=4,
                )

    def _solve_fd(self):
        """Run the finite-difference solution pipeline.

        If any BC is ``'no_outside_loads'``, auto-pads those sides by one
        flexural wavelength (zero load, tapered Te if array), applies
        ``'zero_displacement_zero_slope'`` at the new outer edges, solves on
        the extended domain, then crops ``self.w`` back to the original shape.

        For all other BCs, calls :meth:`_check_warnings_FD`, assembles and
        solves the sparse banded system.  The deflection is available in
        ``self.w`` after the call.
        """
        # ------------------------------------------------------------------
        # Auto-pad any 'no_outside_loads' sides.
        # ------------------------------------------------------------------
        _pad_north = self._bc_north_norm == "no_outside_loads"
        _pad_south = self._bc_south_norm == "no_outside_loads"
        _pad_west  = self._bc_west_norm  == "no_outside_loads"
        _pad_east  = self._bc_east_norm  == "no_outside_loads"
        _qs_inner = _Te_inner = None
        pn = ps = pw = pe = 0

        if _pad_north or _pad_south or _pad_west or _pad_east:
            p = recommended_pad_width(
                self.T_e, min(self.dx, self.dy),
                E=self.E, nu=self.nu, rho_m=self.rho_m,
                rho_fill=self.rho_fill, g=self.g,
            )
            pn = p if _pad_north else 0
            ps = p if _pad_south else 0
            pw = p if _pad_west  else 0
            pe = p if _pad_east  else 0
            _qs_inner = self.qs.copy()
            _Te_inner = self.T_e
            self.qs = np.pad(self.qs, ((pn, ps), (pw, pe)), mode="constant")
            if not np.isscalar(self.T_e):
                self.T_e = _auto_pad_Te_2d(
                    np.asarray(self.T_e, dtype=float), pn, ps, pw, pe
                )
            if _pad_north: self._bc_north_norm = "zero_displacement_zero_slope"
            if _pad_south: self._bc_south_norm = "zero_displacement_zero_slope"
            if _pad_west:  self._bc_west_norm  = "zero_displacement_zero_slope"
            if _pad_east:  self._bc_east_norm  = "zero_displacement_zero_slope"
            sides = []
            if _pad_north: sides.append(f"north +{pn}")
            if _pad_south: sides.append(f"south +{ps}")
            if _pad_west:  sides.append(f"west +{pw}")
            if _pad_east:  sides.append(f"east +{pe}")
            _logger.info(
                "  FD 'no_outside_loads': auto-padding (%s cells).",
                ", ".join(sides),
            )

        self._check_warnings_FD()
        # Only generate coefficient matrix if it is not already provided.
        # In no_check mode the matrix is freed after factorization to save
        # memory; _lu being set is sufficient to skip the rebuild.
        if self.coeff_matrix is not None:
            pass
        elif self.cache_factorization == "no_check" and self._lu is not None:
            pass  # coeff_matrix was freed after factorization; _lu still valid
        else:
            self.elasprep()
            self._build_coefficient_matrix()
        self.fd_solve()

        if pn or ps or pw or pe:
            ny_i, nx_i = _qs_inner.shape
            self.w = self.w[pn : pn + ny_i, pw : pw + nx_i]
            self.qs = _qs_inner
            self.T_e = _Te_inner
            # Restore grid metadata to match the cropped domain.
            self.ny, self.nx = ny_i, nx_i
            if hasattr(self, "T_e_unpadded"):
                self.T_e_unpadded = _Te_inner

    def _solve_fft(self):
        """Spectral (FFT) flexural solution for uniform elastic thickness.

        Applies the analytical 2-D transfer function in the wavenumber domain::

            W(kx, ky) = -Q(kx, ky) / (D(kx²+ky²)²
                        + σ_xx·Te·kx² + σ_yy·Te·ky² + 2σ_xy·Te·kx·ky + Δρg)

        **Boundary conditions and periodicity**

        Each opposite pair of boundaries (west/east, north/south) is handled
        independently:

        * Both sides of a pair set to ``'periodic'`` — that axis is treated as
          genuinely periodic and the load is used as-is along that dimension.

        * Any other value (including ``'no_outside_loads'`` or unset) — that
          axis is zero-padded by ``fft_pad_n_alpha`` × α₂D on each side
          (default 4α₂D), where α₂D = (D/Δρg)^0.25 is the 2-D flexural
          parameter.  Periodic images are separated by at least
          2 × ``fft_pad_n_alpha`` × α₂D of zeros (default 8α₂D).

        Mixing periodic and non-periodic axes is supported: e.g.
        ``bc_west = bc_east = 'periodic'`` with ``bc_north`` and ``bc_south``
        unset gives an x-periodic, y-padded domain.  Setting only *one* side
        of a pair to ``'periodic'`` raises a ``UserWarning`` and the axis
        falls back to zero-padding.

        Requires uniform (scalar) elastic thickness; for variable *Te* use
        the finite-difference method instead.
        """
        if np.isscalar(self.T_e):
            pass
        elif np.all(self.T_e == np.mean(self.T_e)):
            self.T_e = float(np.mean(self.T_e))
        else:
            raise ValueError(
                "The FFT solution requires a scalar (uniform) Te. "
                "For spatially variable Te, use the finite difference method."
            )

        D = self.E * self.T_e**3 / (12.0 * (1.0 - self.nu**2))
        alpha = flexural_wavelengths(
            self.T_e, E=self.E, nu=self.nu,
            rho_m=self.rho_m, rho_fill=self.rho_fill, g=self.g,
        )["alpha_2D"]

        ny, nx = self.qs.shape

        # Per-axis periodicity: each opposite pair is independent.
        periodic_x = (self.bc_west  == "periodic") and (self.bc_east  == "periodic")
        periodic_y = (self.bc_north == "periodic") and (self.bc_south == "periodic")

        # Warn when exactly one side of a pair is periodic (non-physical).
        for (bc_a, bc_b, name_a, name_b) in [
            (self.bc_west,  self.bc_east,  "bc_west",  "bc_east"),
            (self.bc_north, self.bc_south, "bc_north", "bc_south"),
        ]:
            if (bc_a == "periodic") != (bc_b == "periodic"):
                non_periodic = name_a if bc_a != "periodic" else name_b
                periodic_side = name_b if bc_a != "periodic" else name_a
                warnings.warn(
                    f"FFT method: {non_periodic} is not 'periodic' but "
                    f"{periodic_side} is — falling back to no_outside_loads "
                    "zero-padding for this axis.  Set both sides of each pair "
                    "to 'periodic' for exact periodic behavior, or leave both "
                    "unset for no_outside_loads.",
                    UserWarning, stacklevel=4,
                )

        pad_x = 0 if periodic_x else int(np.ceil(self.fft_pad_n_alpha * alpha / self.dx))
        pad_y = 0 if periodic_y else int(np.ceil(self.fft_pad_n_alpha * alpha / self.dy))
        qs_work = np.pad(self.qs, ((pad_y, pad_y), (pad_x, pad_x)), mode="constant")

        ny_work, nx_work = qs_work.shape
        kx = scipy.fft.rfftfreq(nx_work, d=self.dx) * 2.0 * np.pi
        ky = scipy.fft.fftfreq(ny_work, d=self.dy) * 2.0 * np.pi
        Kx, Ky = np.meshgrid(kx, ky)

        Q = scipy.fft.rfft2(qs_work, workers=-1)
        K2 = Kx**2 + Ky**2
        denom = (
            D * K2**2
            + self.sigma_xx * self.T_e * Kx**2
            + self.sigma_yy * self.T_e * Ky**2
            + 2.0 * self.sigma_xy * self.T_e * Kx * Ky
            + self.drho * self.g
        )
        w_work = scipy.fft.irfft2(-Q / denom, s=qs_work.shape, workers=-1)

        self.w = w_work[pad_y : pad_y + ny, pad_x : pad_x + nx]

    def _solve_sas(self):
        """Run the gridded superposition-of-analytical-solutions pipeline."""
        self.spatial_domain_vars_sas()
        self.spatial_domain_gridded()

    def _solve_sas_ng(self):
        """Run the ungridded (non-uniform points) SAS pipeline."""
        self.spatial_domain_vars_sas()
        self.spatial_domain_no_grid()

    ######################################
    ## FUNCTIONS TO SOLVE THE EQUATIONS ##
    ######################################

    ## SPATIAL DOMAIN SUPERPOSITION OF ANALYTICAL SOLUTIONS
    #########################################################

    # SETUP

    def spatial_domain_vars_sas(self):
        """Compute flexural rigidity D, parameter alpha, and Green's-function coefficient for SAS."""
        # Check Te:
        # * If scalar, okay.
        # * If grid, convert to scalar if a singular value
        # * Else, throw an error.
        if np.isscalar(self.T_e):
            pass
        elif np.all(self.T_e == np.mean(self.T_e)):
            self.T_e = np.mean(self.T_e)
        else:
            raise ValueError(
                "The analytical solution requires a scalar (uniform) Te. "
                "For spatially variable Te, use the finite difference method."
            )

        self.D = self.E * self.T_e**3 / (12 * (1 - self.nu**2))  # Flexural rigidity
        self.alpha = (self.D / (self.drho * self.g)) ** 0.25  # 2D flexural parameter
        self.coeff = self.alpha**2 / (2 * np.pi * self.D)

    # GRIDDED

    def spatial_domain_gridded(self):
        """Compute deflection by summing 2D Kelvin-function Green's functions over the load grid."""
        self.nx = self.qs.shape[1]
        self.ny = self.qs.shape[0]

        # Prepare a large grid of solutions beforehand, so we don't have to
        # keep calculating kei (time-intensive!)
        # This pre-prepared solution will be for a unit load
        bigshape = 2 * self.ny + 1, 2 * self.nx + 1  # Tuple shape

        dist_ny = np.arange(bigshape[0]) - self.ny
        dist_nx = np.arange(bigshape[1]) - self.nx

        dist_x, dist_y = np.meshgrid(dist_nx * self.dx, dist_ny * self.dy)

        bigdist = np.sqrt(dist_x**2 + dist_y**2)  # Distances from center
        # Center at [ny,nx]

        biggrid = self.coeff * kei(bigdist / self.alpha)  # Kelvin fcn solution

        # Convolve the load field with the Green's function kernel.
        # biggrid[1:-1, 1:-1] trims the outer ring to give a kernel of shape
        # (2*ny-1, 2*nx-1) centred at [ny-1, nx-1], covering every relative
        # offset that can occur between two cells in the (ny, nx) grid.
        # fftconvolve with mode='same' computes the same sum as the original
        # double loop — each loaded cell's contribution is a phase-shifted copy
        # of the Green's function — but does so in O(N² log N) via the FFT
        # rather than O(N_load × N_grid) Python iterations.
        kernel = biggrid[1:-1, 1:-1]
        self.w = fftconvolve(self.qs * self.dx * self.dy, kernel, mode="same")

    # NO GRID

    def spatial_domain_no_grid(self):
        """Compute deflection by summing 2D Kelvin-function Green's functions at ungridded output points."""
        self.w = np.zeros(self.xw.shape)
        _logger.debug("w = ")
        _logger.debug("%s", self.w.shape)

        # r shape: (N_out, N_load); vectorised over both dimensions at once.
        # Memory cost is O(N_out × N_load) floats — batch if that is too large.
        if self.latlon:
            r = self.greatCircleDistance(
                lat1=self.y[None, :],
                long1=self.x[None, :],
                lat2=self.yw[:, None],
                long2=self.xw[:, None],
                radius=self.planetary_radius,
            )   # (N_out, N_load)
        else:
            r = np.sqrt(
                (self.xw[:, None] - self.x[None, :]) ** 2
                + (self.yw[:, None] - self.y[None, :]) ** 2
            )   # (N_out, N_load)
        self.w = self.coeff * (kei(r / self.alpha) * self.q[None, :]).sum(axis=1)

    ## FINITE DIFFERENCE
    ######################

    def elasprep(self):
        """Precompute dx⁴, dy⁴, dx²dy², and the flexural rigidity array D for the FD solver."""

        if self.method != "sas_ng":
            self.dx4 = self.dx**4
            self.dy4 = self.dy**4
            self.dx2dy2 = self.dx**2 * self.dy**2
        self.D = self.E * self.T_e**3 / (12 * (1 - self.nu**2))

    def _build_coefficient_matrix(self):
        """
        Selects the boundary conditions
        E-W is for inside each panel
        N-S is for the block diagonal matrix ("with fringes")
        Then calls the function to build the diagonal matrix

        The current method of coefficient matrix construction utilizes longer-range
        symmetry in the coefficient matrix to build it block-wise, as opposed to
        the much less computationally efficient row-by-row ("serial") method
        that was previously employed.

        The method is spread across the subroutines here.

        Important to this is the use of np.roll() to properly offset the
        diagonals that end up in the main matrix: spdiags() will put each vector
        on the proper diagonal, but will align them such that their first cell is
        along the first column, instead of using a 45 degrees to matrix corner
        baseline that would stagger them appropriately for this solution method.
        Therefore, np.roll() effectively does this staggering by having the
        appropriate cell in the vector start at the first column.
        """

        # Zeroth, start the timer and print the boundary conditions to the screen
        self.coeff_start_time = time.time()
        _logger.info("  Boundary condition, West: %s", self.bc_west)
        _logger.info("  Boundary condition, East: %s", self.bc_east)
        _logger.info("  Boundary condition, North: %s", self.bc_north)
        _logger.info("  Boundary condition, South: %s", self.bc_south)

        # First, set flexural rigidity boundary conditions to flesh out this padded
        # array
        self._apply_bc_rigidity()

        # Second, build the coefficient arrays -- with the rigidity b.c.'s
        self.get_coeff_values()

        # Third, apply boundary conditions to the coeff_arrays to create the
        # flexural solution
        self._apply_bc_flexure()

        # Fourth, compute RHS correction for prescribed (nonzero) BC values
        self._apply_bc_rhs_inhomogeneous_2d()

        # Fifth, construct the sparse diagonal array
        self.build_diagonals()

        # Finally, compute the total time this process took
        self.coeff_creation_time = time.time() - self.coeff_start_time
        _logger.info(
            "  Time to construct coefficient (operator) array [s]: %.6f",
            self.coeff_creation_time,
        )

    def _apply_bc_rigidity(self):
        """
        Utility function to help implement boundary conditions by specifying
        them for and applying them to the elastic thickness grid
        """

        #########################################
        # FLEXURAL RIGIDITY BOUNDARY CONDITIONS #
        #########################################
        # West
        if self._bc_west_norm == "periodic":
            self.bc_rigidity_west = _RigidityBC.PERIODIC
        elif (
            self._bc_west_norm
            == np.array(["zero_displacement_zero_slope", "zero_moment_zero_shear"])
        ).any():
            self.bc_rigidity_west = _RigidityBC.ZERO_CURVATURE
        elif self._bc_west_norm in ("zero_slope_zero_shear", "zero_displacement_zero_moment"):
            self.bc_rigidity_west = _RigidityBC.MIRROR
        else:
            raise RuntimeError("Invalid Te B.C. case")
        # East
        if self._bc_east_norm == "periodic":
            self.bc_rigidity_east = _RigidityBC.PERIODIC
        elif (
            self._bc_east_norm
            == np.array(["zero_displacement_zero_slope", "zero_moment_zero_shear"])
        ).any():
            self.bc_rigidity_east = _RigidityBC.ZERO_CURVATURE
        elif self._bc_east_norm in ("zero_slope_zero_shear", "zero_displacement_zero_moment"):
            self.bc_rigidity_east = _RigidityBC.MIRROR
        else:
            raise RuntimeError("Invalid Te B.C. case")
        # North
        if self._bc_north_norm == "periodic":
            self.bc_rigidity_north = _RigidityBC.PERIODIC
        elif (
            self._bc_north_norm
            == np.array(["zero_displacement_zero_slope", "zero_moment_zero_shear"])
        ).any():
            self.bc_rigidity_north = _RigidityBC.ZERO_CURVATURE
        elif self._bc_north_norm in ("zero_slope_zero_shear", "zero_displacement_zero_moment"):
            self.bc_rigidity_north = _RigidityBC.MIRROR
        else:
            raise RuntimeError("Invalid Te B.C. case")
        # South
        if self._bc_south_norm == "periodic":
            self.bc_rigidity_south = _RigidityBC.PERIODIC
        elif (
            self._bc_south_norm
            == np.array(["zero_displacement_zero_slope", "zero_moment_zero_shear"])
        ).any():
            self.bc_rigidity_south = _RigidityBC.ZERO_CURVATURE
        elif self._bc_south_norm in ("zero_slope_zero_shear", "zero_displacement_zero_moment"):
            self.bc_rigidity_south = _RigidityBC.MIRROR
        else:
            raise RuntimeError("Invalid Te B.C. case")

        #############
        # PAD ARRAY #
        #############
        if np.isscalar(self.T_e):
            self.D *= np.ones(self.qs.shape)  # And leave Te as a scalar for checks
        else:
            self.T_e_unpadded = self.T_e.copy()
            self.T_e = np.hstack(
                (
                    np.nan * np.zeros((self.T_e.shape[0], 1)),
                    self.T_e,
                    np.nan * np.zeros((self.T_e.shape[0], 1)),
                )
            )
            self.T_e = np.vstack(
                (
                    np.nan * np.zeros(self.T_e.shape[1]),
                    self.T_e,
                    np.nan * np.zeros(self.T_e.shape[1]),
                )
            )
            self.D = np.hstack(
                (
                    np.nan * np.zeros((self.D.shape[0], 1)),
                    self.D,
                    np.nan * np.zeros((self.D.shape[0], 1)),
                )
            )
            self.D = np.vstack(
                (
                    np.nan * np.zeros(self.D.shape[1]),
                    self.D,
                    np.nan * np.zeros(self.D.shape[1]),
                )
            )

        ###############################################################
        # APPLY FLEXURAL RIGIDITY BOUNDARY CONDITIONS TO PADDED ARRAY #
        ###############################################################
        if self.bc_rigidity_west == _RigidityBC.ZERO_CURVATURE:
            self.D[:, 0] = 2 * self.D[:, 1] - self.D[:, 2]
        if self.bc_rigidity_east == _RigidityBC.ZERO_CURVATURE:
            self.D[:, -1] = 2 * self.D[:, -2] - self.D[:, -3]
        if self.bc_rigidity_north == _RigidityBC.ZERO_CURVATURE:
            self.D[0, :] = 2 * self.D[1, :] - self.D[2, :]
        if self.bc_rigidity_south == _RigidityBC.ZERO_CURVATURE:
            self.D[-1, :] = 2 * self.D[-2, :] - self.D[-3, :]

        if self.bc_rigidity_west == _RigidityBC.MIRROR:
            self.D[:, 0] = self.D[:, 2]
        if self.bc_rigidity_east == _RigidityBC.MIRROR:
            self.D[:, -1] = self.D[:, -3]
        if self.bc_rigidity_north == _RigidityBC.MIRROR:
            self.D[0, :] = self.D[
                2, :
            ]  # Yes, will work on corners -- double-reflection
        if self.bc_rigidity_south == _RigidityBC.MIRROR:
            self.D[-1, :] = self.D[-3, :]

        if self.bc_rigidity_west == _RigidityBC.PERIODIC:
            self.D[:, 0] = self.D[:, -2]
        if self.bc_rigidity_east == _RigidityBC.PERIODIC:
            self.D[:, -1] = self.D[:, -3]
        if self.bc_rigidity_north == _RigidityBC.PERIODIC:
            self.D[0, :] = self.D[-2, :]
        if self.bc_rigidity_south == _RigidityBC.PERIODIC:
            self.D[-1, :] = self.D[-3, :]

    def get_coeff_values(self):
        """
        Calculates the matrix of coefficients that is later used via sparse matrix
        solution techniques (scipy.sparse.linalg.spsolve) to compute the flexural
        response to the load. This step need only be performed once, and the
        coefficient matrix can very rapidly compute flexural solutions to any load.
        This makes this particularly good for probelms with time-variable loads or
        that require iteration (e.g., water loading, in which additional water
        causes subsidence, causes additional water detph, etc.).

        These must be linearly combined to solve the equation.

        13 coefficients: 13 matrices of the same size as the load

        NOTATION FOR COEFFICIENT BIULDING MATRICES (e.g., "cj0i_2"):
        c = "coefficient
        j = columns = x-value
        j0 = no offset: look at center of array
        i = rows = y-value
        i_2 = negative 2 offset (i2 = positive 2 offset)
        """

        # don't want to keep typing "self." everwhere!
        D = self.D
        drho = self.drho
        dx4 = self.dx4
        dy4 = self.dy4
        dx2dy2 = self.dx2dy2
        nu = self.nu
        g = self.g

        if np.isscalar(self.T_e):
            # So much simpler with constant D! And symmetrical stencil
            dx2 = self.dx ** 2
            dy2 = self.dy ** 2
            dxdy = self.dx * self.dy
            # x±2, y=0
            self.cj2i0 = D / dx4
            self.cj_2i0 = D / dx4                                                    # Symmetry
            # x=0, y±2
            self.cj0i2 = D / dy4                                                     # Symmetry
            self.cj0i_2 = D / dy4
            # x±1, y=0  (σ_xx acts on x-direction: -σ_xx·Te·∂²w/∂x²)
            self.cj1i0 = -4 * D / dx4 - 4 * D / dx2dy2 - self.sigma_xx * self.T_e / dx2
            self.cj_1i0 = -4 * D / dx4 - 4 * D / dx2dy2 - self.sigma_xx * self.T_e / dx2  # Symmetry
            # x=0, y±1  (σ_yy acts on y-direction: -σ_yy·Te·∂²w/∂y²)
            self.cj0i_1 = -4 * D / dy4 - 4 * D / dx2dy2 - self.sigma_yy * self.T_e / dy2
            self.cj0i1 = (
                -4 * D / dy4 - 4 * D / dx2dy2 - self.sigma_yy * self.T_e / dy2
            )                                                                         # Symmetry
            # x±1, y±1  (σ_xy cross term: -2σ_xy·Te·∂²w/∂x∂y; sign from FD of mixed deriv)
            self.cj1i_1 = 2 * D / dx2dy2 + self.sigma_xy * self.T_e / (2 * dxdy)
            self.cj1i1 = 2 * D / dx2dy2 - self.sigma_xy * self.T_e / (2 * dxdy)
            self.cj_1i_1 = 2 * D / dx2dy2 - self.sigma_xy * self.T_e / (2 * dxdy)   # Symmetry
            self.cj_1i1 = 2 * D / dx2dy2 + self.sigma_xy * self.T_e / (2 * dxdy)    # Symmetry
            # center
            self.cj0i0 = (
                6 * D / dx4
                + 6 * D / dy4
                + 8 * D / dx2dy2
                + drho * g
                + 2 * self.sigma_xx * self.T_e / dx2
                + 2 * self.sigma_yy * self.T_e / dy2
            )
            # Bring up to size
            self.cj2i0 *= np.ones(self.qs.shape)
            self.cj1i_1 *= np.ones(self.qs.shape)
            self.cj1i0 *= np.ones(self.qs.shape)
            self.cj1i1 *= np.ones(self.qs.shape)
            self.cj0i_2 *= np.ones(self.qs.shape)
            self.cj0i_1 *= np.ones(self.qs.shape)
            self.cj0i0 *= np.ones(self.qs.shape)
            self.cj0i1 *= np.ones(self.qs.shape)
            self.cj0i2 *= np.ones(self.qs.shape)
            self.cj_1i_1 *= np.ones(self.qs.shape)
            self.cj_1i0 *= np.ones(self.qs.shape)
            self.cj_1i1 *= np.ones(self.qs.shape)
            self.cj_2i0 *= np.ones(self.qs.shape)
            # Create coefficient arrays to manage boundary conditions
            self.cj2i0_coeff_ij = self.cj2i0.copy()
            self.cj1i_1_coeff_ij = self.cj1i_1.copy()
            self.cj1i0_coeff_ij = self.cj1i0.copy()
            self.cj1i1_coeff_ij = self.cj1i1.copy()
            self.cj0i_2_coeff_ij = self.cj0i_2.copy()
            self.cj0i_1_coeff_ij = self.cj0i_1.copy()
            self.cj0i0_coeff_ij = self.cj0i0.copy()
            self.cj0i1_coeff_ij = self.cj0i1.copy()
            self.cj0i2_coeff_ij = self.cj0i2.copy()
            self.cj_1i_1_coeff_ij = self.cj_1i_1.copy()
            self.cj_1i0_coeff_ij = self.cj_1i0.copy()
            self.cj_1i1_coeff_ij = self.cj_1i1.copy()
            self.cj_2i0_coeff_ij = self.cj_2i0.copy()

        elif isinstance(self.T_e, np.ndarray):
            # All derivatives here, to make reading the equations below easier
            D00 = D[1:-1, 1:-1]
            D10 = D[1:-1, 2:]
            D_10 = D[1:-1, :-2]
            D01 = D[2:, 1:-1]
            D0_1 = D[:-2, 1:-1]
            D11 = D[2:, 2:]
            D_11 = D[2:, :-2]
            D1_1 = D[:-2, 2:]
            D_1_1 = D[:-2, :-2]
            # Derivatives of D -- not including /(dx^a dy^b)
            D0 = D00
            Dx = (-D_10 + D10) / 2.0
            Dy = (-D0_1 + D01) / 2.0
            Dxx = D_10 - 2.0 * D00 + D10
            Dyy = D0_1 - 2.0 * D00 + D01
            Dxy = (D_1_1 - D_11 - D1_1 + D11) / 4.0

            # van Wees and Cloetingh (1994) solution, re-discretized by me
            # using a central difference approx. to 2nd order precision
            # x = -2, y = 0
            self.cj_2i0_coeff_ij = (D0 - Dx) / dx4
            # x = 0, y = -2
            self.cj0i_2_coeff_ij = (D0 - Dy) / dy4
            # x = 0, y = 2
            self.cj0i2_coeff_ij = (D0 + Dy) / dy4
            # x = 2, y = 0
            self.cj2i0_coeff_ij = (D0 + Dx) / dx4
            # x = -1, y = -1
            self.cj_1i_1_coeff_ij = (
                2.0 * D0 - Dx - Dy + Dxy * (1 - nu) / 2.0
            ) / dx2dy2
            # x = -1, y = 1
            self.cj_1i1_coeff_ij = (
                2.0 * D0 - Dx + Dy - Dxy * (1 - nu) / 2.0
            ) / dx2dy2
            # x = 1, y = -1
            self.cj1i_1_coeff_ij = (
                2.0 * D0 + Dx - Dy - Dxy * (1 - nu) / 2.0
            ) / dx2dy2
            # x = 1, y = 1
            self.cj1i1_coeff_ij = (
                2.0 * D0 + Dx + Dy + Dxy * (1 - nu) / 2.0
            ) / dx2dy2
            # x = -1, y = 0
            self.cj_1i0_coeff_ij = (-4.0 * D0 + 2.0 * Dx + Dxx) / dx4 + (
                -4.0 * D0 + 2.0 * Dx + nu * Dyy
            ) / dx2dy2
            # x = 0, y = -1
            self.cj0i_1_coeff_ij = (-4.0 * D0 + 2.0 * Dy + Dyy) / dy4 + (
                -4.0 * D0 + 2.0 * Dy + nu * Dxx
            ) / dx2dy2
            # x = 0, y = 1
            self.cj0i1_coeff_ij = (-4.0 * D0 - 2.0 * Dy + Dyy) / dy4 + (
                -4.0 * D0 - 2.0 * Dy + nu * Dxx
            ) / dx2dy2
            # x = 1, y = 0
            self.cj1i0_coeff_ij = (-4.0 * D0 - 2.0 * Dx + Dxx) / dx4 + (
                -4.0 * D0 - 2.0 * Dx + nu * Dyy
            ) / dx2dy2
            # x = 0, y = 0
            self.cj0i0_coeff_ij = (
                (6.0 * D0 - 2.0 * Dxx) / dx4
                + (6.0 * D0 - 2.0 * Dyy) / dy4
                + (8.0 * D0 - 2.0 * nu * Dxx - 2.0 * nu * Dyy) / dx2dy2
                + drho * g
            )

            ################################################################
            # CREATE COEFFICIENT ARRAYS: PLAIN, WITH NO B.C.'S YET APPLIED #
            ################################################################
            # x = -2, y = 0
            self.cj_2i0 = self.cj_2i0_coeff_ij.copy()
            # x = -1, y = -1
            self.cj_1i_1 = self.cj_1i_1_coeff_ij.copy()
            # x = -1, y = 0
            self.cj_1i0 = self.cj_1i0_coeff_ij.copy()
            # x = -1, y = 1
            self.cj_1i1 = self.cj_1i1_coeff_ij.copy()
            # x = 0, y = -2
            self.cj0i_2 = self.cj0i_2_coeff_ij.copy()
            # x = 0, y = -1
            self.cj0i_1 = self.cj0i_1_coeff_ij.copy()
            # x = 0, y = 0
            self.cj0i0 = self.cj0i0_coeff_ij.copy()
            # x = 0, y = 1
            self.cj0i1 = self.cj0i1_coeff_ij.copy()
            # x = 0, y = 2
            self.cj0i2 = self.cj0i2_coeff_ij.copy()
            # x = 1, y = -1
            self.cj1i_1 = self.cj1i_1_coeff_ij.copy()
            # x = 1, y = 0
            self.cj1i0 = self.cj1i0_coeff_ij.copy()
            # x = 1, y = 1
            self.cj1i1 = self.cj1i1_coeff_ij.copy()
            # x = 2, y = 0
            self.cj2i0 = self.cj2i0_coeff_ij.copy()

        # Provide rows and columns in the 2D input to later functions
        self.ncolsx = self.cj0i0.shape[1]
        self.nrowsy = self.cj0i0.shape[0]

    def _apply_bc_flexure(self):
        """Apply flexural boundary conditions to the 2D coefficient arrays."""
        # The next section of code is split over several functions for the 1D
        # case, but will be all in one function here, at least for now.

        # np.inf is used as the exclusion flag rather than 0 because subsequent
        # stencil additions cannot accidentally un-flag a slot: inf + any finite
        # value = inf. After np.roll along axis=1 (E–W), wrapped ghost-column
        # values land at valid matrix positions; because they remain inf
        # regardless of what is added to them, the isinf→0 pass below zeros
        # them all cleanly before matrix assembly. N–S ghost rows roll off the
        # array ends in C-order and are naturally excluded by spdiags, so no
        # special flag is needed for them.

        #######################################################################
        # DEFINE COEFFICIENTS TO W_j-2 -- W_j+2 WITH B.C.'S APPLIED (x: W, E) #
        #######################################################################

        # Infinitiy is used to flag places where coeff values should be 0,
        # and would otherwise cause boundary condition nan's to appear in the
        # cross-derivatives: infinity is changed into 0 later.

        if self.bc_west == "periodic":
            if self.bc_east == "periodic":
                # For each side, there will be two new diagonals (mostly zeros), and
                # two sets of diagonals that will replace values in current diagonals.
                # This is because of the pattern of fill in the periodic b.c.'s in the
                # x-direction.

                # First, create arrays for the new values.
                # One of the two values here, that from the y -/+ 1, x +/- 1 (E/W)
                # boundary condition, will be in the same location that will be
                # overwritten in the initiating grid by the next perioidic b.c. over
                self.cj_1i1_Periodic_right = np.zeros(self.qs.shape)
                self.cj_2i0_Periodic_right = np.zeros(self.qs.shape)
                j = 0
                self.cj_1i1_Periodic_right[:, j] = self.cj_1i1[:, j]
                self.cj_2i0_Periodic_right[:, j] = self.cj_2i0[:, j]
                j = 1
                self.cj_2i0_Periodic_right[:, j] = self.cj_2i0[:, j]

                # Then, replace existing values with what will be needed to make the
                # periodic boundary condition work.
                j = 0
                # ORDER IS IMPORTANT HERE! Don't change first before it changes other.
                # (We are shuffling down the line)
                self.cj_1i1[:, j] = self.cj_1i0[:, j]
                self.cj_1i0[:, j] = self.cj_1i_1[:, j]

                # And then change remaning off-grid values to np.inf (i.e. those that
                # were not altered to a real value
                # These will be the +/- 2's and the j_1i_1 and the j1i1
                # These are the farthest-out pentadiagonals that can't be reached by
                # the tridiagonals, and the tridiagonals that are farther away on the
                # y (big grid) axis that can't be reached by the single diagonals
                # that are farthest out
                # So 4 diagonals.
                # But ci1j1 is taken care of on -1 end before being rolled forwards
                # (i.e. clockwise, if we are reading from the top of the tread of a
                # tire)
                j = 0
                self.cj_2i0[:, j] += np.inf
                self.cj_1i_1[:, j] += np.inf
                j = 1
                self.cj_2i0[:, j] += np.inf

            else:
                pass  # bc_check() has already warned; proceed with no west periodic stencil
        elif self._bc_west_norm == "zero_displacement_zero_slope":
            # Boundary column (j=0): decouple so c0·w[j=0,i] = 0 → w = 0 exactly.
            j = 0
            self.cj_2i0[:, j] += np.inf
            self.cj_1i_1[:, j] += np.inf
            self.cj_1i0[:, j] += np.inf
            self.cj_1i1[:, j] += np.inf
            self.cj0i_2[:, j] += np.inf
            self.cj0i_1[:, j] += np.inf
            self.cj0i0[:, j] += 0
            self.cj0i1[:, j] += np.inf
            self.cj0i2[:, j] += np.inf
            self.cj1i_1[:, j] += np.inf
            self.cj1i0[:, j] += np.inf
            self.cj1i1[:, j] += np.inf
            self.cj2i0[:, j] += np.inf
            # First interior column (j=1): even reflection w[-1,i]=w[1,i] encodes zero slope.
            j = 1
            self.cj_2i0[:, j] += np.inf
            self.cj0i0[:, j] += self.cj_2i0_coeff_ij[:, j]
        elif self._bc_west_norm == "zero_moment_zero_shear":
            j = 0
            self.cj_2i0[:, j] += np.inf
            self.cj_1i_1[:, j] += np.inf
            self.cj_1i0[:, j] += np.inf
            self.cj_1i1[:, j] += np.inf
            self.cj0i_2[:, j] += 0
            self.cj0i_1[:, j] += 2 * self.cj_1i_1_coeff_ij[:, j]
            self.cj0i0[:, j] += (
                4 * self.cj_2i0_coeff_ij[:, j] + 2 * self.cj_1i0_coeff_ij[:, j]
            )
            self.cj0i1[:, j] += 2 * self.cj_1i1_coeff_ij[:, j]
            self.cj0i2[:, j] += 0
            self.cj1i_1[:, j] += -self.cj_1i_1_coeff_ij[:, j]
            self.cj1i0[:, j] += (
                -4 * self.cj_2i0_coeff_ij[:, j] - self.cj_1i0_coeff_ij[:, j]
            )
            self.cj1i1[:, j] += -self.cj_1i1_coeff_ij[:, j]
            self.cj2i0[:, j] += self.cj_2i0_coeff_ij[:, j]
            # First interior column (j=1): moment ghost w[-1,i]=2w[0,i]-w[1,i]
            j = 1
            self.cj_2i0[:, j] += np.inf
            self.cj_1i0[:, j] += 2 * self.cj_2i0_coeff_ij[:, j]
            self.cj0i0[:, j] -= self.cj_2i0_coeff_ij[:, j]
        elif self._bc_west_norm == "zero_slope_zero_shear":
            j = 0
            self.cj_2i0[:, j] += np.inf
            self.cj_1i_1[:, j] += np.inf
            self.cj_1i0[:, j] += np.inf
            self.cj_1i1[:, j] += np.inf
            self.cj0i_2[:, j] += 0
            self.cj0i_1[:, j] += 0
            self.cj0i0[:, j] += 0
            self.cj0i1[:, j] += 0
            self.cj0i2[:, j] += 0
            self.cj1i_1[:, j] += self.cj_1i_1_coeff_ij[:, j]
            self.cj1i0[:, j] += self.cj_1i0_coeff_ij[:, j]
            self.cj1i1[:, j] += self.cj_1i1_coeff_ij[:, j]
            self.cj2i0[:, j] += self.cj_2i0_coeff_ij[:, j]
            j = 1
            self.cj_2i0[:, j] += np.inf
            self.cj_1i_1[:, j] += 0
            self.cj_1i0[:, j] += 0
            self.cj_1i1[:, j] += 0
            self.cj0i_2[:, j] += 0
            self.cj0i_1[:, j] += 0
            self.cj0i0[:, j] += self.cj_2i0_coeff_ij[:, j]
            self.cj0i1[:, j] += 0
            self.cj0i2[:, j] += 0
            self.cj1i_1[:, j] += 0
            self.cj1i0[:, j] += 0
            self.cj1i1[:, j] += 0
            self.cj2i0[:, j] += 0
        elif self._bc_west_norm == "zero_displacement_zero_moment":
            # Pinned/simply-supported BC: enforce Dirichlet w=0 at the boundary
            # column (j=0) by zeroing all off-diagonal stencil terms there, leaving
            # only c0*w[i,0]=q[i,0].  The moment condition M_x=0 is encoded at
            # j=1 via the odd-reflection ghost: w[j=-1]=-w[j=+1] gives
            # c0[j=1] -= cj_2i0_coeff (the ghost contributes negatively to itself).
            j = 0
            self.cj_2i0[:, j] += np.inf
            self.cj_1i_1[:, j] += np.inf
            self.cj_1i0[:, j] += np.inf
            self.cj_1i1[:, j] += np.inf
            self.cj0i_2[:, j] += np.inf
            self.cj0i_1[:, j] += np.inf
            self.cj0i0[:, j] += 0
            self.cj0i1[:, j] += np.inf
            self.cj0i2[:, j] += np.inf
            self.cj1i_1[:, j] += np.inf
            self.cj1i0[:, j] += np.inf
            self.cj1i1[:, j] += np.inf
            self.cj2i0[:, j] += np.inf
            j = 1
            self.cj_2i0[:, j] += np.inf
            self.cj_1i_1[:, j] += 0
            self.cj_1i0[:, j] += 0
            self.cj_1i1[:, j] += 0
            self.cj0i_2[:, j] += 0
            self.cj0i_1[:, j] += 0
            self.cj0i0[:, j] -= self.cj_2i0_coeff_ij[:, j]
            self.cj0i1[:, j] += 0
            self.cj0i2[:, j] += 0
            self.cj1i_1[:, j] += 0
            self.cj1i0[:, j] += 0
            self.cj1i1[:, j] += 0
            self.cj2i0[:, j] += 0
        else:
            # Possibly redundant safeguard
            raise RuntimeError("Invalid boundary condition")

        if self.bc_east == "periodic":
            # See more extensive comments above (bc_west)

            if self.bc_west == "periodic":
                # New arrays -- new diagonals, but mostly empty. Just corners of blocks
                # (boxes) in block-diagonal matrix
                self.cj1i_1_Periodic_left = np.zeros(self.qs.shape)
                self.cj2i0_Periodic_left = np.zeros(self.qs.shape)
                j = -1
                self.cj1i_1_Periodic_left[:, j] = self.cj1i_1[:, j]
                self.cj2i0_Periodic_left[:, j] = self.cj2i0[:, j]
                j = -2
                self.cj2i0_Periodic_left[:, j] = self.cj2i0[:, j]

                # Then, replace existing values with what will be needed to make the
                # periodic boundary condition work.
                j = -1
                self.cj1i_1[:, j] = self.cj1i0[:, j]
                self.cj1i0[:, j] = self.cj1i1[:, j]

                # And then change remaning off-grid values to np.inf (i.e. those that
                # were not altered to a real value
                j = -1
                self.cj1i1[:, j] += np.inf
                self.cj2i0[:, j] += np.inf
                j = -2
                self.cj2i0[:, j] += np.inf

            else:
                pass  # bc_check() has already warned; proceed with no east periodic stencil

        elif self._bc_east_norm == "zero_displacement_zero_slope":
            # First interior column (j=-2): even reflection w[N,i]=w[N-2,i] encodes zero slope.
            j = -2
            self.cj0i0[:, j] += self.cj2i0_coeff_ij[:, j]
            self.cj2i0[:, j] += np.inf
            # Boundary column (j=-1): decouple so c0·w[j=-1,i] = 0 → w = 0 exactly.
            j = -1
            self.cj_2i0[:, j] += np.inf
            self.cj_1i_1[:, j] += np.inf
            self.cj_1i0[:, j] += np.inf
            self.cj_1i1[:, j] += np.inf
            self.cj0i_2[:, j] += np.inf
            self.cj0i_1[:, j] += np.inf
            self.cj0i0[:, j] += 0
            self.cj0i1[:, j] += np.inf
            self.cj0i2[:, j] += np.inf
            self.cj1i_1[:, j] += np.inf
            self.cj1i0[:, j] += np.inf
            self.cj1i1[:, j] += np.inf
            self.cj2i0[:, j] += np.inf
        elif self._bc_east_norm == "zero_moment_zero_shear":
            j = -1
            self.cj_2i0[:, j] += self.cj2i0_coeff_ij[:, j]
            self.cj_1i_1[:, j] += -self.cj1i_1_coeff_ij[:, j]
            self.cj_1i0[:, j] += (
                -4 * self.cj2i0_coeff_ij[:, j] - self.cj1i0_coeff_ij[:, j]
            )
            self.cj_1i1[:, j] += -self.cj1i1_coeff_ij[:, j]
            self.cj0i_2[:, j] += 0
            self.cj0i_1[:, j] += 2 * self.cj1i_1_coeff_ij[:, j]
            self.cj0i0[:, j] += (
                4 * self.cj2i0_coeff_ij[:, j] + 2 * self.cj1i0_coeff_ij[:, j]
            )
            self.cj0i1[:, j] += 2 * self.cj1i1_coeff_ij[:, j]
            self.cj0i2[:, j] += 0
            self.cj1i_1[:, j] += np.inf
            self.cj1i0[:, j] += np.inf
            self.cj1i1[:, j] += np.inf
            self.cj2i0[:, j] += np.inf
            # First interior column (j=N-2): moment ghost w[N,i]=2w[N-1,i]-w[N-2,i]
            j = -2
            self.cj0i0[:, j] -= self.cj2i0_coeff_ij[:, j]
            self.cj1i0[:, j] += 2 * self.cj2i0_coeff_ij[:, j]
            self.cj2i0[:, j] += np.inf
        elif self._bc_east_norm == "zero_slope_zero_shear":
            j = -1
            self.cj_2i0[:, j] += self.cj2i0_coeff_ij[:, j]
            self.cj_1i_1[:, j] += self.cj1i_1_coeff_ij[:, j]
            self.cj_1i0[:, j] += self.cj1i0_coeff_ij[:, j]
            self.cj_1i1[:, j] += self.cj1i1_coeff_ij[:, j]
            self.cj0i_2[:, j] += 0
            self.cj0i_1[:, j] += 0
            self.cj0i0[:, j] += 0
            self.cj0i1[:, j] += 0
            self.cj0i2[:, j] += 0
            self.cj1i_1[:, j] += np.inf
            self.cj1i0[:, j] += np.inf
            self.cj1i1[:, j] += np.inf
            self.cj2i0[:, j] += np.inf
            j = -2
            self.cj_2i0[:, j] += 0
            self.cj_1i_1[:, j] += 0
            self.cj_1i0[:, j] += 0
            self.cj_1i1[:, j] += 0
            self.cj0i_2[:, j] += 0
            self.cj0i_1[:, j] += 0
            self.cj0i0[:, j] += self.cj2i0_coeff_ij[:, j]
            self.cj0i1[:, j] += 0
            self.cj0i2[:, j] += 0
            self.cj1i_1[:, j] += 0
            self.cj1i0[:, j] += 0
            self.cj1i1[:, j] += 0
            self.cj2i0[:, j] += np.inf
        elif self._bc_east_norm == "zero_displacement_zero_moment":
            # Pinned/simply-supported BC: enforce Dirichlet w=0 at the boundary
            # column (j=-1) by zeroing all off-diagonal stencil terms there.
            # Moment condition M_x=0 is encoded at j=-2 via the odd-reflection
            # ghost w[j=Nx]=-w[j=Nx-2]: c0[j=-2] -= cj2i0_coeff.
            j = -1
            self.cj_2i0[:, j] += np.inf
            self.cj_1i_1[:, j] += np.inf
            self.cj_1i0[:, j] += np.inf
            self.cj_1i1[:, j] += np.inf
            self.cj0i_2[:, j] += np.inf
            self.cj0i_1[:, j] += np.inf
            self.cj0i0[:, j] += 0
            self.cj0i1[:, j] += np.inf
            self.cj0i2[:, j] += np.inf
            self.cj1i_1[:, j] += np.inf
            self.cj1i0[:, j] += np.inf
            self.cj1i1[:, j] += np.inf
            self.cj2i0[:, j] += np.inf
            j = -2
            self.cj_2i0[:, j] += 0
            self.cj_1i_1[:, j] += 0
            self.cj_1i0[:, j] += 0
            self.cj_1i1[:, j] += 0
            self.cj0i_2[:, j] += 0
            self.cj0i_1[:, j] += 0
            self.cj0i0[:, j] -= self.cj2i0_coeff_ij[:, j]
            self.cj0i1[:, j] += 0
            self.cj0i2[:, j] += 0
            self.cj1i_1[:, j] += 0
            self.cj1i0[:, j] += 0
            self.cj1i1[:, j] += 0
            self.cj2i0[:, j] += np.inf
        else:
            # Possibly redundant safeguard
            raise RuntimeError("Invalid boundary condition")

        #######################################################################
        # DEFINE COEFFICIENTS TO W_i-2 -- W_i+2 WITH B.C.'S APPLIED (y: N, S) #
        #######################################################################

        if self.bc_north == "periodic":
            if self.bc_south == "periodic":
                pass  # Will address the N-S (whole-matrix-involving) boundary condition
                # inclusion below, when constructing sparse matrix diagonals
            else:
                pass  # bc_check() has already warned; proceed with no north periodic stencil
        elif self._bc_north_norm == "zero_displacement_zero_slope":
            # Boundary row (i=0): decouple so c0·w[j,i=0] = 0 → w = 0 exactly.
            i = 0
            self.cj_2i0[i, :] += np.inf
            self.cj_1i_1[i, :] += np.inf
            self.cj_1i0[i, :] += np.inf
            self.cj_1i1[i, :] += np.inf
            self.cj0i_2[i, :] += np.inf
            self.cj0i_1[i, :] += np.inf
            self.cj0i0[i, :] += 0
            self.cj0i1[i, :] += np.inf
            self.cj0i2[i, :] += np.inf
            self.cj1i_1[i, :] += np.inf
            self.cj1i0[i, :] += np.inf
            self.cj1i1[i, :] += np.inf
            self.cj2i0[i, :] += np.inf
            # First interior row (i=1): even reflection w[j,-1]=w[j,1] encodes zero slope.
            i = 1
            self.cj0i0[i, :] += self.cj0i_2_coeff_ij[i, :]
        elif self._bc_north_norm == "zero_moment_zero_shear":
            i = 0
            self.cj_2i0[i, :] += 0
            self.cj_1i_1[i, :][self.cj_1i_1[i, :] != np.inf] = np.nan
            self.cj_1i0[i, :] += 2 * self.cj_1i_1_coeff_ij[i, :]
            self.cj_1i1[i, :] += -self.cj_1i_1_coeff_ij[i, :]
            self.cj0i_2[i, :] += 0  # np.nan
            self.cj0i_1[i, :] += 0  # np.nan
            self.cj0i0[i, :] += (
                4 * self.cj0i_2_coeff_ij[i, :] + 2 * self.cj0i_1_coeff_ij[i, :]
            )
            self.cj0i1[i, :] += (
                -4 * self.cj0i_2_coeff_ij[i, :] - self.cj0i_1_coeff_ij[i, :]
            )
            self.cj0i2[i, :] += self.cj0i_2_coeff_ij[i, :]
            self.cj1i_1[i, :][self.cj1i_1[i, :] != np.inf] += 0  # np.nan
            self.cj1i0[i, :] += 2 * self.cj1i_1_coeff_ij[i, :]
            self.cj1i1[i, :] += -self.cj1i_1_coeff_ij[i, :]
            self.cj2i0[i, :] += 0
            # First interior row (i=1): moment ghost w[j,-1]=2w[j,0]-w[j,1]
            i = 1
            self.cj0i_1[i, :] += 2 * self.cj0i_2_coeff_ij[i, :]
            self.cj0i0[i, :] -= self.cj0i_2_coeff_ij[i, :]
        elif self._bc_north_norm == "zero_slope_zero_shear":
            i = 0
            self.cj_2i0[i, :] += 0
            self.cj_1i_1[i, :][self.cj_1i_1[i, :] != np.inf] = np.nan
            self.cj_1i0[i, :] += 0
            self.cj_1i1[i, :] += self.cj_1i_1_coeff_ij[i, :]
            self.cj0i_2[i, :] += 0  # np.nan
            self.cj0i_1[i, :] += 0  # np.nan
            self.cj0i0[i, :] += 0
            self.cj0i1[i, :] += self.cj0i_1_coeff_ij[i, :]
            self.cj0i2[i, :] += self.cj0i_2_coeff_ij[i, :]
            self.cj1i_1[i, :][self.cj1i_1[i, :] != np.inf] += 0  # np.nan
            self.cj1i0[i, :] += 0
            self.cj1i1[i, :] += self.cj1i_1_coeff_ij[i, :]
            self.cj2i0[i, :] += 0
            i = 1
            self.cj_2i0[i, :] += 0
            self.cj_1i_1[i, :] += 0
            self.cj_1i0[i, :] += 0
            self.cj_1i1[i, :] += 0
            self.cj0i_2[i, :] += 0  # np.nan
            self.cj0i_1[i, :] += 0
            self.cj0i0[i, :] += self.cj0i_2_coeff_ij[i, :]
            self.cj0i1[i, :] += 0
            self.cj0i2[i, :] += 0
            self.cj1i_1[i, :] += 0
            self.cj1i0[i, :] += 0
            self.cj1i1[i, :] += 0
            self.cj2i0[i, :] += 0
        elif self._bc_north_norm == "zero_displacement_zero_moment":
            # Pinned/simply-supported BC: enforce Dirichlet w=0 at the boundary
            # row (i=0) by zeroing all off-diagonal stencil terms there.
            # Moment condition M_y=0 is encoded at i=1 via the odd-reflection
            # ghost w[i=-1]=-w[i=+1]: c0[i=1] -= cj0i_2_coeff.
            i = 0
            self.cj_2i0[i, :] += np.inf
            self.cj_1i_1[i, :] += np.inf
            self.cj_1i0[i, :] += np.inf
            self.cj_1i1[i, :] += np.inf
            self.cj0i_2[i, :] += np.inf
            self.cj0i_1[i, :] += np.inf
            self.cj0i0[i, :] += 0
            self.cj0i1[i, :] += np.inf
            self.cj0i2[i, :] += np.inf
            self.cj1i_1[i, :] += np.inf
            self.cj1i0[i, :] += np.inf
            self.cj1i1[i, :] += np.inf
            self.cj2i0[i, :] += np.inf
            i = 1
            self.cj_2i0[i, :] += 0
            self.cj_1i_1[i, :] += 0
            self.cj_1i0[i, :] += 0
            self.cj_1i1[i, :] += 0
            self.cj0i_2[i, :] += 0  # np.nan
            self.cj0i_1[i, :] += 0
            self.cj0i0[i, :] -= self.cj0i_2_coeff_ij[i, :]
            self.cj0i1[i, :] += 0
            self.cj0i2[i, :] += 0
            self.cj1i_1[i, :] += 0
            self.cj1i0[i, :] += 0
            self.cj1i1[i, :] += 0
            self.cj2i0[i, :] += 0
        else:
            # Possibly redundant safeguard
            raise RuntimeError("Invalid boundary condition")

        if self.bc_south == "periodic":
            if self.bc_north == "periodic":
                pass  # Will address the N-S (whole-matrix-involving) boundary condition
                # inclusion below, when constructing sparse matrix diagonals
            else:
                pass  # bc_check() has already warned; proceed with no south periodic stencil
        elif self._bc_south_norm == "zero_displacement_zero_slope":
            # First interior row (i=-2): even reflection w[j,N]=w[j,N-2] encodes zero slope.
            i = -2
            self.cj0i0[i, :] += self.cj0i2_coeff_ij[i, :]
            self.cj0i2[i, :] += np.inf
            # Boundary row (i=-1): decouple so c0·w[j,i=-1] = 0 → w = 0 exactly.
            i = -1
            self.cj_2i0[i, :] += np.inf
            self.cj_1i_1[i, :] += np.inf
            self.cj_1i0[i, :] += np.inf
            self.cj_1i1[i, :] += np.inf
            self.cj0i_2[i, :] += np.inf
            self.cj0i_1[i, :] += np.inf
            self.cj0i0[i, :] += 0
            self.cj0i1[i, :] += np.inf
            self.cj0i2[i, :] += np.inf
            self.cj1i_1[i, :] += np.inf
            self.cj1i0[i, :] += np.inf
            self.cj1i1[i, :] += np.inf
            self.cj2i0[i, :] += np.inf
        elif self._bc_south_norm == "zero_moment_zero_shear":
            # First interior row (i=N-2): moment ghost w[j,N]=2w[j,N-1]-w[j,N-2]
            i = -2
            self.cj0i0[i, :] -= self.cj0i2_coeff_ij[i, :]
            self.cj0i1[i, :] += 2 * self.cj0i2_coeff_ij[i, :]
            i = -1
            self.cj_2i0[i, :] += 0
            self.cj_1i_1[i, :] += -self.cj1i1_coeff_ij[i, :]
            self.cj_1i0[i, :] += 2 * self.cj1i1_coeff_ij[i, :]
            self.cj_1i1[i, :][self.cj_1i1[i, :] != np.inf] += 0  # np.nan
            self.cj0i_2[i, :] += self.cj0i2_coeff_ij[i, :]
            self.cj0i_1[i, :] += (
                -4 * self.cj0i2_coeff_ij[i, :] - self.cj0i1_coeff_ij[i, :]
            )
            self.cj0i0[i, :] += (
                4 * self.cj0i2_coeff_ij[i, :] + 2 * self.cj0i1_coeff_ij[i, :]
            )
            self.cj0i1[i, :] += 0  # np.nan
            self.cj0i2[i, :] += 0  # np.nan
            self.cj1i_1[i, :] += -self.cj_1i1_coeff_ij[i, :]
            self.cj1i0[i, :] += 2 * self.cj_1i1_coeff_ij[i, :]
            self.cj1i1[i, :][self.cj1i1[i, :] != np.inf] += 0  # np.nan
            self.cj2i0[i, :] += 0
        elif self._bc_south_norm == "zero_slope_zero_shear":
            i = -2
            self.cj_2i0[i, :] += 0
            self.cj_1i_1[i, :] += 0
            self.cj_1i0[i, :] += 0
            self.cj_1i1[i, :] += 0
            self.cj0i_2[i, :] += 0
            self.cj0i_1[i, :] += 0
            self.cj0i0[i, :] += self.cj0i2_coeff_ij[i, :]
            self.cj0i1[i, :] += 0
            self.cj0i2[i, :] += 0  # np.nan
            self.cj1i_1[i, :] += 0
            self.cj1i0[i, :] += 0
            self.cj1i1[i, :] += 0
            self.cj2i0[i, :] += 0
            i = -1
            self.cj_2i0[i, :] += 0
            self.cj_1i_1[i, :] += self.cj_1i1_coeff_ij[i, :]
            self.cj_1i0[i, :] += 0
            self.cj_1i1[i, :][self.cj_1i1[i, :] != np.inf] += 0  # np.nan
            self.cj0i_2[i, :] += self.cj0i2_coeff_ij[i, :]
            self.cj0i_1[i, :] += self.cj0i1_coeff_ij[i, :]
            self.cj0i0[i, :] += 0
            self.cj0i1[i, :] += 0  # np.nan
            self.cj0i2[i, :] += 0  # np.nan
            self.cj1i_1[i, :] += self.cj1i1_coeff_ij[i, :]
            self.cj1i0[i, :] += 0
            self.cj1i1[i, :][self.cj1i1[i, :] != np.inf] += 0  # np.nan
            self.cj2i0[i, :] += 0
        elif self._bc_south_norm == "zero_displacement_zero_moment":
            # Pinned/simply-supported BC: enforce Dirichlet w=0 at the boundary
            # row (i=-1) by zeroing all off-diagonal stencil terms there.
            # Moment condition M_y=0 is encoded at i=-2 via the odd-reflection
            # ghost w[i=Ny]=-w[i=Ny-2]: c0[i=-2] -= cj0i2_coeff.
            i = -2
            self.cj_2i0[i, :] += 0
            self.cj_1i_1[i, :] += 0
            self.cj_1i0[i, :] += 0
            self.cj_1i1[i, :] += 0
            self.cj0i_2[i, :] += 0
            self.cj0i_1[i, :] += 0
            self.cj0i0[i, :] -= self.cj0i2_coeff_ij[i, :]
            self.cj0i1[i, :] += 0
            self.cj0i2[i, :] += 0  # np.nan
            self.cj1i_1[i, :] += 0
            self.cj1i0[i, :] += 0
            self.cj1i1[i, :] += 0
            self.cj2i0[i, :] += 0
            i = -1
            self.cj_2i0[i, :] += np.inf
            self.cj_1i_1[i, :] += np.inf
            self.cj_1i0[i, :] += np.inf
            self.cj_1i1[i, :] += np.inf
            self.cj0i_2[i, :] += np.inf
            self.cj0i_1[i, :] += np.inf
            self.cj0i0[i, :] += 0
            self.cj0i1[i, :] += np.inf
            self.cj0i2[i, :] += np.inf
            self.cj1i_1[i, :] += np.inf
            self.cj1i0[i, :] += np.inf
            self.cj1i1[i, :] += np.inf
            self.cj2i0[i, :] += np.inf
        else:
            # Possibly redundant safeguard
            raise RuntimeError("Invalid boundary condition")

        #####################################################
        # CORNERS: INTERFERENCE BETWEEN BOUNDARY CONDITIONS #
        #####################################################

        # In 2D, have to consider diagonals and interference (additive) among
        # boundary conditions

        ############################
        # DIRICHLET -- DO NOTHING. #
        ############################

        # Do nothing.
        # What about combinations?
        # This will mean that dirichlet boundary conditions will implicitly
        # control the corners, so, for examplel, they would be locked all of the
        # way to the edge of the domain instead of becoming free to deflect at the
        # ends.
        # Indeed it is much easier to envision this case than one in which
        # the stationary clamp is released.

        #################
        # 0MOMENT0SHEAR #
        #################
        if self._bc_north_norm == "zero_moment_zero_shear" and self._bc_west_norm == "zero_moment_zero_shear":
            self.cj0i0[0, 0] += 2 * self.cj_1i_1_coeff_ij[0, 0]
            self.cj1i1[0, 0] -= self.cj_1i_1_coeff_ij[0, 0]
        if self._bc_north_norm == "zero_moment_zero_shear" and self._bc_east_norm == "zero_moment_zero_shear":
            self.cj0i0[0, -1] += 2 * self.cj_1i_1_coeff_ij[0, -1]
            self.cj_1i1[0, -1] -= self.cj1i_1_coeff_ij[0, -1]
        if self._bc_south_norm == "zero_moment_zero_shear" and self._bc_west_norm == "zero_moment_zero_shear":
            self.cj0i0[-1, 0] += 2 * self.cj_1i_1_coeff_ij[-1, 0]
            self.cj1i_1[-1, 0] -= self.cj_1i1_coeff_ij[-1, 0]
        if self._bc_south_norm == "zero_moment_zero_shear" and self._bc_east_norm == "zero_moment_zero_shear":
            self.cj0i0[-1, -1] += 2 * self.cj_1i_1_coeff_ij[-1, -1]
            self.cj_1i_1[-1, -1] -= self.cj1i1_coeff_ij[-1, -1]

        ############
        # PERIODIC #
        ############

        # I think that nothing will be needed here.
        # periodic should just take care of all repeating in all directions by
        # its very nature. (I.e. it is embedded in the sparse array structure
        # of diagonals)

        ################
        # COMBINATIONS #
        ################

        ##########
        # MIRROR #
        ##########
        if self._bc_north_norm == "zero_slope_zero_shear" and self._bc_west_norm == "zero_slope_zero_shear":
            self.cj1i1[0, 0] += self.cj_1i_1_coeff_ij[0, 0]
        if self._bc_north_norm == "zero_slope_zero_shear" and self._bc_east_norm == "zero_slope_zero_shear":
            self.cj_1i1[0, -1] += self.cj1i_1_coeff_ij[0, -1]
        if self._bc_south_norm == "zero_slope_zero_shear" and self._bc_west_norm == "zero_slope_zero_shear":
            self.cj1i_1[-1, 0] += self.cj_1i1_coeff_ij[-1, 0]
        if self._bc_south_norm == "zero_slope_zero_shear" and self._bc_east_norm == "zero_slope_zero_shear":
            self.cj_1i_1[-1, -1] += self.cj1i1_coeff_ij[-1, -1]

        ################################
        # 0MOMENT0SHEAR - AND - MIRROR #
        ################################
        # How do multiple types of b.c.'s interfere
        # zero_moment_zero_shear must determine corner conditions in order to be mirrored
        # by the "zero_slope_zero_shear" b.c.
        if (self._bc_north_norm == "zero_slope_zero_shear" and self._bc_west_norm == "zero_moment_zero_shear") or (
            self._bc_west_norm == "zero_slope_zero_shear" and self._bc_north_norm == "zero_moment_zero_shear"
        ):
            self.cj0i0[0, 0] += 2 * self.cj_1i_1_coeff_ij[0, 0]
            self.cj1i1[0, 0] -= self.cj_1i_1_coeff_ij[0, 0]
        if (self._bc_north_norm == "zero_slope_zero_shear" and self._bc_east_norm == "zero_moment_zero_shear") or (
            self._bc_east_norm == "zero_slope_zero_shear" and self._bc_north_norm == "zero_moment_zero_shear"
        ):
            self.cj0i0[0, -1] += 2 * self.cj_1i_1_coeff_ij[0, -1]
            self.cj_1i1[0, -1] -= self.cj_1i_1_coeff_ij[0, -1]
        if (self._bc_south_norm == "zero_slope_zero_shear" and self._bc_west_norm == "zero_moment_zero_shear") or (
            self._bc_west_norm == "zero_slope_zero_shear" and self._bc_south_norm == "zero_moment_zero_shear"
        ):
            self.cj0i0[-1, 0] += 2 * self.cj_1i_1_coeff_ij[-1, 0]
            self.cj1i_1[-1, 0] -= self.cj_1i1_coeff_ij[-1, 0]
        if (self._bc_south_norm == "zero_slope_zero_shear" and self._bc_east_norm == "zero_moment_zero_shear") or (
            self._bc_east_norm == "zero_slope_zero_shear" and self._bc_south_norm == "zero_moment_zero_shear"
        ):
            self.cj0i0[-1, -1] += 2 * self.cj_1i_1_coeff_ij[-1, -1]
            self.cj_1i_1[-1, -1] -= self.cj1i1_coeff_ij[-1, -1]

        ######################
        # 0DISPLACEMENT0MOMENT #
        ######################
        # odd×odd = even (positive), same sign as mirror×mirror.
        # odd×even (mirror) = negative.
        # odd×zero_moment_zero_shear = negative of the mirror×zero_moment_zero_shear correction.
        _D0M = "zero_displacement_zero_moment"

        # NW corner
        if self._bc_north_norm == _D0M and self._bc_west_norm == _D0M:
            self.cj1i1[0, 0] += self.cj_1i_1_coeff_ij[0, 0]
        elif (self._bc_north_norm == _D0M and self._bc_west_norm == "zero_slope_zero_shear") or (
            self._bc_west_norm == _D0M and self._bc_north_norm == "zero_slope_zero_shear"
        ):
            self.cj1i1[0, 0] -= self.cj_1i_1_coeff_ij[0, 0]
        elif (self._bc_north_norm == _D0M and self._bc_west_norm == "zero_moment_zero_shear") or (
            self._bc_west_norm == _D0M and self._bc_north_norm == "zero_moment_zero_shear"
        ):
            self.cj0i0[0, 0] -= 2 * self.cj_1i_1_coeff_ij[0, 0]
            self.cj1i1[0, 0] += self.cj_1i_1_coeff_ij[0, 0]

        # NE corner
        if self._bc_north_norm == _D0M and self._bc_east_norm == _D0M:
            self.cj_1i1[0, -1] += self.cj1i_1_coeff_ij[0, -1]
        elif (self._bc_north_norm == _D0M and self._bc_east_norm == "zero_slope_zero_shear") or (
            self._bc_east_norm == _D0M and self._bc_north_norm == "zero_slope_zero_shear"
        ):
            self.cj_1i1[0, -1] -= self.cj1i_1_coeff_ij[0, -1]
        elif (self._bc_north_norm == _D0M and self._bc_east_norm == "zero_moment_zero_shear") or (
            self._bc_east_norm == _D0M and self._bc_north_norm == "zero_moment_zero_shear"
        ):
            self.cj0i0[0, -1] -= 2 * self.cj_1i_1_coeff_ij[0, -1]
            self.cj_1i1[0, -1] += self.cj1i_1_coeff_ij[0, -1]

        # SW corner
        if self._bc_south_norm == _D0M and self._bc_west_norm == _D0M:
            self.cj1i_1[-1, 0] += self.cj_1i1_coeff_ij[-1, 0]
        elif (self._bc_south_norm == _D0M and self._bc_west_norm == "zero_slope_zero_shear") or (
            self._bc_west_norm == _D0M and self._bc_south_norm == "zero_slope_zero_shear"
        ):
            self.cj1i_1[-1, 0] -= self.cj_1i1_coeff_ij[-1, 0]
        elif (self._bc_south_norm == _D0M and self._bc_west_norm == "zero_moment_zero_shear") or (
            self._bc_west_norm == _D0M and self._bc_south_norm == "zero_moment_zero_shear"
        ):
            self.cj0i0[-1, 0] -= 2 * self.cj_1i_1_coeff_ij[-1, 0]
            self.cj1i_1[-1, 0] += self.cj_1i1_coeff_ij[-1, 0]

        # SE corner
        if self._bc_south_norm == _D0M and self._bc_east_norm == _D0M:
            self.cj_1i_1[-1, -1] += self.cj1i1_coeff_ij[-1, -1]
        elif (self._bc_south_norm == _D0M and self._bc_east_norm == "zero_slope_zero_shear") or (
            self._bc_east_norm == _D0M and self._bc_south_norm == "zero_slope_zero_shear"
        ):
            self.cj_1i_1[-1, -1] -= self.cj1i1_coeff_ij[-1, -1]
        elif (self._bc_south_norm == _D0M and self._bc_east_norm == "zero_moment_zero_shear") or (
            self._bc_east_norm == _D0M and self._bc_south_norm == "zero_moment_zero_shear"
        ):
            self.cj0i0[-1, -1] -= 2 * self.cj_1i_1_coeff_ij[-1, -1]
            self.cj_1i_1[-1, -1] += self.cj1i1_coeff_ij[-1, -1]

        ##############################
        # PERIODIC B.C.'S AND OTHERS #
        ##############################

        # The periodic boundary natively continues the other boundary conditions
        # Nothing to be done here.

    def _apply_bc_rhs_inhomogeneous_2d(self):
        """Compute the RHS correction for prescribed (nonzero) 2D BC values.

        Called after _apply_bc_flexure() and before build_diagonals().  For
        dict-style BCs the ghost nodes carry prescribed moment, shear,
        displacement, or slope values.  The constant parts of those ghost
        expressions move to the RHS as a correction array added to -qs
        before the sparse solve — the same approach used in 1D.

        Cross-derivative stencil terms (cj_1i_1, cj_1i1 and analogues) mean
        that a ghost node at (j−1, k) affects the equations at rows k−1, k,
        and k+1.  The k±1 contributions use np.roll to shift the per-row Δ₁
        vector; for uniform D the roll value equals the interior value so the
        approximation is exact.
        """
        if (self._bc_west_values is None and self._bc_east_values is None
                and self._bc_north_values is None and self._bc_south_values is None):
            self._bc_rhs_correction = None
            return

        ny, nx = self.qs.shape
        correction = np.zeros((ny, nx))
        dx = self.dx
        dy = self.dy

        # D at each boundary edge.
        # Scalar te: self.D is (ny, nx), unpadded.
        # Array te: self.D is (ny+2, nx+2), padded with one nan-row/col each side.
        if np.isscalar(self.T_e):
            D_west  = self.D[:, 0]
            D_east  = self.D[:, -1]
            D_north = self.D[0, :]
            D_south = self.D[-1, :]
        else:
            D_west  = self.D[1:-1, 1]
            D_east  = self.D[1:-1, -2]
            D_north = self.D[1, 1:-1]
            D_south = self.D[-2, 1:-1]

        # --- WEST (j=0) ---
        if self._bc_west_values is not None:
            bv = self._bc_west_values
            if "displacement" in bv:
                w0     = bv["displacement"]
                theta0 = bv["slope"]
                # Decoupled boundary column: enforce w[:, 0] = w0.
                correction[:, 0] += self.cj0i0[:, 0] * w0
                # First interior column (j=1): slope ghost w[-1,k]=w[1,k]-2dx*theta0.
                correction[:, 1] += self.cj_2i0_coeff_ij[:, 1] * 2.0 * dx * theta0
            else:  # moment / shear
                M0 = bv["moment"]
                V0 = bv["shear"]
                Delta1 = M0 * dx**2 / D_west  # constant part of w[-1, k], shape (ny,)
                # Boundary column (j=0): from w[-2,k] and w[-1,k] ghosts.
                correction[:, 0] += (
                    self.cj_2i0_coeff_ij[:, 0] * (2.0 * V0 * dx**3 - 2.0 * M0 * dx**2) / D_west
                    - self.cj_1i0_coeff_ij[:, 0] * Delta1
                    - self.cj_1i_1_coeff_ij[:, 0] * np.roll(Delta1, 1)   # w[-1, k-1] cross term
                    - self.cj_1i1_coeff_ij[:, 0] * np.roll(Delta1, -1)   # w[-1, k+1] cross term
                )
                # First interior column (j=1): cj_2i0 there sees w[-1, k].
                correction[:, 1] -= self.cj_2i0_coeff_ij[:, 1] * Delta1

        # --- EAST (j=nx-1) ---
        if self._bc_east_values is not None:
            bv = self._bc_east_values
            if "displacement" in bv:
                w_e     = bv["displacement"]
                theta_e = bv["slope"]
                correction[:, -1] += self.cj0i0[:, -1] * w_e
                correction[:, -2] -= self.cj2i0_coeff_ij[:, -2] * 2.0 * dx * theta_e
            else:  # moment / shear
                M_e = bv["moment"]
                V_e = bv["shear"]
                Delta1_e = M_e * dx**2 / D_east
                correction[:, -1] += (
                    self.cj2i0_coeff_ij[:, -1] * (-(2.0 * V_e * dx**3 + 2.0 * M_e * dx**2)) / D_east
                    - self.cj1i0_coeff_ij[:, -1] * Delta1_e
                    - self.cj1i_1_coeff_ij[:, -1] * np.roll(Delta1_e, 1)
                    - self.cj1i1_coeff_ij[:, -1] * np.roll(Delta1_e, -1)
                )
                correction[:, -2] -= self.cj2i0_coeff_ij[:, -2] * Delta1_e

        # --- NORTH (i=0) ---
        if self._bc_north_values is not None:
            bv = self._bc_north_values
            if "displacement" in bv:
                w_n     = bv["displacement"]
                theta_n = bv["slope"]
                correction[0, :] += self.cj0i0[0, :] * w_n
                correction[1, :] += self.cj0i_2_coeff_ij[1, :] * 2.0 * dy * theta_n
            else:  # moment / shear
                M_n = bv["moment"]
                V_n = bv["shear"]
                Delta1_n = M_n * dy**2 / D_north
                correction[0, :] += (
                    self.cj0i_2_coeff_ij[0, :] * (2.0 * V_n * dy**3 - 2.0 * M_n * dy**2) / D_north
                    - self.cj0i_1_coeff_ij[0, :] * Delta1_n
                    - self.cj_1i_1_coeff_ij[0, :] * np.roll(Delta1_n, 1)   # w[j-1, -1] cross term
                    - self.cj1i_1_coeff_ij[0, :] * np.roll(Delta1_n, -1)   # w[j+1, -1] cross term
                )
                correction[1, :] -= self.cj0i_2_coeff_ij[1, :] * Delta1_n

        # --- SOUTH (i=ny-1) ---
        if self._bc_south_values is not None:
            bv = self._bc_south_values
            if "displacement" in bv:
                w_s     = bv["displacement"]
                theta_s = bv["slope"]
                correction[-1, :] += self.cj0i0[-1, :] * w_s
                correction[-2, :] -= self.cj0i2_coeff_ij[-2, :] * 2.0 * dy * theta_s
            else:  # moment / shear
                M_s = bv["moment"]
                V_s = bv["shear"]
                Delta1_s = M_s * dy**2 / D_south
                correction[-1, :] += (
                    self.cj0i2_coeff_ij[-1, :] * (-(2.0 * V_s * dy**3 + 2.0 * M_s * dy**2)) / D_south
                    - self.cj0i1_coeff_ij[-1, :] * Delta1_s
                    - self.cj_1i1_coeff_ij[-1, :] * np.roll(Delta1_s, 1)
                    - self.cj1i1_coeff_ij[-1, :] * np.roll(Delta1_s, -1)
                )
                correction[-2, :] -= self.cj0i2_coeff_ij[-2, :] * Delta1_s

        # --- Corner deduplication ---
        # A corner node (e.g. [0,0]) lies on two boundary edges.  When both
        # are displacement Dirichlet, the west loop and north loop both add
        # cj0i0[0,0]*W0 to correction[0,0], yielding 2*W0 instead of W0.
        # Subtract the north/south contribution at each shared corner.
        d_w = self._bc_west_values  is not None and "displacement" in self._bc_west_values
        d_e = self._bc_east_values  is not None and "displacement" in self._bc_east_values
        d_n = self._bc_north_values is not None and "displacement" in self._bc_north_values
        d_s = self._bc_south_values is not None and "displacement" in self._bc_south_values
        if d_n:
            _wn = np.asarray(self._bc_north_values["displacement"])
            if d_w:
                correction[0, 0]   -= self.cj0i0[0, 0]   * float(_wn.flat[0])
            if d_e:
                correction[0, -1]  -= self.cj0i0[0, -1]  * float(_wn.flat[-1])
        if d_s:
            _ws = np.asarray(self._bc_south_values["displacement"])
            if d_w:
                correction[-1, 0]  -= self.cj0i0[-1, 0]  * float(_ws.flat[0])
            if d_e:
                correction[-1, -1] -= self.cj0i0[-1, -1] * float(_ws.flat[-1])

        self._bc_rhs_correction = correction

    def build_diagonals(self):
        """
        Assemble the sparse coefficient matrix from the 2-D stencil arrays.

        Takes the per-cell coefficient arrays populated by
        ``get_coeff_values`` (``cj0i0``, ``cj1i0``, ``cj_1i0``, etc., named
        by column offset *j* and row offset *i*) and shifts each array by the
        appropriate number of rows/columns using :func:`numpy.roll` so that
        off-diagonal contributions land on the correct diagonals when the 2-D
        grid is flattened in row-major (C) order.  Infinite values arising
        from boundary ghost cells are zeroed before assembly.

        The result is stored in ``self.coeff_matrix`` as a
        :class:`scipy.sparse.dia_matrix`, ready for the direct solver called
        by :meth:`F2D._solve_fd`.
        """
        ##########################################################
        # INCORPORATE BOUNDARY CONDITIONS INTO COEFFICIENT ARRAY #
        ##########################################################

        # Roll to keep the proper coefficients at the proper places in the
        # arrays: Python will naturally just do vertical shifts instead of
        # diagonal shifts, so this takes into account the horizontal compoent
        # to ensure that boundary values are at the right place.

        # Roll x
        # ASYMMETRIC RESPONSE HERE -- THIS GETS TOWARDS SOURCE OF PROBLEM!
        self.cj_2i0 = np.roll(self.cj_2i0, -2, 1)
        self.cj_1i0 = np.roll(self.cj_1i0, -1, 1)
        self.cj1i0 = np.roll(self.cj1i0, 1, 1)
        self.cj2i0 = np.roll(self.cj2i0, 2, 1)
        # Roll y
        self.cj0i_2 = np.roll(self.cj0i_2, -2, 0)
        self.cj0i_1 = np.roll(self.cj0i_1, -1, 0)
        self.cj0i1 = np.roll(self.cj0i1, 1, 0)
        self.cj0i2 = np.roll(self.cj0i2, 2, 0)
        # Roll x and y
        self.cj_1i_1 = np.roll(self.cj_1i_1, -1, 1)
        self.cj_1i_1 = np.roll(self.cj_1i_1, -1, 0)
        self.cj_1i1 = np.roll(self.cj_1i1, -1, 1)
        self.cj_1i1 = np.roll(self.cj_1i1, 1, 0)
        self.cj1i_1 = np.roll(self.cj1i_1, 1, 1)
        self.cj1i_1 = np.roll(self.cj1i_1, -1, 0)
        self.cj1i1 = np.roll(self.cj1i1, 1, 1)
        self.cj1i1 = np.roll(self.cj1i1, 1, 0)

        coeff_array_list = [
            self.cj_2i0,
            self.cj_1i0,
            self.cj1i0,
            self.cj2i0,
            self.cj0i_2,
            self.cj0i_1,
            self.cj0i1,
            self.cj0i2,
            self.cj_1i_1,
            self.cj_1i1,
            self.cj1i_1,
            self.cj1i1,
            self.cj0i0,
        ]
        for array in coeff_array_list:
            # np.inf flags excluded stencil slots. Because inf + finite = inf,
            # flagged slots stay flagged regardless of subsequent stencil
            # additions; setting them to 0 here enforces those exclusions
            # before sparse-matrix assembly.
            array[np.isinf(array)] = 0
        # array[np.isnan(array)] = 0 # had been used for testing

        # Reshape to put in solver
        vec_cj_2i0 = np.reshape(self.cj_2i0, -1, order="C")
        vec_cj_1i_1 = np.reshape(self.cj_1i_1, -1, order="C")
        vec_cj_1i0 = np.reshape(self.cj_1i0, -1, order="C")
        vec_cj_1i1 = np.reshape(self.cj_1i1, -1, order="C")
        vec_cj0i_2 = np.reshape(self.cj0i_2, -1, order="C")
        vec_cj0i_1 = np.reshape(self.cj0i_1, -1, order="C")
        vec_cj0i0 = np.reshape(self.cj0i0, -1, order="C")
        vec_cj0i1 = np.reshape(self.cj0i1, -1, order="C")
        vec_cj0i2 = np.reshape(self.cj0i2, -1, order="C")
        vec_cj1i_1 = np.reshape(self.cj1i_1, -1, order="C")
        vec_cj1i0 = np.reshape(self.cj1i0, -1, order="C")
        vec_cj1i1 = np.reshape(self.cj1i1, -1, order="C")
        vec_cj2i0 = np.reshape(self.cj2i0, -1, order="C")

        # Changed this 6 Nov. 2014 in betahaus Berlin to be x-based
        Up2 = vec_cj0i2
        Up1 = np.vstack((vec_cj_1i1, vec_cj0i1, vec_cj1i1))
        Mid = np.vstack((vec_cj_2i0, vec_cj_1i0, vec_cj0i0, vec_cj1i0, vec_cj2i0))
        Dn1 = np.vstack((vec_cj_1i_1, vec_cj0i_1, vec_cj1i_1))
        Dn2 = vec_cj0i_2

        # Number of rows and columns for array size and offsets
        self.ny = self.nrowsy
        self.nx = self.ncolsx

        if (
            self.bc_north == "periodic"
            and self.bc_south == "periodic"
            and self.bc_west == "periodic"
            and self.bc_east == "periodic"
        ):
            # Additional vector creation
            # West
            # Roll
            self.cj_2i0_Periodic_right = np.roll(self.cj_2i0_Periodic_right, -2, 1)
            self.cj_1i1_Periodic_right = np.roll(self.cj_1i1_Periodic_right, -1, 1)
            self.cj_1i1_Periodic_right = np.roll(self.cj_1i1_Periodic_right, 1, 0)
            # Reshape
            vec_cj_2i0_Periodic_right = np.reshape(
                self.cj_2i0_Periodic_right, -1, order="C"
            )
            vec_cj_1i1_Periodic_right = np.reshape(
                self.cj_1i1_Periodic_right, -1, order="C"
            )
            # East
            # Roll
            self.cj1i_1_Periodic_left = np.roll(self.cj1i_1_Periodic_left, 1, 1)
            self.cj1i_1_Periodic_left = np.roll(self.cj1i_1_Periodic_left, -1, 0)
            self.cj2i0_Periodic_left = np.roll(self.cj2i0_Periodic_left, 2, 1)
            # Reshape
            vec_cj1i_1_Periodic_left = np.reshape(
                self.cj1i_1_Periodic_left, -1, order="C"
            )
            vec_cj2i0_Periodic_left = np.reshape(
                self.cj2i0_Periodic_left, -1, order="C"
            )

            # Build diagonals with additional entries
            # I think the fact that everything is rolled will make this work all right
            # without any additional rolling.
            # Checked -- and indeed, what would be in my mind the last value for
            # Mid[3] is the first value in its array. Hooray, patterns!
            self.diags = np.vstack(
                (
                    vec_cj1i_1_Periodic_left,
                    Up1,
                    vec_cj_1i1_Periodic_right,
                    Up2,
                    Dn2,
                    vec_cj1i_1_Periodic_left,
                    Dn1,
                    vec_cj2i0_Periodic_left,
                    Mid,
                    vec_cj_2i0_Periodic_right,
                    Up1,
                    vec_cj_1i1_Periodic_right,
                    Up2,
                    Dn2,
                    vec_cj1i_1_Periodic_left,
                    Dn1,
                    vec_cj_1i1_Periodic_right,
                )
            )
            # Getting too complicated to have everything together
            self.offsets = [
                # New: LL corner of LL box
                -self.ny * self.nx + 1,
                # periodic b.c. tridiag
                self.nx - self.ny * self.nx - 1,
                self.nx - self.ny * self.nx,
                self.nx - self.ny * self.nx + 1,
                # New: UR corner of LL box
                2 * self.nx - self.ny * self.nx - 1,
                # periodic b.c. single diag
                2 * self.nx - self.ny * self.nx,
                -2 * self.nx,
                # New:
                -2 * self.nx + 1,
                # Right term here (-self.nx+1) modified:
                -self.nx - 1,
                -self.nx,
                -self.nx + 1,
                # New:
                -self.nx + 2,
                # -1 and 1 terms here modified:
                -2,
                -1,
                0,
                1,
                2,
                # New:
                self.nx - 2,
                # Left term here (self.nx-1) modified:
                self.nx - 1,
                self.nx,
                self.nx + 1,
                # New:
                2 * self.nx - 1,
                2 * self.nx,
                # periodic b.c. single diag
                self.ny * self.nx - 2 * self.nx,
                # New: LL corner of UR box
                self.ny * self.nx - 2 * self.nx + 1,
                # periodic b.c. tridiag
                self.ny * self.nx - self.nx - 1,
                self.ny * self.nx - self.nx,
                self.ny * self.nx - self.nx + 1,
                # New: UR corner of UR box
                self.ny * self.nx - 1,
            ]

            # create banded sparse matrix
            self.coeff_matrix = scipy.sparse.spdiags(
                self.diags,
                self.offsets,
                self.ny * self.nx,
                self.ny * self.nx,
                format="csr",
            )

            # The double-wrap corner diagonals reuse cj_1i1_Periodic_right and
            # cj1i_1_Periodic_left, which hold the (j-1,i+1) and (j+1,i-1)
            # coefficients respectively.  The two matrix entries that connect
            # cell (0,0)↔(ny-1,nx-1) wrap in BOTH the x and y directions
            # simultaneously, so they need the (j-1,i-1) and (j+1,i+1)
            # coefficients instead.  For sigma_xy=0 these are equal, so the
            # error is silent in the zero-stress case.
            N = self.ny * self.nx
            A_lil = self.coeff_matrix.tolil()
            A_lil[0, N - 1] = self.cj_1i_1_coeff_ij[0, 0]
            A_lil[N - 1, 0] = self.cj1i1_coeff_ij[-1, -1]
            self.coeff_matrix = A_lil.tocsr()

        elif self._bc_west_norm == "periodic" and self._bc_east_norm == "periodic":
            # Additional vector creation
            # West
            # Roll
            self.cj_2i0_Periodic_right = np.roll(self.cj_2i0_Periodic_right, -2, 1)
            self.cj_1i1_Periodic_right = np.roll(self.cj_1i1_Periodic_right, -1, 1)
            self.cj_1i1_Periodic_right = np.roll(self.cj_1i1_Periodic_right, 1, 0)
            # Reshape
            vec_cj_2i0_Periodic_right = np.reshape(
                self.cj_2i0_Periodic_right, -1, order="C"
            )
            vec_cj_1i1_Periodic_right = np.reshape(
                self.cj_1i1_Periodic_right, -1, order="C"
            )
            # East
            # Roll
            self.cj1i_1_Periodic_left = np.roll(self.cj1i_1_Periodic_left, 1, 1)
            self.cj1i_1_Periodic_left = np.roll(self.cj1i_1_Periodic_left, -1, 0)
            self.cj2i0_Periodic_left = np.roll(self.cj2i0_Periodic_left, 2, 1)
            # Reshape
            vec_cj1i_1_Periodic_left = np.reshape(
                self.cj1i_1_Periodic_left, -1, order="C"
            )
            vec_cj2i0_Periodic_left = np.reshape(
                self.cj2i0_Periodic_left, -1, order="C"
            )

            # Build diagonals with additional entries
            self.diags = np.vstack(
                (
                    Dn2,
                    vec_cj1i_1_Periodic_left,
                    Dn1,
                    vec_cj2i0_Periodic_left,
                    Mid,
                    vec_cj_2i0_Periodic_right,
                    Up1,
                    vec_cj_1i1_Periodic_right,
                    Up2,
                )
            )
            # Getting too complicated to have everything together
            self.offsets = [
                -2 * self.nx,
                # New:
                -2 * self.nx + 1,
                # Right term here (-self.nx+1) modified:
                -self.nx - 1,
                -self.nx,
                -self.nx + 1,
                # New:
                -self.nx + 2,
                # -1 and 1 terms here modified:
                -2,
                -1,
                0,
                1,
                2,
                # New:
                self.nx - 2,
                # Left term here (self.nx-1) modified:
                self.nx - 1,
                self.nx,
                self.nx + 1,
                # New:
                2 * self.nx - 1,
                2 * self.nx,
            ]

            # create banded sparse matrix
            self.coeff_matrix = scipy.sparse.spdiags(
                self.diags,
                self.offsets,
                self.ny * self.nx,
                self.ny * self.nx,
                format="csr",
            )

        elif self._bc_north_norm == "periodic" and self._bc_south_norm == "periodic":
            # periodic.
            # If these are periodic, we need to wrap around the ends of the
            # large-scale diagonal structure
            self.diags = np.vstack((Up1, Up2, Dn2, Dn1, Mid, Up1, Up2, Dn2, Dn1))
            # Create banded sparse matrix
            # Rows:
            #      Lower left
            #      Middle
            #      Upper right
            self.coeff_matrix = scipy.sparse.spdiags(
                self.diags,
                [
                    self.nx - self.ny * self.nx - 1,
                    self.nx - self.ny * self.nx,
                    self.nx - self.ny * self.nx + 1,
                    2 * self.nx - self.ny * self.nx,
                    -2 * self.nx,
                    -self.nx - 1,
                    -self.nx,
                    -self.nx + 1,
                    -2,
                    -1,
                    0,
                    1,
                    2,
                    self.nx - 1,
                    self.nx,
                    self.nx + 1,
                    2 * self.nx,
                    self.ny * self.nx - 2 * self.nx,
                    self.ny * self.nx - self.nx - 1,
                    self.ny * self.nx - self.nx,
                    self.ny * self.nx - self.nx + 1,
                ],
                self.ny * self.nx,
                self.ny * self.nx,
                format="csr",
            )

        else:
            # No periodic boundary conditions -- original form of coeff_matrix
            # creator.
            # Arrange in solver
            self.diags = np.vstack((Dn2, Dn1, Mid, Up1, Up2))
            # Create banded sparse matrix
            self.coeff_matrix = scipy.sparse.spdiags(
                self.diags,
                [
                    -2 * self.nx,
                    -self.nx - 1,
                    -self.nx,
                    -self.nx + 1,
                    -2,
                    -1,
                    0,
                    1,
                    2,
                    self.nx - 1,
                    self.nx,
                    self.nx + 1,
                    2 * self.nx,
                ],
                self.ny * self.nx,
                self.ny * self.nx,
                format="csr",
            )  # create banded sparse matrix

    def calc_max_flexural_wavelength(self):
        """
        Returns the approximate maximum flexural wavelength
        This is important when padding of the grid is required: in Flexure (this
        code), grids are padded out to one maximum flexural wavelength, but in any
        case, the flexural wavelength is a good characteristic distance for any
        truncation limit
        """
        if np.isscalar(self.D):
            Dmax = self.D
        else:
            Dmax = self.D.max()
        # This is an approximation if there is fill that evolves with iterations
        # (e.g., water), but should be good enough that this won't do much to it
        alpha = (4 * Dmax / (self.drho * self.g)) ** 0.25  # 2D flexural parameter
        self.max_flexural_wavelength = 2 * np.pi * alpha
        self.maxFlexuralWavelength_ncells_x = int(
            np.ceil(self.max_flexural_wavelength / self.dx)
        )
        self.maxFlexuralWavelength_ncells_y = int(
            np.ceil(self.max_flexural_wavelength / self.dy)
        )

    def fd_solve(self):
        """
        w = fd_solve()
        Sparse flexural response calculation via direct factorization with
        UMFpack.
        Requires the coefficient matrix from "2D.coeff_matrix"
        """

        if self.debug:
            # Will fail if scalar
            with contextlib.suppress(AttributeError):
                _logger.debug("self.T_e %s", self.T_e.shape)
            _logger.debug("self.qs %s", self.qs.shape)
            self.calc_max_flexural_wavelength()
            _logger.debug(
                "maxFlexuralWavelength_ncells: (x, y): %s %s",
                self.maxFlexuralWavelength_ncells_x,
                self.maxFlexuralWavelength_ncells_y,
            )
            _logger.debug("Using direct solution with UMFpack")
        elif self.solver != "direct":
            raise ValueError(
                f"solver={self.solver!r} is not supported; only 'direct' is available "
                "in this release.  An iterative solver may be added in a future version."
            )

        if self.cache_factorization not in (False, True, "no_check"):
            raise ValueError(
                f"cache_factorization must be False, True, or 'no_check'; "
                f"got {self.cache_factorization!r}"
            )

        # qs negative so the plate bends down under a positive (downward) load
        # and up under a negative load (material removed).  The coefficient
        # matrix A encodes D∇⁴ + Δρg, which is positive definite; solving
        # A·w = −q therefore gives w < 0 for q > 0, matching the
        # positive-upward sign convention used throughout gFlex.
        q0vector = -self.qs.reshape(-1, order="C")
        rhs_corr = getattr(self, "_bc_rhs_correction", None)
        if rhs_corr is not None:
            q0vector = q0vector + rhs_corr.reshape(-1, order="C")

        _ls_start = time.perf_counter()
        if self.cache_factorization is False:
            wvector = spsolve(self.coeff_matrix, q0vector, use_umfpack=True)
        elif self.cache_factorization == "no_check":
            if self._lu is None:
                self._lu = factorized(self.coeff_matrix)
                self.coeff_matrix = None  # _lu is sole owner; _solve_fd uses _lu as rebuild-skip signal
            wvector = self._lu(q0vector)
        else:  # True: hash-validated cache
            h = _matrix_hash(self.coeff_matrix)
            if self._lu is None or h != self._lu_matrix_hash:
                self._lu = factorized(self.coeff_matrix)
                self._lu_matrix_hash = h
            wvector = self._lu(q0vector)
        self.linear_solve_time = time.perf_counter() - _ls_start

        # Reshape into grid
        self.w = wvector.reshape(self.qs.shape)
        self.w_padded = self.w.copy()  # for troubleshooting

        # Time to solve used to be here
