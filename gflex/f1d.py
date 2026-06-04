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
import sys
import time
import warnings

import numpy as np
import scipy.fft
from scipy.signal import fftconvolve
from scipy.sparse import spdiags
from scipy.sparse.linalg import factorized, spsolve

from gflex.base import Flexure, _RigidityBC, _matrix_hash


def recommended_pad_width_1d(Te, dx, E=65e9, nu=0.25, rho_m=3300.0, rho_fill=0.0,
                              g=9.8, n_wavelengths=1.0):
    """
    Return the recommended padding width in grid cells for a 1-D variable-Te run.

    The padded domain boundary should be at least one 1-D flexural wavelength
    from the load so that the plate's response is negligible at the boundary.

    The 1-D flexural wavelength is computed from the maximum Te value, giving
    the most conservative (widest) padding estimate.

    Parameters
    ----------
    Te : scalar or 1-D array
        Elastic thickness [m].  The maximum value is used.
    dx : float
        Grid cell size [m].
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

    Examples
    --------
    >>> recommended_pad_width_1d(Te=35e3, dx=5000.)
    94
    """
    drho = rho_m - rho_fill
    D = E * float(np.max(Te)) ** 3 / (12.0 * (1.0 - nu**2))
    alpha_1D = (4.0 * D / (drho * g)) ** 0.25
    lambda_1D = 2.0 * np.pi * alpha_1D
    return int(np.ceil(n_wavelengths * lambda_1D / dx))


def smooth_pad_Te_1d(Te, pad_width, Te_out=None):
    """
    Pad a 1-D elastic thickness array with a smooth linear taper.

    When a spatially variable Te array is padded with a constant value before
    being passed to :class:`F1D`, the abrupt step in flexural rigidity D at
    the inner/outer boundary drives spurious deflections via the D-derivative
    terms in the FD stencil.

    This function eliminates that step by linearly blending the inner-domain
    edge values toward *Te_out* across the padding region.

    The corresponding surface load array should be padded with zeros, e.g.::

        qs_padded = numpy.pad(qs, pad_width, mode='constant')

    Parameters
    ----------
    Te : 1-D array
        Elastic thickness [m] for the inner domain.
    pad_width : int
        Width of the padding on each end in grid cells.  Use
        :func:`recommended_pad_width_1d` to obtain a suitable value.
    Te_out : float, optional
        Te value at the outer edge of the padding.
        Defaults to ``Te.mean()``.

    Returns
    -------
    Te_padded : 1-D array of length ``len(Te) + 2 * pad_width``
        Padded elastic thickness with a smooth linear taper.

    Examples
    --------
    >>> import numpy as np
    >>> Te = 35e3 * np.ones(10)
    >>> Te_pad = smooth_pad_Te_1d(Te, pad_width=4)
    >>> Te_pad.shape
    (18,)

    To run :class:`F1D` with the padded arrays::

        import numpy as np
        from gflex import F1D, smooth_pad_Te_1d, recommended_pad_width_1d
        p = recommended_pad_width_1d(Te, dx=5000.)
        Te_pad = smooth_pad_Te_1d(Te, p)
        qs_pad = np.pad(qs, p, mode='constant')
        flex = F1D()
        flex.te = Te_pad
        flex.qs = qs_pad
        # ... set other parameters and run ...
        w_inner = flex.w[p:-p]   # trim padding from output
    """
    Te = np.asarray(Te, dtype=float)
    if Te.ndim != 1:
        raise ValueError("Te must be a 1-D array")
    if pad_width < 1:
        raise ValueError("pad_width must be >= 1")

    if Te_out is None:
        Te_out = float(Te.mean())
    Te_out = float(Te_out)

    nx = Te.shape[0]
    p = pad_width
    nx_p = nx + 2 * p

    Te_padded = np.full(nx_p, Te_out)
    Te_padded[p:-p] = Te

    for k in range(p):
        alpha = k / p
        Te_padded[k] = (1.0 - alpha) * Te_out + alpha * Te[0]
        Te_padded[nx_p - 1 - k] = (1.0 - alpha) * Te_out + alpha * Te[-1]

    return Te_padded


def pad_domain_1d(Te, qs, dx, n_wavelengths=1.0, Te_out=None,
                  E=65e9, nu=0.25, rho_m=3300.0, rho_fill=0.0, g=9.8):
    """
    Pad both the elastic thickness and surface load arrays for use with F1D.

    Combines :func:`recommended_pad_width_1d` and :func:`smooth_pad_Te_1d`
    into a single call, and zero-pads *qs* to match.  The returned pad width
    *p* can be used to trim the deflection output after the run::

        w_inner = flex.w[p:-p]

    Parameters
    ----------
    Te : 1-D array
        Elastic thickness [m] for the inner domain.
    qs : 1-D array
        Surface load [Pa] for the inner domain.
    dx : float
        Grid cell size [m].
    n_wavelengths : float, optional
        Padding width expressed as a number of flexural wavelengths.
        Default 1.0; use 0.5 for a less conservative (narrower) padding.
    Te_out : float, optional
        Te value at the outer edge of the padding.
        Defaults to ``Te.mean()``.
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
    Te_padded : 1-D array of length ``len(Te) + 2p``
        Smoothly tapered elastic thickness.
    qs_padded : 1-D array of length ``len(qs) + 2p``
        Surface load zero-padded to match.
    p : int
        Pad width in grid cells (same on both ends).

    Examples
    --------
    >>> import numpy as np
    >>> Te = 10e3 * np.ones(5)
    >>> qs = np.zeros(5)
    >>> Te_pad, qs_pad, p = pad_domain_1d(Te, qs, dx=10000.,
    ...     E=65e9, nu=0.25, rho_m=3300., rho_fill=0., g=9.8)
    >>> p
    19
    >>> Te_pad.shape
    (43,)

    To run :class:`F1D` with the padded arrays::

        flex = gflex.F1D()
        flex.te = Te_pad
        flex.qs = qs_pad
        # ... set other parameters and run ...
        w_inner = flex.w[p:-p]   # trim padding from output
    """
    p = recommended_pad_width_1d(
        Te, dx, E=E, nu=nu, rho_m=rho_m, rho_fill=rho_fill,
        g=g, n_wavelengths=n_wavelengths,
    )
    Te_padded = smooth_pad_Te_1d(Te, p, Te_out=Te_out)
    qs_padded = np.pad(qs, p, mode="constant")
    return Te_padded, qs_padded, p


def _sandbox_easter_egg():
    """You found it. The sandbox is not closed — it flexes."""
    import sys
    import warnings

    # ── Physics: 2 m sandbox, 12 mm plywood floor, 4 kg sand castle ─────────
    nx = 200
    dx = 0.01
    qs_load = np.zeros(nx)
    c = nx // 2
    half_load = 10  # 20 cm footprint
    qs_load[c - half_load : c + half_load] = 4.0 * 9.81 / (2 * half_load * dx)

    solver = F1D()
    solver.dx = dx
    solver.qs = qs_load
    solver.te = 0.012      # 12 mm plywood
    solver.E = 10e9
    solver.nu = 0.30
    solver.rho_m = 1.0     # negligible buoyancy (floor on a table)
    solver.rho_fill = 0.0
    solver.g = 9.81
    solver.method = "fd"
    solver.bc_west = "zero_displacement_zero_slope"
    solver.bc_east = "zero_displacement_zero_slope"
    solver.quiet = True
    solver.verbose = False
    solver.debug = False
    solver.initialize()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        solver.run()
    solver.finalize()

    w_sag = -solver.w           # positive = downward
    w_max_mm = float(w_sag.max() * 1000.0)

    # ── Canvas ───────────────────────────────────────────────────────────────
    INT_W = 58
    NROWS = 17
    WALL_L = 1
    WALL_R = WALL_L + INT_W + 1   # = 60

    BOX_INTERIOR_TOP = 8    # first row inside the box at the edges (sag = 0)
    CASTLE_HEIGHT = 9
    SAND_HEIGHT = 3
    MAX_SAG = 3

    idx = np.linspace(0, nx - 1, INT_W).astype(int)
    w_col = w_sag[idx]
    w_col_max = w_col.max()
    if w_col_max > 1e-15:
        sag = np.round(w_col / w_col_max * MAX_SAG).astype(int)
    else:
        sag = np.zeros(INT_W, dtype=int)

    center_sag = int(sag[INT_W // 2])
    castle_base_row = BOX_INTERIOR_TOP + center_sag
    castle_top_row = castle_base_row - (CASTLE_HEIGHT - 1)

    buf = [[' '] * (WALL_R + 2) for _ in range(NROWS)]

    def put(r, c, text):
        for k, ch in enumerate(text):
            if 0 <= r < NROWS and 0 <= c + k < len(buf[0]):
                buf[r][c + k] = ch

    cx = WALL_L + 1 + INT_W // 2

    # Castle: 11 chars wide, 9 rows tall
    castle = [
        "     |     ",
        "    [F]    ",
        " _  |||  _ ",
        "| |     | |",
        "| |  O  | |",
        "| |     | |",
        "| |  |  | |",
        "| |  |  | |",
        "|_|__|__|_|",
    ]
    for i, line in enumerate(castle):
        put(castle_top_row + i, cx - len(line) // 2, line)

    # Box walls (no bottom bar)
    for r in range(BOX_INTERIOR_TOP, NROWS):
        buf[r][WALL_L] = '|'
        buf[r][WALL_R] = '|'

    # Sand grains and floor, column by column
    for i in range(INT_W):
        col = WALL_L + 1 + i
        s = int(sag[i])

        # Sand: alternating dot pattern offset per row → scattered grain look
        for sand_row in range(SAND_HEIGHT):
            r = BOX_INTERIOR_TOP + s + sand_row
            if 0 <= r < NROWS:
                grain = '.' if (col + sand_row) % 2 == 0 else ' '
                if buf[r][col] == ' ':
                    buf[r][col] = grain

        # Floor
        r = BOX_INTERIOR_TOP + SAND_HEIGHT + s
        if 0 <= r < NROWS:
            buf[r][col] = '='

    print()
    for row in buf:
        print(''.join(row))
    print()
    print(f"  gFlex (FD, clamped ends): {w_max_mm:.3f} mm deflection at centre.")
    print(f"  A 4 kg sand castle on 12 mm plywood flexes a 2 m sandbox by {w_max_mm:.3f} mm.")
    print()
    sys.exit(0)


class F1D(Flexure):
    """
    One-dimensional lithospheric flexure solver.

    Computes the deflection *w(x)* of a thin elastic beam overlying an
    inviscid fluid (mantle) given a surface load stress *qs*.  Supports
    spatially variable elastic thickness *Te*.

    Set instance attributes, then call :meth:`initialize`, :meth:`run`, and
    :meth:`finalize` in sequence.  Read ``flex.w`` **before** calling
    :meth:`finalize`; finalize clears all model state including ``w``.

    Attributes
    ----------
    method : str
        Solution method.  ``'fd'`` (finite difference, supports variable
        *Te*), ``'fft'`` (spectral, requires scalar *Te*),
        ``'sas'`` (superposition of analytical solutions, constant *Te*
        only), or ``'sas_ng'`` (SAS on an ungridded point array).
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
    te : float or ndarray of shape (N,)
        Elastic thickness [m].  A scalar is broadcast to the full grid.
    qs : ndarray of shape (N,)
        Surface load stress [Pa].
    dx : float
        Grid spacing [m].
    bc_west, bc_east : str
        Boundary conditions on the west (left) and east (right) ends.
        FD options: ``'zero_displacement_zero_slope'`` (alias ``'clamped'``),
        ``'zero_displacement_zero_moment'``,
        ``'zero_moment_zero_shear'`` (alias ``'free'``),
        ``'mirror'``, ``'periodic'``, ``'sandbox'``.
        ``'zero_slope_zero_shear'`` is a deprecated alias for ``'mirror'``.
        SAS option: ``'no_outside_loads'`` (the default when unset).
    sigma_xx : float, optional
        Normal stress applied at the plate ends [Pa].  FD only.
    quiet : bool
        Suppress timing output.  Default ``False``.
    verbose : bool
        Print progress messages.  Default ``True``.

    Examples
    --------
    Minimal finite-difference run::

        import numpy as np
        from gflex import F1D

        flex = F1D()
        flex.quiet = True
        flex.method = 'fd'
        flex.solver = 'direct'
        flex.g = 9.8
        flex.E = 65e9
        flex.nu = 0.25
        flex.rho_m = 3300.
        flex.rho_fill = 1000.
        flex.te = 30e3
        flex.qs = np.zeros(300)
        flex.qs[100:200] = 1e6      # 100-cell load
        flex.dx = 4000.             # 4 km grid
        flex.bc_west = 'zero_displacement_zero_slope'
        flex.bc_east = 'zero_moment_zero_shear'
        flex.initialize()
        flex.run()
        flex.finalize()
        deflection = flex.w         # (300,) array, negative downward
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
        self.dimension = 1  # Set it here in case it wasn't set for selection before
        super().initialize()
        if self.verbose:
            print("F1D initialized")

    def run(self):
        """
        Execute the flexural solution.

        Selects and runs the method specified by ``self.method``.  The
        deflection array is stored in ``self.w`` on return.  Call
        :meth:`finalize` afterwards to restore any internally modified
        state.
        """
        self.bc_check()
        self.solver_start_time = time.time()
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
            # nonuniform points
            super()._solve_sas_ng()
            self.method_func = self._solve_sas_ng
        else:
            raise ValueError('method must be "fd", "fft", "sas", or "sas_ng"')

        if self.verbose:
            print("F1D run")
        self.method_func()

        self.time_to_solve = time.time() - self.solver_start_time
        if not self.quiet:
            print("Time to solve [s]:", self.time_to_solve)

    def finalize(self):
        """
        Release all model state.

        Restores ``self.te`` to its pre-run value if gFlex padded it
        internally.  Then calls the base ``finalize``, which deletes
        ``w``, ``qs``, and the cached coefficient matrix.  Read ``w``
        before calling this method.
        """
        with contextlib.suppress(AttributeError):
            self.te = self.te_unpadded
        if self.verbose:
            print("F1D finalized")
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
        bc_sides = {"W": self._bc_west_norm, "E": self._bc_east_norm}
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
        loaded = np.nonzero(self.qs)[0]
        if loaded.size == 0:
            return
        nx = self.qs.shape[0]
        Te_arr = (
            self.te
            if isinstance(self.te, np.ndarray)
            else np.full(nx, float(self.te))
        )
        Te_loaded = Te_arr[loaded]
        D_loaded = self.E * Te_loaded**3 / (12 * (1 - self.nu**2))
        alpha_loaded = (4 * D_loaded / (self.drho * self.g)) ** 0.25
        wavelength_loaded = 2 * np.pi * alpha_loaded
        dist_fns = {
            "W": lambda: (loaded + 0.5) * self.dx,
            "E": lambda: (nx - loaded - 0.5) * self.dx,
        }
        for side, bc in bc_sides.items():
            if bc != "zero_displacement_zero_slope":
                continue
            distances = dist_fns[side]()
            frac = distances / wavelength_loaded
            worst = np.argmin(frac)
            if frac[worst] < 1.0:
                i = loaded[worst]
                warnings.warn(
                    f"BC_{side} = 'zero_displacement_zero_slope': nearest loaded cell "
                    f"(index {i}) is {distances[worst]/1e3:.1f} km from the boundary "
                    f"({frac[worst]:.2f} flexural wavelengths; "
                    f"wavelength ≈ {wavelength_loaded[worst]/1e3:.1f} km "
                    f"at Te = {Te_loaded[worst]/1e3:.1f} km). "
                    "The flexural forebulge peaks near one wavelength from the load "
                    "and will be suppressed by this boundary. "
                    "Use pad_domain_1d() to extend the domain before solving, "
                    "then trim w to the original extent.",
                    UserWarning,
                    stacklevel=4,
                )

    def _solve_fd(self):
        """Run the finite-difference solution pipeline.

        Calls :meth:`_check_warnings_FD` first to flag potentially problematic
        boundary conditions, then assembles and solves the sparse banded system.
        After the call, the deflection is available in ``self.w``.
        """
        self._check_warnings_FD()
        self.gridded_x()
        # Only generate coefficient matrix if it is not already provided.
        # In no_check mode the matrix is freed after factorization to save
        # memory; _lu being set is sufficient to skip the rebuild.
        if self.coeff_matrix is not None:
            pass
        elif self.cache_factorization == "no_check" and self._lu is not None:
            pass  # coeff_matrix was freed after factorization; _lu still valid
        else:
            self.elasprepFD()  # define dx4 and D within self
            self._build_coefficient_matrix()
        self.fd_solve()  # Get the deflection, "w"

    def _solve_fft(self):
        """Spectral (FFT) flexural solution for uniform elastic thickness.

        Applies the analytical transfer function in the wavenumber domain::

            W(k) = -Q(k) / (D k⁴ + σ_xx T_e k² + Δρ g)

        **Boundary conditions and periodicity**

        FFT inherently assumes a periodic domain.  Two modes are supported:

        * ``BC_W = BC_E = 'periodic'`` — the load array is used as-is.
          The solution is exact for a load that genuinely repeats with
          period L = N · dx.

        * Any other BC (including ``'no_outside_loads'`` or unset) — the
          load array is zero-padded by four flexural wavelengths (α) on
          each side before the transform, then trimmed back to the original
          length afterwards.  This is mathematically identical to the
          periodic case, but with a padded domain large enough that the
          periodic images of the load are separated by ~8α of zeros and
          therefore do not influence the interior solution.  It is the
          spectral equivalent of the ``'no_outside_loads'`` boundary
          condition used by the SAS solver.

        In both cases the solution is spectral (no spatial discretisation
        error in the transfer function itself); any residual error comes
        from representing a continuous load on a discrete grid.

        Requires uniform (scalar) elastic thickness; for variable *Te* use
        the finite-difference method instead.
        """
        self.gridded_x()

        # Te must be scalar or a uniform array
        if np.isscalar(self.te):
            pass
        elif np.all(self.te == np.mean(self.te)):
            self.te = float(np.mean(self.te))
        else:
            raise ValueError(
                "The FFT solution requires a scalar (uniform) Te. "
                "For spatially variable Te, use the finite difference method."
            )

        D = self.E * self.te**3 / (12.0 * (1.0 - self.nu**2))
        # 1-D flexural parameter α = (4D / Δρg)^0.25
        alpha = (4.0 * D / (self.drho * self.g)) ** 0.25

        periodic = (self.bc_west == "periodic") and (self.bc_east == "periodic")

        if periodic:
            qs_work = self.qs
        else:
            # Zero-pad by 4α on each side (no_outside_loads assumption)
            pad = int(np.ceil(4.0 * alpha / self.dx))
            qs_work = np.pad(self.qs, pad, mode="constant")

        N_work = len(qs_work)
        k = scipy.fft.rfftfreq(N_work, d=self.dx) * 2.0 * np.pi
        Q = scipy.fft.rfft(qs_work, workers=-1)
        denom = D * k**4 + self.sigma_xx * self.te * k**2 + self.drho * self.g
        w_work = scipy.fft.irfft(-Q / denom, n=N_work, workers=-1)

        if periodic:
            self.w = w_work
        else:
            self.w = w_work[pad : pad + self.nx]

    def _solve_sas(self):
        """Run the gridded superposition-of-analytical-solutions pipeline."""
        self.gridded_x()
        self.spatial_domain_vars_sas()
        self.spatial_domain_gridded()

    def _solve_sas_ng(self):
        """Run the ungridded (non-uniform points) SAS pipeline."""
        self.spatial_domain_vars_sas()
        self.spatial_domain_no_grid()

    ######################################
    ## FUNCTIONS TO SOLVE THE EQUATIONS ##
    ######################################

    ## UTILITY
    ############

    def gridded_x(self):
        """Build the x-coordinate array from ``qs`` shape and ``dx``."""
        self.nx = self.qs.shape[0]
        self._x_local = np.arange(0, self.dx * self.nx, self.dx)

    ## SPATIAL DOMAIN SUPERPOSITION OF ANALYTICAL SOLUTIONS
    #########################################################

    # SETUP

    def spatial_domain_vars_sas(self):
        """Compute flexural rigidity D, parameter alpha, and Green's-function coefficient for SAS."""
        # Check Te:
        # * If scalar, okay.
        # * If grid, convert to scalar if a singular value
        # * Else, throw an error.
        if np.isscalar(self.te):
            pass
        elif np.all(self.te == np.mean(self.te)):
            self.te = np.mean(self.te)
        else:
            raise ValueError(
                "The analytical solution requires a scalar (uniform) Te. "
                "For spatially variable Te, use the finite difference method."
            )

        self.D = self.E * self.te**3 / (12 * (1 - self.nu**2))  # Flexural rigidity
        self.alpha = (
            4 * self.D / (self.drho * self.g)
        ) ** 0.25  # 1D flexural parameter
        self.coeff = self.alpha**3 / (8 * self.D)

    # UNIFORM DX ("GRIDDED"): LOADS PROVIDED AS AN ARRAY WITH KNOWN DX TO
    # CONVERT LOAD MAGNITUDE AT A POINT INTO MASS INTEGRATED ACROSS DX

    def spatial_domain_gridded(self):
        """Compute deflection by summing 1D Green's functions over the load grid."""
        # Build the beam Green's function kernel for all relative offsets
        # [-(nx-1), nx-1].  G(r) = exp(-r/α)(cos(r/α) + sin(r/α)) is even,
        # so the kernel is symmetric about its centre.
        r = np.abs(np.arange(-(self.nx - 1), self.nx)) * self.dx
        kernel = np.exp(-r / self.alpha) * (
            np.cos(r / self.alpha) + np.sin(r / self.alpha)
        )
        # fftconvolve is identical to the loop (positive load → negative
        # deflection, hence the minus sign) but O(N log N) instead of O(N²).
        self.w = -self.coeff * fftconvolve(self.qs * self.dx, kernel, mode="same")

    # NONUNIFORM DX (NO GRID): ARBITRARILY-SPACED POINT LOADS
    # So essentially a sum of Green's functions for flexural response

    def spatial_domain_no_grid(self):
        """
        Superposition of analytical solutions without a gridded domain
        """
        if self.debug:
            print("w = ")
            print(self.xw.shape)

        # dist shape: (N_out, N_load); vectorised over both dimensions at once.
        dist = np.abs(self.xw[:, None] - self.x[None, :])
        r = dist / self.alpha
        G = np.exp(-r) * (np.cos(r) + np.sin(r))   # (N_out, N_load)
        self.w = -self.coeff * (G * self.q[None, :]).sum(axis=1)

    ## FINITE DIFFERENCE
    ######################

    def elasprepFD(self):
        """Precompute dx⁴, dx², and the flexural rigidity array D for the FD solver."""
        self.dx4 = self.dx**4
        self.dx2 = self.dx**2  # Needed if horizontal (i.e., tectonic) stresses
        self.D = self.E * self.te**3 / (12 * (1 - self.nu**2))

    def _build_coefficient_matrix(self):
        """
        Selects the boundary conditions
        Then calls the function to build the pentadiagonal matrix to solve
        1D flexure with variable (or constant) elastic thickness
        """

        # Zeroth, start the timer and print the boundary conditions to the screen
        self.coeff_start_time = time.time()
        if self.verbose:
            print("Boundary condition, West:", self.bc_west, type(self.bc_west))
            print("Boundary condition, East:", self.bc_east, type(self.bc_east))

        if self.bc_east == "sandbox" or self.bc_west == "sandbox":
            _sandbox_easter_egg()

        # First, set flexural rigidity boundary conditions to flesh out this padded
        # array
        self._apply_bc_rigidity()

        # Second, build the coefficient arrays -- with the rigidity b.c.'s
        self.get_coeff_values()

        # Third, apply boundary conditions to the coeff_arrays to create the
        # flexural solution
        self._apply_bc_flexure()

        # Fourth, compute RHS correction for prescribed (nonzero) BC values
        self._apply_bc_rhs_inhomogeneous()

        # Fifth, construct the sparse diagonal array
        self.build_diagonals()

        # Finally, compute the total time this process took
        self.coeff_creation_time = time.time() - self.coeff_start_time
        if not self.quiet:
            print(
                "Time to construct coefficient (operator) array [s]:",
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
        elif self._bc_west_norm in ("mirror", "zero_displacement_zero_moment"):
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
        elif self._bc_east_norm in ("mirror", "zero_displacement_zero_moment"):
            self.bc_rigidity_east = _RigidityBC.MIRROR
        else:
            raise RuntimeError("Invalid Te B.C. case")

        #############
        # PAD ARRAY #
        #############
        if np.isscalar(self.te):
            self.D *= np.ones(self.qs.shape)  # And leave Te as a scalar for checks
        else:
            self.te_unpadded = self.te.copy()
        # F2D keeps this inside the "else" and handles this differently,
        # largely because it has different ways of computing the flexural
        # response with variable Te. We'll keep everything simpler here and
        # just pad this array so it can be sent through the same process
        # to create the coefficient arrays.
        self.D = np.hstack([np.nan, self.D, np.nan])

        ###############################################################
        # APPLY FLEXURAL RIGIDITY BOUNDARY CONDITIONS TO PADDED ARRAY #
        ###############################################################
        if self.bc_rigidity_west == _RigidityBC.ZERO_CURVATURE:
            self.D[0] = 2 * self.D[1] - self.D[2]
        if self.bc_rigidity_east == _RigidityBC.ZERO_CURVATURE:
            self.D[-1] = 2 * self.D[-2] - self.D[-3]
        if self.bc_rigidity_west == _RigidityBC.MIRROR:
            self.D[0] = self.D[2]
        if self.bc_rigidity_east == _RigidityBC.MIRROR:
            self.D[-1] = self.D[-3]
        if self.bc_rigidity_west == _RigidityBC.PERIODIC:
            self.D[0] = self.D[-2]
        if self.bc_rigidity_east == _RigidityBC.PERIODIC:
            self.D[-1] = self.D[-3]

    def get_coeff_values(self):
        """Build the five pentadiagonal coefficient arrays for the FD stencil."""
        ##############################
        # BUILD GENERAL COEFFICIENTS #
        ##############################

        # l2 corresponds to top value in solution vector, so to the left (-) side
        # Good reference for how to determine central difference (and other) coefficients is:
        # Fornberg, 1998: Generation of Finite Difference Formulas on Arbitrarily Spaced Grids

        ###################################################
        # DEFINE SUB-ARRAYS FOR DERIVATIVE DISCRETIZATION #
        ###################################################
        Dm1 = self.D[:-2]
        D0 = self.D[1:-1]
        Dp1 = self.D[2:]

        ###########################################################
        # DEFINE COEFFICIENTS TO W_-2 -- W_+2 WITH B.C.'S APPLIED #
        ###########################################################
        self.l2_coeff_i = (Dm1 / 2.0 + D0 - Dp1 / 2.0) / self.dx4
        self.l1_coeff_i = (
            -6.0 * D0 + 2.0 * Dp1
        ) / self.dx4 - self.sigma_xx * self.te / self.dx2
        self.c0_coeff_i = (
            (-2.0 * Dm1 + 10.0 * D0 - 2.0 * Dp1) / self.dx4
            + 2 * self.sigma_xx * self.te / self.dx2
            + self.drho * self.g
        )
        self.r1_coeff_i = (
            2.0 * Dm1 - 6.0 * D0
        ) / self.dx4 - self.sigma_xx * self.te / self.dx2
        self.r2_coeff_i = (-Dm1 / 2.0 + D0 + Dp1 / 2.0) / self.dx4
        # These will be just the 1, -4, 6, -4, 1 for constant Te

        ###################################################################
        # START DIAGONALS AS SIMPLY THE BASE COEFFICIENTS, WITH NO B.C.'S #
        ###################################################################
        self.l2 = self.l2_coeff_i.copy()
        self.l1 = self.l1_coeff_i.copy()
        self.c0 = self.c0_coeff_i.copy()
        self.r1 = self.r1_coeff_i.copy()
        self.r2 = self.r2_coeff_i.copy()

        # Number of columns; equals number of rows too - square coeff matrix
        self.ncolsx = self.c0.shape[0]

        # Either way, the way that Scipy stacks is not the same way that I calculate
        # the rows. It runs offsets down the column instead of across the row. So
        # to simulate this, I need to re-zero everything. To do so, I use
        # numpy.roll. (See self.build_diagonals.)

    def _apply_bc_flexure(self):
        """Apply flexural boundary conditions to the coefficient diagonals."""
        # Some links that helped me teach myself how to set up the boundary conditions
        # in the matrix for the flexure problem:
        #
        # Good explanation of and examples of boundary conditions
        # https://en.wikipedia.org/wiki/Euler%E2%80%93Bernoulli_beam_theory#Boundary_considerations
        #
        # Copy of Fornberg table:
        # https://en.wikipedia.org/wiki/Finite_difference_coefficient
        #
        # Implementing b.c.'s:
        # http://scicomp.stackexchange.com/questions/5355/writing-the-poisson-equation-finite-difference-matrix-with-neumann-boundary-cond
        # http://scicomp.stackexchange.com/questions/7175/trouble-implementing-neumann-boundary-conditions-because-the-ghost-points-cannot

        if self.verbose:
            print("Boundary condition, West:", self.bc_west, type(self.bc_west))
            print("Boundary condition, East:", self.bc_east, type(self.bc_east))

        # In 2D, these are handled inside the function; in 1D, there are separate
        # defined functions. Keeping these due to inertia and fear of cut/paste
        # mistakes
        if self._bc_east_norm == "zero_displacement_zero_slope" or self._bc_west_norm == "zero_displacement_zero_slope":
            self._bc_zero_displacement_zero_slope()
        if self._bc_east_norm == "zero_moment_zero_shear" or self._bc_west_norm == "zero_moment_zero_shear":
            self._bc_zero_moment_zero_shear()
        if self._bc_east_norm == "mirror" or self._bc_west_norm == "mirror":
            self._bc_mirror()
        if self._bc_east_norm == "zero_displacement_zero_moment" or self._bc_west_norm == "zero_displacement_zero_moment":
            self._bc_zero_displacement_zero_moment()
        if self._bc_east_norm == "periodic" and self._bc_west_norm == "periodic":
            self._bc_periodic()

    def build_diagonals(self):
        """
        Builds the diagonals for the coefficient array
        """

        ##########################################################
        # INCORPORATE BOUNDARY CONDITIONS INTO COEFFICIENT ARRAY #
        ##########################################################

        # Roll to keep the proper coefficients at the proper places in the
        # arrays: Python will naturally just do vertical shifts instead of
        # diagonal shifts, so this takes into account the horizontal compoent
        # to ensure that boundary values are at the right place.
        self.l2 = np.roll(self.l2, -2)
        self.l1 = np.roll(self.l1, -1)
        self.r1 = np.roll(self.r1, 1)
        self.r2 = np.roll(self.r2, 2)

        # Then assemble these rows: this is where the periodic boundary condition
        # can matter.
        if self.coeff_matrix is not None:
            pass
        elif self.bc_east == "periodic" and self.bc_west == "periodic":
            # In this case, the boundary-condition-related stacking has already
            # happened inside b.c.-handling function. This is because periodic
            # boundary conditions require extra diagonals to exist on the edges of
            # the solution array
            pass
        else:
            self.diags = np.vstack((self.l2, self.l1, self.c0, self.r1, self.r2))
            self.offsets = np.array([-2, -1, 0, 1, 2])

        # Everybody now (including periodic b.c. cases)
        self.coeff_matrix = spdiags(
            self.diags, self.offsets, self.nx, self.nx, format="csr"
        )

    def _bc_periodic(self):
        """
        periodic boundary conditions: wraparound to the other side.
        """
        if self._bc_east_norm == "periodic" and self._bc_west_norm == "periodic":
            # If both boundaries are periodic, we are good to go (and self-consistent)
            pass  # It is just a shift in the coeff. matrix creation.
        else:
            # If only one boundary is periodic and the other doesn't implicitly
            # involve a periodic boundary, this is illegal!
            # I could allow it, but would have to rewrite the periodic b.c. case,
            # which I don't want to do to allow something that doesn't make
            # physical sense... so if anyone wants to do this for some unforeseen
            # reason, they can just split my function into two pieces themselves.i
            raise ValueError(
                "The boundary opposite a periodic boundary condition must also be "
                "periodic. Having one periodic and one non-periodic boundary on the "
                "same axis is not physically meaningful."
            )
        self.diags = np.vstack(
            (
                self.r1,
                self.r2,
                self.l2,
                self.l1,
                self.c0,
                self.r1,
                self.r2,
                self.l2,
                self.l1,
            )
        )
        self.offsets = np.array(
            [
                1 - self.ncolsx,
                2 - self.ncolsx,
                -2,
                -1,
                0,
                1,
                2,
                self.ncolsx - 2,
                self.ncolsx - 1,
            ]
        )

    def _bc_zero_displacement_zero_slope(self):
        """
        Clamped boundary: zero displacement and zero slope.

        **Boundary node** (i=0 west, i=N-1 east): decoupled from all interior
        nodes so the row reduces to c0·w = 0, enforcing w = 0 exactly for
        zero boundary load.

        **First interior node** (i=1 west, i=N-2 east): the ghost node one
        step outside the domain is eliminated via even reflection
        (w[ghost] = w[interior]), which folds its stencil coefficient into c0.
        This encodes dw/dx = 0 at the boundary.

        Historical note: the original implementation (pre-2026) dropped ghost
        nodes silently rather than reflecting them, and did not decouple the
        boundary row from interior nodes.  The result was a ghost = 0
        truncation that approximated clamped conditions only when the boundary
        was far from any load.  This is now corrected.
        """
        if self._bc_west_norm == "zero_displacement_zero_slope":
            # Boundary node: decouple so c0·w[0] = 0 → w[0] = 0 exactly
            i = 0
            self.l2[i] = np.nan   # ghost w[-2]: out-of-bounds, excluded
            self.l1[i] = np.nan   # ghost w[-1]: out-of-bounds, excluded
            self.c0[i] += 0
            self.r1[i] = 0        # decouple from w[1]
            self.r2[i] = 0        # decouple from w[2]
            # First interior node: even reflection w[-1] = w[1] encodes zero slope
            i = 1
            self.l2[i] = np.nan   # ghost absorbed into c0 below
            self.l1[i] += 0       # w[0] = 0, contributes nothing
            self.c0[i] += self.l2_coeff_i[i]   # fold w[-1] = +w[1] into c0
            self.r1[i] += 0
            self.r2[i] += 0
        if self._bc_east_norm == "zero_displacement_zero_slope":
            # First interior node: even reflection w[N] = w[N-2] encodes zero slope
            i = -2
            self.l2[i] += 0
            self.l1[i] += 0
            self.c0[i] += self.r2_coeff_i[i]   # fold w[N] = +w[N-2] into c0
            self.r1[i] += 0
            self.r2[i] = np.nan   # ghost absorbed into c0 above
            # Boundary node: decouple so c0·w[N-1] = 0 → w[N-1] = 0 exactly
            i = -1
            self.l2[i] = 0        # decouple from w[N-3]
            self.l1[i] = 0        # decouple from w[N-2]
            self.c0[i] += 0
            self.r1[i] = np.nan   # ghost w[N]: out-of-bounds, excluded
            self.r2[i] = np.nan   # ghost w[N+1]: out-of-bounds, excluded

    def _bc_zero_moment_zero_shear(self):
        """
        d2w/dx2 = d3w/dx3 = 0
        (no moment or shear)
        This simulates a free end (broken plate, end of a cantilevered beam:
        think diving board tip)
        It is *not* yet set up to have loads placed on the ends themselves:
        (look up how to do this, thought Wikipedia has some info, but can't find
        it... what I read said something about generalizing)
        """

        # First, just define coefficients for each of the positions in the array
        # These will be added in code instead of being directly combined by
        # the programmer (as I did above (now deleted) for constant Te), which might add
        # rather negligibly to the compute time but save a bunch of possibility
        # for unfortunate typos!

        # Also using 0-curvature boundary condition for D (i.e. Te)
        if self._bc_west_norm == "zero_moment_zero_shear":
            # Boundary node (i=0): eliminate ghosts w[-2] and w[-1] using the
            # moment condition (d²w/dx²=0 → w[-1]=2w[0]-w[1]) and the shear
            # condition (d³w/dx³=0 → w[-2]=4w[0]-4w[1]+w[2]).
            i = 0
            self.l2[i] += np.nan
            self.l1[i] += np.nan
            self.c0[i] += 4 * self.l2_coeff_i[i] + 2 * self.l1_coeff_i[i]
            self.r1[i] += -4 * self.l2_coeff_i[i] - self.l1_coeff_i[i]
            self.r2[i] += self.l2_coeff_i[i]
            # First interior node (i=1): eliminate ghost w[-1] using the same
            # moment condition: w[-1]=2w[0]-w[1] → folds l2*w[-1] into l1 and c0.
            i = 1
            self.l2[i] += np.nan
            self.l1[i] += 2 * self.l2_coeff_i[i]
            self.c0[i] -= self.l2_coeff_i[i]

        if self._bc_east_norm == "zero_moment_zero_shear":
            # First interior node (i=N-2): eliminate ghost w[N] using the moment
            # condition: w[N]=2w[N-1]-w[N-2] → folds r2*w[N] into r1 and c0.
            i = -2
            self.c0[i] -= self.r2_coeff_i[i]
            self.r1[i] += 2 * self.r2_coeff_i[i]
            self.r2[i] += np.nan
            # Boundary node (i=N-1): eliminate ghosts w[N] and w[N+1] using
            # the moment and shear conditions.
            i = -1
            self.l2[i] += self.r2_coeff_i[i]
            self.l1[i] += -4 * self.r2_coeff_i[i] - self.r1_coeff_i[i]
            self.c0[i] += 4 * self.r2_coeff_i[i] + 2 * self.r1_coeff_i[i]
            self.r1[i] += np.nan
            self.r2[i] += np.nan

    def _bc_mirror(self):
        """
        Mirrors qs across the boundary on either the west (left) or east (right)
        side, depending on the selections.

        This can, for example, produce a scenario in which you are observing
        a mountain range up to the range crest (or, more correctly, the halfway
        point across the mountain range).
        """
        if self._bc_west_norm == "mirror":
            i = 0
            # self.l2[i] += np.nan
            # self.l1[i] += np.nan
            self.c0[i] += 0
            self.r1[i] += self.l1_coeff_i[i]
            self.r2[i] += self.l2_coeff_i[i]
            i = 1
            # self.l2[i] += np.nan
            self.l1[i] += 0
            self.c0[i] += self.l2_coeff_i[i]
            self.r1[i] += 0
            self.r2[i] += 0

        if self._bc_east_norm == "mirror":
            i = -2
            self.l2[i] += 0
            self.l1[i] += 0
            self.c0[i] += self.r2_coeff_i[i]
            self.r1[i] += 0
            # self.r2[i] += np.nan
            i = -1
            self.l2[i] += self.r2_coeff_i[i]
            self.l1[i] += self.r1_coeff_i[i]
            self.c0[i] += 0
            # self.r1[i] += np.nan
            # self.r2[i] += np.nan

    def _bc_zero_displacement_zero_moment(self):
        """
        Simply-supported (pinned) BC: zero displacement and zero bending moment
        at the boundary.

        Enforces w = 0 at each boundary node via a Dirichlet condition: all
        off-diagonal stencil entries at that node are set to zero, leaving
        c0·w = q.  For the normal case of no load at the boundary (q = 0),
        this gives w = 0 exactly for any D profile.

        The moment condition (d²w/dx² = 0) is encoded at the first interior
        node via the odd-reflection ghost w[ghost] = -w[interior], which
        contributes -l2_coeff (west) or -r2_coeff (east) to c0 there.

        Off-diagonal entries that after np.roll land at out-of-bounds spdiags
        positions are set to np.nan (dropped automatically); those that land at
        interior matrix positions are set to 0 explicitly.
        """
        if self._bc_west_norm == "zero_displacement_zero_moment":
            # Boundary node: Dirichlet w[0] = 0
            i = 0
            self.l2[i] = np.nan    # ghost at j=-2: out-of-bounds after roll → excluded
            self.l1[i] = np.nan    # ghost at j=-1: out-of-bounds after roll → excluded
            self.c0[i] += 0        # c0·w[0] = q[0]; w[0] = 0 for q[0] = 0
            self.r1[i] = 0         # decouple: lands at interior matrix position (0, 1)
            self.r2[i] = 0         # decouple: lands at interior matrix position (0, 2)
            # First interior node: encode M = 0 via odd-reflection ghost w[-1] = -w[+1]
            i = 1
            self.l2[i] = np.nan    # ghost at j=-1: out-of-bounds after roll → excluded
            self.l1[i] += 0        # coupling to w[0] = 0; contributes 0 to solution
            self.c0[i] -= self.l2_coeff_i[i]   # ghost w[-1] = -w[1] folds into c0
            self.r1[i] += 0
            self.r2[i] += 0

        if self._bc_east_norm == "zero_displacement_zero_moment":
            # First interior node: encode M = 0 via odd-reflection ghost w[N] = -w[N-2]
            i = -2
            self.l2[i] += 0
            self.l1[i] += 0
            self.c0[i] -= self.r2_coeff_i[i]   # ghost w[N] = -w[N-2] folds into c0
            self.r1[i] += 0
            self.r2[i] = np.nan    # ghost at j=N: out-of-bounds after roll → excluded
            # Boundary node: Dirichlet w[N-1] = 0
            i = -1
            self.l2[i] = 0         # decouple: lands at interior matrix position (N-1, N-3)
            self.l1[i] = 0         # decouple: lands at interior matrix position (N-1, N-2)
            self.c0[i] += 0        # c0·w[N-1] = q[N-1]; w[N-1] = 0 for q[N-1] = 0
            self.r1[i] = np.nan    # ghost at j=N: out-of-bounds after roll → excluded
            self.r2[i] = np.nan    # ghost at j=N+1: out-of-bounds after roll → excluded

    def _apply_bc_rhs_inhomogeneous(self):
        """Compute the RHS correction vector for prescribed (nonzero) BC values.

        Called after _apply_bc_flexure() and before build_diagonals().  For
        dict-style BCs the ghost nodes used to eliminate out-of-bounds
        unknowns carry prescribed moment, shear, displacement, or slope
        values.  The constant parts of those ghost expressions move to the
        right-hand side of the linear system as a correction vector added to
        -qs before the sparse solve.
        """
        self._bc_rhs_correction = np.zeros(self.nx)
        dx = self.dx
        D_west = self.D[1]    # padded D array: index 1 = west boundary node
        D_east = self.D[-2]   # padded D array: index -2 = east boundary node

        # --- WEST ---
        if self._bc_west_values is not None:
            bv = self._bc_west_values
            if "displacement" in bv:
                w0     = bv["displacement"]
                theta0 = bv["slope"]
                # Boundary row is decoupled: c0[0]*w[0] = rhs[0] → enforce w[0]=w0
                self._bc_rhs_correction[0] = self.c0[0] * w0
                # First interior: ghost w[-1] = w[1] - 2*dx*theta0 (central diff)
                # The even-reflection term l2*w[-1] leaves a constant l2*(-2dx*theta0)
                self._bc_rhs_correction[1] = (
                    self.l2_coeff_i[1] * 2.0 * dx * theta0
                )
            else:  # "moment" / "shear"
                M0 = bv["moment"]
                V0 = bv["shear"]
                # RHS corrections from non-zero ghost constants.
                # West ghost: w[-2] = 2w[-1] - 2w[1] + w[2] - 2V₀dx³/D + 2M₀dx²/D
                # The V₀ term is subtracted in the ghost → adds to RHS.
                self._bc_rhs_correction[0] = (
                    self.l2_coeff_i[0] * (2.0*V0*dx**3 - 2.0*M0*dx**2) / D_west
                    - self.l1_coeff_i[0] * M0 * dx**2 / D_west
                )
                self._bc_rhs_correction[1] = (
                    -self.l2_coeff_i[1] * M0 * dx**2 / D_west
                )

        # --- EAST ---
        if self._bc_east_values is not None:
            bv = self._bc_east_values
            if "displacement" in bv:
                w_east     = bv["displacement"]
                theta_east = bv["slope"]
                # Boundary row is decoupled: c0[-1]*w[N-1] = rhs[-1] → enforce w[N-1]=w_east
                self._bc_rhs_correction[-1] = self.c0[-1] * w_east
                # First interior: ghost w[N] = w[N-2] + 2*dx*theta_east
                self._bc_rhs_correction[-2] = (
                    self.r2_coeff_i[-2] * 2.0 * dx * theta_east
                )
            else:  # "moment" / "shear"
                M_east = bv["moment"]
                V_east = bv["shear"]
                # RHS corrections
                self._bc_rhs_correction[-1] = (
                    -self.r1_coeff_i[-1] * M_east * dx**2 / D_east
                    - self.r2_coeff_i[-1] * (2.0*V_east*dx**3 + 2.0*M_east*dx**2) / D_east
                )
                self._bc_rhs_correction[-2] += (
                    -self.r2_coeff_i[-2] * M_east * dx**2 / D_east
                )

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
        self.maxFlexuralWavelength_ncells = int(
            np.ceil(self.max_flexural_wavelength / self.dx)
        )

    def fd_solve(self):
        """
        w = fd_solve()
          where coeff is the sparse coefficient matrix output from function
          coeff_matrix and qs is the array of loads (stresses)

        Sparse solver for one-dimensional flexure of an elastic plate
        """

        if self.debug:
            print("qs", self.qs.shape)
            print("Te", self.te.shape)
            self.calc_max_flexural_wavelength()
            print("maxFlexuralWavelength_ncells', self.maxFlexuralWavelength_ncells")

        if self.solver == "direct":
            if self.debug:
                print("Using direct solution with UMFpack")
        else:
            if not self.quiet:
                print("Solution type not understood:")
                print("Defaulting to direct solution with UMFpack")

        if self.cache_factorization not in (False, True, "no_check"):
            raise ValueError(
                f"cache_factorization must be False, True, or 'no_check'; "
                f"got {self.cache_factorization!r}"
            )

        # qs negative so bends down with positive load, bends up with negative load
        # (i.e. material removed)
        rhs = -self.qs + self._bc_rhs_correction
        if self.cache_factorization is False:
            self.w = spsolve(self.coeff_matrix, rhs, use_umfpack=True)
        elif self.cache_factorization == "no_check":
            if self._lu is None:
                self._lu = factorized(self.coeff_matrix)
                self.coeff_matrix = None  # _lu is sole owner; _solve_fd uses _lu as rebuild-skip signal
            self.w = self._lu(rhs)
        else:  # True: hash-validated cache
            h = _matrix_hash(self.coeff_matrix)
            if self._lu is None or h != self._lu_matrix_hash:
                self._lu = factorized(self.coeff_matrix)
                self._lu_matrix_hash = h
            self.w = self._lu(rhs)

        if self.debug:
            print("w.shape:")
            print(self.w.shape)
            print("w:")
            print(self.w)
