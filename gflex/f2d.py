"""
This file is part of gFlex.
gFlex computes lithospheric flexural isostasy with heterogeneous rigidity
Copyright (C) 2010-2020 Andrew D. Wickert

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

import numpy as np
import scipy
from scipy.special import kei
from scipy.sparse.linalg import isolve, spsolve

from gflex.base import Flexure


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
    """
    drho = rho_m - rho_fill
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
        Width of the padding ring in grid cells.  Use
        :func:`recommended_pad_width` to obtain a suitable value.
    Te_out : float, optional
        Te value at the outer edge of the padding ring.
        Defaults to ``Te.mean()``.

    Returns
    -------
    Te_padded : (M + 2*pad_width, N + 2*pad_width) array
        Padded elastic thickness with a smooth linear taper.

    Examples
    --------
    >>> import numpy as np
    >>> from gflex import F2D, smooth_pad_Te, recommended_pad_width
    >>> p = recommended_pad_width(Te, dx=5000.)   # Te is your (M, N) array
    >>> Te_pad = smooth_pad_Te(Te, p)
    >>> qs_pad = np.pad(qs, p, mode='constant')
    >>> flex = F2D()
    >>> flex.Te = Te_pad
    >>> flex.qs = qs_pad
    >>> # ... set other parameters and run ...
    >>> w_inner = flex.w[p:-p, p:-p]              # trim padding from output
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


def pad_domain(Te, qs, dx, dy=None, n_wavelengths=1.0, Te_out=None,
               E=65e9, nu=0.25, rho_m=3300.0, rho_fill=0.0, g=9.8):
    """
    Pad both the elastic thickness and surface load arrays for use with F2D.

    Combines :func:`recommended_pad_width` and :func:`smooth_pad_Te` into a
    single call, and zero-pads *qs* to match.  The returned pad width *p* can
    be used to trim the deflection output after the run::

        w_inner = flex.w[p:-p, p:-p]

    Parameters
    ----------
    Te : (M, N) array
        Elastic thickness [m] for the inner domain.
    qs : (M, N) array
        Surface load [Pa] for the inner domain.
    dx : float
        Grid cell size in the x-direction [m].
    dy : float, optional
        Grid cell size in the y-direction [m].  Defaults to *dx*.
    n_wavelengths : float, optional
        Padding width expressed as a number of flexural wavelengths.
        Default 1.0; use 0.5 for a less conservative (narrower) padding.
    Te_out : float, optional
        Te value at the outer edge of the padding ring.
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
    Te_padded : (M + 2p, N + 2p) array
        Smoothly tapered elastic thickness.
    qs_padded : (M + 2p, N + 2p) array
        Surface load zero-padded to match.
    p : int
        Pad width in grid cells (same on all four sides).

    Examples
    --------
    >>> import numpy as np
    >>> from gflex import F2D, pad_domain
    >>> Te_pad, qs_pad, p = pad_domain(Te, qs, dx=5000.)
    >>> flex = F2D()
    >>> flex.Te = Te_pad
    >>> flex.qs = qs_pad
    >>> # ... set other parameters and run ...
    >>> w_inner = flex.w[p:-p, p:-p]
    """
    if dy is None:
        dy = dx
    p = recommended_pad_width(
        Te, min(dx, dy), E=E, nu=nu, rho_m=rho_m, rho_fill=rho_fill,
        g=g, n_wavelengths=n_wavelengths,
    )
    Te_padded = smooth_pad_Te(Te, p, Te_out=Te_out)
    qs_padded = np.pad(qs, p, mode="constant")
    return Te_padded, qs_padded, p


# class F2D inherits Flexure and overrides __init__ therefore setting up the same
# three parameters as class Isostasy; and it then sets up more parameters specific
# to its own type of simulation.
class F2D(Flexure):
    def initialize(self, filename=None):
        self.dimension = 2  # Set it here in case it wasn't set for selection before
        super().initialize()
        if self.Verbose:
            print("F2D initialized")

    def run(self):
        self.bc_check()
        self.solver_start_time = time.time()

        if self.Method == "FD":
            # Finite difference
            super().FD()
            self.method_func = self.FD
        elif self.Method == "FFT":
            # Fast Fourier transform
            super().FFT()
            self.method_func = self.FFT
        elif self.Method == "SAS":
            # Superposition of analytical solutions
            super().SAS()
            self.method_func = self.SAS
        elif self.Method == "SAS_NG":
            # Superposition of analytical solutions,
            # nonuniform points (no grid)
            super().SAS_NG()
            self.method_func = self.SAS_NG
        else:
            sys.exit('Error: method must be "FD", "FFT", "SAS", or "SAS_NG"')

        if self.Verbose:
            print("F2D run")
        self.method_func()

        self.time_to_solve = time.time() - self.solver_start_time
        if not self.Quiet:
            print("Time to solve [s]:", self.time_to_solve)

    def finalize(self):
        # If elastic thickness has been padded, return it to its original
        # value, so this is not messed up for repeat operations in a
        # model-coupling exercise
        with contextlib.suppress(AttributeError):
            self.Te = self.Te_unpadded
        if self.Verbose:
            print("F2D finalized")
        super().finalize()

    ########################################
    ## FUNCTIONS FOR EACH SOLUTION METHOD ##
    ########################################

    def FD(self):
        # Only generate coefficient matrix if it is not already provided
        if self.coeff_matrix is not None:
            pass
        else:
            self.elasprep()
            self.BC_selector_and_coeff_matrix_creator()
        self.fd_solve()

    def FFT(self):
        sys.exit("The fast Fourier transform solution method is not yet implemented.")

    def SAS(self):
        self.spatialDomainVarsSAS()
        self.spatialDomainGridded()

    def SAS_NG(self):
        self.spatialDomainVarsSAS()
        self.spatialDomainNoGrid()

    ######################################
    ## FUNCTIONS TO SOLVE THE EQUATIONS ##
    ######################################

    ## SPATIAL DOMAIN SUPERPOSITION OF ANALYTICAL SOLUTIONS
    #########################################################

    # SETUP

    def spatialDomainVarsSAS(self):
        # Check Te:
        # * If scalar, okay.
        # * If grid, convert to scalar if a singular value
        # * Else, throw an error.
        if np.isscalar(self.Te):
            pass
        elif np.all(self.Te == np.mean(self.Te)):
            self.Te = np.mean(self.Te)
        else:
            sys.exit(
                "\nINPUT VARIABLE TYPE INCONSISTENT WITH SOLUTION TYPE.\n"
                "The analytical solution requires a scalar Te.\n"
                "(gFlex is smart enough to make this out of a uniform\n"
                "array, but won't know what value you want with a spatially\n"
                "varying array! Try finite difference instead in this case?\n"
                "EXITING."
            )

        self.D = self.E * self.Te**3 / (12 * (1 - self.nu**2))  # Flexural rigidity
        self.alpha = (self.D / (self.drho * self.g)) ** 0.25  # 2D flexural parameter
        self.coeff = self.alpha**2 / (2 * np.pi * self.D)

    # GRIDDED

    def spatialDomainGridded(self):
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

        # Now compute the deflections
        self.w = np.zeros((self.ny, self.nx))  # Deflection array
        for i in range(self.nx):
            for j in range(self.ny):
                # Loop over locations that have loads, and sum
                if self.qs[j, i]:
                    # Solve by summing portions of "biggrid" while moving origin
                    # to location of current cell
                    # Load must be multiplied by grid cell size
                    self.w += (
                        self.qs[j, i]
                        * self.dx
                        * self.dy
                        * biggrid[
                            self.ny - j : 2 * self.ny - j, self.nx - i : 2 * self.nx - i
                        ]
                    )
            # No need to return: w already belongs to "self"

    # NO GRID

    def spatialDomainNoGrid(self):
        self.w = np.zeros(self.xw.shape)
        if self.Debug:
            print("w = ")
            print(self.w.shape)

        if self.latlon:
            for i in range(len(self.x)):
                # More efficient if we have created some 0-load points
                # (e.g., for where we want output)
                if self.q[i] != 0:
                    # Create array of distances from point of load
                    r = self.greatCircleDistance(
                        lat1=self.y[i],
                        long1=self.x[i],
                        lat2=self.yw,
                        long2=self.xw,
                        radius=self.PlanetaryRadius,
                    )
                    self.w += self.q[i] * self.coeff * kei(r / self.alpha)
                    # Compute and sum deflection
                    self.w += self.q[i] * self.coeff * kei(r / self.alpha)
        else:
            for i in range(len(self.x)):
                if self.q[i] != 0:
                    r = ((self.xw - self.x[i]) ** 2 + (self.yw - self.y[i]) ** 2) ** 0.5
                    self.w += self.q[i] * self.coeff * kei(r / self.alpha)

    ## FINITE DIFFERENCE
    ######################

    def elasprep(self):
        """
        dx4, dy4, dx2dy2, D = elasprep(dx,dy,Te,E=1E11,nu=0.25)

        Defines the variables that are required to create the 2D finite
        difference solution coefficient matrix
        """

        if self.Method != "SAS_NG":
            self.dx4 = self.dx**4
            self.dy4 = self.dy**4
            self.dx2dy2 = self.dx**2 * self.dy**2
        self.D = self.E * self.Te**3 / (12 * (1 - self.nu**2))

    def BC_selector_and_coeff_matrix_creator(self):
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
        if self.Verbose:
            print("Boundary condition, West:", self.BC_W, type(self.BC_W))
            print("Boundary condition, East:", self.BC_E, type(self.BC_E))
            print("Boundary condition, North:", self.BC_N, type(self.BC_N))
            print("Boundary condition, South:", self.BC_S, type(self.BC_S))

        # First, set flexural rigidity boundary conditions to flesh out this padded
        # array
        self.BC_Rigidity()

        # Second, build the coefficient arrays -- with the rigidity b.c.'s
        self.get_coeff_values()

        # Third, apply boundary conditions to the coeff_arrays to create the
        # flexural solution
        self.BC_Flexure()

        # Fourth, construct the sparse diagonal array
        self.build_diagonals()

        # Finally, compute the total time this process took
        self.coeff_creation_time = time.time() - self.coeff_start_time
        if not self.Quiet:
            print(
                "Time to construct coefficient (operator) array [s]:",
                self.coeff_creation_time,
            )

    def BC_Rigidity(self):
        """
        Utility function to help implement boundary conditions by specifying
        them for and applying them to the elastic thickness grid
        """

        #########################################
        # FLEXURAL RIGIDITY BOUNDARY CONDITIONS #
        #########################################
        # West
        if self.BC_W == "Periodic":
            self.BC_Rigidity_W = "periodic"
        elif (
            self.BC_W
            == np.array(["0Displacement0Slope", "0Moment0Shear", "0Slope0Shear"])
        ).any():
            self.BC_Rigidity_W = "0 curvature"
        elif self.BC_W == "Mirror":
            self.BC_Rigidity_W = "mirror symmetry"
        else:
            sys.exit("Invalid Te B.C. case")
        # East
        if self.BC_E == "Periodic":
            self.BC_Rigidity_E = "periodic"
        elif (
            self.BC_E
            == np.array(["0Displacement0Slope", "0Moment0Shear", "0Slope0Shear"])
        ).any():
            self.BC_Rigidity_E = "0 curvature"
        elif self.BC_E == "Mirror":
            self.BC_Rigidity_E = "mirror symmetry"
        else:
            sys.exit("Invalid Te B.C. case")
        # North
        if self.BC_N == "Periodic":
            self.BC_Rigidity_N = "periodic"
        elif (
            self.BC_N
            == np.array(["0Displacement0Slope", "0Moment0Shear", "0Slope0Shear"])
        ).any():
            self.BC_Rigidity_N = "0 curvature"
        elif self.BC_N == "Mirror":
            self.BC_Rigidity_N = "mirror symmetry"
        else:
            sys.exit("Invalid Te B.C. case")
        # South
        if self.BC_S == "Periodic":
            self.BC_Rigidity_S = "periodic"
        elif (
            self.BC_S
            == np.array(["0Displacement0Slope", "0Moment0Shear", "0Slope0Shear"])
        ).any():
            self.BC_Rigidity_S = "0 curvature"
        elif self.BC_S == "Mirror":
            self.BC_Rigidity_S = "mirror symmetry"
        else:
            sys.exit("Invalid Te B.C. case")

        #############
        # PAD ARRAY #
        #############
        if np.isscalar(self.Te):
            self.D *= np.ones(self.qs.shape)  # And leave Te as a scalar for checks
        else:
            self.Te_unpadded = self.Te.copy()
            self.Te = np.hstack(
                (
                    np.nan * np.zeros((self.Te.shape[0], 1)),
                    self.Te,
                    np.nan * np.zeros((self.Te.shape[0], 1)),
                )
            )
            self.Te = np.vstack(
                (
                    np.nan * np.zeros(self.Te.shape[1]),
                    self.Te,
                    np.nan * np.zeros(self.Te.shape[1]),
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
        if self.BC_Rigidity_W == "0 curvature":
            self.D[:, 0] = 2 * self.D[:, 1] - self.D[:, 2]
        if self.BC_Rigidity_E == "0 curvature":
            self.D[:, -1] = 2 * self.D[:, -2] - self.D[:, -3]
        if self.BC_Rigidity_N == "0 curvature":
            self.D[0, :] = 2 * self.D[1, :] - self.D[2, :]
        if self.BC_Rigidity_S == "0 curvature":
            self.D[-1, :] = 2 * self.D[-2, :] - self.D[-3, :]

        if self.BC_Rigidity_W == "mirror symmetry":
            self.D[:, 0] = self.D[:, 2]
        if self.BC_Rigidity_E == "mirror symmetry":
            self.D[:, -1] = self.D[:, -3]
        if self.BC_Rigidity_N == "mirror symmetry":
            self.D[0, :] = self.D[
                2, :
            ]  # Yes, will work on corners -- double-reflection
        if self.BC_Rigidity_S == "mirror symmetry":
            self.D[-1, :] = self.D[-3, :]

        if self.BC_Rigidity_W == "periodic":
            self.D[:, 0] = self.D[:, -2]
        if self.BC_Rigidity_E == "periodic":
            self.D[:, -1] = self.D[:, -3]
        if self.BC_Rigidity_N == "periodic":
            self.D[0, :] = self.D[-2, :]
        if self.BC_Rigidity_S == "periodic":
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

        if np.isscalar(self.Te):
            # So much simpler with constant D! And symmetrical stencil
            self.cj2i0 = D / dy4
            self.cj1i_1 = 2 * D / dx2dy2 + 2 * self.sigma_xy * self.Te
            self.cj1i0 = -4 * D / dy4 - 4 * D / dx2dy2 - self.sigma_yy * self.Te
            self.cj1i1 = 2 * D / dx2dy2 - 2 * self.sigma_xy * self.Te
            self.cj0i_2 = D / dx4
            self.cj0i_1 = -4 * D / dx4 - 4 * D / dx2dy2 - self.sigma_xx * self.Te
            self.cj0i0 = (
                6 * D / dx4
                + 6 * D / dy4
                + 8 * D / dx2dy2
                + drho * g
                + 2 * self.sigma_xx * self.Te
                + 2 * self.sigma_yy * self.Te
            )
            self.cj0i1 = (
                -4 * D / dx4 - 4 * D / dx2dy2 - self.sigma_xx * self.Te
            )  # Symmetry
            self.cj0i2 = D / dx4  # Symmetry
            self.cj_1i_1 = 2 * D / dx2dy2 - 2 * self.sigma_xy * self.Te  # Symmetry
            self.cj_1i0 = -4 * D / dy4 - 4 * D / dx2dy2  # Symmetry
            self.cj_1i1 = 2 * D / dx2dy2 + 2 * self.sigma_xy * self.Te  # Symmetry
            self.cj_2i0 = D / dy4  # Symmetry
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

        elif isinstance(self.Te, np.ndarray):
            #######################################################
            # GENERATE COEFFICIENT VALUES FOR EACH SOLUTION TYPE. #
            #    "vWC1994" IS THE BEST: LOOSEST ASSUMPTIONS.      #
            #        OTHERS HERE LARGELY FOR COMPARISON           #
            #######################################################

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

            if self.PlateSolutionType == "vWC1994":
                # van Wees and Cloetingh (1994) solution, re-discretized by me
                # using a central difference approx. to 2nd order precision
                # NEW STENCIL
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

            elif self.PlateSolutionType == "G2009":
                # STENCIL FROM GOVERS ET AL. 2009 -- first-order differences
                # x is j and y is i b/c matrix row/column notation
                # Note that this breaks down with b.c.'s that place too much control
                # on the solution -- harmonic wavetrains
                # x = -2, y = 0
                self.cj_2i0_coeff_ij = D_10 / dx4
                # x = -1, y = -1
                self.cj_1i_1_coeff_ij = (D_10 + D0_1) / dx2dy2
                # x = -1, y = 0
                self.cj_1i0_coeff_ij = -2.0 * (
                    (D0_1 + D00) / dx2dy2 + (D00 + D_10) / dx4
                )
                # x = -1, y = 1
                self.cj_1i1_coeff_ij = (D_10 + D01) / dx2dy2
                # x = 0, y = -2
                self.cj0i_2_coeff_ij = D0_1 / dy4
                # x = 0, y = -1
                self.cj0i_1_coeff_ij = -2.0 * (
                    (D0_1 + D00) / dx2dy2 + (D00 + D0_1) / dy4
                )
                # x = 0, y = 0
                self.cj0i0_coeff_ij = (
                    (D10 + 4.0 * D00 + D_10) / dx4
                    + (D01 + 4.0 * D00 + D0_1) / dy4
                    + (8.0 * D00 / dx2dy2)
                    + drho * g
                )
                # x = 0, y = 1
                self.cj0i1_coeff_ij = -2.0 * ((D01 + D00) / dy4 + (D00 + D01) / dx2dy2)
                # x = 0, y = 2
                self.cj0i2_coeff_ij = D0_1 / dy4
                # x = 1, y = -1
                self.cj1i_1_coeff_ij = (D10 + D0_1) / dx2dy2
                # x = 1, y = 0
                self.cj1i0_coeff_ij = -2.0 * ((D10 + D00) / dx4 + (D10 + D00) / dx2dy2)
                # x = 1, y = 1
                self.cj1i1_coeff_ij = (D10 + D01) / dx2dy2
                # x = 2, y = 0
                self.cj2i0_coeff_ij = D10 / dx4
            else:
                sys.exit(
                    "Not an acceptable plate solution type. Please choose from:\n"
                    + "* vWC1994\n"
                    + "* G2009\n"
                    + ""
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

    def BC_Flexure(self):
        # The next section of code is split over several functions for the 1D
        # case, but will be all in one function here, at least for now.

        # Inf for E-W to separate from nan for N-S. N-S will spill off ends
        # of array (C order, in rows), while E-W will be internal, so I will
        # later change np.inf to 0 to represent where internal boundaries
        # occur.

        #######################################################################
        # DEFINE COEFFICIENTS TO W_j-2 -- W_j+2 WITH B.C.'S APPLIED (x: W, E) #
        #######################################################################

        # Infinitiy is used to flag places where coeff values should be 0,
        # and would otherwise cause boundary condition nan's to appear in the
        # cross-derivatives: infinity is changed into 0 later.

        if self.BC_W == "Periodic":
            if self.BC_E == "Periodic":
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
                self.cj_1i1_Periodic_right[:, j] = self.cj_1i_1[:, j]
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
                sys.exit(
                    "Not physical to have one wrap-around boundary but not its pair."
                )
        elif self.BC_W == "0Displacement0Slope":
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
            self.cj1i_1[:, j] += 0
            self.cj1i0[:, j] += 0
            self.cj1i1[:, j] += 0
            self.cj2i0[:, j] += 0
            j = 1
            self.cj_2i0[:, j] += np.inf
            self.cj_1i_1[:, j] += 0
            self.cj_1i0[:, j] += 0
            self.cj_1i1[:, j] += 0
            self.cj0i_2[:, j] += 0
            self.cj0i_1[:, j] += 0
            self.cj0i0[:, j] += 0
            self.cj0i1[:, j] += 0
            self.cj0i2[:, j] += 0
            self.cj1i_1[:, j] += 0
            self.cj1i0[:, j] += 0
            self.cj1i1[:, j] += 0
            self.cj2i0[:, j] += 0
        elif self.BC_W == "0Moment0Shear":
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
            j = 1
            self.cj_2i0[:, j] += np.inf
            self.cj_1i_1[:, j] += 0
            self.cj_1i0[:, j] += 2 * self.cj_2i0_coeff_ij[:, j]
            self.cj_1i1[:, j] += 0
            self.cj0i_2[:, j] += 0
            self.cj0i_1[:, j] += 0
            self.cj0i0[:, j] += 0
            self.cj0i1[:, j] += 0
            self.cj0i2[:, j] += 0
            self.cj1i_1[:, j] += 0
            self.cj1i0[:, j] += -2 * self.cj_2i0_coeff_ij[:, j]
            self.cj1i1[:, j] += 0
            self.cj2i0[:, j] += self.cj_2i0_coeff_ij[:, j]
        elif self.BC_W == "0Slope0Shear":
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
            self.cj1i1[:, j] += self.cj_1i1_coeff_ij[:, j]  # Interference
            self.cj2i0[:, j] += self.cj_2i0_coeff_ij[:, j]
            j = 1
            self.cj_2i0[:, j] += np.inf
            self.cj_1i_1[:, j] += 0
            self.cj_1i0[:, j] += 0
            self.cj_1i1[:, j] += 0
            self.cj0i_2[:, j] += 0
            self.cj0i_1[:, j] += 0
            self.cj0i0[:, j] += 0
            self.cj0i1[:, j] += 0
            self.cj0i2[:, j] += 0
            self.cj1i_1[:, j] += 0
            self.cj1i0[:, j] += 0
            self.cj1i1[:, j] += 0
            self.cj2i0[:, j] += self.cj_2i0_coeff_ij[:, j]
        elif self.BC_W == "Mirror":
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
        else:
            # Possibly redundant safeguard
            sys.exit("Invalid boundary condition")

        if self.BC_E == "Periodic":
            # See more extensive comments above (BC_W)

            if self.BC_W == "Periodic":
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
                sys.exit(
                    "Not physical to have one wrap-around boundary but not its pair."
                )

        elif self.BC_E == "0Displacement0Slope":
            j = -1
            self.cj_2i0[:, j] += 0
            self.cj_1i_1[:, j] += 0
            self.cj_1i0[:, j] += 0
            self.cj_1i1[:, j] += 0
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
            self.cj0i0[:, j] += 0
            self.cj0i1[:, j] += 0
            self.cj0i2[:, j] += 0
            self.cj1i_1[:, j] += 0
            self.cj1i0[:, j] += 0
            self.cj1i1[:, j] += 0
            self.cj2i0[:, j] += np.inf
        elif self.BC_E == "0Moment0Shear":
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
            j = -2
            self.cj_2i0[:, j] += self.cj2i0_coeff_ij[:, j]
            self.cj_1i_1[:, j] += 0
            self.cj_1i0[:, j] += -2 * self.cj2i0_coeff_ij[:, j]
            self.cj_1i1[:, j] += 0
            self.cj0i_2[:, j] += 0
            self.cj0i_1[:, j] += 0
            self.cj0i0[:, j] += 0
            self.cj0i1[:, j] += 0
            self.cj0i2[:, j] += 0
            self.cj1i_1[:, j] += 0
            self.cj1i0[:, j] += 2 * self.cj2i0_coeff_ij[:, j]
            self.cj1i1[:, j] += 0
            self.cj2i0[:, j] += np.inf
        elif self.BC_E == "0Slope0Shear":
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
            self.cj_2i0[:, j] += self.cj2i0_coeff_ij[:, j]
            self.cj_1i_1[:, j] += 0
            self.cj_1i0[:, j] += 0
            self.cj_1i1[:, j] += 0
            self.cj0i_2[:, j] += 0
            self.cj0i_1[:, j] += 0
            self.cj0i0[:, j] += 0
            self.cj0i1[:, j] += 0
            self.cj0i2[:, j] += 0
            self.cj1i_1[:, j] += 0
            self.cj1i0[:, j] += 0
            self.cj1i1[:, j] += 0
            self.cj2i0[:, j] += np.inf
        elif self.BC_E == "Mirror":
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
        else:
            # Possibly redundant safeguard
            sys.exit("Invalid boundary condition")

        #######################################################################
        # DEFINE COEFFICIENTS TO W_i-2 -- W_i+2 WITH B.C.'S APPLIED (y: N, S) #
        #######################################################################

        if self.BC_N == "Periodic":
            if self.BC_S == "Periodic":
                pass  # Will address the N-S (whole-matrix-involving) boundary condition
                # inclusion below, when constructing sparse matrix diagonals
            else:
                sys.exit(
                    "Not physical to have one wrap-around boundary but not its pair."
                )
        elif self.BC_N == "0Displacement0Slope":
            i = 0
            self.cj_2i0[i, :] += 0
            self.cj_1i_1[i, :][self.cj_1i_1[i, :] != np.inf] = np.nan
            self.cj_1i0[i, :] += 0
            self.cj_1i1[i, :] += 0
            self.cj0i_2[i, :] += 0  # np.nan
            self.cj0i_1[i, :] += 0  # np.nan
            self.cj0i0[i, :] += 0
            self.cj0i1[i, :] += 0
            self.cj0i2[i, :] += 0
            self.cj1i_1[i, :][self.cj1i_1[i, :] != np.inf] += 0  # np.nan
            self.cj1i0[i, :] += 0
            self.cj1i1[i, :] += 0
            self.cj2i0[i, :] += 0
            i = 1
            self.cj_2i0[i, :] += 0
            self.cj_1i_1[i, :] += 0
            self.cj_1i0[i, :] += 0
            self.cj_1i1[i, :] += 0
            self.cj0i_2[i, :] += 0  # np.nan
            self.cj0i_1[i, :] += 0
            self.cj0i0[i, :] += 0
            self.cj0i1[i, :] += 0
            self.cj0i2[i, :] += 0
            self.cj1i_1[i, :] += 0
            self.cj1i0[i, :] += 0
            self.cj1i1[i, :] += 0
            self.cj2i0[i, :] += 0
        elif self.BC_N == "0Moment0Shear":
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
            i = 1
            self.cj_2i0[i, :] += 0
            self.cj_1i_1[i, :] += 0
            self.cj_1i0[i, :] += 0
            self.cj_1i1[i, :] += 0
            self.cj0i_2[i, :] += 0  # np.nan
            self.cj0i_1[i, :] += 2 * self.cj0i_2_coeff_ij[i, :]
            self.cj0i0[i, :] += 0
            self.cj0i1[i, :] += -2 * self.cj0i_2_coeff_ij[i, :]
            self.cj0i2[i, :] += self.cj0i_2_coeff_ij[i, :]
            self.cj1i_1[i, :] += 0
            self.cj1i0[i, :] += 0
            self.cj1i1[i, :] += 0
            self.cj2i0[i, :] += 0
        elif self.BC_N == "0Slope0Shear":
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
            self.cj0i0[i, :] += 0
            self.cj0i1[i, :] += 0
            self.cj0i2[i, :] += self.cj0i_2_coeff_ij[i, :]
            self.cj1i_1[i, :] += 0
            self.cj1i0[i, :] += 0
            self.cj1i1[i, :] += 0
            self.cj2i0[i, :] += 0
        elif self.BC_N == "Mirror":
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
        else:
            # Possibly redundant safeguard
            sys.exit("Invalid boundary condition")

        if self.BC_S == "Periodic":
            if self.BC_N == "Periodic":
                pass  # Will address the N-S (whole-matrix-involving) boundary condition
                # inclusion below, when constructing sparse matrix diagonals
            else:
                sys.exit(
                    "Not physical to have one wrap-around boundary but not its pair."
                )
        elif self.BC_S == "0Displacement0Slope":
            i = -2
            self.cj_2i0[i, :] += 0
            self.cj_1i_1[i, :] += 0
            self.cj_1i0[i, :] += 0
            self.cj_1i1[i, :] += 0
            self.cj0i_2[i, :] += 0
            self.cj0i_1[i, :] += 0
            self.cj0i0[i, :] += 0
            self.cj0i1[i, :] += 0
            self.cj0i2[i, :] += 0  # np.nan
            self.cj1i1[i, :] += 0
            self.cj1i0[i, :] += 0
            self.cj1i_1[i, :] += 0
            self.cj2i0[i, :] += 0
            i = -1
            self.cj_2i0[i, :] += 0
            self.cj_1i_1[i, :] += 0
            self.cj_1i0[i, :] += 0
            self.cj_1i1[i, :][self.cj_1i1[i, :] != np.inf] += 0  # np.nan
            self.cj0i_2[i, :] += 0
            self.cj0i_1[i, :] += 0
            self.cj0i0[i, :] += 0
            self.cj0i1[i, :] += 0  # np.nan
            self.cj0i2[i, :] += 0  # np.nan
            self.cj1i_1[i, :] += 0
            self.cj1i0[i, :] += 0
            self.cj1i1[i, :][self.cj1i1[i, :] != np.inf] += 0  # np.nan
            self.cj2i0[i, :] += 0
        elif self.BC_S == "0Moment0Shear":
            i = -2
            self.cj_2i0[i, :] += 0
            self.cj_1i_1[i, :] += 0
            self.cj_1i0[i, :] += 0
            self.cj_1i1[i, :] += 0
            self.cj0i_2[i, :] += self.cj0i2_coeff_ij[i, :]
            self.cj0i_1[i, :] += -2 * self.cj0i2_coeff_ij[i, :]
            self.cj0i0[i, :] += 0
            self.cj0i1[i, :] += 2 * self.cj0i2_coeff_ij[i, :]
            self.cj0i2[i, :] += 0  # np.nan
            self.cj1i1[i, :] += 0
            self.cj1i0[i, :] += 0
            self.cj1i_1[i, :] += 0
            self.cj2i0[i, :] += 0
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
        elif self.BC_S == "0Slope0Shear":
            i = -2
            self.cj_2i0[i, :] += 0
            self.cj_1i_1[i, :] += 0
            self.cj_1i0[i, :] += 0
            self.cj_1i1[i, :] += 0
            self.cj0i_2[i, :] += self.cj0i2_coeff_ij[i, :]
            self.cj0i_1[i, :] += 0
            self.cj0i0[i, :] += 0
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
        elif self.BC_S == "Mirror":
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
        else:
            # Possibly redundant safeguard
            sys.exit("Invalid boundary condition")

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
        if self.BC_N == "0Moment0Shear" and self.BC_W == "0Moment0Shear":
            self.cj0i0[0, 0] += 2 * self.cj_1i_1_coeff_ij[0, 0]
            self.cj1i1[0, 0] -= self.cj_1i_1_coeff_ij[0, 0]
        if self.BC_N == "0Moment0Shear" and self.BC_E == "0Moment0Shear":
            self.cj0i0[0, -1] += 2 * self.cj_1i_1_coeff_ij[0, -1]
            self.cj_1i1[0, -1] -= self.cj1i_1_coeff_ij[0, -1]
        if self.BC_S == "0Moment0Shear" and self.BC_W == "0Moment0Shear":
            self.cj0i0[-1, 0] += 2 * self.cj_1i_1_coeff_ij[-1, 0]
            self.cj1i_1[-1, 0] -= self.cj_1i1_coeff_ij[-1, 0]
        if self.BC_S == "0Moment0Shear" and self.BC_E == "0Moment0Shear":
            self.cj0i0[-1, -1] += 2 * self.cj_1i_1_coeff_ij[-1, -1]
            self.cj_1i_1[-1, -1] -= self.cj1i1_coeff_ij[-1, -1]

        ############
        # PERIODIC #
        ############

        # I think that nothing will be needed here.
        # Periodic should just take care of all repeating in all directions by
        # its very nature. (I.e. it is embedded in the sparse array structure
        # of diagonals)

        ################
        # COMBINATIONS #
        ################

        ##############################
        # 0SLOPE0SHEAR AND/OR MIRROR #
        ##############################
        # (both end up being the same)
        if (self.BC_N == "0Slope0Shear" or self.BC_N == "Mirror") and (
            self.BC_W == "0Slope0Shear" or self.BC_W == "Mirror"
        ):
            self.cj1i1[0, 0] += self.cj_1i_1_coeff_ij[0, 0]
        if (self.BC_N == "0Slope0Shear" or self.BC_N == "Mirror") and (
            self.BC_E == "0Slope0Shear" or self.BC_E == "Mirror"
        ):
            self.cj_1i1[0, -1] += self.cj1i_1_coeff_ij[0, -1]
        if (self.BC_S == "0Slope0Shear" or self.BC_S == "Mirror") and (
            self.BC_W == "0Slope0Shear" or self.BC_W == "Mirror"
        ):
            self.cj1i_1[-1, 0] += self.cj_1i1_coeff_ij[-1, 0]
        if (self.BC_S == "0Slope0Shear" or self.BC_S == "Mirror") and (
            self.BC_E == "0Slope0Shear" or self.BC_E == "Mirror"
        ):
            self.cj_1i_1[-1, -1] += self.cj1i1_coeff_ij[-1, -1]

        ################################
        # 0MOMENT0SHEAR - AND - MIRROR #
        ################################
        # How do multiple types of b.c.'s interfere
        # 0Moment0Shear must determine corner conditions in order to be mirrored
        # by the "mirror" b.c.
        if (self.BC_N == "Mirror" and self.BC_W == "0Moment0Shear") or (
            self.BC_W == "Mirror" and self.BC_N == "0Moment0Shear"
        ):
            self.cj0i0[0, 0] += 2 * self.cj_1i_1_coeff_ij[0, 0]
            self.cj1i1[0, 0] -= self.cj_1i_1_coeff_ij[0, 0]
        if (self.BC_N == "Mirror" and self.BC_E == "0Moment0Shear") or (
            self.BC_E == "Mirror" and self.BC_N == "0Moment0Shear"
        ):
            self.cj0i0[0, -1] += 2 * self.cj_1i_1_coeff_ij[0, -1]
            self.cj1i1[0, -1] -= self.cj_1i_1_coeff_ij[0, -1]
        if (self.BC_S == "Mirror" and self.BC_W == "0Moment0Shear") or (
            self.BC_W == "Mirror" and self.BC_S == "0Moment0Shear"
        ):
            self.cj0i0[-1, 0] += 2 * self.cj_1i_1_coeff_ij[-1, 0]
            self.cj1i_1[-1, 0] -= self.cj_1i1_coeff_ij[-1, 0]
        if (self.BC_S == "Mirror" and self.BC_E == "0Moment0Shear") or (
            self.BC_E == "Mirror" and self.BC_S == "0Moment0Shear"
        ):
            self.cj0i0[-1, -1] += 2 * self.cj_1i_1_coeff_ij[-1, -1]
            self.cj_1i_1[-1, -1] -= self.cj1i1_coeff_ij[-1, -1]

        ######################################
        # 0MOMENT0SHEAR - AND - 0SLOPE0SHEAR #
        ######################################
        # Just use 0Moment0Shear-style b.c.'s at corners: letting this dominate
        # because it seems to be the more geologically likely b.c.
        if (self.BC_N == "0Slope0Shear" and self.BC_W == "0Moment0Shear") or (
            self.BC_W == "0Slope0Shear" and self.BC_N == "0Moment0Shear"
        ):
            self.cj0i0[0, 0] += 2 * self.cj_1i_1_coeff_ij[0, 0]
            self.cj1i1[0, 0] -= self.cj_1i_1_coeff_ij[0, 0]
        if (self.BC_N == "0Slope0Shear" and self.BC_E == "0Moment0Shear") or (
            self.BC_E == "0Slope0Shear" and self.BC_N == "0Moment0Shear"
        ):
            self.cj0i0[0, -1] += 2 * self.cj_1i_1_coeff_ij[0, -1]
            self.cj1i1[0, -1] -= self.cj_1i_1_coeff_ij[0, -1]
        if (self.BC_S == "0Slope0Shear" and self.BC_W == "0Moment0Shear") or (
            self.BC_W == "0Slope0Shear" and self.BC_S == "0Moment0Shear"
        ):
            self.cj0i0[-1, 0] += 2 * self.cj_1i_1_coeff_ij[-1, 0]
            self.cj1i_1[-1, 0] -= self.cj_1i1_coeff_ij[-1, 0]
        if (self.BC_S == "0Slope0Shear" and self.BC_E == "0Moment0Shear") or (
            self.BC_E == "0Slope0Shear" and self.BC_S == "0Moment0Shear"
        ):
            self.cj0i0[-1, -1] += 2 * self.cj_1i_1_coeff_ij[-1, -1]
            self.cj_1i_1[-1, -1] -= self.cj1i1_coeff_ij[-1, -1]
        # What about 0Moment0SHear on N/S part?

        ##############################
        # PERIODIC B.C.'S AND OTHERS #
        ##############################

        # The Periodic boundary natively continues the other boundary conditions
        # Nothing to be done here.

    def build_diagonals(self):
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
            self.BC_N == "Periodic"
            and self.BC_S == "Periodic"
            and self.BC_W == "Periodic"
            and self.BC_E == "Periodic"
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
                # Periodic b.c. tridiag
                self.nx - self.ny * self.nx - 1,
                self.nx - self.ny * self.nx,
                self.nx - self.ny * self.nx + 1,
                # New: UR corner of LL box
                2 * self.nx - self.ny * self.nx - 1,
                # Periodic b.c. single diag
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
                # Periodic b.c. single diag
                self.ny * self.nx - 2 * self.nx,
                # New: LL corner of UR box
                self.ny * self.nx - 2 * self.nx + 1,
                # Periodic b.c. tridiag
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

        elif self.BC_W == "Periodic" and self.BC_E == "Periodic":
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

        elif self.BC_N == "Periodic" and self.BC_S == "Periodic":
            # Periodic.
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
        self.maxFlexuralWavelength = 2 * np.pi * alpha
        self.maxFlexuralWavelength_ncells_x = int(
            np.ceil(self.maxFlexuralWavelength / self.dx)
        )
        self.maxFlexuralWavelength_ncells_y = int(
            np.ceil(self.maxFlexuralWavelength / self.dy)
        )

    def fd_solve(self):
        """
        w = fd_solve()
        Sparse flexural response calculation.
        Can be performed by direct factorization with UMFpack (defuault)
        or by an iterative minimum residual technique
        These are both the fastest of the standard Scipy builtin techniques in
        their respective classes
        Requires the coefficient matrix from "2D.coeff_matrix"
        """

        if self.Debug:
            # Will fail if scalar
            with contextlib.suppress(AttributeError):
                print("self.Te", self.Te.shape)
            print("self.qs", self.qs.shape)
            self.calc_max_flexural_wavelength()
            print(
                "maxFlexuralWavelength_ncells: (x, y):",
                self.maxFlexuralWavelength_ncells_x,
                self.maxFlexuralWavelength_ncells_y,
            )

        q0vector = self.qs.reshape(-1, order="C")
        if self.Solver == "iterative" or self.Solver == "Iterative":
            if self.Debug:
                print(
                    "Using generalized minimal residual method for iterative solution"
                )
            if self.Verbose:
                print(
                    "Converging to a tolerance of",
                    self.iterative_ConvergenceTolerance,
                    "m between iterations",
                )
            wvector = isolve.lgmres(
                self.coeff_matrix, q0vector
            )  # , tol=1E-10)#,x0=woldvector)#,x0=wvector,tol=1E-15)
            wvector = wvector[0]  # Reach into tuple to get my array back
        else:
            if self.Solver == "direct" or self.Solver == "Direct":
                if self.Debug:
                    print("Using direct solution with UMFpack")
            else:
                if not self.Quiet:
                    print("Solution type not understood:")
                    print("Defaulting to direct solution with UMFpack")
            wvector = spsolve(
                self.coeff_matrix, q0vector, use_umfpack=True
            )

        # Reshape into grid
        self.w = -wvector.reshape(self.qs.shape)
        self.w_padded = self.w.copy()  # for troubleshooting

        # Time to solve used to be here
