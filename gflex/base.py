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

import configparser
import contextlib
import os
import sys
import warnings

import numpy as np
from matplotlib import pyplot as plt

from ._version import __version__


class Utility:

    """
    Generic utility functions
    """

    def configGet(
        self, vartype, category, name, optional=False, specialReturnMessage=None
    ):
        """
        Wraps a try / except and a check for self.filename around ConfigParser
        as it talks to the configuration file.
        Also, checks for existence of configuration file so this won't execute (and fail)
        when no configuration file is provided (e.g., running in coupled mode with CSDMS
        entirely with getters and setters)

        vartype can be 'float', 'str' or 'string' (str and string are the same),
        or 'int' or 'integer' (also the same).

        "Optional" determines whether or not the program will exit if the variable
        fails to load. Set it to "True" if you don't want it to exit. In this case,
        the variable will be set to "None". Otherwise, it defaults to "False".

        "specialReturnMessage" is something that you would like to add at the end
        of a failure to execute message. By default it does not print.
        """

        try:
            if vartype == "float":
                var = self.config.getfloat(category, name)
            elif vartype == "string" or vartype == "str":
                var = self.config.get(category, name)
                if var == "" and not optional:
                    # but "" is acceptable for boundary conditions
                    if name[:17] != "BoundaryCondition":
                        if not self.Quiet:
                            print(
                                "An empty input string here is not an acceptable"
                                " option."
                            )
                            print(name, "is not optional.")
                            print("Program crash likely to occur.")
            elif vartype == "integer" or vartype == "int":
                var = self.config.getint(category, name)
            elif vartype == "boolean" or vartype == "bool":
                var = self.config.getboolean(category, name)
            else:
                print(
                    "Please enter 'float', 'string' (or 'str'), 'integer' (or 'int'),"
                    " or 'boolean (or 'bool') for vartype"
                )
                sys.exit()  # Won't exit, but will lead to exception
            return var
        except:
            if optional:
                # Carry on if the variable is optional
                var = None
                if self.Verbose or self.Debug:
                    if not self.grass:
                        print("")
                        print('No value entered for optional parameter "' + name + '"')
                        print('in category "' + category + '" in configuration file.')
                        print(
                            "No action related to this optional parameter will be taken."
                        )
                        print("")
            else:
                print(
                    "Problem loading "
                    + vartype
                    + ' "'
                    + name
                    + '" in category "'
                    + category
                    + '" from configuration file.'
                )
                if specialReturnMessage:
                    print(specialReturnMessage)
                sys.exit("Exiting.")

    def readyCoeff(self):
        from scipy import sparse

        if sparse.issparse(self.coeff_matrix):
            pass  # Good type
        else:
            try:
                self.coeff_matrix = sparse.dia_matrix(self.coeff_matrix)
            except:
                sys.exit(
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
        This is experimental code that could be used in the spatialDomainNoGrid
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
                    if not self.Quiet:
                        print("dx and dy being overwritten -- supply a full grid")
            else:
                if not self.Quiet:
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
            print("Lat/lon xw and yw must be pre-set: grid will not be square")
            print("and may run into issues with poles, so to ensure the proper")
            print("output points are chosen, the end user should do this.")
            sys.exit()

    def loadFile(self, var, close_on_fail=True):
        """
        A special function to replate a variable name that is a string file path
        with the loaded file.
        var is a string on input
        output is a numpy array or a None-type object (success vs. failure)
        """
        out = None
        try:
            # First see if it is a full path or directly links from the current
            # working directory
            out = np.load(var)
        except:
            try:
                out = np.loadtxt(var)
            except:
                # Then see if it is relative to the location of the configuration file
                try:
                    out = np.load(self.inpath + var)
                except:
                    try:
                        out = np.loadtxt(self.inpath + var)
                    # If failure
                    except:
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
            print(f"Cannot find {var} file")
            print(f"{var} path = {var}")
            print("Looked relative to model python files.")
            print("Also looked relative to configuration file path,")
            print(f"  {self.inpath}")
            print("Exiting.")
            sys.exit()
        elif out is not None and self.Verbose:
            print(f"Loading {var} from {format_name}")

        return out


class Plotting:
    # Plot, if desired
    # 1D all here, 2D in functions
    # Just because there is often more code in 2D plotting functions
    # Also, yes, this portion of the code is NOT efficient or elegant in how it
    # handles functions. But it's just a simple way to visualize results
    # easily! And not too hard to improve with a bit of time. Anyway, the main
    # goal here is the visualization, not the beauty of the code : )
    def plotting(self):
        # try:
        #  self.plotChoice
        # except:
        #  self.plotChoice = None
        if self.plotChoice:
            if self.Verbose:
                print("Starting to plot " + self.plotChoice)
            if self.dimension == 1:
                if self.plotChoice == "q":
                    plt.figure(1)
                    if self.Method == "SAS_NG":
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
                elif self.plotChoice == "w":
                    plt.figure(1)
                    if self.Method == "SAS_NG":
                        plt.plot(self.xw / 1000.0, self.w, "k-")
                    else:
                        plt.plot(self.x / 1000.0, self.w, "k-")
                    plt.ylabel("Deflection [m]", fontsize=12, fontweight="bold")
                    plt.xlabel(
                        "Distance along profile [km]", fontsize=12, fontweight="bold"
                    )
                    plt.tight_layout()
                    plt.show()
                elif self.plotChoice == "both":
                    plt.figure(1, figsize=(6, 9))
                    ax = plt.subplot(212)
                    if self.Method == "SAS_NG":
                        ax.plot(self.xw / 1000.0, self.w, "k-")
                    else:
                        ax.plot(self.x / 1000.0, self.w, "k-")
                    ax.set_ylabel("Deflection [m]", fontsize=12, fontweight="bold")
                    ax.set_xlabel(
                        "Distance along profile [m]", fontsize=12, fontweight="bold"
                    )
                    plt.subplot(211)
                    plt.title("Loads and Lithospheric Deflections", fontsize=16)
                    if self.Method == "SAS_NG":
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
                elif self.plotChoice == "combo":
                    fig = plt.figure(1, figsize=(10, 6))
                    titletext = "Loads and Lithospheric Deflections"
                    ax = fig.add_subplot(1, 1, 1)
                    # Plot undeflected load
                    if self.Method == "SAS_NG":
                        if not self.Quiet:
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
                    if self.Method == "SAS_NG":
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
                    if self.Method == "SAS_NG":
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
                        self.Te
                        if self.Method == "FD":
                            if type(self.Te) is np.ndarray:
                                if (self.Te != (self.Te).mean()).any():
                                    plt.title(titletext, fontsize=16)
                                else:
                                    plt.title(
                                        titletext
                                        + ", $T_e$ = "
                                        + str((self.Te / 1000).mean())
                                        + " km",
                                        fontsize=16,
                                    )
                            else:
                                plt.title(
                                    titletext
                                    + ", $T_e$ = "
                                    + str(self.Te / 1000)
                                    + " km",
                                    fontsize=16,
                                )
                        else:
                            plt.title(
                                titletext + ", $T_e$ = " + str(self.Te / 1000) + " km",
                                fontsize=16,
                            )
                    except:
                        plt.title(titletext, fontsize=16)
                    # x and y labels
                    plt.ylabel("Loads and flexural response [m]", fontsize=16)
                    plt.xlabel("Distance along profile [km]", fontsize=16)
                    # legend -- based on lables
                    plt.legend(loc=0, numpoints=1, fancybox=True)
                    plt.tight_layout()
                    plt.show()
                else:
                    if not self.Quiet:
                        print(
                            'Incorrect plotChoice input, "'
                            + self.plotChoice
                            + '" provided.'
                        )
                        print(
                            "Possible input strings are: q, w, both, and (for 1D) combo"
                        )
                        print("Unable to produce plot.")
            elif self.dimension == 2:
                if self.plotChoice == "q":
                    fig = plt.figure(1, figsize=(8, 6))
                    if self.Method != "SAS_NG":
                        self.surfplot(
                            self.qs / (self.rho_m * self.g),
                            "Load thickness, mantle equivalent [m]",
                        )
                        plt.show()
                    else:
                        self.xyzinterp(
                            self.x,
                            self.y,
                            self.q,
                            "Load volume, mantle equivalent [m$^3$]",
                        )
                    plt.tight_layout()
                    plt.show()
                elif self.plotChoice == "w":
                    fig = plt.figure(1, figsize=(8, 6))
                    if self.Method != "SAS_NG":
                        self.surfplot(self.w, "Deflection [m]")
                        plt.show()
                    else:
                        self.xyzinterp(self.xw, self.yw, self.w, "Deflection [m]")
                    plt.tight_layout()
                    plt.show()
                elif self.plotChoice == "both":
                    plt.figure(1, figsize=(6, 9))
                    if self.Method != "SAS_NG":
                        self.twoSurfplots()
                        plt.show()
                    else:
                        plt.subplot(211)
                        self.xyzinterp(
                            self.x,
                            self.y,
                            self.q,
                            "Load volume, mantle equivalent [m$^3$]",
                        )
                        plt.subplot(212)
                        self.xyzinterp(self.xw, self.yw, self.w, "Deflection [m]")
                        plt.tight_layout()
                        plt.show()
                else:
                    if not self.Quiet:
                        print(
                            'Incorrect plotChoice input, "'
                            + self.plotChoice
                            + '" provided.'
                        )
                        print(
                            "Possible input strings are: q, w, both, and (for 1D) combo"
                        )
                        print("Unable to produce plot.")

    def surfplot(self, z, titletext):
        """
        Plot if you want to - for troubleshooting - 1 figure
        """
        if self.latlon:
            plt.imshow(
                z, extent=(0, self.dx * z.shape[0], self.dy * z.shape[1], 0)
            )  # ,interpolation='nearest'
            plt.xlabel("longitude [deg E]", fontsize=12, fontweight="bold")
            plt.ylabel("latitude [deg N]", fontsize=12, fontweight="bold")
        else:
            plt.imshow(
                z,
                extent=(
                    0,
                    self.dx / 1000.0 * z.shape[0],
                    self.dy / 1000.0 * z.shape[1],
                    0,
                ),
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

        plt.subplot(211)
        plt.title("Load thickness, mantle equivalent [m]", fontsize=16)
        if self.latlon:
            plt.imshow(
                self.qs / (self.rho_m * self.g),
                extent=(0, self.dx * self.qs.shape[0], self.dy * self.qs.shape[1], 0),
            )
            plt.xlabel("longitude [deg E]", fontsize=12, fontweight="bold")
            plt.ylabel("latitude [deg N]", fontsize=12, fontweight="bold")
        else:
            plt.imshow(
                self.qs / (self.rho_m * self.g),
                extent=(
                    0,
                    self.dx / 1000.0 * self.qs.shape[0],
                    self.dy / 1000.0 * self.qs.shape[1],
                    0,
                ),
            )
            plt.xlabel("x [km]", fontsize=12, fontweight="bold")
            plt.ylabel("y [km]", fontsize=12, fontweight="bold")
        plt.colorbar()

        plt.subplot(212)
        plt.title("Deflection [m]")
        if self.latlon:
            plt.imshow(
                self.w,
                extent=(0, self.dx * self.w.shape[0], self.dy * self.w.shape[1], 0),
            )
            plt.xlabel("longitude [deg E]", fontsize=12, fontweight="bold")
            plt.ylabel("latitude [deg N]", fontsize=12, fontweight="bold")
        else:
            plt.imshow(
                self.w,
                extent=(
                    0,
                    self.dx / 1000.0 * self.w.shape[0],
                    self.dy / 1000.0 * self.w.shape[1],
                    0,
                ),
            )
            plt.xlabel("x [km]", fontsize=12, fontweight="bold")
            plt.ylabel("y [km]", fontsize=12, fontweight="bold")
        plt.colorbar()

    def xyzinterp(self, x, y, z, titletext):
        """
        Interpolates and plots ungridded model outputs from SAS_NG solution
        """
        # Help from http://wiki.scipy.org/Cookbook/Matplotlib/Gridding_irregularly_spaced_data

        if self.Verbose:
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
        if self.latlon:
            plt.contourf(xi, yi, zi, 100, cmap=plt.cm.jet)
        else:
            plt.contourf(xi / 1000.0, yi / 1000.0, zi, 100, cmap=plt.cm.jet)
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
    def __init__(self, filename=None):
        """
        WhichModel is a copy of initialization features inside the main class
        """
        self.filename = filename
        if self.filename:
            try:
                # only let this function imoprt things once
                self.whichModel_AlreadyRun
            except AttributeError:
                # Open parser and get what kind of model
                _fileisvalid = self.config = configparser.ConfigParser()
                _fileisvalid = len(_fileisvalid)
                if _fileisvalid:
                    try:
                        self.config.read(filename)
                        # Need to change this and all slashes to be Windows compatible
                        self.inpath = os.path.dirname(os.path.realpath(filename)) + "/"
                        # Need to have these guys inside "try" to make sure it is set up OK
                        # (at least for them)
                        self.dimension = self.configGet("integer", "mode", "dimension")
                        self.whichModel_AlreadyRun = True
                    except:
                        sys.exit(
                            ">>>> Error: cannot locate specified configuration file. <<<<"
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
        # 17 Nov 2014: Splitting out initialize from __init__ to allow space
        # to use getters and setters to define values

        # Use standard routine to pull out values
        # If no filename provided, will not initialize configuration file.
        self.filename = filename

        # DEFAULT VERBOSITY
        # Set default "quiet" to False, unless set by setter or overwritten by
        # the configuration file.
        self.Quiet = False
        # And also set default verbosity
        self.Verbose = True
        self.Debug = False

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

        # Default values for lat/lon usage -- defaulting not to use it
        try:
            self.latlon
        except AttributeError:
            self.latlon = False
        try:
            self.PlanetaryRadius
        except AttributeError:
            self.PlanetaryRadius = None

    def initialize(self, filename=None):
        # Values from configuration file

        # If a filename is provided here, overwrite any prior value
        if filename:
            if self.filename:
                pass  # Don't overwrite if filename is None-type
                # "Debug" not yet defined.
                # if self.Debug:
                #  print("Overwriting filename from '__init__' step with that from\n"+\
                #        "initialize step."
            else:
                # Update to new filename
                self.filename = filename

        if self.filename:
            # Set up ConfigParser
            self.config = configparser.ConfigParser()
            try:
                self.config.read(self.filename)
                self.inpath = os.path.dirname(os.path.realpath(self.filename)) + "/"
                # Need to have these guys inside "try" to make sure it is set up OK
                # (at least for them)
                self.dimension = self.configGet("integer", "mode", "dimension")
                self.whichModel_AlreadyRun = True
            except:
                sys.exit(
                    "No configuration file at specified path, or configuration file"
                    " configured incorrectly"
                )

            # Set verbosity for model run
            # Default is "verbose" with no debug or quiet
            # Verbose
            try:
                self.Verbose = self.configGet(
                    "bool", "verbosity", "Verbose", optional=False
                )
            except:
                pass
            # Deebug means that whole arrays, etc., can be printed
            try:
                self.Debug = self.configGet(
                    "bool", "verbosity", "Debug", optional=False
                )
            except:
                pass
            # Deebug means that whole arrays, etc., can be printed
            try:
                self.Quiet = self.configGet(
                    "bool", "verbosity", "Quiet", optional=False
                )
            except:
                pass
        # Quiet overrides all others
        if self.Quiet:
            self.Debug = False
            self.Verbose = False

        # Introduce model
        # After configuration file can define "Quiet", and getter/setter should be done
        # by this point if we are going that way.
        if not self.Quiet:
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

            self.Method = self.configGet("string", "mode", "method")
            # Boundary conditions
            # This used to be nested inside an "if self.Method == 'FD'", but it seems
            # better to define these to ensure there aren't mistaken impressions
            # about what they do for the SAS case
            # Not optional: flexural solutions can be very sensitive to b.c.'s
            self.BC_E = self.configGet(
                "string", "numerical", "BoundaryCondition_East", optional=False
            )
            self.BC_W = self.configGet(
                "string", "numerical", "BoundaryCondition_West", optional=False
            )
            if self.dimension == 2:
                self.BC_N = self.configGet(
                    "string", "numerical2D", "BoundaryCondition_North", optional=False
                )
                self.BC_S = self.configGet(
                    "string", "numerical2D", "BoundaryCondition_South", optional=False
                )

            # Parameters
            self.g = self.configGet("float", "parameter", "GravAccel")
            self.rho_m = self.configGet("float", "parameter", "MantleDensity")
            self.rho_fill = self.configGet(
                "float", "parameter", "InfillMaterialDensity"
            )

            # Grid spacing
            if self.Method != "SAS_NG":
                # No meaning for ungridded superimposed analytical solutions
                # From configuration file
                self.dx = self.configGet("float", "numerical", "GridSpacing_x")
                if self.dimension == 2:
                    self.dy = self.configGet("float", "numerical2D", "GridSpacing_y")

            # Mode: solution method and type of plate solution (if applicable)
            if self.filename:
                self.Method = self.configGet("string", "mode", "method")
                if self.dimension == 2:
                    self.PlateSolutionType = self.configGet(
                        "string", "mode", "PlateSolutionType"
                    )

            # Loading grid
            # q0 is either a load array or an x,y,q array.
            # Therefore q_0, initial q, before figuring out what it really is
            # for grid, q0 could also be written as $q_\sigma$ or q/(dx*(dy))
            # it is a surface normal stress that is h_load * rho_load * g
            # it later is combined with dx and (if 2D) dy for FD cases
            # for point loads, need mass: q0 should be written as [x, (y), force])
            self.q0 = self.configGet("string", "input", "Loads")

        # Parameters -- rho_m and rho_fill defined, so this outside
        # of if-statement (to work with getters/setters as well)
        self.drho = self.rho_m - self.rho_fill
        if self.filename:
            self.E = self.configGet("float", "parameter", "YoungsModulus")
            self.nu = self.configGet("float", "parameter", "PoissonsRatio")

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
                    sys.exit(
                        "Must define q0, q, or qs by this stage in the initialization"
                        " step from either configuration file (string) or direct array"
                        " import"
                    )
        else:
            # Stop program if q0 is None-type
            if self.q0 is None:  # if is None type, just be patient
                sys.exit(
                    "Must define non-None-type q0 by this stage in the initialization"
                    " step from either configuration file (string) or direct array"
                    " import"
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
            if self.Method != "SAS_NG":
                if self.q0.ndim != self.dimension:
                    print("Number of dimensions in loads file is inconsistent with")
                    print("number of dimensions in solution technique.")
                    print("Loads", self.q0.ndim)
                    print("Dimensions", self.dimension)
                    print(self.q0)
                    print("Exiting.")
                    sys.exit()

        # Plotting selection
        self.plotChoice = self.configGet("string", "output", "Plot", optional=True)

        # Ensure that Te is of floating-point type to avoid integer math
        # and floor division
        try:
            self.Te = self.Te.astype(float)  # array
        except:
            # Integer scalar Te does not seem to be a problem, but taking this step
            # anyway for consistency
            try:
                self.Te = float(self.Te)  # integer
            except:
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
            if self.Method != "FD":
                warnings.warn(
                    "End loads have been set but will not be implemented because the"
                    " solution method is not finite difference",
                    category=RuntimeWarning,
                    stacklevel=2,
                )
        try:
            self.sigma_xy
        except AttributeError:
            self.sigma_xy = 0
        else:
            if self.Method != "FD":
                warnings.warn(
                    "End loads have been set but will not be implemented because the"
                    " solution method is not finite difference",
                    category=RuntimeWarning,
                    stacklevel=2,
                )
        try:
            self.sigma_yy
        except AttributeError:
            self.sigma_yy = 0
        else:
            if self.Method != "FD":
                warnings.warn(
                    "End loads have been set but will not be implemented because the"
                    " solution method is not finite difference",
                    category=RuntimeWarning,
                    stacklevel=2,
                )

    # Finalize
    def finalize(self):
        # Can include an option for this later, but for the moment, this will
        # clear the coefficient array so it doens't cause problems for model runs
        # searching for the proper rigidity
        with contextlib.suppress(AttributeError):
            del self.coeff_matrix
        if not self.Quiet:
            print("")

    # SAVING TO FILE AND PLOTTING STEPS

    # Output: One of the functions run by isostasy.py; not part of IRF
    # (for standalone model use)
    def output(self):
        if self.Verbose:
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
            # If wOutFile exists, has already been set by a setter
            self.wOutFile
        except AttributeError:
            # Otherwise, it needs to be set by an configuration file
            try:
                self.wOutFile = self.configGet(
                    "string", "output", "DeflectionOut", optional=True
                )
            except:
                # if there is no parsable output string, do not generate output;
                # this allows the user to leave the line blank and produce no output
                if self.Debug:
                    print("No output filename provided:")
                    print("  not writing any deflection output to file")
        else:
            if self.Verbose:
                print("Output filename provided.")
        if self.wOutFile:
            if self.wOutFile[-4:] == ".npy":
                from numpy import save

                save(self.wOutFile, self.w)
            else:
                from numpy import savetxt

                # Shouldn't need more than mm precision, at very most
                savetxt(self.wOutFile, self.w, fmt="%.3f")
                if self.Verbose:
                    print("Saving deflections --> " + self.wOutFile)

    def bc_check(self):
        # Check that boundary conditions are acceptable with code implementation
        # Acceptable b.c.'s
        if self.Method == "FD":
            # Check if a coefficient array has been defined
            # It would only be by a getter or setter;
            # no way to do I/O with this with present configuration files
            # Define as None for use later.
            try:
                self.coeff_matrix
            except AttributeError:
                self.coeff_matrix = None
            # No need to create a coeff_matrix if one already exists
            if self.coeff_matrix is None:
                # Acceptable boundary conditions
                self.bc1D = np.array(
                    [
                        "0Displacement0Slope",
                        "Periodic",
                        "Mirror",
                        "0Moment0Shear",
                        "0Slope0Shear",
                    ]
                )
                self.bc2D = np.array(
                    [
                        "0Displacement0Slope",
                        "Periodic",
                        "Mirror",
                        "0Moment0Shear",
                        "0Slope0Shear",
                    ]
                )
                # Boundary conditions should be defined by this point -- whether via
                # the configuration file or the getters and setters
                self.bclist = [self.BC_E, self.BC_W]
                if self.dimension == 2:
                    self.bclist += [self.BC_N, self.BC_S]
                # Now check that these are valid boundary conditions
                for bc in self.bclist:
                    if self.dimension == 1:
                        if bc not in self.bc1D:
                            sys.exit(
                                f"{bc!r} is not an acceptable 1D finite difference"
                                " boundary condition and/or is not yet implement in"
                                " the code. Acceptable boundary conditions are:"
                                f" {', '.join(repr(bc) for bc in self.bc1D)}\n"
                                "Exiting."
                            )
                    elif self.dimension == 2:
                        if bc not in self.bc2D:
                            sys.exit(
                                f"{bc!r} is not an acceptable 2D finite difference"
                                " boundary condition and/or is not yet implement in"
                                " the code. Acceptable boundary conditions are:"
                                f" {', '.join(repr(bc) for bc in self.bc2D)}\n"
                                "Exiting."
                            )
                    else:
                        sys.exit(
                            "For a flexural solution, grid must be 1D or 2D. Exiting."
                        )
        else:
            # Analytical solution boundary conditions
            # If they aren't set, it is because no input file has been used
            # Just set them to an empty string (like input file would do)
            try:
                self.BC_E
            except AttributeError:
                self.BC_E = ""
            try:
                self.BC_W
            except AttributeError:
                self.BC_W = ""
            if self.dimension == 2:
                try:
                    self.BC_S
                except AttributeError:
                    self.BC_S = ""
                try:
                    self.BC_N
                except AttributeError:
                    self.BC_N = ""
            else:
                # Simplifies flow control a few lines down to define these as None-type
                self.BC_S = None
                self.BC_N = None
            if (
                self.BC_E == "NoOutsideLoads"
                or self.BC_E == ""
                and self.BC_W == "NoOutsideLoads"
                or self.BC_W == ""
            ) and (
                self.dimension != 2
                or (
                    self.BC_E == "NoOutsideLoads"
                    or self.BC_E == ""
                    and self.BC_W == "NoOutsideLoads"
                    or self.BC_W == ""
                )
            ):
                if (
                    self.BC_E == ""
                    or self.BC_W == ""
                    or self.BC_S == ""
                    or self.BC_N == ""
                ):
                    if self.Verbose:
                        print(
                            "Assuming NoOutsideLoads boundary condition, as this is"
                            " implicit in the superposition-based analytical solution"
                        )
            else:
                if not self.Quiet:
                    print("")
                    print(">>> BOUNDARY CONDITIONS IMPROPERLY DEFINED <<<")
                    print("")
                    print("For analytical solutions the boundaries must be either:")
                    print("")
                    print("* NoOutsideLoads (explicitly)")
                    print("* <left blank>")
                    print("")
                    print(
                        "The latter is to implictly indicate a desire to use the only"
                    )
                    print("boundary condition available for the superposition-based")
                    print("analytical solutions.")
                    print(
                        "This check is in place to ensure that the user does not apply"
                    )
                    print("boundary conditions for finite difference solutions to the")
                    print("analytical solutions and expect them to work.")
                    print("")
                    sys.exit()

    def coeffArraySizeCheck(self):
        """
        Make sure that q0 and coefficient array are the right size compared to
        each other (for finite difference if loading a pre-build coefficient
        array). Otherwise, exit.
        """
        if np.prod(self.coeff_matrix.shape) != np.long(
            np.prod(np.array(self.qs.shape, dtype=np.int64) + 2) ** 2
        ):
            print("Inconsistent size of q0 array and coefficient mattrix")
            print("Exiting.")
            sys.exit()

    def TeArraySizeCheck(self):
        """
        Checks that Te and q0 array sizes are compatible
        For finite difference solution.
        """
        # Only if they are both defined and are arrays
        # Both being arrays is a possible bug in this check routine that I have
        # intentionally introduced
        if isinstance(self.Te, np.ndarray) and isinstance(self.qs, np.ndarray):
            # Doesn't touch non-arrays or 1D arrays
            if type(self.Te) is np.ndarray:
                if (np.array(self.Te.shape) != np.array(self.qs.shape)).any():
                    sys.exit("q0 and Te arrays have incompatible shapes. Exiting.")
            else:
                if self.Debug:
                    print("Te and qs array sizes pass consistency check")

    ### need to determine its interface, it is best to have a uniform interface
    ### no matter it is 1D or 2D; but if it can't be that way, we can set up a
    ### variable-length arguments, which is the way how Python overloads functions.

    def FD(self):
        """
        Set-up for the finite difference solution method
        """
        if self.Verbose:
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
        # Is there a solver defined
        try:
            self.Solver  # See if it exists already
        except AttributeError:
            # Well, will fail if it doesn't see this, maybe not the most reasonable
            # error message.
            if self.filename:
                self.Solver = self.configGet("string", "numerical", "Solver")
            else:
                sys.exit("No solver defined!")
        # Check consistency of size if coeff array was loaded
        if self.filename:
            # In the case that it is iterative, find the convergence criterion
            self.iterative_ConvergenceTolerance = self.configGet(
                "float", "numerical", "ConvergenceTolerance"
            )
            # Try to import Te grid or scalar for the finite difference solution
            try:
                self.Te = self.configGet(
                    "float", "input", "ElasticThickness", optional=False
                )
                if self.Te is None:
                    Tepath = self.configGet(
                        "string", "input", "ElasticThickness", optional=False
                    )
                    self.Te = Tepath
                else:
                    Tepath = None
            except:
                Tepath = self.configGet(
                    "string", "input", "ElasticThickness", optional=False
                )
                self.Te = Tepath
            if self.Te is None:
                if self.coeff_matrix is not None:
                    pass
                else:
                    # Have to bring this out here in case it was discovered in the
                    # try statement that there is no value given
                    sys.exit(
                        "No input elastic thickness or coefficient matrix supplied."
                    )
        # or if getter/setter
        if isinstance(self.Te, str):
            # Try to import Te grid or scalar for the finite difference solution
            Tepath = self.Te
        else:
            Tepath = None  # in case no self.filename present (like for GRASS GIS)
        # If there is a Tepath, import Te
        # Assume that even if a coeff_matrix is defined
        # That the user wants Te if they gave the path
        if Tepath:
            self.Te = self.loadFile(self.Te, close_on_fail=False)
            if self.Te is None:
                print("Requested Te file is provided but cannot be located.")
                print("No scalar elastic thickness is provided in configuration file")
                print("(Typo in path to input Te grid?)")
                if self.coeff_matrix is not None:
                    print("But a coefficient matrix has been found.")
                    print("Calculations will be carried forward using it.")
                else:
                    print("Exiting.")
                    sys.exit()

            # Check that Te is the proper size if it was loaded
            # Will be array if it was loaded
            if self.Te.any():
                self.TeArraySizeCheck()

    ### need work
    def FFT(self):
        pass

    # SAS and SAS_NG are the exact same here; leaving separate just for symmetry
    # with other functions

    def SAS(self):
        """
        Set-up for the rectangularly-gridded superposition of analytical solutions
        method for solving flexure
        """
        if self.x is None:
            self.x = np.arange(self.dx / 2.0, self.dx * self.qs.shape[0], self.dx)
        if self.filename:
            # Define the (scalar) elastic thickness
            self.Te = self.configGet("float", "input", "ElasticThickness")
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

    def SAS_NG(self):
        """
        Set-up for the ungridded superposition of analytical solutions
        method for solving flexure
        """
        if self.filename:
            # Define the (scalar) elastic thickness
            self.Te = self.configGet("float", "input", "ElasticThickness")
            # See if it wants to be run in lat/lon
            # Could put under in 2D if-statement, but could imagine an eventual desire
            # to change this and have 1D lat/lon profiles as well.
            # So while the options will be under "numerical2D", this place here will
            # remain held for an eventual future.
            self.latlon = self.configGet(
                "string", "numerical2D", "latlon", optional=True
            )
            self.PlanetaryRadius = self.configGet(
                "float", "numerical2D", "PlanetaryRadius", optional=True
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
                    sys.exit(
                        "For 1D (ungridded) SAS_NG configuration file, need [x,w]"
                        f" array. Your dimensions are: {self.q0.shape}"
                    )
        else:
            try:
                # If these have already been set, e.g., by getters/setters, great!
                self.x
                self.u
                self.q
            except AttributeError:
                # Using [x, y, w] configuration file
                if self.q0.shape[1] == 3:
                    self.x = self.q0[:, 0]
                    self.y = self.q0[:, 1]
                    self.q = self.q0[:, 2]
                else:
                    sys.exit(
                        "For 2D (ungridded) SAS_NG configuration file, need [x,y,w]"
                        f" array. Your dimensions are: {self.q0.shape}"
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
            try:
                self.xw = self.configGet("string", "input", "xw", optional=True)
                if self.xw == "":
                    self.xw = None
            except:
                self.xw = None
        # If strings, load arrays
        if isinstance(self.xw, str):
            self.xw = self.loadFile(self.xw)
        if self.dimension == 2:
            try:
                # already set by setter?
                self.yw
            except AttributeError:
                try:
                    self.yw = self.configGet("string", "input", "yw", optional=True)
                    if self.yw == "":
                        self.yw = None
                except:
                    self.yw = None
            # At this point, can check if we have both None or both defined
            if (self.xw is not None and self.yw is None) or (
                self.xw is None and self.yw is not None
            ):
                sys.exit(
                    "SAS_NG output at specified points requires both xw and yw to be defined"
                )
            # All right, now just finish defining
            if isinstance(self.yw, str):
                self.yw = self.loadFile(self.yw)
            elif self.yw is None:
                self.yw = self.y.copy()
        if self.xw is None:
            self.xw = self.x.copy()
