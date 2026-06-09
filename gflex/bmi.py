"""
CSDMS Basic Model Interface (BMI) v2 implementation for gFlex.

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

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import NDArray

try:
    from bmipy.bmi import Bmi as _BmiBase
except ImportError as _err:
    _BmiBase = object  # type: ignore[assignment,misc]
    _bmipy_import_error: ImportError | None = _err
else:
    _bmipy_import_error = None

from gflex.base import WhichModel
from gflex.f1d import F1D
from gflex.f2d import F2D


class BmiGflex(_BmiBase):
    """BMI wrapper for gFlex lithospheric flexure.

    Implements the CSDMS Basic Model Interface v2 specification.
    Supports 1-D and 2-D gridded flexure solutions (FD, FFT, SAS methods).
    The SAS_NG point-load method is not suited to the BMI grid model.

    Grid
    ----
    All variables share grid identifier 0, a uniform rectilinear grid.
    Grid shape follows NumPy / C order: (nrows,) in 1-D or (nrows, ncols)
    in 2-D, with spacing (dy, dx) and origin at (0, 0).

    Time
    ----
    gFlex solves instantaneous elastic equilibrium. Time is therefore
    nominal: start=0, step=1, end=inf. Each call to update() applies
    the current load and computes deflection.

    Variables
    ---------
    Input:  ``lithosphere__load_pressure`` [Pa]  — surface normal stress q0
    Output: ``lithosphere__vertical_displacement`` [m] — downward deflection w
    """

    _name = "gFlex Lithospheric Flexure"
    _input_var_names = ("lithosphere__load_pressure",)
    _output_var_names = ("lithosphere__vertical_displacement",)

    _var_units = {
        "lithosphere__load_pressure": "Pa",
        "lithosphere__vertical_displacement": "m",
    }
    _var_grids = {
        "lithosphere__load_pressure": 0,
        "lithosphere__vertical_displacement": 0,
    }
    _var_loc = {
        "lithosphere__load_pressure": "node",
        "lithosphere__vertical_displacement": "node",
    }

    def __init__(self) -> None:
        """Initialize internal BMI state arrays.

        Requires the optional ``bmipy`` dependency (``pip install gflex[bmi]``).
        Call :meth:`initialize` with a configuration file before calling
        :meth:`update`.
        """
        if _bmipy_import_error is not None:
            raise ImportError(
                "bmipy is required to use BmiGflex. "
                "Install it with: pip install gflex[bmi]"
            ) from _bmipy_import_error
        self._model = None
        self._load = None
        self._w = None
        self._values: dict[str, NDArray[Any]] = {}
        self._shape: tuple[int, ...] = ()
        self._spacing: tuple[float, ...] = ()
        self._origin: tuple[float, ...] = ()
        self._current_time = 0.0

    # ------------------------------------------------------------------
    # Control functions
    # ------------------------------------------------------------------

    def initialize(self, config_file: str) -> None:
        """Initialize gFlex from a configuration file.

        Parameters
        ----------
        config_file : str
            Path to a gFlex YAML configuration file.
        """
        obj = WhichModel(config_file)
        if obj.dimension == 1:
            self._model = F1D(config_file)
        elif obj.dimension == 2:
            self._model = F2D(config_file)
        else:
            raise ValueError(f"Unsupported dimension: {obj.dimension}")

        self._model.initialize(config_file)

        # Own arrays for the two BMI-exposed variables.  The model's internal
        # q0 is consumed (renamed to qs, then deleted) during run(), so we
        # keep a separate copy that survives across update() calls.
        self._load = self._model.q0.copy()
        self._w = np.zeros(self._load.shape)

        if self._model.dimension == 1:
            self._spacing = (float(self._model.dx),)
        else:
            self._spacing = (float(self._model.dy), float(self._model.dx))
        self._shape = self._load.shape
        self._origin = (0.0,) * self._model.dimension
        self._current_time = 0.0

        self._values = {
            "lithosphere__load_pressure": self._load,
            "lithosphere__vertical_displacement": self._w,
        }

    def update(self) -> None:
        """Compute flexural deflection for the current load.

        Writes the current ``lithosphere__load_pressure`` array into the
        model's internal ``qs`` field, runs the solver, then copies the
        result into ``lithosphere__vertical_displacement``.
        """
        # Sync the BMI-owned load array into the model before each run.
        # Setting qs directly bypasses the q0 → qs copy that run() would
        # otherwise do, ensuring set_value() changes are always honoured.
        self._model.qs = self._load.copy()
        self._model.run()
        self._w[:] = self._model.w
        self._current_time += 1.0

    def update_until(self, time: float) -> None:
        """Advance model time to *time* by repeated calls to update().

        Parameters
        ----------
        time : float
            Target model time (must be >= current time).
        """
        while self._current_time < time:
            self.update()

    def finalize(self) -> None:
        """Tear down the model and release resources."""
        if self._model is not None:
            self._model.finalize()
            self._model = None

    # ------------------------------------------------------------------
    # Info functions
    # ------------------------------------------------------------------

    def get_component_name(self) -> str:
        """Return the human-readable name of this BMI component."""
        return self._name

    def get_input_item_count(self) -> int:
        """Return the number of input variables."""
        return len(self._input_var_names)

    def get_output_item_count(self) -> int:
        """Return the number of output variables."""
        return len(self._output_var_names)

    def get_input_var_names(self) -> tuple[str, ...]:
        """Return CSDMS Standard Names for all input variables."""
        return self._input_var_names

    def get_output_var_names(self) -> tuple[str, ...]:
        """Return CSDMS Standard Names for all output variables."""
        return self._output_var_names

    # ------------------------------------------------------------------
    # Variable info functions
    # ------------------------------------------------------------------

    def get_var_grid(self, name: str) -> int:
        """Return the grid identifier for variable *name*."""
        return self._var_grids[name]

    def get_var_type(self, name: str) -> str:
        """Return the NumPy dtype string for variable *name*."""
        return str(self.get_value_ptr(name).dtype)

    def get_var_units(self, name: str) -> str:
        """Return the UDUNITS-compatible unit string for variable *name*."""
        return self._var_units[name]

    def get_var_itemsize(self, name: str) -> int:
        """Return the size in bytes of one element of variable *name*."""
        return self.get_value_ptr(name).itemsize

    def get_var_nbytes(self, name: str) -> int:
        """Return the total number of bytes used by variable *name*."""
        return self.get_value_ptr(name).nbytes

    def get_var_location(self, name: str) -> str:
        """Return the grid location ('node', 'edge', or 'face') of variable *name*."""
        return self._var_loc[name]

    # ------------------------------------------------------------------
    # Time functions
    # ------------------------------------------------------------------

    def get_start_time(self) -> float:
        """Return the model start time (always 0.0)."""
        return 0.0

    def get_end_time(self) -> float:
        """Return the model end time (unbounded; returns ``inf``)."""
        return float("inf")

    def get_current_time(self) -> float:
        """Return the current model time (incremented by 1 each update)."""
        return self._current_time

    def get_time_step(self) -> float:
        """Return the model time step (always 1.0)."""
        return 1.0

    def get_time_units(self) -> str:
        """Return the time-unit string (``'s'``)."""
        return "s"

    # ------------------------------------------------------------------
    # Getters and setters
    # ------------------------------------------------------------------

    def get_value(self, name: str, dest: NDArray[Any]) -> NDArray[Any]:
        """Copy the flattened values of variable *name* into *dest* and return it."""
        dest[:] = self.get_value_ptr(name).flat
        return dest

    def get_value_ptr(self, name: str) -> NDArray[Any]:
        """Return a live reference to the internal array for variable *name*."""
        return self._values[name]

    def get_value_at_indices(
        self,
        name: str,
        dest: NDArray[Any],
        inds: NDArray[np.intp],
    ) -> NDArray[Any]:
        """Copy selected flat-indexed elements of variable *name* into *dest*."""
        dest[:] = self.get_value_ptr(name).flat[inds]
        return dest

    def set_value(self, name: str, src: NDArray[Any]) -> None:
        """Overwrite the entire array for variable *name* with values from *src*."""
        self.get_value_ptr(name).flat[:] = src

    def set_value_at_indices(
        self,
        name: str,
        inds: NDArray[np.intp],
        src: NDArray[Any],
    ) -> None:
        """Set selected flat-indexed elements of variable *name* from *src*."""
        self.get_value_ptr(name).flat[inds] = src

    # ------------------------------------------------------------------
    # Grid functions — uniform rectilinear
    # ------------------------------------------------------------------

    def get_grid_rank(self, grid: int) -> int:
        """Return the number of dimensions of grid *grid*."""
        return len(self._shape)

    def get_grid_size(self, grid: int) -> int:
        """Return the total number of nodes in grid *grid*."""
        return int(np.prod(self._shape))

    def get_grid_type(self, grid: int) -> str:
        """Return the grid type string (``'uniform_rectilinear'``)."""
        return "uniform_rectilinear"

    def get_grid_shape(
        self, grid: int, shape: NDArray[np.intp]
    ) -> NDArray[np.intp]:
        """Fill *shape* with the grid dimensions and return it."""
        shape[:] = self._shape
        return shape

    def get_grid_spacing(
        self, grid: int, spacing: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """Fill *spacing* with the grid cell spacings [m] and return it."""
        spacing[:] = self._spacing
        return spacing

    def get_grid_origin(
        self, grid: int, origin: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """Fill *origin* with the grid origin coordinates [m] and return it."""
        origin[:] = self._origin
        return origin

    # ------------------------------------------------------------------
    # Grid functions — not applicable for uniform rectilinear
    # ------------------------------------------------------------------

    def get_grid_x(
        self, grid: int, x: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """Not implemented — uniform rectilinear grids have no unstructured node coordinates."""
        raise NotImplementedError("get_grid_x")

    def get_grid_y(
        self, grid: int, y: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """Not implemented — uniform rectilinear grids have no unstructured node coordinates."""
        raise NotImplementedError("get_grid_y")

    def get_grid_z(
        self, grid: int, z: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """Not implemented — uniform rectilinear grids have no unstructured node coordinates."""
        raise NotImplementedError("get_grid_z")

    def get_grid_node_count(self, grid: int) -> int:
        """Not implemented — use :meth:`get_grid_size` for uniform rectilinear grids."""
        raise NotImplementedError("get_grid_node_count")

    def get_grid_edge_count(self, grid: int) -> int:
        """Not implemented — uniform rectilinear grids have no explicit edge topology."""
        raise NotImplementedError("get_grid_edge_count")

    def get_grid_face_count(self, grid: int) -> int:
        """Not implemented — uniform rectilinear grids have no explicit face topology."""
        raise NotImplementedError("get_grid_face_count")

    def get_grid_edge_nodes(
        self, grid: int, edge_nodes: NDArray[np.intp]
    ) -> NDArray[np.intp]:
        """Not implemented — uniform rectilinear grids have no unstructured edge topology."""
        raise NotImplementedError("get_grid_edge_nodes")

    def get_grid_face_edges(
        self, grid: int, face_edges: NDArray[np.intp]
    ) -> NDArray[np.intp]:
        """Not implemented — uniform rectilinear grids have no unstructured face topology."""
        raise NotImplementedError("get_grid_face_edges")

    def get_grid_face_nodes(
        self, grid: int, face_nodes: NDArray[np.intp]
    ) -> NDArray[np.intp]:
        """Not implemented — uniform rectilinear grids have no unstructured face topology."""
        raise NotImplementedError("get_grid_face_nodes")

    def get_grid_nodes_per_face(
        self, grid: int, nodes_per_face: NDArray[np.intp]
    ) -> NDArray[np.intp]:
        """Not implemented — uniform rectilinear grids have no unstructured face topology."""
        raise NotImplementedError("get_grid_nodes_per_face")
