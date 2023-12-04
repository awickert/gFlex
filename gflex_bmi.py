#! /usr/bin/env python
import numpy as np

# Imports for gFlex
from base import WhichModel
from f1d import F1D
from f2d import F2D
from prattairy import PrattAiry


class BmiGflex:
    _name = "Isostasy and Lithospheric Flexure"
    _input_var_names = ("earth_material_load__mass",)
    _output_var_names = (
        "lithosphere__vertical_displacement",
        "earth_material_load__mass",
    )
    _var_units = {
        "earth_material_load__mass": "kg",
        "lithosphere__vertical_displacement": "m",
    }

    def __init__(self):
        self._model = None
        self._values = {}
        self._shape = ()
        self._spacing = ()
        self._origin = ()
        self._coords = ()

    def initialize(self, config_file=None):
        ##############################
        # Deleted some of Eric's stuff here b/c I internally manage an input file.
        # Eric -- could you tell me if you have a better way for consistent
        # input files you would like to see CSDMS-compliant models employ?
        ##############################
        if config_file is None:
            pass
        else:
            self._model = WhichModel(
                config_file
            )  # This line should work outside if-statement as well.
            if self._model.model == "flexure":  # Really need to rename self.model!!!
                if self._model.dimension == 1:
                    self._model = F1D(config_file)
                elif self._model.dimension == 2:
                    self._model = F2D(config_file)
            elif self._model.model == "PrattAiry":
                self._model = PrattAiry(config_file)

        self._model.initialize(config_file)  # Does nothing

        if self._model.dimension == 1:
            self._spacing = (self._model.dx,)
            self._coords = (np.arange(self._model.q0.shape[0]) * self._model.dx,)
        elif self._model.dimension == 2:
            self._spacing = (self._model.dy, self._model.dx)
            self._coords = (
                np.arange(self._model.q0.shape[0]) * self._model.dy,
                np.arange(self._model.q0.shape[1]) * self._model.dx,
            )
        self._shape = self._model.q0.shape
        self._origin = (0.0,) * self._model.dimension

        self._w = np.empty_like(self._model.q0)

        # PROBABLY SHOULD RENAME "self.model"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # can remove plotting and file output options -- not meant to be part of
        # BMI interface!!!!!!!!!
        self._values = {
            "earth_material_load__mass": self._model.q0,
            "lithosphere__vertical_displacement": self._w,
        }

    def update(self):
        self._model.run()
        self._w[:] = self._model.w

    def update_frac(self, time_frac):
        self.update()

    def update_until(self, then):
        self.update()

    def finalize(self):
        self._model.finalize()
        self._model = None

    def get_var_type(self, var_name):
        return str(self.get_value_ptr(var_name).dtype)

    def get_var_units(self, var_name):
        return self._var_units[var_name]

    def get_var_rank(self, var_name):
        return self.get_value_ptr(var_name).ndim

    def get_var_size(self, var_name):
        return self.get_value_ptr(var_name).size

    def get_var_nbytes(self, var_name):
        return self.get_value_ptr(var_name).nbytes

    def get_value_ptr(self, var_name):
        return self._values[var_name]

    def get_value(self, var_name):
        return self.get_value_ptr(var_name).copy()

    def get_value_at_indices(self, var_name, indices):
        return self.get_value_ptr(var_name).take(indices)

    def set_value(self, var_name, src):
        val = self.get_value_ptr(var_name)
        val[:] = src

    def set_value_at_indices(self, var_name, src, indices):
        val = self.get_value_ptr(var_name)
        val.flat[indices] = src

    def get_component_name(self):
        return self._name

    def get_input_var_names(self):
        return self._input_var_names

    def get_output_var_names(self):
        return self._output_var_names

    def get_grid_shape(self, var_name):
        return self.get_value_ptr(var_name).shape

    def get_grid_spacing(self, var_name):
        if var_name in self._values:
            return self._spacing

    def get_grid_origin(self, var_name):
        if var_name in self._values:
            return self._origin

    def get_grid_type(self, var_name):
        if var_name in self._values:
            return "uniform_rectilinear"
        else:
            raise KeyError(var_name)

    def get_grid_x(self, var_name):
        if var_name in self._values:
            return self._coords[-1]
        else:
            raise KeyError(var_name)

    def get_grid_y(self, var_name):
        if var_name in self._values:
            return self._coords[-2]
        else:
            raise KeyError(var_name)

    def get_start_time(self):
        raise NotImplementedError("get_start_time")

    def get_end_time(self):
        raise NotImplementedError("get_end_time")

    def get_current_time(self):
        raise NotImplementedError("get_current_time")

    def get_time_step(self):
        raise NotImplementedError("get_time_step")
