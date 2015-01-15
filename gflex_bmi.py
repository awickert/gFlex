#! /usr/bin/env python

# Imports for BMI
import types
import numpy as np
from bmi import Bmi, BmiGridType

# Imports for gFlex 
from base import *
from f1d import *
from f2d import *
from prattairy import *
import time


class BmiHeat(Bmi):
  _name = 'Isostasy and Lithospheric Flexure'
  ##############################
  # Not sure if I should include non-physical input variables that tell
  # the model how to behave.
  # I haven't included variables to set varbosity...
  ##############################
  _input_var_names = ['model__solution_type',
                      'model__dimensions',
                      'model__solution_method',
                      'model__elastic_plate_solution_type',
                      'lithosphere__young_modulus',
                      'lithosphere__poisson_ratio',
                      'standard_gravity_constant',
                      'mantle__mass-per-volume_density',
                      'infill_earth_material__mass-per-volume_density',
                      'lithosphere__elastic_thickness',
                      'earth_material_load__pressure',
                      'model__plot_output',
                      'model_grid_cell__x_length',
                      'model_grid_cell__y_length',
                      'model__solver_type',
                      'model__iterative_convergence_tolerance',
                      'model__finite_difference_coefficient_array',
                      'model__boundary_condition_west',
                      'model__boundary_condition_east',
                      'model__boundary_condition_north',
                      'model__boundary_condition_south',
                      'model__output_file_name']
  _output_var_names = ['lithosphere__vertical_displacement']
  _var_units = {'model__solution_type' : '',
                'model__dimensions' : '',
                'model__solution_method' : '',
                'model__elastic_plate_solution_type' : '',
                'lithosphere__young_modulus' : 'Pa',
                'lithosphere__poisson_ratio' : '',
                'standard_gravity_constant' : 'm s**(-2)',
                'mantle__mass-per-volume_density' : 'kg m**3',
                'infill_earth_material__mass-per-volume_density' : 'kg m**3',
                'lithosphere__elastic_thickness' : 'm',
                'earth_material_load__pressure' : 'Pa',
                'model__plot_output' : '',
                'model_grid_cell__x_length' : 'm',
                'model_grid_cell__y_length' : 'm',
                'model__solver_type' : '',
                'model__iterative_convergence_tolerance' : 'm',
                'model__finite_difference_coefficient_array' : '',
                'model__boundary_condition_west' : '',
                'model__boundary_condition_east' : '',
                'model__boundary_condition_north' : '',
                'model__boundary_condition_south' : '',
                'lithosphere__vertical_displacement' : 'm',
                'model__output_file_name' : ''}

  def __init__(self):
    self._model = None
    self._values = {}

  def initialize(self, config_file=None):
    ##############################
    # Deleted some of Eric's stuff here b/c I internally manage an input file.
    # Eric -- could you tell me if you have a better way for consistent
    # input files you would like to see CSDMS-compliant models employ?
    ##############################
    if config_file is None:
      pass
    else:
      self._model = WhichModel(config_file) # This line should work outside if-statement as well.
      if self._model.model == 'flexure': # Really need to rename self.model!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if self._model.dimension == 1:
          self._model = F1D(config_file)
        elif obj.dimension == 2:
          self._model = F2D(config_file)
      elif obj.model == 'PrattAiry':
        self._model = PrattAiry(config_file)

    obj.initialize(config_file) # Does nothing

    self._values = {
        'plate_surface__temperature': self._model.z,
    }
    
    # PROBABLY SHOULD RENAME "self.model"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # can remove plotting and file output options -- not meant to be part of BMI interface!!!!!!!!!
    self._values = {'model__solution_type' : self.model,
                    'model__dimensions' : self.dimension,
                    'model__solution_method' : self.method,
                    'model__elastic_plate_solution_type' : self.PlateSolutionType,
                    'lithosphere__young_modulus' : self.E,
                    'lithosphere__poisson_ratio' : self.nu,
                    'standard_gravity_constant' : self.g,
                    'mantle__mass-per-volume_density' : self.rho_m,
                    'infill_earth_material__mass-per-volume_density' : self.rho_fill,
                    'lithosphere__elastic_thickness' : self.Te,
                    'earth_material_load__mass' : self.q0,
                    'model__plot_output' : self.plotChoice,
                    'model_grid_cell__x_length' : self.dx,
                    'model_grid_cell__y_length' : self.dy,
                    'model__solver_type' : self.solver,
                    'model__iterative_convergence_tolerance' : 'm',
                    'model__finite_difference_coefficient_array' : '',
                    'model__boundary_condition_west' : self.BC_W,
                    'model__boundary_condition_east' : self.BC_E,
                    'model__boundary_condition_north' : self.BC_N,
                    'model__boundary_condition_south' : self.BC_S,
                    'lithosphere__vertical_displacement' : self.w,
                    'model__output_file_name' : self.wOutFile}

  def update(self):
    # NOT TIME-STEPPING!!! ONLY DOES THIS ONCE!
    self.run()

  def update_frac(self, time_frac):
    sys.exit("Non-time-evolving")

  def update_until(self, then):
    sys.exit("Non-time-evolving")

  def finalize(self):
    obj.finalize()
    self._model = None

  def get_var_type (self, var_name):
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

  def get_grid_shape (self, var_name):
    return self.get_value_ptr(var_name).shape

  def get_grid_spacing(self, var_name):
    if var_name in self._values:
      return self._model.spacing

  def get_grid_origin(self, var_name):
    if var_name in self._values:
      return self._model.origin

  #################################
  # HOW DO I DO RECTILINEAR??????????????????????????????????????????????????????????
  #######################
  def get_grid_type(self, var_name):
    if var_name in self._values:
      return BmiGridType.UNIFORM
    else:
      return BmiGridType.UNKNOWN

  def get_start_time (self):
    sys.exit("Non-time-evolving")

  def get_end_time (self):
    sys.exit("Non-time-evolving")

  def get_current_time (self):
    sys.exit("Non-time-evolving")

  def get_time_step (self):
    sys.exit("Non-time-evolving")

