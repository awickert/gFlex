#!/usr/bin/env python

import numpy as np

# LANDLAB
from landlab import RasterModelGrid, Component

# GFLEX
from base import *
from f1d import *
from f2d import *
from prattairy import *

# PYTHON
import numpy as np
import time


_VALID_METHODS = set(['airy', 'flexure'])


def assert_method_is_valid(method):
    if method not in _VALID_METHODS:
        raise ValueError('%s: Invalid method name' % method)


class gFlex(Component):

    _name = 'gFlex'

    _input_var_names = set(['model__solution_type',
                            'model__dimensions',
                            'model__solution_method',
                            'model__elastic_plate_solution_type',
                            'lithosphere__young_modulus',
                            'lithosphere__poisson_ratio',
                            'standard_gravity_constant',
                            'mantle__mass-per-volume_density',
                            'infill_earth_material__mass-per-volume_density',
                            'lithosphere__elastic_thickness',
                            'earth_material_load__mass',
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
                            'model__output_file_name'])

    _output_var_names = set(['lithosphere__vertical_displacement'])

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
                  'earth_material_load__mass' : 'kg',
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
    
    # THIS CURRENTLY EXISTS AS A PLACEHOLDER UNTIL I HEAR BACK FROM DAN ON A 
    # UNIFORM WAY OF INCLUDING MODULES
    # ALSO, IS THERE ANY WAY OF CHECKING THAT IT IS A RECTANGULAR GRID AND/OR INTERPOLATING???

    def __init__(self, grid, **kwds):
        self._eet = kwds.pop('eet', 65000.)
        self._youngs = kwds.pop('youngs', 7e10)
        self._method = kwds.pop('method', 'airy')

        super(FlexureComponent, self).__init__(grid, **kwds)

        for name in self._input_var_names:
            if not name in self.grid.at_node:
                self.grid.add_zeros('node', name, units=self._var_units[name])

        for name in self._output_var_names:
            if not name in self.grid.at_node:
                self.grid.add_zeros('node', name, units=self._var_units[name])

        self._last_load = self.grid.field_values('node', 'lithosphere__overlying_pressure').copy()

        self._nodal_values = self.grid['node']

    def update(self, n_procs=1):
        elevation = self._nodal_values['lithosphere__elevation']
        load = self._nodal_values['lithosphere__overlying_pressure']
        deflection = self._nodal_values['lithosphere__elevation_increment']
        deposition = self._nodal_values['planet_surface_sediment__deposition_increment']

        new_load = ((load - self._last_load) +
                    (deposition * 2650. * 9.81).flat)

        self._last_load = load.copy()

        deflection.fill(0.)

        if self._method == 'airy':
            deflection[:] = new_load / (3300. * 9.81)
        else:
            self.subside_loads(new_load, self.coords, deflection=deflection,
                               n_procs=n_procs)

        elevation -= deflection

    def subside_loads(self, loads, locs, deflection=None, n_procs=1):
        if deflection is None:
            deflection = np.empty(self.shape, dtype=np.float)

        subside_point_loads(loads, locs, self.coords, self._eet, self._youngs,
                            deflection=deflection, n_procs=n_procs)

        return deflection

    def subside_load(self, load, loc, deflection=None):
        subside_point_load(
            load, loc, self.coords, self._eet, self._youngs,
            deflection=self._nodal_values['lithosphere__elevation_increment'])

        return deflection


if __name__ == "__main__":
    import doctest
    doctest.testmod()
