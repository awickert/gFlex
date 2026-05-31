#! /usr/bin/env python
"""Tests for BmiGflex."""

import pathlib

import numpy as np
import pytest

bmipy = pytest.importorskip("bmipy")

from gflex.bmi import BmiGflex

INPUT_DIR = pathlib.Path(__file__).parent.parent / "input"
CONFIG_1D = str(INPUT_DIR / "input_f1d.yaml")


@pytest.fixture
def bmi_1d():
    """Initialized 1-D BmiGflex; finalized after the test."""
    bmi = BmiGflex()
    bmi.initialize(CONFIG_1D)
    yield bmi
    bmi.finalize()


# ------------------------------------------------------------------
# Smoke test
# ------------------------------------------------------------------


def test_initialize_update_finalize():
    bmi = BmiGflex()
    bmi.initialize(CONFIG_1D)
    bmi.update()
    bmi.finalize()


# ------------------------------------------------------------------
# Variable introspection
# ------------------------------------------------------------------


def test_get_input_var_names(bmi_1d):
    assert bmi_1d.get_input_var_names() == ("lithosphere__load_pressure",)


def test_get_output_var_names(bmi_1d):
    assert bmi_1d.get_output_var_names() == ("lithosphere__vertical_displacement",)


def test_get_input_item_count(bmi_1d):
    assert bmi_1d.get_input_item_count() == 1


def test_get_output_item_count(bmi_1d):
    assert bmi_1d.get_output_item_count() == 1


def test_get_var_units_load(bmi_1d):
    assert bmi_1d.get_var_units("lithosphere__load_pressure") == "Pa"


def test_get_var_units_deflection(bmi_1d):
    assert bmi_1d.get_var_units("lithosphere__vertical_displacement") == "m"


def test_get_var_location(bmi_1d):
    assert bmi_1d.get_var_location("lithosphere__load_pressure") == "node"
    assert bmi_1d.get_var_location("lithosphere__vertical_displacement") == "node"


def test_get_var_grid(bmi_1d):
    assert bmi_1d.get_var_grid("lithosphere__load_pressure") == 0
    assert bmi_1d.get_var_grid("lithosphere__vertical_displacement") == 0


def test_get_component_name(bmi_1d):
    assert "gFlex" in bmi_1d.get_component_name()


# ------------------------------------------------------------------
# Grid introspection
# ------------------------------------------------------------------


def test_get_grid_type(bmi_1d):
    assert bmi_1d.get_grid_type(0) == "uniform_rectilinear"


def test_get_grid_rank_1d(bmi_1d):
    assert bmi_1d.get_grid_rank(0) == 1


def test_get_grid_size(bmi_1d):
    n = bmi_1d.get_grid_size(0)
    assert n > 0
    shape = np.zeros(1, dtype=np.intp)
    bmi_1d.get_grid_shape(0, shape)
    assert n == shape[0]


def test_get_grid_spacing_positive(bmi_1d):
    spacing = np.zeros(1)
    bmi_1d.get_grid_spacing(0, spacing)
    assert spacing[0] > 0.0


def test_get_grid_origin_zero(bmi_1d):
    origin = np.zeros(1)
    bmi_1d.get_grid_origin(0, origin)
    assert origin[0] == pytest.approx(0.0)


# ------------------------------------------------------------------
# Time functions
# ------------------------------------------------------------------


def test_start_time(bmi_1d):
    assert bmi_1d.get_start_time() == pytest.approx(0.0)


def test_end_time_infinite(bmi_1d):
    assert bmi_1d.get_end_time() == float("inf")


def test_time_step(bmi_1d):
    assert bmi_1d.get_time_step() == pytest.approx(1.0)


def test_current_time_advances(bmi_1d):
    t0 = bmi_1d.get_current_time()
    bmi_1d.update()
    assert bmi_1d.get_current_time() == pytest.approx(t0 + 1.0)


def test_update_until(bmi_1d):
    bmi_1d.update_until(3.0)
    assert bmi_1d.get_current_time() == pytest.approx(3.0)


# ------------------------------------------------------------------
# get_value / set_value round-trip
# ------------------------------------------------------------------


def test_deflection_nonzero_after_update(bmi_1d):
    """Non-zero load must produce non-zero deflection after update()."""
    n = bmi_1d.get_grid_size(0)
    load = np.zeros(n)
    bmi_1d.get_value("lithosphere__load_pressure", load)
    # The fixture config has a non-zero load from the file.
    bmi_1d.update()
    w = np.zeros(n)
    bmi_1d.get_value("lithosphere__vertical_displacement", w)
    assert np.any(w != 0.0)


def test_set_value_affects_deflection(bmi_1d):
    """Setting load to zero then updating must yield zero deflection."""
    n = bmi_1d.get_grid_size(0)
    bmi_1d.set_value("lithosphere__load_pressure", np.zeros(n))
    bmi_1d.update()
    w = np.zeros(n)
    bmi_1d.get_value("lithosphere__vertical_displacement", w)
    assert np.allclose(w, 0.0)


def test_get_value_ptr_is_live(bmi_1d):
    """get_value_ptr() must return the same object after update()."""
    ptr_before = bmi_1d.get_value_ptr("lithosphere__vertical_displacement")
    bmi_1d.update()
    ptr_after = bmi_1d.get_value_ptr("lithosphere__vertical_displacement")
    assert ptr_before is ptr_after


def test_set_value_at_indices(bmi_1d):
    n = bmi_1d.get_grid_size(0)
    bmi_1d.set_value("lithosphere__load_pressure", np.zeros(n))
    inds = np.array([n // 2], dtype=np.intp)
    bmi_1d.set_value_at_indices("lithosphere__load_pressure", inds, np.array([1e9]))
    src = np.zeros(1)
    bmi_1d.get_value_at_indices("lithosphere__load_pressure", src, inds)
    assert src[0] == pytest.approx(1e9)
