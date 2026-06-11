#! /usr/bin/env python
"""Tests for BmiGflex."""

import pathlib

import numpy as np
import pytest

bmipy = pytest.importorskip("bmipy")

from gflex.bmi import BmiGflex

INPUT_DIR = pathlib.Path(__file__).parent.parent / "input"
CONFIG_1D = str(INPUT_DIR / "input_f1d.yaml")
CONFIG_2D = str(INPUT_DIR / "input_f2d.yaml")


@pytest.fixture
def bmi_1d():
    """Initialized 1-D BmiGflex; finalized after the test."""
    bmi = BmiGflex()
    bmi.initialize(CONFIG_1D)
    yield bmi
    bmi.finalize()


@pytest.fixture
def bmi_2d():
    """Initialized 2-D BmiGflex; finalized after the test."""
    bmi = BmiGflex()
    bmi.initialize(CONFIG_2D)
    yield bmi
    bmi.finalize()


# ------------------------------------------------------------------
# Smoke tests
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
    assert bmi_1d.get_input_var_names() == (
        "load__normal_component_of_stress",
        "lithosphere__elastic_thickness",
        "lithosphere__young_modulus",
        "lithosphere__poisson_ratio",
        "mantle__mass-per-volume_density",
        "infill_material__mass-per-volume_density",
        "planet_surface__gravitational_acceleration",
    )


def test_get_output_var_names(bmi_1d):
    assert bmi_1d.get_output_var_names() == ("lithosphere__vertical_displacement",)


def test_get_input_item_count(bmi_1d):
    assert bmi_1d.get_input_item_count() == 7


def test_get_output_item_count(bmi_1d):
    assert bmi_1d.get_output_item_count() == 1


def test_get_var_units_load(bmi_1d):
    assert bmi_1d.get_var_units("load__normal_component_of_stress") == "Pa"


def test_get_var_units_te(bmi_1d):
    assert bmi_1d.get_var_units("lithosphere__elastic_thickness") == "m"


def test_get_var_units_deflection(bmi_1d):
    assert bmi_1d.get_var_units("lithosphere__vertical_displacement") == "m"


def test_get_var_units_constants(bmi_1d):
    assert bmi_1d.get_var_units("lithosphere__young_modulus") == "Pa"
    assert bmi_1d.get_var_units("lithosphere__poisson_ratio") == "1"
    assert bmi_1d.get_var_units("mantle__mass-per-volume_density") == "kg m-3"
    assert bmi_1d.get_var_units("infill_material__mass-per-volume_density") == "kg m-3"
    assert bmi_1d.get_var_units("planet_surface__gravitational_acceleration") == "m s-2"


def test_get_var_location(bmi_1d):
    assert bmi_1d.get_var_location("load__normal_component_of_stress") == "node"
    assert bmi_1d.get_var_location("lithosphere__elastic_thickness") == "node"
    assert bmi_1d.get_var_location("lithosphere__vertical_displacement") == "node"
    assert bmi_1d.get_var_location("lithosphere__young_modulus") == "node"
    assert bmi_1d.get_var_location("planet_surface__gravitational_acceleration") == "node"


def test_get_var_grid(bmi_1d):
    assert bmi_1d.get_var_grid("load__normal_component_of_stress") == 0
    assert bmi_1d.get_var_grid("lithosphere__elastic_thickness") == 0
    assert bmi_1d.get_var_grid("lithosphere__vertical_displacement") == 0
    assert bmi_1d.get_var_grid("lithosphere__young_modulus") == 1
    assert bmi_1d.get_var_grid("lithosphere__poisson_ratio") == 1
    assert bmi_1d.get_var_grid("mantle__mass-per-volume_density") == 1
    assert bmi_1d.get_var_grid("infill_material__mass-per-volume_density") == 1
    assert bmi_1d.get_var_grid("planet_surface__gravitational_acceleration") == 1


# ------------------------------------------------------------------
# Grid 1 (scalar constants)
# ------------------------------------------------------------------


def test_get_grid_rank_scalar(bmi_1d):
    assert bmi_1d.get_grid_rank(1) == 1


def test_get_grid_size_scalar(bmi_1d):
    assert bmi_1d.get_grid_size(1) == 1


def test_get_grid_shape_scalar(bmi_1d):
    shape = np.zeros(1, dtype=np.intp)
    bmi_1d.get_grid_shape(1, shape)
    assert shape[0] == 1


def test_get_grid_spacing_scalar(bmi_1d):
    spacing = np.zeros(1)
    bmi_1d.get_grid_spacing(1, spacing)
    assert spacing[0] == pytest.approx(1.0)


def test_get_grid_origin_scalar(bmi_1d):
    origin = np.zeros(1)
    bmi_1d.get_grid_origin(1, origin)
    assert origin[0] == pytest.approx(0.0)


# ------------------------------------------------------------------
# Scalar constant get/set round-trips
# ------------------------------------------------------------------


def test_get_constants_nonzero(bmi_1d):
    """All scalar physical constants must be positive after initialize."""
    for name in (
        "lithosphere__young_modulus",
        "lithosphere__poisson_ratio",
        "mantle__mass-per-volume_density",
        "planet_surface__gravitational_acceleration",
    ):
        dest = np.zeros(1)
        bmi_1d.get_value(name, dest)
        assert dest[0] > 0, f"{name} should be positive"


def test_get_infill_density_nonnegative(bmi_1d):
    dest = np.zeros(1)
    bmi_1d.get_value("infill_material__mass-per-volume_density", dest)
    assert dest[0] >= 0


def test_set_rho_fill_changes_deflection(bmi_1d):
    """Higher infill density reduces restoring buoyancy, increasing deflection.

    Δρ = ρ_m − ρ_fill, so rho_fill=1030 (ocean) < rho_fill=0 (air) restoring
    force → larger-magnitude deflection with ocean infill.
    """
    n = bmi_1d.get_grid_size(0)

    bmi_1d.set_value("infill_material__mass-per-volume_density", np.array([0.0]))
    bmi_1d.update()
    w_air = np.zeros(n)
    bmi_1d.get_value("lithosphere__vertical_displacement", w_air)

    bmi_1d.set_value("infill_material__mass-per-volume_density", np.array([1030.0]))
    bmi_1d.update()
    w_ocean = np.zeros(n)
    bmi_1d.get_value("lithosphere__vertical_displacement", w_ocean)

    assert np.abs(w_ocean).max() > np.abs(w_air).max(), (
        "Ocean infill (rho_fill=1030) reduces Δρ, so peak deflection must exceed air (rho_fill=0)"
    )


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


def test_update_until_zero_does_not_advance(bmi_1d):
    """update_until(0.0) at t=0 must not call update() at all."""
    assert bmi_1d.get_current_time() == pytest.approx(0.0)
    bmi_1d.update_until(0.0)
    assert bmi_1d.get_current_time() == pytest.approx(0.0)


def test_update_until_past_time_is_no_op(bmi_1d):
    """update_until with a time <= current_time must not advance the model."""
    bmi_1d.update()                          # advance to t = 1.0
    t_before = bmi_1d.get_current_time()
    bmi_1d.update_until(0.5)                 # target in the past
    assert bmi_1d.get_current_time() == pytest.approx(t_before)


def test_update_until_2d(bmi_2d):
    """update_until works for 2-D BMI across multiple steps."""
    bmi_2d.update_until(3.0)
    assert bmi_2d.get_current_time() == pytest.approx(3.0)


# ------------------------------------------------------------------
# get_value / set_value round-trip
# ------------------------------------------------------------------


def test_deflection_nonzero_after_update(bmi_1d):
    """Non-zero load must produce non-zero deflection after update()."""
    n = bmi_1d.get_grid_size(0)
    load = np.zeros(n)
    bmi_1d.get_value("load__normal_component_of_stress", load)
    # The fixture config has a non-zero load from the file.
    bmi_1d.update()
    w = np.zeros(n)
    bmi_1d.get_value("lithosphere__vertical_displacement", w)
    assert np.any(w != 0.0)


def test_get_te_nonzero(bmi_1d):
    """lithosphere__elastic_thickness must be positive after initialize."""
    n = bmi_1d.get_grid_size(0)
    te = np.zeros(n)
    bmi_1d.get_value("lithosphere__elastic_thickness", te)
    assert np.all(te > 0)


def test_set_te_changes_deflection(bmi_1d):
    """A stiffer Te produces smaller-magnitude deflection under the same load."""
    n = bmi_1d.get_grid_size(0)
    te = np.zeros(n)
    bmi_1d.get_value("lithosphere__elastic_thickness", te)
    te_original = te.copy()

    bmi_1d.update()
    w_soft = np.zeros(n)
    bmi_1d.get_value("lithosphere__vertical_displacement", w_soft)

    # Double elastic thickness → much stiffer plate → smaller deflection
    bmi_1d.set_value("lithosphere__elastic_thickness", te_original * 2.0)
    bmi_1d.update()
    w_stiff = np.zeros(n)
    bmi_1d.get_value("lithosphere__vertical_displacement", w_stiff)

    assert np.abs(w_stiff).max() < np.abs(w_soft).max(), (
        "Doubling Te must reduce peak deflection magnitude"
    )


def test_set_value_affects_deflection(bmi_1d):
    """Setting load to zero then updating must yield zero deflection."""
    n = bmi_1d.get_grid_size(0)
    bmi_1d.set_value("load__normal_component_of_stress", np.zeros(n))
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
    bmi_1d.set_value("load__normal_component_of_stress", np.zeros(n))
    inds = np.array([n // 2], dtype=np.intp)
    bmi_1d.set_value_at_indices("load__normal_component_of_stress", inds, np.array([1e9]))
    src = np.zeros(1)
    bmi_1d.get_value_at_indices("load__normal_component_of_stress", src, inds)
    assert src[0] == pytest.approx(1e9)


# ==================================================================
# 2-D BMI tests
# ==================================================================


def test_2d_initialize_update_finalize():
    bmi = BmiGflex()
    bmi.initialize(CONFIG_2D)
    bmi.update()
    bmi.finalize()


# ------------------------------------------------------------------
# 2-D grid introspection
# ------------------------------------------------------------------


def test_get_grid_rank_2d(bmi_2d):
    assert bmi_2d.get_grid_rank(0) == 2


def test_get_grid_shape_2d(bmi_2d):
    shape = np.zeros(2, dtype=np.intp)
    bmi_2d.get_grid_shape(0, shape)
    assert len(shape) == 2
    assert np.all(shape > 0)


def test_get_grid_size_2d(bmi_2d):
    n = bmi_2d.get_grid_size(0)
    shape = np.zeros(2, dtype=np.intp)
    bmi_2d.get_grid_shape(0, shape)
    assert n == shape[0] * shape[1]


def test_get_grid_spacing_2d(bmi_2d):
    spacing = np.zeros(2)
    bmi_2d.get_grid_spacing(0, spacing)
    assert np.all(spacing > 0.0)


def test_get_grid_origin_2d(bmi_2d):
    origin = np.zeros(2)
    bmi_2d.get_grid_origin(0, origin)
    assert origin[0] == pytest.approx(0.0)
    assert origin[1] == pytest.approx(0.0)


# ------------------------------------------------------------------
# 2-D get_value / set_value round-trip
# ------------------------------------------------------------------


def test_2d_deflection_nonzero_after_update(bmi_2d):
    """Non-zero load must produce non-zero deflection after update()."""
    bmi_2d.update()
    n = bmi_2d.get_grid_size(0)
    w = np.zeros(n)
    bmi_2d.get_value("lithosphere__vertical_displacement", w)
    assert np.any(w != 0.0)


def test_2d_set_value_affects_deflection(bmi_2d):
    """Setting load to zero then updating must yield zero deflection."""
    n = bmi_2d.get_grid_size(0)
    bmi_2d.set_value("load__normal_component_of_stress", np.zeros(n))
    bmi_2d.update()
    w = np.zeros(n)
    bmi_2d.get_value("lithosphere__vertical_displacement", w)
    assert np.allclose(w, 0.0)
