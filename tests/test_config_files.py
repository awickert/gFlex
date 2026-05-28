#! /usr/bin/env python
"""Tests for INI and YAML configuration file loading.

Verifies that both formats parse correctly and produce identical physics.
"""

import pathlib

import numpy as np
import pytest

from gflex.base import WhichModel
from gflex.f1d import F1D
from gflex.f2d import F2D

INPUT_DIR = pathlib.Path(__file__).parent.parent / "input"


def _run_from_config(path):
    """Run gFlex from a config file; return the flex object after finalize()."""
    path = str(path)
    obj = WhichModel(path)
    if obj.dimension == 1:
        flex = F1D(path)
    elif obj.dimension == 2:
        flex = F2D(path)
    else:
        raise ValueError(f"Unexpected dimension: {obj.dimension}")
    flex.initialize(path)
    flex.run()
    flex.finalize()
    return flex


# ---------------------------------------------------------------------------
# Smoke tests — each format loads and runs without error
# ---------------------------------------------------------------------------

def test_yaml_1d_runs():
    flex = _run_from_config(INPUT_DIR / "input_f1d.yaml")
    assert flex.w.ndim == 1
    assert flex.w.size > 0
    assert not np.any(np.isnan(flex.w))
    assert np.any(flex.w != 0)


def test_yaml_2d_runs():
    flex = _run_from_config(INPUT_DIR / "input_f2d.yaml")
    assert flex.w.ndim == 2
    assert flex.w.size > 0
    assert not np.any(np.isnan(flex.w))
    assert np.any(flex.w != 0)


def test_ini_1d_runs():
    flex = _run_from_config(INPUT_DIR / "input_f1d")
    assert flex.w.ndim == 1
    assert not np.any(np.isnan(flex.w))
    assert np.any(flex.w != 0)


def test_ini_2d_runs():
    flex = _run_from_config(INPUT_DIR / "input_f2d")
    assert flex.w.ndim == 2
    assert not np.any(np.isnan(flex.w))
    assert np.any(flex.w != 0)


# ---------------------------------------------------------------------------
# Equivalence tests — YAML and INI for the same run must give identical w
# ---------------------------------------------------------------------------

def test_yaml_ini_equivalence_1d():
    yaml_flex = _run_from_config(INPUT_DIR / "input_f1d.yaml")
    ini_flex = _run_from_config(INPUT_DIR / "input_f1d")
    np.testing.assert_array_equal(yaml_flex.w, ini_flex.w)


def test_yaml_ini_equivalence_2d():
    yaml_flex = _run_from_config(INPUT_DIR / "input_f2d.yaml")
    ini_flex = _run_from_config(INPUT_DIR / "input_f2d")
    np.testing.assert_array_equal(yaml_flex.w, ini_flex.w)
