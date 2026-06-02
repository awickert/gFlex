#! /usr/bin/env python
"""Tests for YAML configuration file loading."""

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


def test_missing_required_key_raises():
    """A required config key missing from the config raises ValueError, not SystemExit."""
    flex = F1D()
    flex.config = {"input": {}}
    # ElasticThickness is absent — configGet with optional=False must raise ValueError
    with pytest.raises(ValueError, match="ElasticThickness"):
        flex.configGet("float", "input", "ElasticThickness", optional=False)
