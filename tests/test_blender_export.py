#! /usr/bin/env python
"""Tests for export_for_blender()."""

import numpy as np
import pytest

from gflex import export_for_blender


def _run(tmp_path, w, dx=5000.0, dy=None, **kwargs):
    """Write mesh to tmp_path and exec() it; return the resulting namespace."""
    path = str(tmp_path / "mesh.py")
    export_for_blender(w, dx=dx, dy=dy, path=path, **kwargs)
    ns = {}
    with open(path) as fh:
        exec(fh.read(), ns)
    return ns


@pytest.fixture
def simple_w():
    """4×6 deflection array with a central depression."""
    w = np.zeros((4, 6))
    w[1:3, 2:4] = -500.0
    return w


def test_output_is_valid_python(tmp_path, simple_w):
    """The written file must exec() without error."""
    _run(tmp_path, simple_w)


def test_grid_dimensions(tmp_path, simple_w):
    ny, nx = simple_w.shape
    ns = _run(tmp_path, simple_w)
    assert ns["nrows"] == ny
    assert ns["ncols"] == nx


def test_vertex_and_face_counts(tmp_path, simple_w):
    ny, nx = simple_w.shape
    ns = _run(tmp_path, simple_w)
    assert len(ns["vertices"]) == ny * nx
    assert len(ns["faces"]) == (ny - 1) * (nx - 1)


def test_w_norm_length(tmp_path, simple_w):
    ny, nx = simple_w.shape
    ns = _run(tmp_path, simple_w)
    assert len(ns["w_norm"]) == ny * nx


def test_deflection_extremes_sign(tmp_path, simple_w):
    """w_min_bu should be negative (downward deflection present)."""
    ns = _run(tmp_path, simple_w)
    assert ns["w_min_bu"] < 0.0
    assert ns["w_max_bu"] == pytest.approx(0.0)


def test_qs_attributes_written(tmp_path, simple_w):
    ny, nx = simple_w.shape
    qs = np.zeros_like(simple_w)
    qs[1:3, 2:4] = 1e6
    ns = _run(tmp_path, simple_w, qs=qs)
    assert "qs_verts" in ns
    assert "qs_norm" in ns
    assert len(ns["qs_verts"]) == ny * nx
    assert max(ns["qs_norm"]) == pytest.approx(1.0)


def test_qs_absent_by_default(tmp_path, simple_w):
    ns = _run(tmp_path, simple_w)
    assert "qs_verts" not in ns
    assert "qs_norm" not in ns


def test_te_attributes_written_scalar(tmp_path, simple_w):
    ns = _run(tmp_path, simple_w, Te=35e3)
    assert "te_verts" in ns
    assert ns["te_min"] == pytest.approx(35e3)
    assert ns["te_max"] == pytest.approx(35e3)


def test_te_attributes_written_array(tmp_path, simple_w):
    Te = np.full_like(simple_w, 30e3)
    Te[0, :] = 40e3
    ns = _run(tmp_path, simple_w, Te=Te)
    assert ns["te_min"] == pytest.approx(30e3)
    assert ns["te_max"] == pytest.approx(40e3)


def test_te_absent_by_default(tmp_path, simple_w):
    ns = _run(tmp_path, simple_w)
    assert "te_verts" not in ns


def test_dy_defaults_to_dx(tmp_path, simple_w):
    dx = 4000.0
    ns_square = _run(tmp_path, simple_w, dx=dx)
    ns_explicit = _run(tmp_path, simple_w, dx=dx, dy=dx)
    assert ns_square["vertices"] == ns_explicit["vertices"]
