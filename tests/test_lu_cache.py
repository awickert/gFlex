"""Tests for the LU factorization cache (cache_factorization attribute)."""

import numpy as np
import pytest

from gflex.base import _matrix_hash
from gflex.f1d import F1D
from gflex.f2d import F2D


def _make_f1d(cache=False):
    flex = F1D()
    flex.quiet = True
    flex.method = "fd"
    flex.bc_west = "zero_moment_zero_shear"
    flex.bc_east = "zero_moment_zero_shear"
    flex.E = 65e9
    flex.nu = 0.25
    flex.rho_m = 3300.0
    flex.rho_fill = 0.0
    flex.g = 9.8
    flex.T_e = 35e3
    flex.dx = 10e3
    nx = 100
    flex.qs = np.zeros(nx)
    flex.qs[nx // 2] = 1e6
    flex.cache_factorization = cache
    flex.initialize()
    return flex


def _make_f2d(cache=False):
    flex = F2D()
    flex.quiet = True
    flex.method = "fd"
    flex.bc_west = "zero_moment_zero_shear"
    flex.bc_east = "zero_moment_zero_shear"
    flex.bc_north = "zero_moment_zero_shear"
    flex.bc_south = "zero_moment_zero_shear"
    flex.E = 65e9
    flex.nu = 0.25
    flex.rho_m = 3300.0
    flex.rho_fill = 0.0
    flex.g = 9.8
    flex.T_e = 35e3
    flex.dx = 10e3
    flex.dy = 10e3
    ny, nx = 40, 40
    flex.qs = np.zeros((ny, nx))
    flex.qs[ny // 2, nx // 2] = 1e6
    flex.cache_factorization = cache
    flex.initialize()
    return flex


# ── correctness: cached result matches uncached ───────────────────────────────
# factorized() and spsolve() use different UMFpack entry points, so allow
# a small tolerance (~1e-9 relative) beyond strict floating-point identity.

def test_1d_cache_true_matches_uncached():
    ref = _make_f1d(cache=False)
    ref.run()
    cached = _make_f1d(cache=True)
    cached.run()
    np.testing.assert_allclose(cached.w, ref.w, rtol=1e-9)


def test_1d_cache_no_check_matches_uncached():
    ref = _make_f1d(cache=False)
    ref.run()
    cached = _make_f1d(cache="no_check")
    cached.run()
    np.testing.assert_allclose(cached.w, ref.w, rtol=1e-9)


def test_2d_cache_true_matches_uncached():
    ref = _make_f2d(cache=False)
    ref.run()
    cached = _make_f2d(cache=True)
    cached.run()
    np.testing.assert_allclose(cached.w, ref.w, rtol=1e-9)


def test_2d_cache_no_check_matches_uncached():
    ref = _make_f2d(cache=False)
    ref.run()
    cached = _make_f2d(cache="no_check")
    cached.run()
    np.testing.assert_allclose(cached.w, ref.w, rtol=1e-9)


# ── reuse: second run() with same matrix reuses the cached factorization ──────

def test_1d_cache_reused_on_second_run():
    flex = _make_f1d(cache=True)
    flex.run()
    lu_first = flex._lu
    flex.run()
    assert flex._lu is lu_first, "LU callable should be reused when matrix is unchanged"


def test_1d_no_check_reused_on_second_run():
    flex = _make_f1d(cache="no_check")
    flex.run()
    lu_first = flex._lu
    flex.run()
    assert flex._lu is lu_first


# ── invalidation: changing the load (qs only) does NOT refactorize ────────────

def test_1d_cache_load_change_reuses_lu():
    """Changing qs does not affect the coefficient matrix — LU must be reused."""
    flex = _make_f1d(cache=True)
    flex.run()
    w_first = flex.w.copy()
    lu_first = flex._lu

    flex.qs[len(flex.qs) // 4] += 5e5
    flex.run()

    assert flex._lu is lu_first, "LU should be reused; only qs changed"
    assert not np.allclose(flex.w, w_first), "deflection must change when load changes"


# ── invalidation: stale hash triggers refactorization ────────────────────────

def test_1d_cache_stale_hash_refactorizes():
    """Corrupting the stored hash forces a refactorization on the next run()."""
    flex = _make_f1d(cache=True)
    flex.run()
    lu_first = flex._lu

    # Simulate a stale cache (e.g. Te changed and matrix was rebuilt externally)
    flex._lu_matrix_hash = b"\x00" * 16

    flex.run()

    assert flex._lu is not lu_first, "stale hash must trigger refactorization"
    assert flex._lu_matrix_hash == _matrix_hash(flex.coeff_matrix), \
        "hash must be updated to the current matrix after refactorization"


# ── invalid value raises ──────────────────────────────────────────────────────

def test_invalid_cache_factorization_raises():
    flex = _make_f1d(cache="bad_value")
    with pytest.raises(ValueError, match="cache_factorization"):
        flex.run()


# ── finalize clears the cache ─────────────────────────────────────────────────

def test_finalize_clears_lu_cache():
    flex = _make_f1d(cache=True)
    flex.run()
    assert flex._lu is not None
    flex.finalize()
    assert not hasattr(flex, "_lu")
    assert not hasattr(flex, "_lu_matrix_hash")


# ── no_check free-solo: coeff_matrix freed after factorization ────────────────

def test_1d_no_check_coeff_matrix_freed_after_factorization():
    """coeff_matrix is released to None immediately after the LU is built."""
    flex = _make_f1d(cache="no_check")
    flex.run()
    assert flex.coeff_matrix is None, "coeff_matrix should be freed after factorization"
    assert flex._lu is not None, "_lu must be retained"


def test_2d_no_check_coeff_matrix_freed_after_factorization():
    flex = _make_f2d(cache="no_check")
    flex.run()
    assert flex.coeff_matrix is None
    assert flex._lu is not None


def test_no_check_second_run_reuses_lu_and_keeps_coeff_matrix_none():
    """Second run reuses _lu without rebuilding the matrix."""
    flex = _make_f1d(cache="no_check")
    flex.run()
    lu_first = flex._lu

    flex.qs[len(flex.qs) // 4] += 5e5
    flex.run()

    assert flex._lu is lu_first, "_lu must be reused on second run"
    assert flex.coeff_matrix is None, "coeff_matrix must remain None after second run"


def test_no_check_te_change_triggers_full_rebuild():
    """Changing te invalidates the freed coeff_matrix and forces a full rebuild."""
    flex = _make_f1d(cache="no_check")
    flex.run()
    w_first = flex.w.copy()

    flex.T_e = 20e3  # different Te → smart invalidation fires
    assert flex._lu is None, "_lu must be cleared when T_e changes"
    flex.run()

    assert flex.coeff_matrix is None, "coeff_matrix freed again after rebuild"
    assert flex._lu is not None, "_lu must be repopulated after rebuild"
    assert not np.allclose(flex.w, w_first), "deflection must differ for different Te"
