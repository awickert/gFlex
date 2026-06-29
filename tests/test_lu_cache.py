"""Tests for the LU factorization cache (``cache_factorization`` attribute).

After the cache redesign:

* ``cache_factorization=False`` (default) — no LU cache; ``spsolve`` each run.
* ``cache_factorization=True`` — cache the LU factorisation and reuse it,
  trusting setter-based invalidation (matrix-determining inputs invalidate the
  cache on reassignment; array inputs are read-only).  The coefficient matrix
  is freed once the LU is built.  There is no per-run matrix hash.
* ``cache_factorization="no_check"`` — deprecated alias for ``True``.
"""

import numpy as np
import pytest

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


def test_2d_cache_true_matches_uncached():
    ref = _make_f2d(cache=False)
    ref.run()
    cached = _make_f2d(cache=True)
    cached.run()
    np.testing.assert_allclose(cached.w, ref.w, rtol=1e-9)


# ── reuse: second run() with same matrix reuses the cached factorization ──────

def test_1d_cache_reused_on_second_run():
    flex = _make_f1d(cache=True)
    flex.run()
    lu_first = flex._lu
    flex.run()
    assert flex._lu is lu_first, "LU callable should be reused when matrix is unchanged"


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


# ── cached mode frees the coefficient matrix once the LU is built ─────────────

def test_1d_cache_true_frees_coeff_matrix():
    flex = _make_f1d(cache=True)
    flex.run()
    assert flex.coeff_matrix is None, "coeff_matrix should be freed after factorization"
    assert flex._lu is not None, "_lu must be retained"


def test_2d_cache_true_frees_coeff_matrix():
    flex = _make_f2d(cache=True)
    flex.run()
    assert flex.coeff_matrix is None
    assert flex._lu is not None


def test_cache_true_second_run_keeps_coeff_matrix_none():
    flex = _make_f1d(cache=True)
    flex.run()
    lu_first = flex._lu
    flex.qs[len(flex.qs) // 4] += 5e5
    flex.run()
    assert flex._lu is lu_first, "_lu must be reused on second run"
    assert flex.coeff_matrix is None, "coeff_matrix must remain None after second run"


# ── invalidation: a Te change clears the cache and forces a rebuild ───────────

def test_cache_te_change_triggers_rebuild():
    flex = _make_f1d(cache=True)
    flex.run()
    w_first = flex.w.copy()

    flex.T_e = 20e3  # different Te → setter invalidation fires
    assert flex._lu is None, "_lu must be cleared when T_e changes"
    flex.run()

    assert flex.coeff_matrix is None, "coeff_matrix freed again after rebuild"
    assert flex._lu is not None, "_lu must be repopulated after rebuild"
    assert not np.allclose(flex.w, w_first), "deflection must differ for different Te"


# ── 'no_check' is a deprecated alias for True ─────────────────────────────────

def test_no_check_is_deprecated_alias_for_true():
    flex = _make_f1d(cache="no_check")
    with pytest.warns(DeprecationWarning, match="no_check"):
        flex.run()
    assert flex.cache_factorization is True
    assert flex.coeff_matrix is None  # behaves as True: matrix freed
    assert flex._lu is not None


def test_no_check_result_matches_uncached():
    ref = _make_f1d(cache=False)
    ref.run()
    cached = _make_f1d(cache="no_check")
    with pytest.warns(DeprecationWarning):
        cached.run()
    np.testing.assert_allclose(cached.w, ref.w, rtol=1e-9)


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
