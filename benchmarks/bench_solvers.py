"""
Solver benchmarks for gFlex.

Run with:
    python benchmarks/bench_solvers.py

Prints timing tables to stdout; redirect to a file to keep a record.

FD benchmarks cover:
  - direct sparse LU vs ILU-preconditioned LGMRES (iterative)
  - three Te profiles: constant, sinusoidally varying, abrupt step

This makes it easy to measure the benefit of an FFT preconditioner
(constant-Te FFT solve used as M for the iterative solve) relative to
the ILU baseline.

Grid-size notes:
  - SAS:  O(N²) in 1D, O(N⁴) in 2D — kept small to avoid long runs.
  - FFT:  O(N log N) — handles large grids easily.
  - FD direct:  roughly O(N^1.5) for the sparse LU factorisation.
  - FD iterative:  problem-dependent; may be slower than direct for small N.
"""

import time

import numpy as np

from gflex.f1d import F1D
from gflex.f2d import F2D


# ── physical constants ────────────────────────────────────────────────────────

_COMMON = dict(g=9.8, E=65e9, nu=0.25, rho_m=3300.0, rho_fill=0.0)
_TE_REF = 35000.0     # reference elastic thickness [m]


# ── Te profiles ───────────────────────────────────────────────────────────────

def _te_1d(n, profile, dx=5000.0):
    """Return a Te scalar or array for a 1-D grid of n points."""
    if profile == "constant":
        return _TE_REF
    x = np.arange(n) * dx
    L = x[-1]
    if profile == "sinusoidal":
        # One full cycle over the domain; amplitude ±50 % of reference
        return _TE_REF * (1.0 + 0.5 * np.sin(2.0 * np.pi * x / L))
    if profile == "abrupt":
        # Step: west half thin (½ × ref), east half thick (1½ × ref)
        Te = np.full(n, 0.5 * _TE_REF)
        Te[n // 2:] = 1.5 * _TE_REF
        return Te
    raise ValueError(f"unknown Te profile: {profile!r}")


def _te_2d(ny, nx, profile, dx=5000.0, dy=5000.0):
    """Return a Te scalar or array for a 2-D grid of (ny, nx) points."""
    if profile == "constant":
        return _TE_REF
    x = np.arange(nx) * dx
    y = np.arange(ny) * dy
    xx, yy = np.meshgrid(x, y)
    Lx, Ly = x[-1], y[-1]
    if profile == "sinusoidal":
        # One full cycle in each direction; amplitude ±50 % of reference
        return _TE_REF * (1.0 + 0.5 * np.sin(2.0 * np.pi * xx / Lx)
                                      * np.sin(2.0 * np.pi * yy / Ly))
    if profile == "abrupt":
        # NW quadrant thin, SE quadrant thick, others at reference
        Te = np.full((ny, nx), _TE_REF)
        Te[:ny // 2, :nx // 2] = 0.5 * _TE_REF
        Te[ny // 2:, nx // 2:] = 1.5 * _TE_REF
        return Te
    raise ValueError(f"unknown Te profile: {profile!r}")


# ── object builders ───────────────────────────────────────────────────────────

def _make_f1d(n, method, te=_TE_REF, solver="direct", bc="0Displacement0Slope"):
    flex = F1D()
    flex.Quiet = True
    flex.Method = method
    flex.Solver = solver
    flex.dx = 5000.0
    flex.Te = te
    flex.qs = np.zeros(n)
    flex.qs[n // 4 : 3 * n // 4] = 1e6
    flex.BC_W = bc
    flex.BC_E = bc
    for k, v in _COMMON.items():
        setattr(flex, k, v)
    flex.initialize()
    return flex


def _make_f2d(nx, ny, method, te=_TE_REF, solver="direct", bc="0Displacement0Slope"):
    flex = F2D()
    flex.Quiet = True
    flex.Method = method
    flex.Solver = solver
    flex.dx = 5000.0
    flex.dy = 5000.0
    flex.Te = te
    flex.qs = np.zeros((ny, nx))
    flex.qs[ny // 4 : 3 * ny // 4, nx // 4 : 3 * nx // 4] = 1e6
    flex.BC_W = bc
    flex.BC_E = bc
    flex.BC_N = bc
    flex.BC_S = bc
    for k, v in _COMMON.items():
        setattr(flex, k, v)
    flex.initialize()
    return flex


# ── formatting helpers ────────────────────────────────────────────────────────

def _tick():
    return time.perf_counter()


def _hdr(cols):
    print("  ".join(f"{c:>{w}}" for c, w in cols))
    print("  ".join("-" * w for _, w in cols))


# ── FD benchmarks: direct vs iterative, three Te profiles ────────────────────

_TE_PROFILES = ("constant", "sinusoidal", "abrupt")


def bench_1d_fd(sizes, iter_max=2000, profiles=_TE_PROFILES):
    """Benchmark 1-D FD.

    Iterative solve is skipped for n > iter_max: the ILU-preconditioned
    LGMRES converges reliably up to ~2000 nodes but offers little
    advantage over the direct sparse LU at that scale.
    """
    print("\n1D FD  (direct vs iterative, three Te profiles)")
    cols = [
        ("n", 7), ("Te profile", 12),
        ("assemble(s)", 12), ("direct(s)", 10), ("iterative(s)", 13), ("iter/direct", 11),
    ]
    _hdr(cols)
    for n in sizes:
        for prof in profiles:
            te = _te_1d(n, prof)
            flex = _make_f1d(n, "FD", te)
            flex.bc_check()
            flex.gridded_x()

            t0 = _tick()
            flex.elasprepFD()
            flex.BC_selector_and_coeff_matrix_creator()
            t_asm = _tick() - t0

            flex.Solver = "direct"
            t0 = _tick()
            flex.fd_solve()
            t_direct = _tick() - t0

            if n <= iter_max:
                flex.Solver = "iterative"
                t0 = _tick()
                flex.fd_solve()
                t_iter = _tick() - t0
                ratio = t_iter / t_direct if t_direct > 1e-9 else float("nan")
                iter_str = f"{t_iter:>13.4f}"
                ratio_str = f"{ratio:>11.2f}"
            else:
                iter_str = f"{'—':>13}"
                ratio_str = f"{'—':>11}"

            print(f"  {n:>7}  {prof:>12}  {t_asm:>12.4f}"
                  f"  {t_direct:>10.4f}  {iter_str}  {ratio_str}")
        if n != sizes[-1]:
            print()


def bench_2d_fd(sizes, iter_max=100, profiles=_TE_PROFILES):
    """Benchmark 2-D FD.

    Iterative solve is skipped for n > iter_max: the biharmonic
    operator's condition number grows as O(N⁴), so ILU-preconditioned
    LGMRES can be very slow or non-convergent for large grids.  The cap
    keeps the benchmark run to a reasonable wall time; raise iter_max
    to explore larger sizes once a better preconditioner is in place.
    """
    print("\n2D FD  (direct vs iterative, three Te profiles)")
    cols = [
        ("n×n", 9), ("Te profile", 12),
        ("assemble(s)", 12), ("direct(s)", 10), ("iterative(s)", 13), ("iter/direct", 11),
    ]
    _hdr(cols)
    for n in sizes:
        label = f"{n}×{n}"
        for prof in profiles:
            te = _te_2d(n, n, prof)
            flex = _make_f2d(n, n, "FD", te)
            flex.bc_check()

            t0 = _tick()
            flex.elasprep()
            flex.BC_selector_and_coeff_matrix_creator()
            t_asm = _tick() - t0

            flex.Solver = "direct"
            t0 = _tick()
            flex.fd_solve()
            t_direct = _tick() - t0

            if n <= iter_max:
                flex.Solver = "iterative"
                t0 = _tick()
                flex.fd_solve()
                t_iter = _tick() - t0
                ratio = t_iter / t_direct if t_direct > 1e-9 else float("nan")
                iter_str = f"{t_iter:>13.4f}"
                ratio_str = f"{ratio:>11.2f}"
            else:
                iter_str = f"{'—':>13}"
                ratio_str = f"{'—':>11}"

            print(f"  {label:>9}  {prof:>12}  {t_asm:>12.4f}"
                  f"  {t_direct:>10.4f}  {iter_str}  {ratio_str}")
        if n != sizes[-1]:
            print()


# ── FFT and SAS benchmarks (constant Te only) ─────────────────────────────────

def bench_1d_fft(sizes):
    print("\n1D FFT")
    cols = [("n", 8), ("total (s)", 12)]
    _hdr(cols)
    for n in sizes:
        flex = _make_f1d(n, "FFT")
        t0 = _tick()
        flex.run()
        print("  ".join([f"{n:>8}", f"{_tick() - t0:>12.4f}"]))


def bench_1d_sas(sizes):
    print("\n1D SAS")
    cols = [("n", 8), ("total (s)", 12)]
    _hdr(cols)
    for n in sizes:
        flex = _make_f1d(n, "SAS")
        t0 = _tick()
        flex.run()
        print("  ".join([f"{n:>8}", f"{_tick() - t0:>12.4f}"]))


def bench_2d_fft(sizes):
    print("\n2D FFT")
    cols = [("n×n", 9), ("total (s)", 12)]
    _hdr(cols)
    for n in sizes:
        label = f"{n}×{n}"
        flex = _make_f2d(n, n, "FFT")
        t0 = _tick()
        flex.run()
        print("  ".join([f"{label:>9}", f"{_tick() - t0:>12.4f}"]))


def bench_2d_sas(sizes):
    print("\n2D SAS  (O(N⁴) Python loops — kept small)")
    cols = [("n×n", 9), ("total (s)", 12)]
    _hdr(cols)
    for n in sizes:
        label = f"{n}×{n}"
        flex = _make_f2d(n, n, "SAS")
        t0 = _tick()
        flex.run()
        print("  ".join([f"{label:>9}", f"{_tick() - t0:>12.4f}"]))


# ── entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("gFlex solver benchmarks")
    print("=" * 70)

    print("\n--- 1D ---")
    bench_1d_fd(sizes=[100, 500, 2000, 5000])
    bench_1d_fft(sizes=[100, 500, 2000, 5000, 20000])
    bench_1d_sas(sizes=[100, 500, 2000, 5000])

    print("\n--- 2D ---")
    bench_2d_fd(sizes=[50, 100, 200, 400])
    bench_2d_fft(sizes=[50, 100, 500, 1000])
    bench_2d_sas(sizes=[10, 25, 50, 100])
