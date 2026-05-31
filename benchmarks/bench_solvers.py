"""
Solver benchmarks for gFlex.

Run with:
    python benchmarks/bench_solvers.py

Each section prints a table with timing in seconds.  Run before and after
optimisation work to measure improvement.  Results are not saved
automatically; redirect stdout if you want a record.

Grid sizes are chosen so each solver finishes in a reasonable wall time:
- SAS is O(N^2) in 1D and O(N^4) in 2D (Python loops), so kept small.
- FFT scales as O(N log N) and can handle large grids quickly.
- FD scales roughly as O(N^1.5) for the direct sparse solve.
"""

import time

import numpy as np

from gflex.f1d import F1D
from gflex.f2d import F2D


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_COMMON = dict(
    g=9.8,
    E=65e9,
    nu=0.25,
    rho_m=3300.0,
    rho_fill=0.0,
)


def _make_f1d(n, method, bc="0Displacement0Slope"):
    flex = F1D()
    flex.Quiet = True
    flex.Method = method
    flex.Solver = "direct"
    flex.dx = 5000.0
    flex.Te = 35000.0
    flex.qs = np.zeros(n)
    flex.qs[n // 4 : 3 * n // 4] = 1e6
    flex.BC_W = bc
    flex.BC_E = bc
    for k, v in _COMMON.items():
        setattr(flex, k, v)
    flex.initialize()
    return flex


def _make_f2d(nx, ny, method, bc="0Displacement0Slope"):
    flex = F2D()
    flex.Quiet = True
    flex.Method = method
    flex.Solver = "direct"
    flex.dx = 5000.0
    flex.dy = 5000.0
    flex.Te = 35000.0
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


def _tick():
    return time.perf_counter()


def _hdr(cols):
    print("  ".join(f"{c:>{w}}" for c, w in cols))
    print("  ".join("-" * w for _, w in cols))


# ---------------------------------------------------------------------------
# 1D benchmarks
# ---------------------------------------------------------------------------

def bench_1d_fd(sizes):
    print("\n1D FD")
    cols = [("n", 8), ("assembly (s)", 14), ("solve (s)", 12), ("total (s)", 12)]
    _hdr(cols)
    for n in sizes:
        flex = _make_f1d(n, "FD")
        flex.bc_check()   # initialises coeff_matrix and BC defaults
        flex.gridded_x()

        t0 = _tick()
        flex.elasprepFD()
        flex.BC_selector_and_coeff_matrix_creator()
        t_assemble = _tick() - t0

        t0 = _tick()
        flex.fd_solve()
        t_solve = _tick() - t0

        row = [n, t_assemble, t_solve, t_assemble + t_solve]
        print("  ".join(f"{v:>{w}}" if i == 0 else f"{v:>{w}.4f}"
                        for i, ((_, w), v) in enumerate(zip(cols, row))))


def bench_1d_fft(sizes):
    print("\n1D FFT")
    cols = [("n", 8), ("total (s)", 12)]
    _hdr(cols)
    for n in sizes:
        flex = _make_f1d(n, "FFT")
        t0 = _tick()
        flex.run()
        t = _tick() - t0
        print("  ".join([f"{n:>8}", f"{t:>12.4f}"]))


def bench_1d_sas(sizes):
    print("\n1D SAS")
    cols = [("n", 8), ("total (s)", 12)]
    _hdr(cols)
    for n in sizes:
        flex = _make_f1d(n, "SAS")
        t0 = _tick()
        flex.run()
        t = _tick() - t0
        print("  ".join([f"{n:>8}", f"{t:>12.4f}"]))


# ---------------------------------------------------------------------------
# 2D benchmarks
# ---------------------------------------------------------------------------

def bench_2d_fd(sizes):
    print("\n2D FD")
    cols = [("n×n", 10), ("assembly (s)", 14), ("solve (s)", 12), ("total (s)", 12)]
    _hdr(cols)
    for n in sizes:
        flex = _make_f2d(n, n, "FD")
        flex.bc_check()   # initialises coeff_matrix and BC defaults

        t0 = _tick()
        flex.elasprep()
        flex.BC_selector_and_coeff_matrix_creator()
        t_assemble = _tick() - t0

        t0 = _tick()
        flex.fd_solve()
        t_solve = _tick() - t0

        label = f"{n}×{n}"
        print(f"  {label:>10}  {t_assemble:>14.4f}  {t_solve:>12.4f}"
              f"  {t_assemble + t_solve:>12.4f}")


def bench_2d_fft(sizes):
    print("\n2D FFT")
    cols = [("n×n", 10), ("total (s)", 12)]
    _hdr(cols)
    for n in sizes:
        flex = _make_f2d(n, n, "FFT")
        label = f"{n}×{n}"
        t0 = _tick()
        flex.run()
        t = _tick() - t0
        print("  ".join([f"{label:>10}", f"{t:>12.4f}"]))


def bench_2d_sas(sizes):
    print("\n2D SAS  (O(N^4) Python loops — kept small)")
    cols = [("n×n", 10), ("total (s)", 12)]
    _hdr(cols)
    for n in sizes:
        flex = _make_f2d(n, n, "SAS")
        label = f"{n}×{n}"
        t0 = _tick()
        flex.run()
        t = _tick() - t0
        print("  ".join([f"{label:>10}", f"{t:>12.4f}"]))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("gFlex solver benchmarks")
    print("=" * 56)

    print("\n--- 1D ---")
    bench_1d_fd(sizes=[100, 500, 2000, 5000])
    bench_1d_fft(sizes=[100, 500, 2000, 5000, 20000])
    bench_1d_sas(sizes=[100, 500, 2000, 5000])

    print("\n--- 2D ---")
    bench_2d_fd(sizes=[50, 100, 200, 400])
    bench_2d_fft(sizes=[50, 100, 500, 1000])
    bench_2d_sas(sizes=[10, 25, 50, 100])
