"""
Solver benchmarks for gFlex.

Run with:
    python benchmarks/bench_solvers.py

Results are printed to stdout and saved to benchmarks/results/ as a text
file named {short_commit_hash}_{UTC_timestamp}.txt.  The file header
records the full commit hash, branch, whether the working tree is clean,
and key hardware/software details so runs are reproducible and comparable.

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

import os
import platform
import signal
import subprocess
import sys
import time

import numpy as np
import scipy
import scipy.sparse
from scipy.sparse.linalg import lgmres, spilu, LinearOperator

from gflex.f1d import F1D
from gflex.f2d import F2D


# ── provenance: git and system info ───────────────────────────────────────────

def _git_info():
    """Return dict with commit hash, branch, and dirty-flag."""
    def _run(*args):
        try:
            return subprocess.check_output(args, text=True,
                                           stderr=subprocess.DEVNULL).strip()
        except Exception:
            return "unknown"

    full  = _run("git", "rev-parse", "HEAD")
    short = _run("git", "rev-parse", "--short", "HEAD")
    branch = _run("git", "rev-parse", "--abbrev-ref", "HEAD")
    dirty = bool(_run("git", "status", "--porcelain"))
    return {"full": full, "short": short, "branch": branch, "dirty": dirty}


def _cpu_model():
    """Best-effort CPU model string, Linux and macOS."""
    try:
        with open("/proc/cpuinfo") as f:
            for line in f:
                if line.startswith("model name"):
                    return line.split(":", 1)[1].strip()
    except OSError:
        pass
    try:
        return subprocess.check_output(
            ["sysctl", "-n", "machdep.cpu.brand_string"],
            text=True, stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        pass
    return platform.processor() or "unknown"


def _ram_gb():
    """Total RAM in GiB, Linux and macOS."""
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemTotal:"):
                    return int(line.split()[1]) / (1024 ** 2)
    except OSError:
        pass
    try:
        b = subprocess.check_output(
            ["sysctl", "-n", "hw.memsize"],
            text=True, stderr=subprocess.DEVNULL,
        ).strip()
        return int(b) / (1024 ** 3)
    except Exception:
        pass
    return None


def _header(git, timestamp):
    """Format a provenance block for the top of the results file."""
    dirty = " (dirty)" if git["dirty"] else ""
    ram = _ram_gb()
    ram_str = f"{ram:.1f} GiB" if ram is not None else "unknown"
    lines = [
        "gFlex solver benchmarks",
        "=" * 70,
        f"timestamp : {timestamp}",
        f"commit    : {git['full']}{dirty}",
        f"branch    : {git['branch']}",
        f"cpu       : {_cpu_model()}",
        f"cpu cores : {os.cpu_count()}",
        f"ram       : {ram_str}",
        f"os        : {platform.platform()}",
        f"python    : {sys.version.split()[0]}",
        f"numpy     : {np.__version__}",
        f"scipy     : {scipy.__version__}",
    ]
    return "\n".join(lines)


def _push_result(local_path, hostname):
    """Upload result file to awickert/gflex-benchmarks via gh CLI.

    Path in the repo: results/{hostname}/{filename}.
    Requires `gh` to be installed and authenticated; fails gracefully if not.
    """
    import base64

    filename    = os.path.basename(local_path)
    remote_path = f"results/{hostname}/{filename}"
    with open(local_path, "rb") as f:
        content = base64.b64encode(f.read()).decode()

    result = subprocess.run(
        [
            "gh", "api",
            f"repos/awickert/gflex-benchmarks/contents/{remote_path}",
            "--method", "PUT",
            "--field", f"message=Add benchmark result from {hostname}",
            "--field", f"content={content}",
        ],
        capture_output=True, text=True,
    )
    if result.returncode == 0:
        print(
            f"Pushed to github.com/awickert/gflex-benchmarks/{remote_path}",
            file=sys.__stdout__,
        )
    else:
        print(
            f"Warning: push to gflex-benchmarks failed:\n{result.stderr.strip()}",
            file=sys.__stdout__,
        )


class _Tee:
    """Mirror writes to stdout and a file simultaneously."""

    def __init__(self, path):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        self._f = open(path, "w")
        self._stdout = sys.stdout
        sys.stdout = self

    def write(self, s):
        self._stdout.write(s)
        self._f.write(s)

    def flush(self):
        self._stdout.flush()
        self._f.flush()

    def close(self):
        sys.stdout = self._stdout
        self._f.close()


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
    if profile == "gradient":
        # Linear ramp: ½ × ref at west edge, 1½ × ref at east edge
        return _TE_REF * (0.5 + np.linspace(0.0, 1.0, n))
    if profile == "wide_range":
        # 10:1 ramp: 10 km west, 100 km east (craton–ocean range)
        return np.linspace(10000.0, 100000.0, n)
    if profile == "noisy_mild":
        # White noise ±20 % of reference, seeded for reproducibility
        rng = np.random.default_rng(42)
        return _TE_REF * (1.0 + 0.2 * rng.uniform(-1.0, 1.0, n))
    if profile == "noisy_strong":
        # White noise ±50 % of reference
        rng = np.random.default_rng(42)
        return _TE_REF * (1.0 + 0.5 * rng.uniform(-1.0, 1.0, n))
    raise ValueError(f"unknown Te profile: {profile!r}")


def _te_2d(ny, nx, profile, dx=5000.0, dy=5000.0):
    """Return a Te scalar or array for a 2-D grid of (ny, nx) points.

    Profiles ending in ``_45`` vary along the diagonal direction
    (x + y) / (Lx + Ly), which maximises the cross-derivative coupling
    term (1 - ν) ∂²D/∂x∂y · ∂²w/∂x∂y in the variable-D 2-D equation.
    """
    if profile == "constant":
        return _TE_REF
    x = np.arange(nx) * dx
    y = np.arange(ny) * dy
    xx, yy = np.meshgrid(x, y)
    Lx, Ly = x[-1], y[-1]

    # Axis-aligned profiles
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
    if profile == "gradient":
        # Linear ramp along x: ½ × ref at west, 1½ × ref at east
        return _TE_REF * (0.5 + xx / Lx)

    # 45-degree profiles: variation along normalised diagonal r ∈ [0, 1]
    r = (xx / Lx + yy / Ly) / 2.0   # 0 at (0,0), 1 at (Lx,Ly)
    if profile == "sinusoidal_45":
        return _TE_REF * (1.0 + 0.5 * np.sin(2.0 * np.pi * r))
    if profile == "abrupt_45":
        # Step across the main diagonal
        Te = np.where(r < 0.5, 0.5 * _TE_REF, 1.5 * _TE_REF)
        return Te
    if profile == "gradient_45":
        # Linear ramp along diagonal: ½ × ref at (0,0), 1½ × ref at (Lx,Ly)
        return _TE_REF * (0.5 + r)

    # Wide dynamic range
    if profile == "wide_range":
        # 10:1 ramp along x: 10 km west, 100 km east
        return np.linspace(10000.0, 100000.0, nx)[np.newaxis, :] * np.ones((ny, 1))
    if profile == "wide_range_45":
        return np.linspace(10000.0, 100000.0, 1)[0] + r * (100000.0 - 10000.0)

    # Noisy Te (spatially uncorrelated white noise, seeded for reproducibility)
    rng = np.random.default_rng(42)
    if profile == "noisy_mild":
        return _TE_REF * (1.0 + 0.2 * rng.uniform(-1.0, 1.0, (ny, nx)))
    if profile == "noisy_strong":
        return _TE_REF * (1.0 + 0.5 * rng.uniform(-1.0, 1.0, (ny, nx)))

    # Circular inclusion centred in domain (radius = Lx / 5)
    cx, cy = Lx / 2.0, Ly / 2.0
    dist = np.sqrt((xx - cx) ** 2 + (yy - cy) ** 2)
    mask = dist < Lx / 5.0
    if profile == "disk_thin":
        # Thin disk (0.3 × ref) in reference background
        Te = np.full((ny, nx), _TE_REF)
        Te[mask] = 0.3 * _TE_REF
        return Te
    if profile == "disk_thick":
        # Thick disk (3 × ref) — rigid craton or seamount
        Te = np.full((ny, nx), _TE_REF)
        Te[mask] = 3.0 * _TE_REF
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


# ── iterative solver with diagnostics ────────────────────────────────────────

def _iter_solve(matrix, rhs, tol=1e-3, maxiter=200):
    """ILU-preconditioned LGMRES with iteration counting.

    Returns (elapsed_s, n_iters, w, converged).
    Mirrors the preconditioner logic in gFlex's fd_solve() so timings
    are directly comparable.

    maxiter caps the iteration count so that poorly-conditioned large
    problems (the 2D biharmonic condition number grows as O(N⁴)) do not
    run indefinitely.  Non-convergence is flagged with a "!" in the table.
    """
    iters = [0]

    def _cb(_r):
        iters[0] += 1

    try:
        ilu = spilu(matrix.tocsc(), fill_factor=20, drop_tol=1e-4)
        M = LinearOperator(matrix.shape, ilu.solve)
    except RuntimeError:
        M = scipy.sparse.diags(1.0 / matrix.diagonal())

    t0 = _tick()
    w, info = lgmres(matrix, rhs, M=M, rtol=tol, maxiter=maxiter, callback=_cb)
    t = _tick() - t0
    return t, iters[0], w, info == 0


class _SolverTimeout(Exception):
    pass


def _timeout_iter_solve(matrix, rhs, tol=1e-3, maxiter=10, timeout_s=10):
    """Wrap _iter_solve with a wall-time cap via SIGALRM (Linux/macOS only).

    Returns (elapsed_s, n_iters, w, converged) on success, or None if the
    timeout fires before the solve completes (covers both the ILU
    factorisation and the lgmres iterations).
    """
    def _handler(sig, frame):
        raise _SolverTimeout()

    old_handler = signal.signal(signal.SIGALRM, _handler)
    signal.alarm(timeout_s)
    try:
        result = _iter_solve(matrix, rhs, tol=tol, maxiter=maxiter)
    except _SolverTimeout:
        result = None
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old_handler)
    return result


# ── FD benchmarks: direct vs iterative, three Te profiles ────────────────────

_TE_PROFILES_1D = (
    "constant",
    "sinusoidal", "abrupt", "gradient",
    "wide_range",
    "noisy_mild", "noisy_strong",
)
_TE_PROFILES_2D = (
    "constant",
    "sinusoidal", "sinusoidal_45",
    "abrupt",     "abrupt_45",
    "gradient",   "gradient_45",
    "wide_range", "wide_range_45",
    "noisy_mild", "noisy_strong",
    "disk_thin",  "disk_thick",
)
# Subset used for non-square domain benchmarks: one smooth, one abrupt diagonal,
# one high-dynamic-range, one disk — enough to reveal anisotropic effects
_TE_PROFILES_2D_NONSQ = (
    "constant", "abrupt_45", "wide_range", "disk_thick",
)


def bench_1d_fd(sizes, iter_max=2000, profiles=_TE_PROFILES_1D):
    """Benchmark 1-D FD.

    Iterative solve is skipped for n > iter_max: the ILU-preconditioned
    LGMRES converges reliably up to ~2000 nodes but offers little
    advantage over the direct sparse LU at that scale.
    """
    print("\n1D FD  (direct vs iterative)")
    cols = [
        ("n", 7), ("Te profile", 14),
        ("assemble(s)", 12), ("direct(s)", 10),
        ("iter(s)", 9), ("iters", 6), ("rel_err", 9), ("iter/dir", 9),
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
            w_direct = flex.w.copy()

            if n <= iter_max:
                rhs = -flex.qs
                t_iter, n_iter, w_iter, ok = _iter_solve(flex.coeff_matrix, rhs)
                rel_err = (np.linalg.norm(w_iter - w_direct)
                           / np.linalg.norm(w_direct))
                ratio = t_iter / t_direct if t_direct > 1e-9 else float("nan")
                sfx = "" if ok else "!"   # "!" flags non-convergence
                print(f"  {n:>7}  {prof:>14}  {t_asm:>12.4f}"
                      f"  {t_direct:>10.4f}"
                      f"  {t_iter:>9.4f}  {n_iter:>5}{sfx}"
                      f"  {rel_err:>9.2e}  {ratio:>9.2f}")
            else:
                print(f"  {n:>7}  {prof:>14}  {t_asm:>12.4f}"
                      f"  {t_direct:>10.4f}"
                      f"  {'—':>9}  {'—':>6}  {'—':>9}  {'—':>9}")
        if n != sizes[-1]:
            print()


def bench_2d_fd(sizes, iter_timeout=60, profiles=_TE_PROFILES_2D):
    """Benchmark 2-D FD.

    The iterative solve is attempted for every grid size but is abandoned
    after iter_timeout seconds (wall time).  This cap covers both the ILU
    factorisation and the lgmres iterations: the biharmonic operator's
    condition number grows as O(N⁴), so ILU preconditioning becomes
    prohibitively expensive for large grids.  Timed-out rows are marked
    "T/O"; raise iter_timeout to explore larger sizes once a better
    preconditioner (e.g. FFT-based) is in place.
    """
    print("\n2D FD  (direct vs iterative)")
    cols = [
        ("n×n", 9), ("Te profile", 14),
        ("assemble(s)", 12), ("direct(s)", 10),
        ("iter(s)", 9), ("iters", 6), ("rel_err", 9), ("iter/dir", 9),
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
            w_direct = flex.w.flatten()

            rhs = -flex.qs.reshape(-1, order="C")
            result = _timeout_iter_solve(flex.coeff_matrix, rhs,
                                         timeout_s=iter_timeout)
            if result is not None:
                t_iter, n_iter, w_iter, ok = result
                rel_err = (np.linalg.norm(w_iter - w_direct)
                           / np.linalg.norm(w_direct))
                ratio = t_iter / t_direct if t_direct > 1e-9 else float("nan")
                sfx = "" if ok else "!"
                print(f"  {label:>9}  {prof:>14}  {t_asm:>12.4f}"
                      f"  {t_direct:>10.4f}"
                      f"  {t_iter:>9.4f}  {n_iter:>5}{sfx}"
                      f"  {rel_err:>9.2e}  {ratio:>9.2f}")
            else:
                print(f"  {label:>9}  {prof:>14}  {t_asm:>12.4f}"
                      f"  {t_direct:>10.4f}"
                      f"  {'T/O':>9}  {'—':>6}  {'—':>9}  {'—':>9}")
        if n != sizes[-1]:
            print()


# ── 2D FD non-square domains ─────────────────────────────────────────────────

def bench_2d_fd_nonsquare(shapes, iter_max=100, profiles=_TE_PROFILES_2D_NONSQ):
    """FD benchmark on non-square grids (nx ≠ ny), direct solve only.

    Aspect ratios up to 4:1 test anisotropic stencil behaviour.
    A representative subset of Te profiles is used to keep run time manageable.
    """
    print("\n2D FD  non-square domains  (direct only, subset of Te profiles)")
    cols = [
        ("nx×ny", 11), ("Te profile", 14),
        ("assemble(s)", 12), ("direct(s)", 10),
    ]
    _hdr(cols)
    for (nx, ny) in shapes:
        label = f"{nx}×{ny}"
        for prof in profiles:
            te = _te_2d(ny, nx, prof)
            flex = _make_f2d(nx, ny, "FD", te)
            flex.bc_check()

            t0 = _tick()
            flex.elasprep()
            flex.BC_selector_and_coeff_matrix_creator()
            t_asm = _tick() - t0

            flex.Solver = "direct"
            t0 = _tick()
            flex.fd_solve()
            t_direct = _tick() - t0

            print(f"  {label:>11}  {prof:>14}"
                  f"  {t_asm:>12.4f}  {t_direct:>10.4f}")
        if shapes.index((nx, ny)) != len(shapes) - 1:
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
    from datetime import datetime, timezone

    git = _git_info()
    ts  = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    _RESULTS_DIR = os.path.join(os.path.dirname(__file__), "results")
    outpath = os.path.join(_RESULTS_DIR, f"{git['short']}_{ts}.txt")

    tee = _Tee(outpath)
    try:
        print(_header(git, ts))

        print("\n--- 1D ---")
        bench_1d_fd(sizes=[100, 500, 2000, 5000])
        bench_1d_fft(sizes=[100, 500, 2000, 5000, 20000])
        bench_1d_sas(sizes=[100, 500, 2000, 5000])

        print("\n--- 2D ---")
        bench_2d_fd(sizes=[50, 100, 200, 400], iter_timeout=10)
        bench_2d_fd_nonsquare(shapes=[(200, 50), (400, 100), (200, 25)])
        bench_2d_fft(sizes=[50, 100, 500, 1000])
        bench_2d_sas(sizes=[10, 25, 50, 100])
    finally:
        tee.close()

    hostname = platform.node()
    print(f"\nResults saved to {outpath}", file=sys.__stdout__)
    _push_result(outpath, hostname)
