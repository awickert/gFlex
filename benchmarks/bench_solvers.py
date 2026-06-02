"""
Solver benchmarks for gFlex.

Run with:
    python benchmarks/bench_solvers.py

Results are printed to stdout and saved to benchmarks/results/ as a text
file named {short_commit_hash}_{UTC_timestamp}.txt.  The file header
records the full commit hash, branch, whether the working tree is clean,
and key hardware/software details so runs are reproducible and comparable.

FD benchmarks cover:
  - direct sparse LU across a range of Te profiles and grid sizes
  - square and non-square (aspect-ratio) domains

Grid-size notes:
  - SAS:  O(N²) in 1D; O(N² log N) in 2D (fftconvolve, load-pattern-independent).
  - FFT:  O(N log N) — handles large grids easily.
  - FD direct:  roughly O(N^1.5) for the sparse LU factorisation.
"""

import os
import platform
import subprocess
import sys
import time

import numpy as np
import scipy
import scipy.sparse
from scipy.ndimage import gaussian_filter, gaussian_filter1d

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

# Flexural wavelength for reference Te — used to set correlation lengths and
# tanh transition widths so they scale with the physical problem.
_D_REF = _COMMON["E"] * _TE_REF**3 / (12.0 * (1.0 - _COMMON["nu"]**2))
_LAMBDA_FLEX = (
    2.0 * np.pi
    * (_D_REF / ((_COMMON["rho_m"] - _COMMON["rho_fill"]) * _COMMON["g"])) ** 0.25
)   # ≈ 330 km for the defaults above


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
    if profile == "tanh_step":
        # Smooth step: same contrast as abrupt (½ → 1½ × ref) but with a
        # transition width of one flexural wavelength — more geologically realistic
        return _TE_REF * (1.0 + 0.5 * np.tanh((x - L / 2.0) / _LAMBDA_FLEX))
    if profile == "corr_mild":
        # Spatially correlated noise ±20 % of reference; Gaussian-smoothed with
        # sigma = ½ flexural wavelength, representative of real Te heterogeneity
        rng = np.random.default_rng(42)
        raw = rng.uniform(-1.0, 1.0, n)
        sigma = 0.5 * _LAMBDA_FLEX / dx
        smoothed = gaussian_filter1d(raw, sigma=sigma)
        smoothed /= np.abs(smoothed).max()
        return _TE_REF * (1.0 + 0.2 * smoothed)
    if profile == "corr_strong":
        # Spatially correlated noise ±50 % of reference; same smoothing as corr_mild
        rng = np.random.default_rng(42)
        raw = rng.uniform(-1.0, 1.0, n)
        sigma = 0.5 * _LAMBDA_FLEX / dx
        smoothed = gaussian_filter1d(raw, sigma=sigma)
        smoothed /= np.abs(smoothed).max()
        return _TE_REF * (1.0 + 0.5 * smoothed)
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

    # Smooth step profiles: same contrast as abrupt but tanh transition of width λ_flex
    if profile == "tanh_step":
        return _TE_REF * (1.0 + 0.5 * np.tanh((xx - Lx / 2.0) / _LAMBDA_FLEX))
    if profile == "tanh_step_45":
        scale = 0.5 * (Lx + Ly)
        return _TE_REF * (1.0 + 0.5 * np.tanh((r - 0.5) * scale / _LAMBDA_FLEX))

    # Spatially correlated noise: Gaussian-smoothed with sigma = ½ flexural wavelength
    sigma = 0.5 * _LAMBDA_FLEX / dx
    rng2 = np.random.default_rng(42)
    if profile == "corr_mild":
        raw = rng2.uniform(-1.0, 1.0, (ny, nx))
        smoothed = gaussian_filter(raw, sigma=sigma)
        smoothed /= np.abs(smoothed).max()
        return _TE_REF * (1.0 + 0.2 * smoothed)
    if profile == "corr_strong":
        raw = rng2.uniform(-1.0, 1.0, (ny, nx))
        smoothed = gaussian_filter(raw, sigma=sigma)
        smoothed /= np.abs(smoothed).max()
        return _TE_REF * (1.0 + 0.5 * smoothed)

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

def _make_f1d(n, method, te=_TE_REF, solver="direct",
              bc="zero_displacement_zero_slope"):
    flex = F1D()
    flex.quiet = True
    flex.method = method
    flex.solver = solver
    flex.dx = 5000.0
    flex.te = te
    flex.qs = np.zeros(n)
    flex.qs[n // 4 : 3 * n // 4] = 1e6
    flex.bc_west = bc
    flex.bc_east = bc
    for k, v in _COMMON.items():
        setattr(flex, k, v)
    flex.initialize()
    return flex


def _make_f2d(nx, ny, method, te=_TE_REF, solver="direct",
              bc="zero_displacement_zero_slope", qs=None):
    flex = F2D()
    flex.quiet = True
    flex.method = method
    flex.solver = solver
    flex.dx = 5000.0
    flex.dy = 5000.0
    flex.te = te
    if qs is None:
        flex.qs = np.zeros((ny, nx))
        # Central quarter-area load (central 50 % of each axis = 25 % of domain
        # area).  SAS timing is now load-pattern-independent (fftconvolve).
        flex.qs[ny // 4 : 3 * ny // 4, nx // 4 : 3 * nx // 4] = 1e6
    else:
        flex.qs = qs.copy()
    flex.bc_west = bc
    flex.bc_east = bc
    flex.bc_north = bc
    flex.bc_south = bc
    for k, v in _COMMON.items():
        setattr(flex, k, v)
    flex.initialize()
    return flex


def _make_qs_2d(ny, nx, pattern):
    """Return a (ny, nx) load array for the named pattern.

    ``'small'``   : 3×3 central patch (N_load = 9)
    ``'quarter'`` : central 50 % of each axis (N_load ≈ nx·ny / 4)
    ``'full'``    : entire domain (N_load = nx·ny)
    """
    qs = np.zeros((ny, nx))
    if pattern == "small":
        cy, cx = ny // 2, nx // 2
        qs[cy - 1 : cy + 2, cx - 1 : cx + 2] = 1e6
    elif pattern == "quarter":
        qs[ny // 4 : 3 * ny // 4, nx // 4 : 3 * nx // 4] = 1e6
    elif pattern == "full":
        qs[:] = 1e6
    else:
        raise ValueError(f"unknown load pattern: {pattern!r}")
    return qs


# ── formatting helpers ────────────────────────────────────────────────────────

def _tick():
    return time.perf_counter()


def _hdr(cols):
    print("  ".join(f"{c:>{w}}" for c, w in cols))
    print("  ".join("-" * w for _, w in cols))


# ── FD benchmarks ────────────────────────────────────────────────────────────

_TE_PROFILES_1D = (
    "constant",
    "sinusoidal", "abrupt", "gradient",
    "wide_range",
    "tanh_step",
    "noisy_mild", "noisy_strong",
    "corr_mild",  "corr_strong",
)
_TE_PROFILES_2D = (
    "constant",
    "sinusoidal",  "sinusoidal_45",
    "abrupt",      "abrupt_45",
    "gradient",    "gradient_45",
    "wide_range",  "wide_range_45",
    "tanh_step",   "tanh_step_45",
    "noisy_mild",  "noisy_strong",
    "corr_mild",   "corr_strong",
    "disk_thin",   "disk_thick",
)
# Subset used for non-square domain benchmarks: representative mix of smooth,
# transitional, correlated, and high-dynamic-range Te structures
_TE_PROFILES_2D_NONSQ = (
    "constant", "tanh_step", "corr_mild", "wide_range", "abrupt_45", "disk_thick",
)


def bench_1d_fd(sizes, profiles=_TE_PROFILES_1D):
    """Benchmark 1-D FD direct solver across Te profiles and grid sizes."""
    print("\n1D FD  (direct sparse LU)")
    cols = [
        ("n", 7), ("Te profile", 14),
        ("assemble(s)", 12), ("direct(s)", 10),
    ]
    _hdr(cols)
    for n in sizes:
        for prof in profiles:
            te = _te_1d(n, prof)
            flex = _make_f1d(n, "fd", te)
            flex.bc_check()
            flex.gridded_x()

            t0 = _tick()
            flex.elasprepFD()
            flex._build_coefficient_matrix()
            t_asm = _tick() - t0

            t0 = _tick()
            flex.fd_solve()
            t_direct = _tick() - t0

            print(f"  {n:>7}  {prof:>14}  {t_asm:>12.4f}  {t_direct:>10.4f}")
        if n != sizes[-1]:
            print()


def bench_2d_fd(sizes, profiles=_TE_PROFILES_2D):
    """Benchmark 2-D FD direct solver across Te profiles and grid sizes."""
    print("\n2D FD  (direct sparse LU)")
    cols = [
        ("n×n", 9), ("Te profile", 14),
        ("assemble(s)", 12), ("direct(s)", 10),
    ]
    _hdr(cols)
    for n in sizes:
        label = f"{n}×{n}"
        for prof in profiles:
            te = _te_2d(n, n, prof)
            flex = _make_f2d(n, n, "fd", te)
            flex.bc_check()

            t0 = _tick()
            flex.elasprep()
            flex._build_coefficient_matrix()
            t_asm = _tick() - t0

            t0 = _tick()
            flex.fd_solve()
            t_direct = _tick() - t0

            print(f"  {label:>9}  {prof:>14}  {t_asm:>12.4f}  {t_direct:>10.4f}")
        if n != sizes[-1]:
            print()


# ── 2D FD non-square domains ─────────────────────────────────────────────────

def bench_2d_fd_nonsquare(shapes, profiles=_TE_PROFILES_2D_NONSQ):
    """FD benchmark on non-square grids (nx ≠ ny): direct solver.

    Aspect ratios up to 4:1 test anisotropic stencil behaviour and the
    effect of grid shape on sparse-LU factorisation cost.
    A representative subset of Te profiles is used to keep run time manageable.
    """
    print("\n2D FD  non-square domains  (direct sparse LU)")
    cols = [
        ("nx×ny", 11), ("Te profile", 14),
        ("assemble(s)", 12), ("direct(s)", 10),
    ]
    _hdr(cols)
    for (nx, ny) in shapes:
        label = f"{nx}×{ny}"
        for prof in profiles:
            te = _te_2d(ny, nx, prof)
            flex = _make_f2d(nx, ny, "fd", te)
            flex.bc_check()

            t0 = _tick()
            flex.elasprep()
            flex._build_coefficient_matrix()
            t_asm = _tick() - t0

            t0 = _tick()
            flex.fd_solve()
            t_direct = _tick() - t0

            print(f"  {label:>11}  {prof:>14}  {t_asm:>12.4f}  {t_direct:>10.4f}")
        if shapes.index((nx, ny)) != len(shapes) - 1:
            print()


# ── FFT and SAS benchmarks (constant Te only) ─────────────────────────────────

def bench_1d_fft(sizes):
    print("\n1D FFT")
    cols = [("n", 8), ("total (s)", 12)]
    _hdr(cols)
    for n in sizes:
        flex = _make_f1d(n, "fft")
        t0 = _tick()
        flex.run()
        print("  ".join([f"{n:>8}", f"{_tick() - t0:>12.4f}"]))


def bench_1d_sas(sizes):
    print("\n1D SAS")
    cols = [("n", 8), ("total (s)", 12)]
    _hdr(cols)
    for n in sizes:
        flex = _make_f1d(n, "sas")
        t0 = _tick()
        flex.run()
        print("  ".join([f"{n:>8}", f"{_tick() - t0:>12.4f}"]))


def bench_2d_fft(sizes):
    print("\n2D FFT")
    cols = [("n×n", 9), ("total (s)", 12)]
    _hdr(cols)
    for n in sizes:
        label = f"{n}×{n}"
        flex = _make_f2d(n, n, "fft")
        t0 = _tick()
        flex.run()
        print("  ".join([f"{label:>9}", f"{_tick() - t0:>12.4f}"]))


def bench_2d_sas(sizes):
    print("\n2D SAS  (O(N² log N) fftconvolve — load-pattern-independent)")
    cols = [("n×n", 9), ("total (s)", 12)]
    _hdr(cols)
    for n in sizes:
        label = f"{n}×{n}"
        flex = _make_f2d(n, n, "sas")
        t0 = _tick()
        flex.run()
        print("  ".join([f"{label:>9}", f"{_tick() - t0:>12.4f}"]))


# Maximum grid side length (cells) at which SAS is run per load pattern.
# SAS now uses fftconvolve (O(N² log N), load-pattern-independent), so the
# caps are conservative; adjust upward if runtimes are acceptable.
_SAS_CAP = {"small": 200, "quarter": 100, "full": 50}


def bench_2d_load_geometry(sizes_fd):
    """Compare FD, FFT, and SAS timing across load-pattern geometries.

    Three load patterns are tested on square grids (constant Te):
    - ``'small'``   : 3×3 central patch (N_load = 9)
    - ``'quarter'`` : central 50 % per axis (N_load ≈ N²/4)
    - ``'full'``    : entire domain (N_load = N²)

    All three methods are load-pattern-independent: FD and FFT because the
    stiffness matrix and spectral operator depend only on the grid and Te;
    SAS because the fftconvolve rewrite is O(N² log N) regardless of how
    many cells are loaded.  The SAS column should be flat across patterns
    at a given grid size.

    A ``'--'`` in the SAS column means n exceeds the per-pattern cap
    (small ≤ 200, quarter ≤ 100, full ≤ 50); caps are conservative and
    can be raised now that SAS no longer scales with N_load.

    Note: ``'full'`` with FD (0Displacement0Slope BCs) produces a
    bowl-shaped deflection — the plate is clamped at its edges, so
    deflection drops to zero there.  The timing is a valid scaling
    measurement regardless.
    """
    print("\n2D load geometry  (FD direct + FFT + SAS, constant Te, square grids)")
    print("  All three methods are load-pattern-independent (SAS uses fftconvolve).")
    print("  '--' in SAS column = n exceeds per-pattern cap (conservative; can raise).")
    print("  'full' with FD: 0Displacement0Slope BCs give bowl-shaped deflection.")
    cols = [
        ("n×n",   9), ("pattern",  9), ("N_load",  9),
        ("FD(s)", 9), ("FFT(s)",   9), ("SAS(s)",  9),
    ]
    _hdr(cols)

    patterns = ("small", "quarter", "full")
    for n in sizes_fd:
        label = f"{n}×{n}"
        for i, pat in enumerate(patterns):
            qs = _make_qs_2d(n, n, pat)
            n_load = int(np.count_nonzero(qs))

            flex_fd = _make_f2d(n, n, "fd", qs=qs)
            t0 = _tick()
            flex_fd.run()
            t_fd = _tick() - t0

            flex_fft = _make_f2d(n, n, "fft", qs=qs)
            t0 = _tick()
            flex_fft.run()
            t_fft = _tick() - t0

            if n <= _SAS_CAP[pat]:
                flex_sas = _make_f2d(n, n, "sas", qs=qs)
                t0 = _tick()
                flex_sas.run()
                t_sas = f"{_tick() - t0:9.4f}"
            else:
                t_sas = f"{'--':>9}"

            row_label = label if i == 0 else ""
            print(f"  {row_label:>9}  {pat:>9}  {n_load:>9}  "
                  f"{t_fd:9.4f}  {t_fft:9.4f}  {t_sas}")
        if n != sizes_fd[-1]:
            print()


# ── LU cache benchmark ───────────────────────────────────────────────────────

def bench_lu_cache(sizes_1d, sizes_2d, n_solves=10):
    """Benchmark the LU factorization cache for repeated solves on a fixed domain.

    Compares three cache modes across n_solves back-to-back ``run()`` calls
    with only the load (``qs``) changing between solves:

    ``False``       — ``spsolve()`` on every call (full re-factorization)
    ``True``        — hash-validated ``factorized()`` reuse
    ``"no_check"``  — unconditional ``factorized()`` reuse (user guarantees
                      matrix stability)

    The reported speedup is ``t_false / t_no_check``.
    """
    import gflex.f1d as _f1d_mod
    import gflex.f2d as _f2d_mod

    rng = np.random.default_rng(42)

    print("\nLU cache  (repeated solves, load-only changes)")
    print(f"  n_solves = {n_solves}")
    cols = [
        ("grid", 9), ("t_False(s)", 11), ("t_True(s)", 10),
        ("t_no_check(s)", 14), ("speedup", 8),
    ]
    _hdr(cols)

    # ── 1D ────────────────────────────────────────────────────────────────────
    for n in sizes_1d:
        # Generate n_solves random sparse point loads
        loads = []
        for _ in range(n_solves):
            qs = np.zeros(n)
            idx = rng.integers(0, n, size=3)
            qs[idx] = rng.uniform(1e5, 1e7, size=3)
            loads.append(qs)

        times = {}
        for mode in (False, True, "no_check"):
            def _setup1d(n=n):
                f = _make_f1d(n, "fd")
                f.bc_check()
                f.gridded_x()
                f.elasprepFD()
                f._build_coefficient_matrix()
                f.cache_factorization = mode
                f._lu = None
                return f

            flex = _setup1d()
            t0 = _tick()
            for qs in loads:
                flex.qs = qs
                flex.fd_solve()
            times[str(mode)] = _tick() - t0

        speedup = times["False"] / times["no_check"]
        print(f"  {'1D-' + str(n):>9}  {times['False']:>11.4f}  {times['True']:>10.4f}"
              f"  {times['no_check']:>14.4f}  {speedup:>8.2f}x")

    # ── 2D ────────────────────────────────────────────────────────────────────
    for n in sizes_2d:
        loads = []
        for _ in range(n_solves):
            qs = np.zeros((n, n))
            rows = rng.integers(0, n, size=3)
            cols_ = rng.integers(0, n, size=3)
            qs[rows, cols_] = rng.uniform(1e5, 1e7, size=3)
            loads.append(qs)

        times = {}
        for mode in (False, True, "no_check"):
            def _setup2d(n=n):
                f = _make_f2d(n, n, "fd")
                f.bc_check()
                f.elasprep()
                f._build_coefficient_matrix()
                f.cache_factorization = mode
                f._lu = None
                return f

            flex = _setup2d()
            t0 = _tick()
            for qs in loads:
                flex.qs = qs
                flex.fd_solve()
            times[str(mode)] = _tick() - t0

        speedup = times["False"] / times["no_check"]
        label = f"{n}×{n}"
        print(f"  {label:>9}  {times['False']:>11.4f}  {times['True']:>10.4f}"
              f"  {times['no_check']:>14.4f}  {speedup:>8.2f}x")


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
        bench_2d_fd(sizes=[50, 100, 200, 400])
        bench_2d_fd_nonsquare(shapes=[
            # 2:1 aspect ratio — moderate elongation, both orientations
            (100, 50),   (50, 100),    # ~5 k cells
            (200, 100),  (100, 200),   # ~20 k cells
            # 4:1 aspect ratio — strong elongation, both orientations
            (200, 50),   (50, 200),    # ~10 k cells
            (400, 100),  (100, 400),   # ~40 k cells
            # 8:1 — very elongated (stress-tests thin-strip performance)
            (200, 25),   (400, 50),    # ~5 k / ~20 k cells
        ])
        bench_2d_fft(sizes=[50, 100, 500, 1000])
        bench_2d_sas(sizes=[10, 25, 50, 100])
        bench_2d_load_geometry(sizes_fd=[25, 50, 100, 200])

        print("\n--- LU cache ---")
        bench_lu_cache(sizes_1d=[500, 2000, 5000], sizes_2d=[50, 100, 200],
                       n_solves=10)
    finally:
        tee.close()

    hostname = platform.node()
    print(f"\nResults saved to {outpath}", file=sys.__stdout__)
    _push_result(outpath, hostname)
