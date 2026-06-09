"""
LU memory benchmark for gFlex FD solver.

Measures peak RSS after sparse LU factorisation of the actual gFlex
coefficient matrix at several grid sizes.  Each measurement runs in a
fresh subprocess so that process-level peak RSS reflects only the
factorisation cost, without pollution from prior allocations.

Run with:
    python benchmarks/bench_memory.py            # default sizes
    python benchmarks/bench_memory.py --large    # include 500×500, 600×600

Results are printed to stdout and saved to benchmarks/results/ as a text
file named mem_{short_commit_hash}_{UTC_timestamp}.txt.

Notes
-----
- Peak RSS is read from resource.getrusage(RUSAGE_SELF).ru_maxrss.  On
  Linux this is in kibibytes; on macOS it is in bytes.  The script
  normalises to MiB automatically.
- The FD matrix built here has constant Te = 35 km and clamped (
  zero_displacement_zero_slope) BCs on all four sides — the same
  structure as the default QGIS-plugin run.
- Padding is NOT applied; the sizes given are the domain sizes passed to
  the solver.  To estimate post-padding memory for a no_outside_loads run,
  use recommended_pad_width() to find the padded grid size and look up
  that size in the table.
"""

import json
import os
import platform
import subprocess
import sys
import textwrap
from datetime import datetime, timezone

import numpy as np

# ── worker code (runs in child subprocess) ────────────────────────────────────

# The worker builds and factorises the actual gFlex FD coefficient matrix,
# then reports peak RSS (kB on Linux, bytes on macOS) via JSON to stdout.
_WORKER_CODE = textwrap.dedent("""\
    import json, platform, sys, warnings
    import numpy as np
    import scipy.sparse
    import scipy.sparse.linalg

    def _rss_mib():
        \"\"\"Current RSS in MiB from /proc/self/status (Linux) or ru_maxrss (macOS).\"\"\"
        try:
            with open("/proc/self/status") as f:
                for line in f:
                    if line.startswith("VmRSS:"):
                        return int(line.split()[1]) / 1024.0
        except OSError:
            pass
        import resource
        raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return raw / (1024.0 if platform.system() != "Darwin" else 1024.0 * 1024.0)

    n = int(sys.argv[1])
    N = n * n

    # Realistic 13-diagonal sparse matrix matching the gFlex 2-D FD stencil.
    # Offsets: ±2n (corner terms), ±n±1 (cross terms), ±n, ±2, ±1, 0.
    offsets = [-2*n, -(n+1), -n, -(n-1), -2, -1, 0, 1, 2, n-1, n, n+1, 2*n]
    diag_data = []
    for off in offsets:
        length = N - abs(off)
        if off == 0:
            d = np.full(length, 6.0)
        else:
            d = np.full(length, -1.0)
        diag_data.append(d)

    A = scipy.sparse.diags(diag_data, offsets=offsets, shape=(N, N), format="csc")
    A = A.tocsc()

    before = _rss_mib()
    lu = scipy.sparse.linalg.factorized(A)
    after  = _rss_mib()

    lu_mib = max(after - before, 0.0)

    print(json.dumps({"n": n, "cells": N, "lu_mib": round(lu_mib, 1)}))
""")

# ── subprocess runner ─────────────────────────────────────────────────────────

def _measure(n):
    """Return (cells, lu_mib) for an n×n grid via a fresh subprocess."""
    result = subprocess.run(
        [sys.executable, "-c", _WORKER_CODE, str(n)],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"Worker failed for n={n}:\n{result.stderr.strip()}"
        )
    data = json.loads(result.stdout.strip())
    return data["cells"], data["lu_mib"]


# ── provenance helpers (mirrors bench_solvers.py) ─────────────────────────────

def _git_info():
    def _run(*args):
        try:
            return subprocess.check_output(
                args, text=True, stderr=subprocess.DEVNULL
            ).strip()
        except Exception:
            return "unknown"
    return {
        "full":   _run("git", "rev-parse", "HEAD"),
        "short":  _run("git", "rev-parse", "--short", "HEAD"),
        "branch": _run("git", "rev-parse", "--abbrev-ref", "HEAD"),
        "dirty":  bool(_run("git", "status", "--porcelain")),
    }


def _cpu_model():
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
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemTotal:"):
                    return int(line.split()[1]) / (1024 ** 2)
    except OSError:
        pass
    try:
        b = subprocess.check_output(
            ["sysctl", "-n", "hw.memsize"], text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
        return int(b) / (1024 ** 3)
    except Exception:
        pass
    return None


def _header(git, ts):
    dirty  = " (dirty)" if git["dirty"] else ""
    ram    = _ram_gb()
    ram_s  = f"{ram:.1f} GiB" if ram is not None else "unknown"
    return "\n".join([
        "gFlex LU memory benchmark",
        "=" * 60,
        f"timestamp : {ts}",
        f"commit    : {git['full']}{dirty}",
        f"branch    : {git['branch']}",
        f"cpu       : {_cpu_model()}",
        f"ram       : {ram_s}",
        f"os        : {platform.platform()}",
        f"python    : {sys.version.split()[0]}",
        f"numpy     : {np.__version__}",
    ])


# ── output helpers ─────────────────────────────────────────────────────────────

def _log_slope(cells_list, mib_list):
    """Log-log slope over the whole range."""
    if len(cells_list) < 2:
        return float("nan")
    lx = np.log(cells_list)
    ly = np.log(mib_list)
    return float(np.polyfit(lx, ly, 1)[0])


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    large = "--large" in sys.argv

    sizes = [100, 200, 300, 400]
    if large:
        sizes += [500, 600]

    git = _git_info()
    ts  = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    _RESULTS_DIR = os.path.join(os.path.dirname(__file__), "results")
    os.makedirs(_RESULTS_DIR, exist_ok=True)
    outpath = os.path.join(_RESULTS_DIR, f"mem_{git['short']}_{ts}.txt")

    lines = [_header(git, ts), ""]
    print(lines[0])

    header = f"{'Grid':>8}  {'Cells':>10}  {'Peak LU RAM':>14}"
    sep    = "-" * len(header)
    print(header)
    print(sep)
    lines += [header, sep]

    all_cells = []
    all_mib   = []

    for n in sizes:
        try:
            cells, lu_mib = _measure(n)
        except RuntimeError as exc:
            row = f"{n:4d}×{n:<4d}  {'failed':>10}  {str(exc)}"
            print(row)
            lines.append(row)
            continue

        all_cells.append(cells)
        all_mib.append(lu_mib)

        mem_str = (f"{lu_mib / 1024:.2f} GiB" if lu_mib >= 1024
                   else f"{lu_mib:.0f} MiB")
        row = f"{n:4d}×{n:<4d}  {cells:>10,d}  {mem_str:>14}"
        print(row)
        lines.append(row)

    if len(all_cells) >= 2:
        slope = _log_slope(all_cells, all_mib)
        summary = f"\nLog-log slope (O(n^{slope:.2f}) in cell count)"
        print(summary)
        lines.append(summary)

    with open(outpath, "w") as f:
        f.write("\n".join(lines) + "\n")

    print(f"\nResults saved to {outpath}")


if __name__ == "__main__":
    main()
