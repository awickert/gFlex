#!/usr/bin/env python3
"""
Generate benchmark figures from a gFlex bench_solvers results file.

Usage::

    python benchmarks/plot_results.py [results_file] [output_dir]

``results_file``
    Path to a ``*.txt`` file produced by ``bench_solvers.py``.  Defaults to
    the most recent file in ``benchmarks/results/``.

``output_dir``
    Directory in which to write the output PNGs.  Defaults to
    ``docs/_static/``.  Three files are written:

    * ``bench_scaling.png``   — solver scaling: time vs. grid size
    * ``bench_lu_cache.png``  — LU cache speedup via the run() path
    * ``bench_te_sweep.png``  — load-only vs. Te-change cost per solve
"""

import glob
import os
import re
import sys
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── paths ─────────────────────────────────────────────────────────────────────

_HERE        = os.path.dirname(os.path.abspath(__file__))
_RESULTS_DIR = os.path.join(_HERE, "results")
_STATIC_DIR  = os.path.join(_HERE, "..", "docs", "_static")


def _latest_results(results_dir=_RESULTS_DIR):
    files = glob.glob(os.path.join(results_dir, "*.txt"))
    if not files:
        return None
    # Sort by ISO timestamp embedded in the filename (_{YYYY-MM-DDTHH:MM:SSZ}.txt)
    def _ts(f):
        m = re.search(r"_(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z)", os.path.basename(f))
        return m.group(1) if m else ""
    return sorted(files, key=_ts)[-1]


# ── results-file parser ───────────────────────────────────────────────────────

def _parse_file(path):
    """Parse a bench_solvers results file into a dict of row lists.

    Each value is a list of rows; each row is a list of stripped string fields
    split on 2-or-more whitespace characters.

    Keys returned
    -------------
    fd1d, fft1d, sas1d     — 1-D method sections
    fd2d, fft2d, sas2d     — 2-D method sections
    lu_run                 — LU cache, run() path
    te_sweep               — Te-sweep section
    """
    with open(path) as f:
        lines = [l.rstrip() for l in f]

    sections = {}
    title     = None
    past_sep  = False
    rows      = []

    for line in lines:
        s = line.strip()

        # Skip blank lines, file-header block (=====), and --- X --- top-level
        # section markers.  The marker format is "--- text ---" (exactly 3
        # leading dashes then a space); column separators like
        # "-------  ------  ----------" must NOT be skipped here — they are
        # detected inside the subsection loop below.
        if not s:
            continue
        if re.match(r"^={5,}", s):
            continue
        if s.startswith("--- ") and s.endswith(" ---"):
            continue

        # Non-indented, non-separator line → new subsection title
        if not line[0].isspace() and not re.match(r"^[-\s]+$", s):
            if title and rows:
                sections[title] = rows
            title    = s
            past_sep = False
            rows     = []
            continue

        if title is None:
            continue

        # Column separator: a line of dashes (possibly with spaces)
        if not past_sep:
            if "-" in s and re.match(r"^[-\s]+$", s):
                past_sep = True
            continue

        # Data row
        parts = re.split(r"\s{2,}", s)
        if parts and parts[0]:
            rows.append(parts)

    if title and rows:
        sections[title] = rows

    # Classify sections by title
    data = {}
    for t, r in sections.items():
        if "1D FD" in t and "direct" in t and "non-square" not in t:
            data["fd1d"] = r
        elif t.startswith("1D FFT"):
            data["fft1d"] = r
        elif t.startswith("1D SAS"):
            data["sas1d"] = r
        elif "2D FD" in t and "direct" in t and "non-square" not in t:
            data["fd2d"] = r
        elif t.startswith("2D FFT"):
            data["fft2d"] = r
        elif t.startswith("2D SAS"):
            data["sas2d"] = r
        elif "LU cache" in t and "run()" in t:
            data["lu_run"] = r
        elif t.startswith("Te sweep"):
            data["te_sweep"] = r

    return data


# ── data-extraction helpers ───────────────────────────────────────────────────

def _first_int(s):
    """Return the first integer embedded in string s (handles '50×50', '500')."""
    m = re.search(r"\d+", s)
    return int(m.group()) if m else None


def _grid_n(grid):
    """Extract the primary dimension N from a grid label.

    '1D-500'  → 500   (1-D: number after the dash)
    '50×50'   → 50    (2-D: first number)
    """
    if grid.startswith("1D-"):
        return int(grid.split("-", 1)[1])
    return _first_int(grid)


def _fd_scaling(rows, n_col, t_asm_col, t_direct_col):
    """Aggregate FD rows → (ns, mean_total, min_total, max_total) arrays."""
    by_n = defaultdict(list)
    for row in rows:
        n = _first_int(row[n_col])
        t = float(row[t_asm_col]) + float(row[t_direct_col])
        by_n[n].append(t)
    ns = sorted(by_n)
    return (np.array(ns),
            np.array([np.mean(by_n[n]) for n in ns]) * 1000,   # → ms
            np.array([np.min(by_n[n])  for n in ns]) * 1000,
            np.array([np.max(by_n[n])  for n in ns]) * 1000)


def _simple_scaling(rows, n_col=0, t_col=1):
    """Extract (ns, times_ms) from FFT/SAS-style rows with one time column."""
    ns, ts = [], []
    for row in rows:
        n = _first_int(row[n_col])
        if n is not None:
            ns.append(n)
            ts.append(float(row[t_col]) * 1000)
    idx = np.argsort(ns)
    return np.array(ns)[idx], np.array(ts)[idx]


def _lu_speedup(rows):
    """Return aggregated speedup data for 1-D and 2-D grids.

    Each returned tuple: (ns, mean_nc, min_nc, max_nc, mean_true)
    where _nc means no_check speedup and _true means True-mode speedup.
    """
    d_nc_1d   = defaultdict(list)
    d_true_1d = defaultdict(list)
    d_nc_2d   = defaultdict(list)
    d_true_2d = defaultdict(list)

    for row in rows:
        grid  = row[0]
        t_f   = float(row[2])
        t_t   = float(row[3])
        t_nc  = float(row[4])
        if t_f <= 0 or t_nc <= 0:
            continue
        sp_nc   = t_f / t_nc
        sp_true = t_f / t_t if t_t > 0 else None
        n = _grid_n(grid)
        # Distinguish 1-D ("1D-500") from 2-D ("50×50") by prefix
        if grid.startswith("1D"):
            d_nc_1d[n].append(sp_nc)
            if sp_true:
                d_true_1d[n].append(sp_true)
        else:
            d_nc_2d[n].append(sp_nc)
            if sp_true:
                d_true_2d[n].append(sp_true)

    def _agg(nc, true):
        ns = sorted(nc)
        return (np.array(ns),
                np.array([np.mean(nc[n])   for n in ns]),
                np.array([np.min(nc[n])    for n in ns]),
                np.array([np.max(nc[n])    for n in ns]),
                np.array([np.mean(true[n]) for n in ns]))

    return _agg(d_nc_1d, d_true_1d), _agg(d_nc_2d, d_true_2d)


def _sweep_vs_load(lu_rows, sw_rows, n_solves=10):
    """Return (1d_tuple, 2d_tuple) where each is (ns, t_load_ms, t_sweep_ms)."""
    lu_1d = defaultdict(list)
    lu_2d = defaultdict(list)
    for row in lu_rows:
        grid = row[0]
        t_nc = float(row[4])
        n    = _grid_n(grid)
        if grid.startswith("1D"):
            lu_1d[n].append(t_nc)
        else:
            lu_2d[n].append(t_nc)

    sw_1d = {}
    sw_2d = {}
    for row in sw_rows:
        grid = row[0]
        t_nc = float(row[3])
        n    = _grid_n(grid)
        if grid.startswith("1D"):
            sw_1d[n] = t_nc
        else:
            sw_2d[n] = t_nc

    def _combine(lu, sw):
        ns = sorted(set(lu) & set(sw))
        t_load  = np.array([np.mean(lu[n]) / n_solves * 1000 for n in ns])
        t_sweep = np.array([sw[n]          / n_solves * 1000 for n in ns])
        return np.array(ns), t_load, t_sweep

    return _combine(lu_1d, sw_1d), _combine(lu_2d, sw_2d)


# ── figure style ──────────────────────────────────────────────────────────────

_C_FD   = "#2166ac"   # blue
_C_FFT  = "#d6604d"   # red
_C_SAS  = "#4dac26"   # green
_C_NC   = "#2166ac"   # no_check — blue
_C_TRUE = "#d6604d"   # True mode — red
_C_LOAD = "#2166ac"   # load-only — blue
_C_SWEP = "#d6604d"   # Te sweep — red

plt.rcParams.update({
    "figure.facecolor": "white",
    "axes.facecolor":   "white",
    "axes.edgecolor":   "#444",
    "axes.grid":        True,
    "grid.color":       "#ccc",
    "grid.linestyle":   "-",
    "grid.linewidth":   0.6,
    "axes.spines.top":  False,
    "axes.spines.right": False,
    "font.size":        10,
    "axes.labelsize":   10,
    "axes.titlesize":   10,
    "legend.fontsize":  8,
    "legend.framealpha": 0.8,
    "xtick.direction":  "in",
    "ytick.direction":  "in",
})


# ── Figure 1: solver scaling ──────────────────────────────────────────────────

def plot_scaling(data, outpath):
    """Solve time vs. grid size for FD, FFT, and SAS (1-D and 2-D)."""
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    fig.subplots_adjust(wspace=0.38)

    # ── 1-D ──────────────────────────────────────────────────────────────────
    ax = axes[0]
    ns_fd, t_fd, t_fd_lo, t_fd_hi = _fd_scaling(data["fd1d"], 0, 2, 3)
    ax.fill_between(ns_fd, t_fd_lo, t_fd_hi, alpha=0.18, color=_C_FD)
    ax.loglog(ns_fd, t_fd, "o-", color=_C_FD, lw=1.6, ms=5,
              label="FD (mean, shaded=range)")

    ns, ts = _simple_scaling(data["fft1d"])
    ax.loglog(ns, ts, "s--", color=_C_FFT, lw=1.6, ms=5, label="FFT")

    ns, ts = _simple_scaling(data["sas1d"])
    ax.loglog(ns, ts, "^:", color=_C_SAS, lw=1.6, ms=5, label="SAS")

    ax.set_xlabel("Grid cells $N$")
    ax.set_ylabel("Solve time (ms)")
    ax.set_title("1-D solvers")
    ax.legend()

    # ── 2-D ──────────────────────────────────────────────────────────────────
    ax = axes[1]
    ns_fd, t_fd, t_fd_lo, t_fd_hi = _fd_scaling(data["fd2d"], 0, 2, 3)
    cells_fd = ns_fd ** 2
    ax.fill_between(cells_fd, t_fd_lo, t_fd_hi, alpha=0.18, color=_C_FD)
    ax.loglog(cells_fd, t_fd, "o-", color=_C_FD, lw=1.6, ms=5,
              label="FD (mean, shaded=range)")

    ns, ts = _simple_scaling(data["fft2d"])
    ax.loglog(ns**2, ts, "s--", color=_C_FFT, lw=1.6, ms=5, label="FFT")

    ns, ts = _simple_scaling(data["sas2d"])
    ax.loglog(ns**2, ts, "^:", color=_C_SAS, lw=1.6, ms=5, label="SAS")

    ax.set_xlabel("Grid cells $N^2$")
    ax.set_ylabel("Solve time (ms)")
    ax.set_title("2-D solvers")
    ax.legend()

    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {outpath}")


# ── Figure 2: LU cache speedup ────────────────────────────────────────────────

def plot_lu_cache(data, outpath):
    """Speedup vs. uncached for no_check and True cache modes (run() path)."""
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    fig.subplots_adjust(wspace=0.38)

    lu1d, lu2d = _lu_speedup(data["lu_run"])

    panels = [
        (axes[0], lu1d, "1-D LU cache speedup  (run() path)", "Grid cells $N$"),
        (axes[1], lu2d, "2-D LU cache speedup  (run() path)", "Grid cells $N^2$"),
    ]
    for ax, (ns, sp_nc, sp_lo, sp_hi, sp_t), title, xlabel in panels:
        x = ns if ax is axes[0] else ns**2
        ax.fill_between(x, sp_lo, sp_hi, alpha=0.18, color=_C_NC)
        ax.semilogx(x, sp_nc, "o-", color=_C_NC, lw=1.6, ms=5,
                    label="no_check (mean, shaded=range)")
        ax.semilogx(x, sp_t,  "s--", color=_C_TRUE, lw=1.6, ms=5,
                    label="True (mean)")
        ax.axhline(1.0, ls=":", color="#666", lw=1.0)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Speedup vs uncached ($t_{\\rm False} / t_{\\rm mode}$)")
        ax.set_title(title)
        ax.legend()
        # Minor x-ticks at each measured point
        ax.set_xticks(x)
        ax.set_xticklabels([str(int(v)) for v in x], rotation=30, ha="right",
                           fontsize=8)

    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {outpath}")


# ── Figure 3: Te sweep vs load-only ───────────────────────────────────────────

def plot_te_sweep(data, outpath):
    """Per-solve cost: load-only (no_check) vs Te changes every call."""
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    fig.subplots_adjust(wspace=0.38)

    d1d, d2d = _sweep_vs_load(data["lu_run"], data["te_sweep"])

    panels = [
        (axes[0], d1d, "1-D: load-only vs. $T_e$ change", "Grid cells $N$"),
        (axes[1], d2d, "2-D: load-only vs. $T_e$ change", "Grid cells $N^2$"),
    ]
    for ax, (ns, t_load, t_sweep), title, xlabel in panels:
        x = ns if ax is axes[0] else ns**2
        ax.semilogy(x, t_load,  "o-",  color=_C_LOAD, lw=1.6, ms=5,
                    label="Load-only  (no_check)")
        ax.semilogy(x, t_sweep, "s--", color=_C_SWEP, lw=1.6, ms=5,
                    label="$T_e$ change (any mode)")
        ax.set_xscale("log")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Time per solve (ms)")
        ax.set_title(title)
        ax.legend()
        ax.set_xticks(x)
        ax.set_xticklabels([str(int(v)) for v in x], rotation=30, ha="right",
                           fontsize=8)

    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {outpath}")


# ── entry point ───────────────────────────────────────────────────────────────

def main(results_file=None, output_dir=None):
    if results_file is None:
        results_file = _latest_results()
        if results_file is None:
            print("No results files found in benchmarks/results/", file=sys.stderr)
            sys.exit(1)
    if output_dir is None:
        output_dir = _STATIC_DIR

    os.makedirs(output_dir, exist_ok=True)
    print(f"Parsing {results_file}")
    data = _parse_file(results_file)

    missing = [k for k in ("fd1d", "fft1d", "sas1d", "fd2d", "fft2d", "sas2d",
                            "lu_run", "te_sweep") if k not in data]
    if missing:
        print(f"Warning: sections not found in results file: {missing}",
              file=sys.stderr)

    plot_scaling(data,   os.path.join(output_dir, "bench_scaling.png"))
    plot_lu_cache(data,  os.path.join(output_dir, "bench_lu_cache.png"))
    plot_te_sweep(data,  os.path.join(output_dir, "bench_te_sweep.png"))
    print("Done.")


if __name__ == "__main__":
    args = sys.argv[1:]
    main(
        results_file=args[0] if len(args) > 0 else None,
        output_dir=args[1]   if len(args) > 1 else None,
    )
