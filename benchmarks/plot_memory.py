#!/usr/bin/env python3
"""
Generate memory-scaling figures for the gFlex FD LU factorisation.

Reads the most recent benchmarks/results/mem_*.txt file produced by
bench_memory.py and writes one figure to docs/_static/:

    bench_memory.png

Two panels:
  Left:  LU RAM vs. cell count (log-log): empirical points, power-law fit,
         and reference O(N^1.0) / O(N^1.5) scaling lines.
  Right: No-outside-loads padding overhead — LU memory for the original vs.
         all-sides-padded domain (T_e = 35 km, dx = 5 km, 67 cells/side).

Usage::

    python benchmarks/plot_memory.py [results_file] [output_dir]

``results_file``
    Path to a mem_*.txt file from bench_memory.py.  Defaults to the most
    recent file in benchmarks/results/.

``output_dir``
    Directory for the output PNG.  Defaults to docs/_static/.
"""

import glob
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


_HERE        = os.path.dirname(os.path.abspath(__file__))
_RESULTS_DIR = os.path.join(_HERE, "results")
_STATIC_DIR  = os.path.join(_HERE, "..", "docs", "_static")


# ── results parser ─────────────────────────────────────────────────────────────

def _latest_mem_results(results_dir=_RESULTS_DIR):
    files = glob.glob(os.path.join(results_dir, "mem_*.txt"))
    if not files:
        return None
    def _ts(f):
        m = re.search(r"_(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z)", os.path.basename(f))
        return m.group(1) if m else ""
    return sorted(files, key=_ts)[-1]


def _parse_mem_file(path):
    """Return list of (cells, lu_mib) from a bench_memory results file."""
    rows = []
    in_data = False
    with open(path) as f:
        for line in f:
            s = line.strip()
            if re.match(r"^-{3,}", s):
                in_data = True
                continue
            if not in_data:
                continue
            if not s or s.lower().startswith("log"):
                continue
            # "100×100       10,000          26 MiB"  or  "500×500  250,000  1.76 GiB"
            m = re.match(
                r"\d+×\d+\s+([\d,]+)\s+([\d.]+)\s+(MiB|GiB)", s
            )
            if m:
                cells = int(m.group(1).replace(",", ""))
                val   = float(m.group(2))
                unit  = m.group(3)
                mib   = val if unit == "MiB" else val * 1024
                rows.append((cells, mib))
    return rows


# ── figure style (matches plot_results.py) ────────────────────────────────────

_C_FD   = "#2166ac"   # blue  — unpadded / empirical
_C_PAD  = "#d6604d"   # red   — padded domain
_C_REF  = "#888888"   # gray  — reference scaling lines

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


# ── figure ─────────────────────────────────────────────────────────────────────

def plot_memory(rows, outpath):
    """Two-panel memory figure: scaling (log-log) + padding overhead."""
    cells = np.array([r[0] for r in rows], dtype=float)
    mib   = np.array([r[1] for r in rows], dtype=float)

    # Power-law fit in log-log space
    slope, intercept = np.polyfit(np.log(cells), np.log(mib), 1)
    C = np.exp(intercept)

    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    fig.subplots_adjust(wspace=0.42)

    # ── left panel: scaling (log-log) ─────────────────────────────────────────
    ax = axes[0]

    N_line = np.logspace(np.log10(cells.min()), np.log10(cells.max() * 3))
    anchor = C * cells.min() ** slope

    ax.loglog(N_line, anchor * (N_line / cells.min()) ** 1.0,
              "--", color=_C_REF, lw=1.0, label="$O(N^{1.0})$")
    ax.loglog(N_line, anchor * (N_line / cells.min()) ** 1.5,
              ":",  color=_C_REF, lw=1.0, label="$O(N^{1.5})$")
    ax.loglog(N_line, C * N_line ** slope,
              "-", color=_C_FD, lw=1.4, alpha=0.55,
              label=f"Fit $O(N^{{{slope:.2f}}})$")
    ax.loglog(cells, mib, "o", color=_C_FD, ms=7, zorder=5,
              label="Measured (SuperLU, COLAMD)")

    ax.set_xlabel("Grid cells $N$")
    ax.set_ylabel("Peak LU RAM (MiB)")
    ax.set_title("FD LU memory scaling (2-D)")
    ax.legend()

    # ── right panel: padding overhead ─────────────────────────────────────────
    ax = axes[1]

    # recommended_pad_width(Te=35e3, dx=5e3) → 67 cells per padded side
    pad_cells_side = 67
    ns_orig   = np.arange(50, 651, 25)
    c_orig    = ns_orig ** 2
    c_padded  = (ns_orig + 2 * pad_cells_side) ** 2

    m_orig   = C * c_orig   ** slope
    m_padded = C * c_padded ** slope

    ax.semilogy(ns_orig, m_orig   / 1024, "-",  color=_C_FD,  lw=1.6,
                label="Original domain")
    ax.semilogy(ns_orig, m_padded / 1024, "--", color=_C_PAD, lw=1.6,
                label="All-sides no_outside_loads\n"
                      r"($T_e = 35$ km, $dx = 5$ km)")

    ax.set_xlabel("Original domain (cells per side $n$)")
    ax.set_ylabel("LU RAM (GiB)")
    ax.set_title("no_outside_loads padding overhead")
    ax.legend()

    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {outpath}")


# ── entry point ────────────────────────────────────────────────────────────────

def main(results_file=None, output_dir=None):
    if results_file is None:
        results_file = _latest_mem_results()
        if results_file is None:
            print("No mem_*.txt files in benchmarks/results/", file=sys.stderr)
            sys.exit(1)
    if output_dir is None:
        output_dir = _STATIC_DIR

    os.makedirs(output_dir, exist_ok=True)
    print(f"Parsing {results_file}")
    rows = _parse_mem_file(results_file)
    if not rows:
        print("No data rows parsed — check results file format.", file=sys.stderr)
        sys.exit(1)
    print(f"  {len(rows)} data points: cells = {[r[0] for r in rows]}")

    outpath = os.path.join(output_dir, "bench_memory.png")
    plot_memory(rows, outpath)
    print("Done.")


if __name__ == "__main__":
    args = sys.argv[1:]
    main(
        results_file=args[0] if len(args) > 0 else None,
        output_dir=args[1]   if len(args) > 1 else None,
    )
