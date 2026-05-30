#!/usr/bin/env python3
"""
Generate ball-and-stick finite-difference boundary condition diagrams.

Each boundary condition gets its own PNG saved to docs/_static/.
Nodes sit on a schematic deflection profile so the physical meaning of each
reflection is immediately visible.  Ghost nodes are shown at their
mathematically correct reflected positions.

Usage
-----
    python docs/scripts/generate_bc_diagrams.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from scipy.interpolate import CubicSpline

HERE   = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, "..", "_static")
os.makedirs(OUTDIR, exist_ok=True)

# ── shared style ───────────────────────────────────────────────────────────
R  = 0.10          # node radius (data units)

C_REAL  = "#1a6faf"   # blue  — physical domain node
C_GHOST = "white"     # open  — ghost node (reflected)
C_EXCL  = "#d0d0d0"   # grey  — excluded ghost
C_EDGE  = "#1a6faf"
C_EDGX  = "#999999"
C_EQN   = "#c03a2b"   # red   — reflection equation on label
C_CURVE = "#aaaaaa"   # grey  — schematic deflection curve
C_BDY   = "#222222"   # black — boundary dashed line

FIG_W, FIG_H = 7.0, 2.8
DPI = 150


# ── drawing primitives ────────────────────────────────────────────────────

def _circle(ax, x, y, kind):
    fc = {"real": C_REAL, "pin": C_REAL, "ghost": C_GHOST, "excl": C_EXCL}[kind]
    ec = C_EDGX if kind == "excl" else C_EDGE
    lw = 1.2 if kind == "excl" else 1.6
    ls = "--" if kind == "excl" else "-"
    ax.add_patch(Circle((x, y), R, fc=fc, ec=ec, lw=lw, ls=ls,
                         zorder=6, clip_on=False))


def _pin(ax, x, y):
    """Triangle-and-hatch pin symbol below node (x, y)."""
    yt = y - R
    yb = yt - 0.30
    hw = 0.18
    ax.plot([x, x],          [yt, yb],    "k-", lw=1.4, zorder=7)
    ax.plot([x-hw, x+hw],    [yb, yb],    "k-", lw=1.4, zorder=7)
    for xi in np.linspace(x - hw, x + hw, 6):
        ax.plot([xi, xi - 0.06], [yb, yb - 0.09], "k-", lw=0.9, zorder=7)


def _curve(ax, xs, ys):
    """Smooth cubic-spline curve through node positions."""
    order = np.argsort(xs)
    xs_s, ys_s = xs[order], ys[order]
    if len(xs_s) >= 4:
        cs   = CubicSpline(xs_s, ys_s)
        xf   = np.linspace(xs_s[0], xs_s[-1], 300)
        ax.plot(xf, cs(xf), color=C_CURVE, lw=0.9, zorder=2, alpha=0.7)
    else:
        ax.plot(xs_s, ys_s, color=C_CURVE, lw=0.9, zorder=2, alpha=0.7)


# ── main drawing function ──────────────────────────────────────────────────

def draw_bc(*, name, title, subtitle,
            bdy_x,
            real_xy,          # [(x, y), ...] — physical domain nodes, left→right
            ghost_xy,         # [(x, y), ...] — ghost nodes, inner→outer
            ghost_kinds,      # ['ghost'|'excl', ...] matching ghost_xy
            ghost_equations=(),  # one string per ghost node (or None); shown in red
                                 # below the node label, e.g. r"$= -w_1$"
            pins=(),          # indices into real_xy that get a pin symbol
            stencil_note=None,   # optional small text (Mirror/0S0S distinction)
            ):
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))

    rx = np.array([p[0] for p in real_xy])
    ry = np.array([p[1] for p in real_xy])
    gx = np.array([p[0] for p in ghost_xy])
    gy = np.array([p[1] for p in ghost_xy])

    all_x = np.concatenate([gx, rx])
    all_y = np.concatenate([gy, ry])

    # axes limits with padding — extra bottom room for two-line ghost labels
    xlo = all_x.min() - 0.70
    xhi = all_x.max() + 0.70
    ylo = all_y.min() - 0.95
    yhi = all_y.max() + 0.72
    ax.set_xlim(xlo, xhi)
    ax.set_ylim(ylo, yhi)
    ax.set_aspect("equal")
    ax.axis("off")

    # equilibrium baseline
    ax.axhline(0, color="#c8c8c8", lw=0.8, zorder=1)

    # schematic deflection curve through ghost + real nodes
    # (exclude excluded ghosts from the curve so it doesn't dip to 0 spuriously)
    curve_mask = [k != "excl" for k in ghost_kinds]
    curve_x = np.concatenate([gx[curve_mask], rx])
    curve_y = np.concatenate([gy[curve_mask], ry])
    _curve(ax, curve_x, curve_y)

    # boundary line
    ax.axvline(bdy_x, color=C_BDY, lw=1.3, ls="--", zorder=3, alpha=0.85)
    ax.text(bdy_x, ylo + 0.04, "boundary",
            ha="center", va="bottom", fontsize=7.5, color="#555555")

    # ghost nodes
    for (x, y), kind in zip(ghost_xy, ghost_kinds):
        _circle(ax, x, y, kind)

    # real nodes and pins
    for i, (x, y) in enumerate(real_xy):
        _circle(ax, x, y, "real")
        if i in pins:
            _pin(ax, x, y)

    # node labels (and optional red reflection equation) below the panel baseline
    label_y = ylo + 0.12
    n_ghost = len(ghost_xy)
    eqs = list(ghost_equations) + [None] * n_ghost   # pad to at least n_ghost
    for i, (x, _) in enumerate(ghost_xy):
        idx = i - n_ghost          # -2, -1
        lbl = rf"$w_{{{idx}}}$"
        col = C_EDGX if ghost_kinds[i] == "excl" else C_EDGE
        ax.text(x, label_y, lbl, ha="center", va="top", fontsize=9, color=col)
        if eqs[i] is not None:
            ax.text(x, label_y - 0.20, eqs[i],
                    ha="center", va="top", fontsize=8, color=C_EQN)
    for i, (x, _) in enumerate(real_xy):
        lbl = rf"$w_{i}$" if i > 0 else r"$w_0$"
        ax.text(x, label_y, lbl, ha="center", va="top", fontsize=9, color=C_EDGE)

    # title
    ax.text(xlo + 0.05, yhi - 0.02, title,
            ha="left", va="top", fontsize=10, fontweight="bold")
    ax.text(xlo + 0.05, yhi - 0.22, subtitle,
            ha="left", va="top", fontsize=8.5, color="#444444", style="italic")

    if stencil_note:
        ax.text(xhi - 0.05, yhi - 0.02, stencil_note,
                ha="right", va="top", fontsize=7.5, color="#666666", style="italic")

    fig.tight_layout(pad=0.3)
    outpath = os.path.join(OUTDIR, f"bc_diagram_{name}.png")
    fig.savefig(outpath, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  saved {outpath}")


# ── BC definitions ─────────────────────────────────────────────────────────
# Node x-positions: ghosts at -2, -1; domain at 0, 1, 2, 3.
# Node y-positions: schematic deflection (negative = downward, typical for
# lithospheric loading, but here shown positive for visual clarity).

def run():

    # ------------------------------------------------------------------
    # 0Displacement0Slope — clamped end
    # w = 0 AND dw/dx = 0 at boundary node.
    # Ghost cells excluded from stencil (stencil truncated at boundary).
    # ------------------------------------------------------------------
    real = [(0, 0.00), (1, 0.18), (2, 0.32), (3, 0.40)]
    draw_bc(
        name         = "0Displacement0Slope",
        title        = "0Displacement0Slope",
        subtitle     = "clamped end:  w = 0,  dw/dx = 0",
        bdy_x        = -0.5,
        real_xy      = real,
        ghost_xy     = [(-2, 0.00), (-1, 0.00)],
        ghost_kinds  = ["excl", "excl"],
        ghost_equations = [None, None],
        pins         = [0],
        stencil_note = "ghost cells excluded\n(stencil truncated)",
    )

    # ------------------------------------------------------------------
    # 0Displacement0Moment — simply-supported (pinned) end
    # w = 0 at boundary node (Dirichlet); M = 0 via odd-reflection ghost.
    # Boundary IS node 0; ghost w[-1] = -w[+1].
    # ------------------------------------------------------------------
    ry = [0.00, 0.44, 0.60, 0.50]
    real = list(zip([0, 1, 2, 3], ry))
    draw_bc(
        name        = "0Displacement0Moment",
        title       = "0Displacement0Moment",
        subtitle    = "simply-supported end:  w = 0,  M = 0",
        bdy_x       = -0.5,
        real_xy     = real,
        ghost_xy    = [(-2, -ry[2]), (-1, -ry[1])],   # odd reflection
        ghost_kinds = ["ghost", "ghost"],
        ghost_equations = [r"$= -w_2$", r"$= -w_1$"],
        pins        = [0],
    )

    # ------------------------------------------------------------------
    # Mirror — even reflection; symmetry boundary
    # Boundary at cell face between ghost and node 0.
    # Ghost w[-1] = +w[+1] (even), w[-2] = +w[+2].
    # Node 0 is the first domain node; may have nonzero deflection.
    # ------------------------------------------------------------------
    ry = [0.58, 0.72, 0.68, 0.52]
    real = list(zip([0, 1, 2, 3], ry))
    draw_bc(
        name        = "Mirror",
        title       = "Mirror",
        subtitle    = "symmetry (even reflection):  dw/dx = 0 at boundary face",
        bdy_x       = -0.5,
        real_xy     = real,
        ghost_xy    = [(-2, ry[2]), (-1, ry[1])],   # even reflection
        ghost_kinds = ["ghost", "ghost"],
        ghost_equations = [r"$= +w_2$", r"$= +w_1$"],
        pins        = [],
        stencil_note = "stencil at w₁: ghost\nexcluded from i−1 slot",
    )

    # ------------------------------------------------------------------
    # 0Slope0Shear — level, shear-free boundary
    # Same ghost values as Mirror (even reflection), but the stencil at
    # node 1 carries an additional ghost contribution to the w[3] slot.
    # ------------------------------------------------------------------
    ry = [0.58, 0.72, 0.68, 0.52]   # same profile as Mirror for comparison
    real = list(zip([0, 1, 2, 3], ry))
    draw_bc(
        name        = "0Slope0Shear",
        title       = "0Slope0Shear",
        subtitle    = "level, shear-free:  dw/dx = 0,  V = 0  (free vertical displacement)",
        bdy_x       = -0.5,
        real_xy     = real,
        ghost_xy    = [(-2, ry[2]), (-1, ry[1])],
        ghost_kinds = ["ghost", "ghost"],
        ghost_equations = [r"$= +w_2$", r"$= +w_1$"],
        pins        = [],
        stencil_note = "stencil at w₁: ghost also\nadds to w₃ slot\n(≠ Mirror)",
    )

    # ------------------------------------------------------------------
    # 0Moment0Shear — free (broken) end
    # M = 0: w[-1] = 2w[0] - w[1]  (zero curvature extrapolation)
    # V = 0: w[-2] = 4w[0] - 4w[1] + w[2]
    # ------------------------------------------------------------------
    ry = [0.72, 0.54, 0.38, 0.24]
    real = list(zip([0, 1, 2, 3], ry))
    gy1 = 2*ry[0] - ry[1]              # = 0.90
    gy2 = 4*ry[0] - 4*ry[1] + ry[2]   # = 1.14
    draw_bc(
        name        = "0Moment0Shear",
        title       = "0Moment0Shear",
        subtitle    = "free (broken) end:  M = 0,  V = 0",
        bdy_x       = -0.5,
        real_xy     = real,
        ghost_xy    = [(-2, gy2), (-1, gy1)],
        ghost_kinds = ["ghost", "ghost"],
        ghost_equations = [r"$= 4w_0 - 4w_1 + w_2$", r"$= 2w_0 - w_1$"],
        pins        = [],
    )

    # ------------------------------------------------------------------
    # Periodic — wrap-around
    # Ghost cells come from the far end of the domain.
    # Shown schematically: ghost positions labelled w_{N-2}, w_{N-1}.
    # ------------------------------------------------------------------
    ry = [0.32, 0.55, 0.65, 0.52]
    gy1 = 0.18   # ≈ w[N-2] from right end of a sine-like profile
    gy2 = 0.10   # ≈ w[N-3]
    real = list(zip([0, 1, 2, 3], ry))
    draw_bc(
        name        = "Periodic",
        title       = "Periodic",
        subtitle    = "wrap-around:  domain tiles infinitely",
        bdy_x       = -0.5,
        real_xy     = real,
        ghost_xy    = [(-2, gy2), (-1, gy1)],
        ghost_kinds = ["ghost", "ghost"],
        ghost_equations = [r"$= w_{N-3}$", r"$= w_{N-2}$"],
        pins        = [],
    )


if __name__ == "__main__":
    print("Generating BC diagrams …")
    run()
    print("Done.")
