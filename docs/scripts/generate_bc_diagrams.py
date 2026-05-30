#!/usr/bin/env python3
"""
Generate ball-and-stick finite-difference boundary condition diagrams.

Each boundary condition gets its own SVG saved to docs/_static/.
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
C_ARC   = "#e07000"   # orange — stencil coupling arc
C_CURVE = "#555555"   # dark grey — schematic deflection curve
C_BDY   = "#222222"   # black — boundary dashed line

FIG_W, FIG_H = 7.0, 2.8
DPI = 150


# ── drawing primitives ────────────────────────────────────────────────────

def _circle(ax, x, y, kind):
    fc = {"real": C_REAL, "ghost": C_GHOST, "excl": C_EXCL}[kind]
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
        ax.plot(xf, cs(xf), color=C_CURVE, lw=1.5, zorder=2, alpha=0.8)
    else:
        ax.plot(xs_s, ys_s, color=C_CURVE, lw=1.5, zorder=2, alpha=0.8)


def _stencil_arc(ax, x0, y0, x1, y1, label="", color=C_ARC):
    """Curved arrow from (x0, y0+R) to (x1, y1+R) showing extra stencil coupling."""
    ax.annotate("",
                xy=(x1, y1 + R + 0.04),
                xytext=(x0, y0 + R + 0.04),
                arrowprops=dict(
                    arrowstyle="-|>",
                    color=color,
                    lw=1.1,
                    connectionstyle="arc3,rad=0.38",
                    mutation_scale=8,
                ),
                zorder=8, clip_on=False)
    if label:
        mx = (x0 + x1) / 2
        my = max(y0, y1) + R + 0.28
        ax.text(mx, my, label, ha="center", va="bottom",
                fontsize=7.5, color=color)


# ── main drawing function ──────────────────────────────────────────────────

def draw_bc(*, name, title, subtitle,
            bdy_x,
            real_xy,          # [(x, y), ...] — physical domain nodes, left→right
            ghost_xy,         # [(x, y), ...] — ghost nodes, inner→outer
            ghost_kinds,      # ['ghost'|'excl', ...] matching ghost_xy
            ghost_labels=None,    # custom labels for ghost nodes (default w_{-2}, w_{-1})
            ghost_equations=(),  # one string per ghost node (or None); shown in red
            pins=(),          # indices into real_xy that get a pin symbol
            stencil_note=None,   # optional small text, upper-right corner
            stencil_arcs=(),     # list of dicts {x0,y0,x1,y1,label,color}
            gap_xy=(),           # [(x,y,label,kind), ...] — far-end nodes shown with
                                 # a dashed gap from last real node, then solid between
                                 # each consecutive pair; use 'real' kind for blue balls
            gap_label="· · ·",   # text drawn in the middle of the dashed gap segment
            ):
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))

    rx = np.array([p[0] for p in real_xy])
    ry = np.array([p[1] for p in real_xy])
    gx = np.array([p[0] for p in ghost_xy])
    gy = np.array([p[1] for p in ghost_xy])

    gapx = np.array([p[0] for p in gap_xy]) if gap_xy else np.array([])
    gapy = np.array([p[1] for p in gap_xy]) if gap_xy else np.array([])

    all_x = np.concatenate([gx, rx, gapx]) if len(gapx) else np.concatenate([gx, rx])
    all_y = np.concatenate([gy, ry, gapy]) if len(gapy) else np.concatenate([gy, ry])

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

    # Deflection curve through active ghost nodes + all real nodes (spline).
    # Excluded ghost nodes (0Displacement0Slope) are drawn as a flat horizontal
    # line to w₁ instead: the clamped BC enforces zero slope at the boundary.
    excl_mask = np.array([k == "excl" for k in ghost_kinds])
    curve_mask = ~excl_mask
    curve_x = np.concatenate([gx[curve_mask], rx])
    curve_y = np.concatenate([gy[curve_mask], ry])
    _curve(ax, curve_x, curve_y)
    if excl_mask.any():
        x_start = gx[excl_mask].min()
        ax.plot([x_start, rx[0]], [ry[0], ry[0]],
                color=C_CURVE, lw=1.5, zorder=2, alpha=0.8)

    # Far-end gap nodes: straight dashed line from last real node → first gap node,
    # then solid segments between consecutive gap nodes.
    if gap_xy:
        x_last, y_last = rx[-1], ry[-1]
        x_g0, y_g0 = gap_xy[0][0], gap_xy[0][1]
        ax.plot([x_last, x_g0], [y_last, y_g0],
                color=C_CURVE, lw=1.2, ls=(0, (4, 4)), zorder=2, alpha=0.65)
        if gap_label:
            mx = (x_last + x_g0) / 2
            my = (y_last + y_g0) / 2 + 0.08
            ax.text(mx, my, gap_label, ha="center", va="bottom",
                    fontsize=9, color=C_CURVE, alpha=0.65)
        for i in range(len(gap_xy) - 1):
            ax.plot([gap_xy[i][0], gap_xy[i+1][0]],
                    [gap_xy[i][1], gap_xy[i+1][1]],
                    color=C_CURVE, lw=1.5, zorder=2, alpha=0.8)

    # boundary line
    ax.axvline(bdy_x, color=C_BDY, lw=1.3, ls="--", zorder=3, alpha=0.85)
    ax.text(bdy_x, ylo + 0.04, "boundary",
            ha="center", va="bottom", fontsize=7.5, color="#555555")

    # ghost nodes
    for (x, y), kind in zip(ghost_xy, ghost_kinds):
        _circle(ax, x, y, kind)

    # gap nodes
    for item in gap_xy:
        x, y, _lbl, kind = item
        _circle(ax, x, y, kind)

    # real nodes and pins
    for i, (x, y) in enumerate(real_xy):
        _circle(ax, x, y, "real")
        if i in pins:
            _pin(ax, x, y)

    # stencil coupling arcs (drawn above nodes, before labels)
    for arc in stencil_arcs:
        _stencil_arc(ax,
                     arc['x0'], arc['y0'],
                     arc['x1'], arc['y1'],
                     label=arc.get('label', ''),
                     color=arc.get('color', C_ARC))

    # node labels (and optional red reflection equation) below the panel baseline
    label_y = ylo + 0.12
    n_ghost = len(ghost_xy)
    eqs = list(ghost_equations) + [None] * n_ghost   # pad to at least n_ghost
    for i, (x, _) in enumerate(ghost_xy):
        idx = i - n_ghost          # -2, -1
        lbl = (ghost_labels[i] if ghost_labels else rf"$w_{{{idx}}}$")
        col = C_EDGX if ghost_kinds[i] == "excl" else C_EDGE
        ax.text(x, label_y, lbl, ha="center", va="top", fontsize=9, color=col)
        if eqs[i] is not None:
            ax.text(x, label_y - 0.20, eqs[i],
                    ha="center", va="top", fontsize=8, color=C_EQN)
    for i, (x, _) in enumerate(real_xy):
        lbl = rf"$w_{{{i+1}}}$"
        ax.text(x, label_y, lbl, ha="center", va="top", fontsize=9, color=C_EDGE)
    for item in gap_xy:
        x, _y, lbl, _kind = item
        ax.text(x, label_y, lbl, ha="center", va="top", fontsize=9, color=C_EDGE)

    # title and subtitle — always anchored at left edge; "\n" inserts a line break
    title_x = xlo + 0.05
    ty = yhi - 0.02
    for line in title.split("\n"):
        ax.text(title_x, ty, line,
                ha="left", va="top", fontsize=10, fontweight="bold")
        ty -= 0.20
    for line in subtitle.split("\n"):
        ax.text(title_x, ty, line,
                ha="left", va="top", fontsize=8.5, color="#444444", style="italic")
        ty -= 0.18

    if stencil_note:
        ax.text(xhi - 0.05, yhi - 0.02, stencil_note,
                ha="right", va="top", fontsize=7.5, color="#666666", style="italic")

    fig.tight_layout(pad=0.3)
    outpath = os.path.join(OUTDIR, f"bc_diagram_{name}.svg")
    fig.savefig(outpath, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  saved {outpath}")


# ── BC definitions ─────────────────────────────────────────────────────────
# Node x-positions: ghosts at -2, -1; domain at 0, 1, 2, 3.
# Node y-positions: schematic deflection (positive = upward, for visual clarity).

def run():

    # ------------------------------------------------------------------
    # 0Displacement0Slope — clamped end
    # w = 0 AND dw/dx = 0 at boundary node.
    # Ghost columns are absent from the stiffness matrix entirely
    # (stencil truncated at boundary).
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
        stencil_note = "ghost columns absent\nfrom stiffness matrix",
    )

    # ------------------------------------------------------------------
    # 0Displacement0Moment — simply-supported (pinned) end
    # w = 0 at boundary node (Dirichlet); M = 0 via odd-reflection ghost.
    # Boundary IS node 1; ghost w₋₁ = −w₂.
    # ------------------------------------------------------------------
    ry = [0.00, 0.44, 0.60, 0.50]
    real = list(zip([0, 1, 2, 3], ry))
    draw_bc(
        name        = "0Displacement0Moment",
        title       = "0Displacement0Moment",
        subtitle     = "simply-supported end:\n  w = 0,  M = 0",
        bdy_x       = -0.5,
        real_xy     = real,
        ghost_xy    = [(-2, -ry[2]), (-1, -ry[1])],   # odd reflection
        ghost_kinds = ["ghost", "ghost"],
        ghost_equations = [r"$= -w_3$", r"$= -w_2$"],
        pins        = [0],
    )

    # ------------------------------------------------------------------
    # Mirror — even reflection; symmetry boundary
    # Boundary at cell face between ghost and node 1.
    # Ghost w₋₁ = +w₂ (even), w₋₂ = +w₃.
    # At node w₁: the i−2 ghost slot is dropped; only the i−1 ghost enters.
    # ------------------------------------------------------------------
    ry = [0.58, 0.72, 0.68, 0.52]
    real = list(zip([0, 1, 2, 3], ry))
    draw_bc(
        name        = "Mirror",
        title       = "Mirror",
        subtitle     = "symmetry (even reflection):\n  dw/dx = 0 at boundary face",
        bdy_x       = -0.5,
        real_xy     = real,
        ghost_xy    = [(-2, ry[2]), (-1, ry[1])],   # even reflection
        ghost_kinds = ["ghost", "ghost"],
        ghost_equations = [r"$= +w_3$", r"$= +w_2$"],
        pins        = [],
        stencil_note = "at w₁: i−2 ghost slot\ndropped from stencil",
    )

    # ------------------------------------------------------------------
    # 0Slope0Shear — level, shear-free boundary
    # Same ghost values as Mirror (even reflection), but the reflected ghost
    # w₋₁ also contributes an extra term to the w₃ slot of the stencil at w₁.
    # Orange arc shows this extra coupling (absent in Mirror).
    # ------------------------------------------------------------------
    ry = [0.58, 0.72, 0.68, 0.52]   # same profile as Mirror for comparison
    real = list(zip([0, 1, 2, 3], ry))
    draw_bc(
        name        = "0Slope0Shear",
        title       = "0Slope0Shear",
        subtitle     = "level, shear-free:  dw/dx = 0,  V = 0",
        bdy_x       = -0.5,
        real_xy     = real,
        ghost_xy    = [(-2, ry[2]), (-1, ry[1])],
        ghost_kinds = ["ghost", "ghost"],
        ghost_equations = [r"$= +w_3$", r"$= +w_2$"],
        pins        = [],
        stencil_note = "ghost also adds to\nw₃ slot  (≠ Mirror)",
        stencil_arcs = [{'x0': -1, 'y0': ry[1], 'x1': 2, 'y1': ry[2],
                         'label': 'extra', 'color': C_ARC}],
    )

    # ------------------------------------------------------------------
    # 0Moment0Shear — free (broken) end
    # M = 0: w₋₁ = 2w₁ − w₂     (zero curvature extrapolation)
    # V = 0: w₋₂ = 4w₁ − 4w₂ + w₃
    # Profile ascends into the domain; ghosts extrapolate downward (below
    # equilibrium), producing the downward-curving free-end shape.
    # ------------------------------------------------------------------
    ry = [0.15, 0.40, 0.60, 0.72]
    real = list(zip([0, 1, 2, 3], ry))
    gy1 = 2*ry[0] - ry[1]              # = −0.10
    gy2 = 4*ry[0] - 4*ry[1] + ry[2]   # = −0.40
    draw_bc(
        name        = "0Moment0Shear",
        title       = "0Moment0Shear",
        subtitle    = "free (broken) end:  M = 0,  V = 0",
        bdy_x       = -0.5,
        real_xy     = real,
        ghost_xy    = [(-2, gy2), (-1, gy1)],
        ghost_kinds = ["ghost", "ghost"],
        ghost_equations = [r"$= 4w_1 - 4w_2 + w_3$", r"$= 2w_1 - w_2$"],
        pins        = [],
    )

    # ------------------------------------------------------------------
    # Periodic — wrap-around
    # Ghost cells come from the far end of the domain.
    # Four domain nodes (w₁–w₄) shown near the left boundary, plus the
    # two far-end source nodes w_{N-1} and w_N shown in blue; a dashed
    # gap segment indicates the (unspecified) domain interior distance.
    # The wider figure relative to the other BCs reflects that Periodic
    # is categorically different: no physical wall, just wrap-around.
    #   w₋₁ = w_N   (last node wraps to inner ghost)
    #   w₋₂ = w_{N-1}  (second-to-last wraps to outer ghost)
    # ------------------------------------------------------------------
    ry    = [0.32, 0.55, 0.65, 0.52]
    y_Nm1 = 0.10   # w_{N-1} value (= ghost w₋₂)
    y_N   = 0.18   # w_N value   (= ghost w₋₁)
    real  = list(zip([0, 1, 2, 3], ry))
    draw_bc(
        name        = "Periodic",
        title       = "Periodic",
        subtitle    = "wrap-around:  domain tiles infinitely",
        bdy_x       = -0.5,
        real_xy     = real,
        ghost_xy    = [(-2, y_Nm1), (-1, y_N)],   # outer=w_{N-1}, inner=w_N
        ghost_kinds = ["ghost", "ghost"],
        ghost_equations = [r"$= w_{N-1}$", r"$= w_N$"],
        pins        = [],
        gap_xy      = [(4, y_Nm1, r"$w_{N-1}$", "real"),
                       (5, y_N,   r"$w_N$",      "real")],
        gap_label   = "· · ·",
    )


if __name__ == "__main__":
    print("Generating BC diagrams …")
    run()
    print("Done.")
