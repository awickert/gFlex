#!/usr/bin/env python3
"""
Generate finite-difference stencil diagrams for gFlex documentation.

Produces three SVGs in docs/_static/:
  stencil_1d_interior.svg      — 1-D biharmonic 5-point interior stencil
  stencil_2d_interior.svg      — 2-D biharmonic 13-point interior stencil
                                 (uniform D, dx = dy)
  stencil_mirror_vs_0d0m.svg   — even vs odd ghost-node reflection,
                                 side-by-side comparison

Run from the repo root:
    python docs/scripts/generate_stencil_diagrams.py
"""

import os

import matplotlib
import numpy as np
from matplotlib.patches import Circle
from scipy.interpolate import CubicSpline

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE   = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, "..", "_static")
os.makedirs(OUTDIR, exist_ok=True)

# ── style constants (match generate_bc_diagrams.py) ──────────────────────────
R       = 0.10          # node radius (data units)
C_REAL  = "#1a6faf"     # blue  — interior node
C_CTR   = "#c03a2b"     # red   — stencil centre
C_GHOST = "white"       # open  — ghost node
C_FAINT = "#d8d8d8"     # grey  — out-of-stencil grid node (2-D only)
C_EDGE  = "#1a6faf"
C_EDGF  = "#aaaaaa"     # faint edge for out-of-stencil nodes
C_COEFF = "#e07000"     # orange — coefficient labels
C_CURVE = "#555555"
C_BDY   = "#222222"
C_EQN   = "#c03a2b"
SCALE   = 1.0
MARGIN  = 0.35


# ── primitives ────────────────────────────────────────────────────────────────

def _circle(ax, x, y, kind="real"):
    styles = {
        "real":   dict(fc=C_REAL,  ec=C_EDGE, lw=1.6),
        "center": dict(fc=C_CTR,   ec=C_CTR,  lw=1.6),
        "ghost":  dict(fc="white", ec=C_EDGE, lw=1.6),
        "faint":  dict(fc=C_FAINT, ec=C_EDGF, lw=0.9),
    }[kind]
    ax.add_patch(Circle((x, y), R, zorder=6, clip_on=False, **styles))


def _curve(ax, xs, ys):
    order = np.argsort(xs)
    xs_s, ys_s = xs[order], ys[order]
    if len(xs_s) >= 4:
        cs = CubicSpline(xs_s, ys_s)
        xf = np.linspace(xs_s[0], xs_s[-1], 300)
        ax.plot(xf, cs(xf), color=C_CURVE, lw=1.5, zorder=2, alpha=0.8)
    else:
        ax.plot(xs_s, ys_s, color=C_CURVE, lw=1.5, zorder=2, alpha=0.8)


def _pin(ax, x, y):
    yt = y - R
    yb = yt - 0.28
    hw = 0.16
    ax.plot([x, x],       [yt, yb],  "k-", lw=1.4, zorder=7)
    ax.plot([x-hw, x+hw], [yb, yb],  "k-", lw=1.4, zorder=7)
    for xi in np.linspace(x - hw, x + hw, 5):
        ax.plot([xi, xi - 0.05], [yb, yb - 0.08], "k-", lw=0.9, zorder=7)


# ── 1-D interior stencil ─────────────────────────────────────────────────────

def draw_1d_stencil():
    """
    Five-point biharmonic stencil at an interior node.
    Coefficients (uniform D, normalised by dx⁴/D): +1, −4, +6, −4, +1.
    """
    xs = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    y0 = 0.0
    ys = np.full_like(xs, y0)
    center = 2

    coeffs = ["+1", "−4", "+6", "−4", "+1"]
    labels = [
        r"$w_{i-2}$", r"$w_{i-1}$", r"$w_i$", r"$w_{i+1}$", r"$w_{i+2}$",
    ]

    xlo, xhi = -0.70, 4.70
    ylo, yhi = -0.52, 0.90
    fig_w = (xhi - xlo) * SCALE + MARGIN
    fig_h = (yhi - ylo) * SCALE + MARGIN
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_xlim(xlo, xhi)
    ax.set_ylim(ylo, yhi)
    ax.set_aspect("equal")
    ax.axis("off")

    # Baseline through all nodes
    ax.plot([xs[0] - 0.3, xs[-1] + 0.3], [y0, y0],
            color="#c8c8c8", lw=0.8, zorder=1)

    # Lines from centre to all other stencil nodes
    for i in range(len(xs)):
        if i != center:
            ax.plot([xs[center], xs[i]], [y0, y0],
                    color=C_COEFF, lw=1.0, zorder=3, alpha=0.50)

    # Nodes
    for i, (x, y) in enumerate(zip(xs, ys)):
        _circle(ax, x, y, "center" if i == center else "real")

    # Coefficient labels (just above each node)
    for i, (x, y, coeff) in enumerate(zip(xs, ys, coeffs)):
        color = C_CTR if i == center else C_COEFF
        ax.text(x, y + R + 0.12, coeff,
                ha="center", va="bottom", fontsize=9.5,
                color=color, fontweight="bold")

    # Node labels (bottom strip)
    label_y = ylo + 0.12
    for x, lbl in zip(xs, labels):
        ax.text(x, label_y, lbl, ha="center", va="top", fontsize=9, color=C_EDGE)

    # Title and equation
    ax.text(xlo + 0.05, yhi - 0.02, "1-D interior stencil",
            ha="left", va="top", fontsize=10, fontweight="bold")
    ax.text(
        xlo + 0.05, yhi - 0.22,
        r"$\frac{D}{\Delta x^4}(+1,\,-4,\,+6,\,-4,\,+1)"
        r"\cdot(w_{i-2},\,\ldots,\,w_{i+2})$"
        "   (uniform $D$)",
        ha="left", va="top", fontsize=7.8, color="#444444",
    )

    fig.tight_layout(pad=0.3)
    out = os.path.join(OUTDIR, "stencil_1d_interior.svg")
    fig.savefig(out, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  saved {out}")


# ── 2-D interior stencil ─────────────────────────────────────────────────────

def draw_2d_stencil():
    """
    Thirteen-point 2-D biharmonic stencil for uniform D and dx = dy = h.

    Coefficients (×h⁴/D):
        centre  (0, 0) :  +20
        (±1, 0) and (0, ±1) :  −8   (4 nodes)
        (±2, 0) and (0, ±2) :  +1   (4 nodes)
        (±1, ±1)            :  +2   (4 nodes)
    The 12 remaining nodes of the surrounding 5×5 grid are shown faint.
    """
    # Stencil: (di, dj) → coefficient
    stencil = {
        ( 0,  0): 20,
        ( 1,  0): -8, (-1,  0): -8,
        ( 0,  1): -8, ( 0, -1): -8,
        ( 2,  0):  1, (-2,  0):  1,
        ( 0,  2):  1, ( 0, -2):  1,
        ( 1,  1):  2, (-1,  1):  2,
        ( 1, -1):  2, (-1, -1):  2,
    }
    all_pts   = [(di, dj) for di in range(-2, 3) for dj in range(-2, 3)]
    faint_pts = [p for p in all_pts if p not in stencil]

    margin = 0.90
    xlo, xhi = -2 - margin, 2 + margin
    ylo, yhi = -2 - margin, 2 + margin

    fig_w = (xhi - xlo) * SCALE + MARGIN
    fig_h = (yhi - ylo) * SCALE + MARGIN
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_xlim(xlo, xhi)
    ax.set_ylim(ylo, yhi)
    ax.set_aspect("equal")
    ax.axis("off")

    # Light grid lines
    for v in range(-2, 3):
        ax.axvline(v, color="#ececec", lw=0.7, zorder=0)
        ax.axhline(v, color="#ececec", lw=0.7, zorder=0)

    # Lines from centre to each stencil node
    for (di, dj) in stencil:
        if (di, dj) == (0, 0):
            continue
        ax.plot([0, di], [0, dj], color=C_COEFF, lw=1.0, zorder=3, alpha=0.50)

    # Out-of-stencil nodes
    for (di, dj) in faint_pts:
        _circle(ax, di, dj, "faint")

    # Stencil nodes
    for (di, dj) in stencil:
        _circle(ax, di, dj, "center" if (di, dj) == (0, 0) else "real")

    # Coefficient labels — offset away from the centre so they don't overlap
    offsets = {
        ( 0,  0): ( 0.00,  0.00),   # inside circle
        ( 1,  0): ( 0.24,  0.00),
        (-1,  0): (-0.24,  0.00),
        ( 0,  1): ( 0.00,  0.24),
        ( 0, -1): ( 0.00, -0.24),
        ( 2,  0): ( 0.00,  0.24),
        (-2,  0): ( 0.00,  0.24),
        ( 0,  2): ( 0.24,  0.00),
        ( 0, -2): ( 0.24,  0.00),
        ( 1,  1): ( 0.20,  0.16),
        (-1,  1): (-0.20,  0.16),
        ( 1, -1): ( 0.20, -0.16),
        (-1, -1): (-0.20, -0.16),
    }
    for (di, dj), coeff in stencil.items():
        sign  = "+" if coeff > 0 else ""
        ox, oy = offsets[(di, dj)]
        color = C_CTR if (di, dj) == (0, 0) else C_COEFF
        fs    = 9.5 if (di, dj) == (0, 0) else 8.5
        ax.text(di + ox, dj + oy, f"{sign}{coeff}",
                ha="center", va="center", fontsize=fs,
                color=color, fontweight="bold", zorder=7)

    # Axis direction arrows
    for (dx, dy, lbl) in [(0.38, 0, "i"), (0, 0.38, "j")]:
        ax.annotate("",
                    xy=(xhi - 0.12 + dx * 0, yhi - 0.12 + dy * 0)
                       if dx else (0 + dx * 0, yhi - 0.12),
                    xytext=(xhi - 0.50, yhi - 0.50) if dx else (0, yhi - 0.50),
                    xycoords="data", textcoords="data",
                    arrowprops=dict(arrowstyle="-|>", color="#777777",
                                    lw=0.9, mutation_scale=7))
    ax.annotate("", xy=(xhi - 0.14, -0.0), xytext=(xhi - 0.52, -0.0),
                arrowprops=dict(arrowstyle="-|>", color="#777777",
                                lw=0.9, mutation_scale=7))
    ax.text(xhi - 0.08, 0.0, "$i$", ha="left", va="center",
            fontsize=9, color="#777777")
    ax.annotate("", xy=(0.0, yhi - 0.14), xytext=(0.0, yhi - 0.52),
                arrowprops=dict(arrowstyle="-|>", color="#777777",
                                lw=0.9, mutation_scale=7))
    ax.text(0.0, yhi - 0.07, "$j$", ha="center", va="bottom",
            fontsize=9, color="#777777")

    # Title
    ax.text(xlo + 0.05, yhi - 0.02, "2-D interior stencil",
            ha="left", va="top", fontsize=10, fontweight="bold")
    ax.text(xlo + 0.05, yhi - 0.22,
            r"uniform $D$, $\Delta x = \Delta y = h$"
            r";  coefficients $\times\,h^4/D$",
            ha="left", va="top", fontsize=8, color="#444444")

    fig.tight_layout(pad=0.3)
    out = os.path.join(OUTDIR, "stencil_2d_interior.svg")
    fig.savefig(out, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  saved {out}")


# ── Mirror vs 0Displacement0Moment ghost comparison ──────────────────────────

def draw_ghost_comparison():
    """
    Side-by-side: even reflection (Mirror) and odd reflection (0D0M).
    Both panels use the same real-node deflection profile so that the only
    visual difference is the sign of the ghost values.
    """
    # Shared real-node profile (w₁ … w₄)
    rx  = [0, 1, 2, 3]
    ry  = [0.18, 0.48, 0.64, 0.54]
    bdy_x = -0.5

    panels = [
        dict(
            title    = "Mirror",
            subtitle = "even reflection",
            ghost_ys = [+ry[2], +ry[1]],     # w₀ = +w₃, w₋₁ = +w₂  (ghost_x[0]=-2, [1]=-1)
            ghost_eqs= [r"$w_0 = +w_3$",
                        r"$w_{-1} = +w_2$"],
            pin      = False,
        ),
        dict(
            title    = "0Displacement0Moment",
            subtitle = "odd reflection",
            ghost_ys = [-ry[2], -ry[1]],     # w₀ = −w₃, w₋₁ = −w₂
            ghost_eqs= [r"$w_0 = -w_3$",
                        r"$w_{-1} = -w_2$"],
            pin      = True,                  # w₁ = 0 enforced (Dirichlet)
        ),
    ]
    ghost_x = [-1, -2]   # inner ghost first, outer ghost second

    # ── figure: two fixed-size axes side by side ──────────────────────────
    xlo, xhi = -2.65, 3.65
    ylo, yhi = -1.15, 1.25

    pw = (xhi - xlo) * SCALE        # panel width  (inches)
    ph = (yhi - ylo) * SCALE        # panel height (inches)
    gap  = 0.55                     # gap between panels (inches)
    pad  = MARGIN
    fig_w = 2 * pw + gap + 2 * pad
    fig_h = ph + 2 * pad
    fig = plt.figure(figsize=(fig_w, fig_h))

    for p_idx, panel in enumerate(panels):
        left   = (pad + p_idx * (pw + gap)) / fig_w
        bottom = pad / fig_h
        ax = fig.add_axes([left, bottom, pw / fig_w, ph / fig_h])
        ax.set_xlim(xlo, xhi)
        ax.set_ylim(ylo, yhi)
        ax.set_aspect("equal")
        ax.axis("off")

        # Baseline and boundary
        ax.axhline(0, color="#c8c8c8", lw=0.8, zorder=1)
        ax.axvline(bdy_x, color=C_BDY, lw=1.3, ls="--", zorder=3, alpha=0.85)
        ax.text(bdy_x - 0.03, ylo + 0.04, "boundary",
                ha="right", va="bottom", fontsize=7.5, color="#555555",
                rotation=90)

        # Deflection curve through ghost nodes + real nodes
        gys = panel["ghost_ys"]
        all_x = np.array(ghost_x + rx, dtype=float)
        all_y = np.array(gys + ry, dtype=float)
        _curve(ax, all_x, all_y)

        # Ghost nodes
        for x, y in zip(ghost_x, gys):
            _circle(ax, x, y, "ghost")

        # Real nodes (+ pin for 0D0M w₁ = 0)
        for i, (x, y) in enumerate(zip(rx, ry)):
            _circle(ax, x, y, "real")
            if i == 0 and panel["pin"]:
                _pin(ax, x, y)

        # Node labels (bottom strip)
        label_y = ylo + 0.10
        for x, lbl in zip(ghost_x, [r"$w_0$", r"$w_{-1}$"]):
            ax.text(x, label_y, lbl, ha="center", va="top",
                    fontsize=9, color=C_EDGE)
        for i, x in enumerate(rx):
            ax.text(x, label_y, rf"$w_{{{i+1}}}$",
                    ha="center", va="top", fontsize=9, color=C_EDGE)

        # Ghost equations in red (below labels)
        eq_y = label_y - 0.22
        for x, eq in zip(ghost_x, panel["ghost_eqs"]):
            ax.text(x, eq_y, eq, ha="center", va="top",
                    fontsize=8, color=C_EQN)

        # Panel title and subtitle
        ax.text(xlo + 0.05, yhi - 0.02, panel["title"],
                ha="left", va="top", fontsize=10, fontweight="bold")
        ax.text(xlo + 0.05, yhi - 0.22, panel["subtitle"],
                ha="left", va="top", fontsize=8.5, color="#444444",
                style="italic")

    out = os.path.join(OUTDIR, "stencil_mirror_vs_0d0m.svg")
    fig.savefig(out, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  saved {out}")


# ── entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("Generating stencil diagrams …")
    draw_1d_stencil()
    draw_2d_stencil()
    draw_ghost_comparison()
    print("Done.")
