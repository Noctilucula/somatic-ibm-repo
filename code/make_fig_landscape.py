#!/usr/bin/env python3
"""Fig. 2 (Section IV): the continuum of module integrity.

Repository entry point: a pure conceptual schematic with no data
dependency; run from anywhere with
    python code/make_fig_landscape.py
and it writes Fig_LandscapeContinuum.{svg,pdf,png} into ../figures.

Two rows of four square panels sharing one integrity gradient:

  a  attractor landscapes as sums of Gaussian basins on the unit state
     plane; integrity declines left to right (five shallow basins ->
     three -> two -> a single deep basin, with rising barriers)
  b  topological abstraction of the same landscapes: each basin becomes
     a node, drawn as a solid coloured core with ONE dashed ring whose
     radius equals the basin-width sigma of the corresponding basin in
     panel a (so the a <-> b mapping is geometrically honest), and an
     orange dot at the energy minimum.  Edges are drift pathways:
     gently curved when drift between equivalent states is free
     (uniform curvature, high integrity), straight when constrained,
     dotted and crossed in red when the return path is topologically
     blocked; at very low integrity only the single deep basin remains.

Layout (reference-style redesign of panel b):
  * both rows are framed as the SAME unit state plane -- identical
    light-grey frames at identical column positions -- so the
    basin -> node mapping of panel b reads as a spatial map of panel a;
  * the axis frame follows the reference grammar: a shared horizontal
    "Declining module integrity" axis on top (arrow in the reading
    direction, matching the High -> Very low column labels) and rotated
    row captions ("Attractor landscape", "State network") on the left;
    axis metadata uses one muted-violet accent colour;
  * a glyph legend beneath panel b decodes the drawing vocabulary
    (state node, basin width, energy minimum, drift path, blocked
    path), mirroring the legend convention of the reference schematic.

Pure conceptual schematic: no simulation output, no data files.
On-chart text is identity-only; all interpretation lives in the figure
caption (figures/Fig_LandscapeContinuum_caption.tex).

Delivery follows the project convention: editable vector SVG
(svg.fonttype "none", text as text) and PDF (Type 42), PNG raster
preview only.  The figure is saved at exactly 7.2 x 4.05 in (183 mm,
double column); bbox_inches="tight" is deliberately NOT used so the
physical size is preserved.
"""

import os

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch

_HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(_HERE, os.pardir, "figures")
os.makedirs(OUT, exist_ok=True)

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 6,
    "figure.dpi": 300,
    "savefig.dpi": 600,
    "svg.fonttype": "none",   # SVG: text stays editable text
    "pdf.fonttype": 42,       # PDF: TrueType, text stays editable
})

GREY_LIGHT, GREY_MID, GREY_TEXT = "#bdbdbd", "#969696", "#636363"
ORANGE, RED = "#e6550d", "#cb181d"
NODE_FILL = ["#9ecae1", "#6baed6", "#4292c6", "#08519c"]
METADATA = "#7e6ba8"          # axis-metadata accent (reference style)
LABELS = ["High", "Moderate", "Low", "Very low"]

# one config per column: (basin positions, depth, basin-width sigma)
CONFIGS = [
    ([(0.30, 0.35), (0.70, 0.32), (0.50, 0.55), (0.32, 0.72), (0.70, 0.72)],
     0.35, 0.11),
    ([(0.30, 0.42), (0.70, 0.42), (0.50, 0.70)], 0.65, 0.13),
    ([(0.32, 0.50), (0.68, 0.50)], 1.05, 0.15),
    ([(0.50, 0.50)], 1.70, 0.17),
]
EDGES = [
    [(0, 1), (0, 2), (0, 3), (1, 2), (1, 4), (2, 3), (2, 4), (3, 4)],
    [(0, 1), (0, 2), (1, 2)],
    [(0, 1)],
    [],
]

FIGW, FIGH = 7.2, 4.05        # 183 mm double column
XL, XR = 0.075, 0.945         # panel span (figure fractions)
GAPX = 0.028
PW = (XR - XL - 3.0 * GAPX) / 4.0          # 0.1965
PH = PW * FIGW / FIGH                       # square panels
XS = [XL + i * (PW + GAPX) for i in range(4)]
ROW_A_Y, ROW_B_Y = 0.578, 0.150
CORE_K = 0.45                  # node core radius, in units of basin sigma
LEVELS = np.linspace(-1.8, 0.0, 22)

fig = plt.figure(figsize=(FIGW, FIGH))


def fx(inches):
    return inches / FIGW


def fy(inches):
    return inches / FIGH


def fig_artist(a):
    """Add a patch/line placed in raw figure coordinates."""
    a.set_transform(fig.transFigure)
    fig.add_artist(a)
    return a


def frame(ax):
    """Unit state plane: identical light-grey box frame, no ticks."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(True)
        sp.set_edgecolor(GREY_LIGHT)
        sp.set_linewidth(0.4)


# ---------------------------------------------------------------- panel a
for i, (basins, depth, sigma) in enumerate(CONFIGS):
    ax = fig.add_axes([XS[i], ROW_A_Y, PW, PH])
    frame(ax)
    x = y = np.linspace(0.0, 1.0, 120)
    X, Y = np.meshgrid(x, y)
    V = -sum(depth * np.exp(-((X - bx) ** 2 + (Y - by) ** 2)
                            / (2.0 * sigma ** 2)) for bx, by in basins)
    ax.contourf(X, Y, V, levels=LEVELS, cmap="Blues_r")
    ax.contour(X, Y, V, levels=LEVELS, colors="#2c7fb8",
               linewidths=0.2, alpha=0.5)
    ax.plot([bx for bx, _ in basins], [by for _, by in basins], "o",
            color=ORANGE, ms=2.8, mew=0, zorder=5)

# column labels sit between the two rows and label both of them
LBL_Y = (ROW_A_Y + ROW_B_Y + PH) / 2.0
for i, lab in enumerate(LABELS):
    fig.text(XS[i] + PW / 2.0, LBL_Y, lab, ha="center", va="center",
             fontsize=5, color="#333333")

# ---------------------------------------------------------------- panel b
for i, (basins, depth, sigma) in enumerate(CONFIGS):
    ax = fig.add_axes([XS[i], ROW_B_Y, PW, PH])
    frame(ax)

    # drift pathways
    for a, b in EDGES[i]:
        (x1, y1), (x2, y2) = basins[a], basins[b]
        if i == 0:
            # high integrity: free drift, uniform gentle curvature
            ax.add_patch(FancyArrowPatch(
                (x1, y1), (x2, y2), arrowstyle="-",
                connectionstyle="arc3,rad=0.18", color=GREY_MID,
                lw=0.6, alpha=0.8, zorder=1))
        elif i == 1:
            # constrained drift: straight, fewer paths
            ax.plot([x1, x2], [y1, y2], color=GREY_MID, lw=0.65,
                    alpha=0.85, zorder=1)
        else:
            # blocked return path: red dotted
            ax.plot([x1, x2], [y1, y2], color=RED, lw=0.75,
                    ls=(0, (1.2, 1.2)), alpha=0.9, zorder=1)

    # blocking marker
    if i == 2:
        (x1, y1), (x2, y2) = basins
        mx, my = 0.5 * (x1 + x2), 0.5 * (y1 + y2)
        ax.add_patch(Circle((mx, my), 0.062, facecolor="white",
                            edgecolor="none", zorder=3))
        ax.plot(mx, my, "x", color=RED, ms=5.5, mew=1.2, zorder=4)

    # nodes: dashed ring at the basin-width sigma (matches panel a),
    # solid core scaled with sigma, orange dot at the minimum
    for (px, py) in basins:
        ax.add_patch(Circle((px, py), sigma, facecolor="none",
                            edgecolor=NODE_FILL[i], lw=0.55,
                            ls=(0, (2.0, 1.6)), alpha=0.85, zorder=2))
        ax.add_patch(Circle((px, py), CORE_K * sigma,
                            facecolor=NODE_FILL[i], edgecolor="white",
                            lw=0.6, zorder=3))
        ax.plot(px, py, "o", color=ORANGE, ms=2.2, mew=0, zorder=5)

# ------------------------------------------------- panel letters (a, b)
fig.text(0.013, ROW_A_Y + PH + 0.006, "a", fontsize=8, fontweight="bold",
         va="bottom", ha="left")
fig.text(0.013, ROW_B_Y + PH + 0.006, "b", fontsize=8, fontweight="bold",
         va="bottom", ha="left")

# ------------------------------------------- rotated row captions (left)
for txt, ycen in [("Attractor landscape", ROW_A_Y + PH / 2.0),
                  ("State network", ROW_B_Y + PH / 2.0)]:
    fig.text(0.042, ycen, txt, rotation=90, ha="center", va="center",
             fontsize=6, color=METADATA)

# ------------------------------------------------- shared integrity axis
axis_y = 0.950
fig_artist(Line2D([XL, XR], [axis_y, axis_y], color=METADATA, lw=0.6))
for dy in (0.004, -0.004):     # drawn arrowhead (label text stays plain)
    fig_artist(Line2D([XR - 0.006, XR], [axis_y + dy, axis_y],
                      color=METADATA, lw=0.6))
fig.text((XL + XR) / 2.0, 0.958, "Declining module integrity",
         ha="center", va="bottom", fontsize=5.5, color=METADATA)

# ---------------------------------------------------------------- legend
# Reference-style glyph legend under panel b; every entry decodes one
# element of the drawing vocabulary with an identity-only label.
LY_IN = 0.30                   # glyph centre line, inches from bottom
RING_C = "#4292c6"
LEG = [                        # (kind, glyph half-width in, label)
    ("node", 0.125, "state node"),
    ("ring", 0.105, "basin width"),
    ("dot", 0.055, "energy minimum"),
    ("path", 0.150, "drift path"),
    ("block", 0.145, "blocked path"),
]
PAD, GAP_IN, CW = 0.055, 0.20, 0.042   # glyph->text pad, item gap, char width
widths = [2.0 * g + PAD + CW * len(t) for _, g, t in LEG]
cx = (FIGW - (sum(widths) + GAP_IN * (len(LEG) - 1))) / 2.0
for kind, g, txt in LEG:
    gcx = cx + g
    if kind == "node":
        fig_artist(Ellipse((fx(gcx), fy(LY_IN)), fx(0.19), fy(0.19),
                           fill=False, edgecolor=RING_C, lw=0.55,
                           ls=(0, (2.0, 1.6))))
        fig_artist(Ellipse((fx(gcx), fy(LY_IN)), fx(0.10), fy(0.10),
                           facecolor=RING_C, edgecolor="white", lw=0.5))
        fig_artist(Line2D([fx(gcx)], [fy(LY_IN)], ls="none", marker="o",
                          ms=2.4, mfc=ORANGE, mec="none"))
    elif kind == "ring":
        fig_artist(Ellipse((fx(gcx), fy(LY_IN)), fx(0.15), fy(0.15),
                           fill=False, edgecolor=RING_C, lw=0.6,
                           ls=(0, (2.0, 1.6))))
    elif kind == "dot":
        fig_artist(Line2D([fx(gcx)], [fy(LY_IN)], ls="none", marker="o",
                          ms=3.4, mfc=ORANGE, mec="none"))
    elif kind == "path":
        fig_artist(FancyArrowPatch(
            (fx(gcx - 0.13), fy(LY_IN - 0.045)),
            (fx(gcx + 0.13), fy(LY_IN + 0.045)),
            arrowstyle="-", connectionstyle="arc3,rad=0.25",
            color=GREY_MID, lw=0.7))
    elif kind == "block":
        fig_artist(Line2D([fx(gcx - 0.12), fx(gcx + 0.12)],
                          [fy(LY_IN), fy(LY_IN)], color=RED, lw=0.75,
                          ls=(0, (1.2, 1.2))))
        fig_artist(Ellipse((fx(gcx), fy(LY_IN)), fx(0.075), fy(0.075),
                           facecolor="white", edgecolor="none"))
        fig_artist(Line2D([fx(gcx)], [fy(LY_IN)], ls="none", marker="x",
                          ms=4.5, mew=1.0, color=RED))
    fig.text(fx(cx + 2.0 * g + PAD), fy(LY_IN), txt, fontsize=5.5,
             ha="left", va="center", color="#333333")
    cx += 2.0 * g + PAD + CW * len(txt) + GAP_IN

# ---------------------------------------------------------------- export
for ext in ("png", "pdf", "svg"):
    fig.savefig(os.path.join(OUT, "Fig_LandscapeContinuum.%s" % ext))
plt.close(fig)
print("saved Fig_LandscapeContinuum.{png,pdf,svg} to", OUT)
print("panel geometry: PW=%.4f PH=%.4f (in: %.3f x %.3f)"
      % (PW, PH, PW * FIGW, PH * FIGH))
