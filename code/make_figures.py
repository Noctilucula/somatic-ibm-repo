#!/usr/bin/env python3
"""Publication figures for the somatic-channel IBM results.

Repository entry point: reads the archived simulation output from
../data/ibm_results and writes the manuscript figures (SVG + PDF + a
raster PNG preview) and somatic_ibm_summary.json into ../figures.
Run from anywhere:  python code/make_figures.py

Main figure (5 panels, labelled a-e):
  a  model schematic, Muller-style: two tissue strips (fragmented C0 vs
     mixed C2-C3) with event markers (mutation = blue dot, nuclear
     migration through hyphal fusion = yellow dot)
  b  spatial structure at a common checkpoint: representative-replicate
     compartment grids for all four architectures (top) and the
     per-replicate sectoring coefficient F across the connectivity
     gradient (bottom), column-aligned small multiples
  c  fixation probability vs deleterious effect, per connectivity level,
     with the rescaled-Kimura analytic overlay
  d  erosion time of the core gene under three selective regimes
  e  paired-arm disconnection experiment (E3-style rate shift)

Design rule: on-chart text carries identity only (which line / which
event); exact numbers and statistical detail live in the figure caption
(figures/Fig_SomaticIBM_main_caption.tex). Panel labels are lowercase
bold letters (a-e), matching the caption convention of the manuscript
("Fig.1 | ..." with panels a, b, c, ...).

Vector delivery: both figures are exported as editable vector files --
SVG (text kept as editable text objects, svg.fonttype "none") and PDF
(Type 42 TrueType fonts, text preserved as text). The PNG is a raster
preview only (git-ignored).

The snapshot family in ../data/ibm_results (snap_*.npz, all four
levels) was generated with the released deterministic seeding, so the
snapshot grids and every panel-b statistic are exactly reproducible; the
exported figure files differ between runs only in embedded timestamps
and internal element IDs (the dot jitter uses a fixed-seed generator).

Supplementary figure: snapshot time series per level (median replicate).

English labels: the figures are destined for insertion into the English
manuscript.
"""

import json
import os

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(_HERE, os.pardir, "data", "ibm_results")
OUT = os.path.join(_HERE, os.pardir, "figures")
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "font.sans-serif": ["DejaVu Sans"],
    "axes.unicode_minus": False,
    "figure.facecolor": "#FFFFFF",
    "axes.facecolor": "#FFFFFF",
    "axes.edgecolor": "#E5E7EB",
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": False,
    "xtick.major.size": 0,
    "ytick.major.size": 0,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "axes.labelsize": 10.5,
    "axes.titlesize": 11,
    "axes.titleweight": "bold",
    "axes.titlepad": 8,
    "legend.frameon": False,
    "legend.fontsize": 9,
    "figure.dpi": 200,
    "savefig.dpi": 300,
    "savefig.facecolor": "#FFFFFF",
    "svg.fonttype": "none",   # SVG: keep text as editable text objects
    "pdf.fonttype": 42,       # PDF: TrueType, text stays editable text
    "ps.fonttype": 42,
})

CB = {"C0": "#CC3311", "C1": "#EE7733", "C2": "#0077BB", "C3": "#009988"}
G900, G700, G500, G400, G300 = ("#111827", "#374151", "#6B7280",
                                "#9CA3AF", "#D1D5DB")
LVL = ["C0", "C1", "C2", "C3"]
LNAME = {"C0": "Fragmented (C0)", "C1": "Local fusion (C1)",
         "C2": "Local + rare long-range (C2)",
         "C3": "Strong long-range (C3)"}

# schematic palette (Panel A)
WT_BG = "#F1F5FA"       # wild-type shading
MUT_FILL = "#3182BD"    # mutant shading (Blues family, matches panel B)
MUT_DOT = "#12457E"     # mutation event marker
MIG_DOT = "#E6A817"     # migration (hyphal fusion) event marker


def clean(ax, grid=True):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid:
        ax.yaxis.grid(True, alpha=0.10, color=G300)
        ax.set_axisbelow(True)


def panel_letter(ax, s, x=None):
    """Panel letter at ONE constant absolute offset outside the panel's
    top-left corner (0.30 in left, 0.06 in above).  With x fixed across
    panels the letters form a single straight column; called after the
    constrained layout has been finalised."""
    p = ax.get_position()
    if x is None:
        x = p.x0 - 0.30 / FIGW_IN
    fig.text(x, p.y1 + 0.06 / FIGH_IN, s, fontsize=14,
             fontweight="bold", color=G900, ha="left", va="bottom")


# ---------------------------------------------------------------- data
calib = {lv: dict(np.load(os.path.join(RES, "calib_%s.npz" % lv)))
         for lv in LVL}
selresp = {lv: dict(np.load(os.path.join(RES, "sel_%s.npz" % lv)))
           for lv in LVL}
fix = {}
for lv in LVL:
    for S in [0.0, 0.01, 0.02, 0.05, 0.1]:
        z = dict(np.load(os.path.join(RES, "fix_%s_S%.4g.npz" % (lv, S))))
        fix[(lv, S)] = z
ero = {}
for tag, S, w in [("costly", 0.02, 0.0), ("neutral", 0.0, 0.0),
                  ("engine", 0.0, 0.02)]:
    for lv in LVL:
        z = dict(np.load(os.path.join(
            RES, "ero_%s_S%.3g_w%.3g.npz" % (lv, S, w))))
        ero[(tag, lv)] = z
disc = dict(np.load(os.path.join(RES, "disc_S0.02.npz")))

GAMMA = float(calib["C0"]["gamma"])
NTOT = int(calib["C0"]["N_tot"])
NE = {lv: float(calib[lv]["Ne_var"]) for lv in LVL}
FEFF = {lv: float(selresp[lv]["F_eff"]) for lv in LVL}
FSTRUCT = {lv: float(calib[lv]["F_struct"]) for lv in LVL}
MU = 3e-4


def u_kimura(N, s, x0):
    """Haploid diffusion fixation probability, genic selection s."""
    a = 2.0 * N * s
    a = np.asarray(a, dtype=float)
    with np.errstate(over="ignore", invalid="ignore"):
        num = -np.expm1(-a * x0)          # 1 - e^{-a x0}
        den = -np.expm1(-a)               # 1 - e^{-a}
    out = np.where(np.abs(a) < 1e-12, x0, num / den)
    return np.where(np.isfinite(out), out, x0)


# ---------------------------------------------------------------- figure
FIGW_IN, FIGH_IN = 13.8, 15.1
fig = plt.figure(figsize=(FIGW_IN, FIGH_IN), constrained_layout=True)
# slightly larger inter-panel padding than the CL default: keeps the
# b -> c/d and c/d -> e row transitions evenly breathing
fig.get_layout_engine().set(h_pad=0.09, w_pad=0.07)
gs = gridspec.GridSpec(4, 2, figure=fig, wspace=0.16, hspace=0.12,
                       width_ratios=[1, 1],
                       height_ratios=[2.85, 3.60, 3.70, 2.95])

# ── Panel A: model schematic (Muller-style strips) ────────────────────
gA = gs[0, :].subgridspec(2, 1, height_ratios=[1.0, 0.17])
gA_strips = gA[0].subgridspec(1, 2, wspace=0.10)
axA0 = fig.add_subplot(gA_strips[0, 0])
axA1 = fig.add_subplot(gA_strips[0, 1])
axAleg = fig.add_subplot(gA[1])
axAleg.axis("off")


def row_curve(keys, phase=0.0, amp=0.016, n=500):
    """Piecewise-linear story curve with a light deterministic wobble."""
    t = np.linspace(0.0, 1.0, n)
    w = np.interp(t, [k[0] for k in keys], [k[1] for k in keys])
    wob = (amp * np.sin(2 * np.pi * (2.3 * t + phase))
           + 0.6 * amp * np.sin(2 * np.pi * (5.1 * t + 2.3 * phase)))
    mask = (w > 0.02) & (w < 0.98)
    return t, np.clip(w + wob * mask, 0.0, 1.0)


def draw_strip(ax, rows, mut_dots, mig_dots, title, ylabel=True):
    """One tissue strip: 3 compartment rows, Muller-style shading."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 3)
    ax.invert_yaxis()
    ax.set_facecolor(WT_BG)
    for i, (keys, phase, amp) in enumerate(rows):
        t, w = row_curve(keys, phase, amp)
        ax.fill_between(t, i, i + w, color=MUT_FILL, lw=0, zorder=2)
    for yb in (1, 2):                       # white separators between rows
        ax.axhline(yb, color="white", lw=1.8, zorder=3)
    for (tm, row) in mut_dots:              # blue = mutation
        ax.plot([tm], [row + 0.45], "o", ms=7.5, mfc=MUT_DOT,
                mec="white", mew=1.3, zorder=6)
    for (tm, bnd) in mig_dots:              # yellow = migration
        ax.plot([tm], [bnd], "o", ms=7.5, mfc=MIG_DOT,
                mec="white", mew=1.3, zorder=6)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(True); sp.set_color(G400); sp.set_linewidth(1.0)
    ax.set_title(title, loc="left", pad=12, fontsize=10, color=G900)
    if ylabel:
        ax.set_ylabel("Compartments", fontsize=9.5, color=G500)
    # time arrow, top-right (reference-style)
    ax.annotate("", xy=(0.985, 1.045), xytext=(0.84, 1.045),
                xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color=G700, lw=1.3,
                                shrinkA=0, shrinkB=0),
                annotation_clip=False)
    ax.text(0.815, 1.045, "Time", transform=ax.transAxes, ha="right",
            va="center", fontsize=9, color=G700)


# fragmented: mutation confined to its compartment (sectors, F ~ 1)
rows_frag = [
    ([(0.0, 0.0), (1.0, 0.0)], 0.0, 0.0),                       # row 1: WT
    ([(0.0, 0.0), (0.13, 0.0), (0.24, 0.30), (0.38, 0.60),
      (0.52, 0.74), (0.66, 0.70), (0.82, 0.78), (1.0, 0.76)],
     1.7, 0.020),                                               # row 2
    ([(0.0, 0.0), (0.52, 0.0), (0.64, 0.16), (0.78, 0.30),
      (0.92, 0.34), (1.0, 0.32)], 3.1, 0.014),                  # row 3
]
# mixed: migration homogenizes rows, tissue-wide erosion
rows_mix = [
    ([(0.0, 0.0), (0.32, 0.02), (0.46, 0.18), (0.62, 0.42),
      (0.76, 0.62), (0.90, 0.80), (1.0, 0.92)], 0.4, 0.015),
    ([(0.0, 0.0), (0.13, 0.0), (0.24, 0.22), (0.34, 0.40),
      (0.48, 0.52), (0.62, 0.66), (0.76, 0.78), (0.90, 0.88),
      (1.0, 0.95)], 1.7, 0.015),
    ([(0.0, 0.0), (0.32, 0.02), (0.46, 0.15), (0.62, 0.38),
      (0.76, 0.58), (0.90, 0.78), (1.0, 0.90)], 2.9, 0.015),
]
draw_strip(axA0, rows_frag,
           mut_dots=[(0.13, 1), (0.52, 2)], mig_dots=[],
           title="Fragmented tissue (C0): sectors stay confined")
draw_strip(axA1, rows_mix,
           mut_dots=[(0.13, 1)],
           mig_dots=[(0.34, 1), (0.34, 2), (0.62, 1), (0.62, 2)],
           title="Mixed tissue (C2\u2013C3): fusion spreads the lineage",
           ylabel=False)
legA = axAleg.legend(
    handles=[
        Patch(facecolor=WT_BG, edgecolor=G400, label="wild-type shading"),
        Patch(facecolor=MUT_FILL, edgecolor="none",
              label="core-gene-loss shading"),
        Line2D([], [], marker="o", ls="none", mfc=MUT_DOT, mec="white",
               ms=8.5, mew=1.3, label="core-gene loss (mutation, $\\mu$)"),
        Line2D([], [], marker="o", ls="none", mfc=MIG_DOT, mec="white",
               ms=8.5, mew=1.3,
               label="nuclear migration (hyphal fusion)"),
    ],
    loc="center", ncol=4, frameon=False, fontsize=8.8,
    handlelength=1.25, columnspacing=1.1, handletextpad=0.5)

# ── Panel B: spatial structure at a common checkpoint ──────────────────
# Column-wise small multiples: each architecture gets a representative-
# replicate compartment grid (top; the replicate whose whole-body mutant
# frequency is closest to the ensemble mean at t = 1800) above a per-
# replicate sectoring panel (bottom: one dot per polymorphic replicate,
# 0.02 < xbar < 0.98, median marked), so the two rows share columns.
# Dashed line: F = 1/K, the value attained when compartments sample
# independently from the tissue mean.
T_SHOW = 1800
gB = gs[1, :].subgridspec(2, 5, width_ratios=[1, 1, 1, 1, 0.045],
                          height_ratios=[1.90, 1.0])
axesB = [fig.add_subplot(gB[0, i]) for i in range(4)]
axesF = [fig.add_subplot(gB[1, i]) for i in range(4)]
cax = fig.add_subplot(gB[0, 4])
snap_meta = {}
jitter_rng = np.random.default_rng(20260317)
for k, lv in enumerate(LVL):
    z = np.load(os.path.join(RES, "snap_%s.npz" % lv))
    grids = z["g%d" % T_SHOW].astype(float) / 8.0
    R_snap = int(z["R"])
    ens = float(z["g%d" % T_SHOW].sum()) / (NTOT * R_snap)
    xb_r = grids.mean(axis=(1, 2))
    F_r = grids.var(axis=(1, 2)) / np.maximum(xb_r * (1.0 - xb_r), 1e-12)
    poly = (xb_r > 0.02) & (xb_r < 0.98)
    Fs = F_r[poly]
    ri = int(np.argmin(np.abs(xb_r - ens)))
    grid = grids[ri]
    xb = float(grid.mean())
    F = float(grid.var() / (xb * (1.0 - xb)))
    rc = np.bincount(np.rint(grid * 8.0).astype(int).ravel(), minlength=9)
    snap_meta[lv] = dict(
        t=T_SHOW, xbar_ens=ens, R=R_snap, rep=int(ri), xbar_rep=xb,
        F_rep=F, n_poly=int(poly.sum()),
        F_med=float(np.median(Fs)), F_q25=float(np.percentile(Fs, 25)),
        F_q75=float(np.percentile(Fs, 75)), F_min=float(Fs.min()),
        F_max=float(Fs.max()),
        rep_frac0=float(rc[0] / rc.sum()), rep_frac8=float(rc[8] / rc.sum()),
        rep_frac_int=float((rc.sum() - rc[0] - rc[8]) / rc.sum()))
    im = axesB[k].imshow(grid, cmap="Blues", vmin=0.0, vmax=1.0,
                         interpolation="nearest")
    axesB[k].set_title(LNAME[lv], fontsize=10)
    axesB[k].set_xticks([]); axesB[k].set_yticks([])
    for sp in axesB[k].spines.values():
        sp.set_visible(True); sp.set_color(G300)
    af = axesF[k]
    # symmetric jitter: centred by construction, small extra noise
    n_d = Fs.size
    jit = (np.linspace(-1.0, 1.0, n_d) * 0.115
           if n_d > 1 else np.zeros(1))
    jit = jit + jitter_rng.uniform(-0.02, 0.02, size=n_d)
    af.scatter(jit, Fs, s=17, color=CB[lv], alpha=0.55, linewidths=0,
               zorder=3)
    af.hlines(float(np.median(Fs)), -0.20, 0.20, color=CB[lv], lw=2.6,
              zorder=4)
    af.axhline(1.0 / 8.0, color=G400, ls="--", lw=1.0, zorder=1)
    af.set_xlim(-0.5, 0.5)
    af.set_ylim(0.0, 1.06)
    af.set_xticks([])
    af.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
    if k:
        af.tick_params(labelleft=False)
    clean(af)
axesF[0].set_ylabel("sectoring coefficient  $F$")
axesF[0].set_title("Sectoring across the connectivity gradient", loc="left")
axesF[3].text(0.46, 0.145, "$1/K$", fontsize=8.5, color=G500, ha="right",
              va="bottom")
axesF[3].text(1.0, 1.04, "$t=%d$" % T_SHOW, transform=axesF[3].transAxes,
              ha="right", va="bottom", fontsize=9, color=G500)
cb = fig.colorbar(im, cax=cax, shrink=0.85, pad=0.04)
cb.outline.set_visible(False)
cb.ax.tick_params(labelsize=8)
cb.set_label("mutant fraction", fontsize=8)

# ── Panel C: fixation probability ─────────────────────────────────────
axC = fig.add_subplot(gs[2, 0])
Sgrid = [0.0, 0.01, 0.02, 0.05, 0.1]
Sx = np.logspace(-2.6, -0.75, 200)
for lv in LVL:
    c = CB[lv]
    # analytic overlay: s_eff = -gamma * F_eff * S, Ne = Ne_var
    s_eff = -GAMMA * FEFF[lv] * Sx
    u = u_kimura(NE[lv], s_eff, 1.0 / NTOT)
    axC.plot(Sx, u, color=c, lw=1.6, zorder=2)
    xs, ys, lo, hi, zero = [], [], [], [], []
    for S in Sgrid:
        z = fix[(lv, S)]
        k, n = int(z["fixed"]), int(z["R"])
        p = k / n
        if k == 0:
            zero.append(S)
            continue
        # Wilson score interval for the binomial proportion (good
        # coverage at the small event counts of the lower points)
        z = 1.959963985
        den = 1.0 + z * z / n
        ctr = (p + z * z / (2.0 * n)) / den
        half = z * np.sqrt(p * (1.0 - p) / n + z * z / (4.0 * n * n)) / den
        xs.append(S); ys.append(p)
        lo.append(max(ctr - half, 1e-7)); hi.append(ctr + half)
    axC.errorbar(xs, ys, yerr=[np.array(ys) - np.array(lo),
                               np.array(hi) - np.array(ys)],
                 fmt="o", color=c, ms=5.5, lw=1.1, capsize=2.5, zorder=4,
                 label=None)
    if zero:
        axC.scatter(zero, [2.0 / fix[(lv, s)]["R"] for s in zero],
                    marker="v", color=c, s=42, zorder=4)
axC.axhline(1.0 / NTOT, color=G400, ls="--", lw=1.1, zorder=1)
# classical single-level reference (census size, unattenuated selection)
u_cl = u_kimura(NTOT, -Sx, 1.0 / NTOT)
axC.plot(Sx, u_cl, color=G400, ls=":", lw=1.3, zorder=1)
axC.set_yscale("log")
axC.set_xlim(0, 0.107)
axC.set_ylim(3e-6, 3e-3)
axC.set_xlabel("deleterious effect of core-gene loss, $S$")
axC.set_ylabel("fixation probability  $P_{\\mathrm{fix}}$")
axC.set_title("Selection collapse with network mixing", loc="left")
clean(axC)
hC = [plt.Line2D([], [], color=CB[lv], lw=1.8) for lv in LVL]
lC = ["%s" % LNAME[lv] for lv in LVL]
axC.legend(hC, lC, loc="upper right", bbox_to_anchor=(1.0, 1.0),
           fontsize=8.2, handlelength=1.6, labelspacing=0.35, ncol=2,
           frameon=True, framealpha=1.0, facecolor="white",
           edgecolor="none", borderaxespad=0.3)

# ── Panel D: erosion time under three regimes ─────────────────────────
axD = fig.add_subplot(gs[2, 1])
reg_style = {
    "costly": ("#0077BB", "o", "$S=0.02$, $w=0$ (costly core)"),
    "neutral": (G400, "s", "$S=0$, $w=0$ (short-circuited)"),
    "engine": ("#CC3311", "^", "$S=0$, $w=0.02$ (selfish saving)"),
}
for tag, (c, m, lab) in reg_style.items():
    meds, cens = [], []
    for lv in LVL:
        z = ero[(tag, lv)]
        t = z["t_ero"]; er = t >= 0
        meds.append(np.median(t[er]) if er.any() else np.nan)
        cens.append(1.0 - er.mean())
    axD.plot(range(4), meds, color=c, marker=m, ms=7, lw=1.6, zorder=3,
             label=lab)
axD.axhline(1.0 / MU, color=G300, ls="--", lw=1.0)
axD.set_yscale("log")
axD.set_xticks(range(4))
axD.set_xticklabels(["C0", "C1", "C2", "C3"])
axD.set_ylim(300, 30000)
axD.set_xlabel("network connectivity level")
axD.set_ylabel("median time to core-gene loss")
axD.set_title("Erosion regimes across connectivity", loc="left")
clean(axD)
axD.legend(loc="upper right", bbox_to_anchor=(1.0, 1.0), fontsize=8.2,
           handlelength=1.6, labelspacing=0.35)

# ── Panel E: disconnection experiment ─────────────────────────────────
axE = fig.add_subplot(gs[3, :])
trec = disc["t_rec"] + disc["t_star"]
axE.plot(trec, disc["xA"], color="#0077BB", lw=2.2, zorder=3,
         label="arm A, connected ($m$ unchanged)")
axE.plot(trec, disc["xB"], color="#CC3311", lw=2.2, zorder=3,
         label="arm B, disconnected at $t^*$ ($m=0$)")
axE.plot(trec, disc["F_B"], color=G500, ls=":", lw=2.0, zorder=2,
         label="sectoring coefficient $F$, arm B")
axE.axvline(disc["t_star"], color=G400, ls="--", lw=1.1)
axE.text(disc["t_star"] + 280, 0.035, "$t^{*}$", fontsize=10,
         color=G500, va="bottom")
tA = disc["tA"]; tB = disc["tB"]
mA = np.median(tA[tA >= 0]); mB = np.median(tB[tB >= 0])
mA_post = mA - disc["t_star"]; mB_post = mB - disc["t_star"]
axE.set_xlabel("time (nuclear generations)")
axE.set_ylabel("mutant frequency $\\bar{x}$ / sectoring $F$")
axE.set_ylim(0, 1.02)
axE.set_xlim(0, trec[-1] * 1.02)
axE.set_title("Same lineage, one architectural variable: fusion on/off",
              loc="left")
clean(axE)
axE.legend(loc="lower right", bbox_to_anchor=(1.0, 0.0), fontsize=8.6)

# ── unified panel letters ─────────────────────────────────────────
# One constant absolute offset for every panel; a, b, c and e share
# the figure's left margin so their letters form one straight column,
# d (right-hand column) keeps the same offset from its own panel.
fig.canvas.draw()               # finalise constrained-layout positions
XL = min(axA0.get_position().x0, axesB[0].get_position().x0,
         axC.get_position().x0, axE.get_position().x0) - 0.30 / FIGW_IN
panel_letter(axA0, "a", XL)
panel_letter(axesB[0], "b", XL)
panel_letter(axC, "c", XL)
panel_letter(axD, "d")
panel_letter(axE, "e", XL)

fig.savefig(os.path.join(OUT, "Fig_SomaticIBM_main.png"), dpi=300)
fig.savefig(os.path.join(OUT, "Fig_SomaticIBM_main.pdf"))
fig.savefig(os.path.join(OUT, "Fig_SomaticIBM_main.svg"))
print("panel b grid axes (fig fractions):",
      [tuple(round(v, 3) for v in ax.get_position().bounds)
       for ax in axesB])
print("panel b F axes (fig fractions):",
      [tuple(round(v, 3) for v in ax.get_position().bounds)
       for ax in axesF])
plt.close(fig)
print("main figure saved; median erosion A=%.0f B=%.0f" % (mA, mB))

# ---------------------------------------------------------------- fig S1
fig2 = plt.figure(figsize=(13.6, 5.6), constrained_layout=True)
gs2 = gridspec.GridSpec(1, 4, figure=fig2,
                        width_ratios=[1, 1, 1, 0.04], wspace=0.06)
axes = [fig2.add_subplot(gs2[0, i]) for i in range(3)]
cax2 = fig2.add_subplot(gs2[0, 3])
show_times = [1200, 1800, 2400, 3600, 6000]
for k, lv in enumerate(["C0", "C2", "C3"]):
    z = np.load(os.path.join(RES, "snap_%s.npz" % lv))
    times = list(z["times"])
    for j, tt in enumerate(show_times):
        ti = times.index(tt)
        grids = z["g%d" % tt].astype(float) / 8.0
        xs = grids.sum(axis=(1, 2)) / 144.0
        ri = int(np.argsort(np.abs(xs - np.median(xs)))[0])
        grid = grids[ri]
        sub = fig2.add_subplot(gs2[0, k].subgridspec(1, len(show_times))[0, j])
        im = sub.imshow(grid, cmap="Blues", vmin=0, vmax=1,
                        interpolation="nearest")
        sub.set_xticks([]); sub.set_yticks([])
        for sp in sub.spines.values():
            sp.set_visible(True); sp.set_color("#E5E7EB")
        if k == 0:
            sub.set_title("$t=%d$" % tt, fontsize=8.5, color=G500)
        if j == 0:
            sub.set_ylabel(LNAME[lv].split(" (")[0], fontsize=10)
cb2 = fig2.colorbar(im, cax=cax2, shrink=0.8, pad=0.05)
cb2.outline.set_visible(False)
cb2.ax.tick_params(labelsize=8)
cb2.set_label("mutant fraction per compartment", fontsize=8)
for ax in axes:
    ax.set_visible(False)
# (no in-figure suptitle: interpretive text lives in the caption,
#  consistent with the identity-only design rule of the main figure)
fig2.savefig(os.path.join(OUT, "Fig_SomaticIBM_snapshots.png"), dpi=300)
fig2.savefig(os.path.join(OUT, "Fig_SomaticIBM_snapshots.pdf"))
fig2.savefig(os.path.join(OUT, "Fig_SomaticIBM_snapshots.svg"))
plt.close(fig2)
print("snapshot figure saved")

# ---------------------------------------------------------------- summary
summary = dict(
    gamma=GAMMA, N_tot=NTOT, mu=MU,
    levels={lv: dict(Ne_var=NE[lv], F_eff=FEFF[lv],
                     F_struct=FSTRUCT[lv]) for lv in LVL},
    P_fix={("%s|S=%g" % (lv, S)): dict(
        p=float(fix[(lv, S)]["P_fix"]), n_fix=int(fix[(lv, S)]["fixed"]),
        R=int(fix[(lv, S)]["R"])) for lv in LVL for S in Sgrid},
    erosion={("%s|%s" % (tag, lv)): dict(
        median=float(np.median(ero[(tag, lv)]["t_ero"][
            ero[(tag, lv)]["t_ero"] >= 0])),
        eroded=int((ero[(tag, lv)]["t_ero"] >= 0).sum()),
        R=int(ero[(tag, lv)]["R"])) for tag in reg_style for lv in LVL},
    disconnect=dict(t_star=int(disc["t_star"]),
                    thr=float(disc["thr"]),
                    xbar_at_tstar=float(disc["xA"][0]),
                    median_A=float(mA), median_B=float(mB),
                    median_A_post=float(mA_post),
                    median_B_post=float(mB_post),
                    nA=int((tA >= 0).sum()), nB=int((tB >= 0).sum())),
    panelB=snap_meta,
)
with open(os.path.join(OUT, "somatic_ibm_summary.json"), "w") as f:
    json.dump(summary, f, indent=1)
print(json.dumps(summary["levels"], indent=1))
print("panel B stats:")
print(json.dumps(summary["panelB"], indent=1))
print("summary saved")
