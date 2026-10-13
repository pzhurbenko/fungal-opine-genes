#!/usr/bin/env python3
"""Supplementary figure: distributions behind every threshold of the paper, with the threshold as a vertical line.
Style as paper/figures/make_fig1.py. Inputs: analysis/thresholds/{margins_116.tsv, backbone/, chs/, pkerr/}.
Output: analysis/thresholds/thresholds_distributions.{png,svg}. Env: PlantCad2.
"""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
T = f"{GMO_ROOT}/analysis/thresholds/"
GREY, VERM, BLUE = "#8c8c8c", "#D55E00", "#0072B2"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "font.size": 6.5, "mathtext.fontset": "custom", "mathtext.rm": "sans", "mathtext.it": "sans:italic", "mathtext.bf": "sans:bold", "pdf.fonttype": 42, "svg.fonttype": "none"})
m = pd.read_csv(T + "margins_116.tsv", sep="\t")
A = m[m.base_fate.isin(["B", "C", "single"])]
nc = m[m.base_fate != "contamination"]


def panel(ax, groups, thr, xlabel, title, bins, stable=None, logx=False):
    """groups: list of (values, label, color); thr: threshold; stable: (lo, hi) range with unchanged result."""
    if stable is not None: ax.axvspan(*stable, color="#e6e6e6", lw=0, zorder=0)
    ax.hist([np.asarray(v, float) for v, _, _ in groups], bins=bins, stacked=True, color=[c for _, _, c in groups],
            label=[f"{l} ({len(v)})" for v, l, _ in groups], edgecolor="white", lw=0.3)
    ax.axvline(thr, color="black", lw=0.9, ls="--")
    if logx: ax.set_xscale("log")
    ax.set_xlabel(xlabel); ax.set_ylabel("count"); ax.set_title(title, loc="left", fontsize=7, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False); ax.legend(frameon=False, fontsize=5.5, loc="upper left", handlelength=1)


PANELS = []
st = pd.read_csv(T + "funnel_stability.tsv", sep="\t").set_index("param")
S = lambda k: (st.loc[k, "stable_from"], st.loc[k, "stable_to"])
PANELS.append(lambda ax: panel(ax, [(m.o_bact_pid[m.base_fate != "contamination"], "other proteins", GREY), (m.o_bact_pid[m.base_fate == "contamination"], "probable contamination", VERM)],
                               90, "best bacterial hit, % identity", "a  probable contamination:\n≥90% identity", np.arange(20, 101, 2.5), S("contam")))
PANELS.append(lambda ax: panel(ax, [(nc.o_top50_bact[nc.base_fate == "marginal"], "single bacterial hit", GREY), (A.o_top50_bact, "bacterial-type", BLUE)],
                               35, "non-fungal among 50 best hits", "b  bacterial-type:\n≥35/50 non-fungal hits", np.arange(0, 52, 2), S("A_hits")))
PANELS.append(lambda ax: panel(ax, [(nc.o_h[nc.base_fate == "marginal"], "single bacterial hit", GREY), (A.o_h, "bacterial-type", BLUE)],
                               50, "bit-score margin vs fungal hit", "c  bacterial-type:\nmargin ≥50 bits", np.logspace(0, 3.4, 30), S("A_margin"), logx=True))
PANELS.append(lambda ax: panel(ax, [(A.pident_seed[A.base_fate != "B"], "bacterial-type, other", GREY), (A.pident_seed[A.base_fate == "B"], "T-DNA candidate", VERM)],
                               35, "identity to T-DNA seed, %", "d  T-DNA candidate:\nseed identity ≥35%", np.arange(20, 66, 1.5), S("B_ident")))
PANELS.append(lambda ax: panel(ax, [(A.cov_seed[A.base_fate != "B"], "bacterial-type, other", GREY), (A.cov_seed[A.base_fate == "B"], "T-DNA candidate", VERM)],
                               0.4, "coverage of T-DNA seed", "e  T-DNA candidate:\ncoverage ≥0.4", np.arange(0.2, 1.02, 0.04), S("B_cov")))
PANELS.append(lambda ax: panel(ax, [((A.pident_seed / A.o_bact_pid)[A.base_fate != "B"], "bacterial-type, other", GREY), ((A.pident_seed / A.o_bact_pid)[A.base_fate == "B"], "T-DNA candidate", VERM)],
                               0.9, "seed id. / best bacterial id.", "f  T-DNA candidate:\nratio ≥0.9", np.arange(0.3, 1.4, 0.04), S("B_ratio")))
Ar = A[A.rhizobial_best_bact]
PANELS.append(lambda ax: panel(ax, [(Ar.o_bact_pid[Ar.base_fate != "C"], "bacterial-type, other", GREY), (Ar.o_bact_pid[Ar.base_fate == "C"], "non-T-DNA candidate", VERM)],
                               50, "best Rhizobiales hit, % identity", "g  non-T-DNA candidate:\nRhizobiales hit ≥50%", np.arange(25, 82, 2), S("C_ident")))
EXTRA = T + "extra_panels.py"  # backbone / chs / pkerr panels (written once those analyses are done)
try: exec(open(EXTRA).read())
except FileNotFoundError: pass
n = len(PANELS); nc_ = 4; nr = int(np.ceil(n / nc_))
fig, axs = plt.subplots(nr, nc_, figsize=(7.2, 1.9 * nr), constrained_layout=True)
for ax, f in zip(axs.flat, PANELS): f(ax)
for ax in list(axs.flat)[n:]: ax.axis("off")
for ext in ["png", "svg", "pdf"]: fig.savefig(T + "thresholds_distributions." + ext, dpi=300)
import os
FO = f"{GMO_ROOT}/paper/supplementary/figures/"
os.makedirs(FO, exist_ok=True)
for ext in ["png", "svg", "pdf"]: fig.savefig(FO + "FigureS1_thresholds_distributions." + ext, dpi=300)
print("panels:", n)
