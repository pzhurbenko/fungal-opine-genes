"""Figure 2 — protein-level screen by gene (dot plot; formerly panel b of Figure 1; 2026-10-12aa).
Output: figures/Fig2_screen.{png,svg,pdf}. Data loading as in make_fig1.py.bak_20261012aa (same counting base as Table 1).
"""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import os
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
D = f"{GMO_ROOT}/analysis/refset/screen/seeds_ext/"
OUT = f"{GMO_ROOT}/paper/figures/"
GREY, VERM, BLUE = "#8c8c8c", "#D55E00", "#0072B2"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "font.size": 7.5, "mathtext.fontset": "custom", "mathtext.rm": "sans", "mathtext.it": "sans:italic", "mathtext.bf": "sans:bold",
                     "pdf.fonttype": 42, "svg.fonttype": "none"})
best = pd.read_csv(D + "seeds_ext_best.tsv", sep="\t")
best = best[best.dup_dropped == 0]
rec = pd.read_csv(D + "reciprocal/reciprocal.tsv", sep="\t")
rec = rec.assign(gene=rec.gene.str.split(",")).explode("gene")
rec["gene"] = rec.gene.str.strip()
rec = rec[rec.o_h > 0]
GENE_ORDER = ["sus", "cus", "chs", "acs", "acs_pKerr14a", "nos", "ocs", "vis", "ags", "mas1", "mas2", "mis"]
rows = [(g, "acs (pKerr14a)" if g == "acs_pKerr14a" else g) for g in GENE_ORDER]
labels_wanted = {"cus": ["Mycena", "Monosporozyma"], "chs": ["Fusarium"], "sus": ["Melanogaster"], "acs": ["Boletus", "Butyriboletus"]}

fig = plt.figure(figsize=(6.7, 5.4))
axb = fig.add_axes([0.20, 0.17, 0.66, 0.80])
rng = np.random.default_rng(1)
ypos = {}; y = 0; ticks = []; tl = []
for r in rows:
    ypos[r[0]] = y; ticks.append(y); tl.append(r[1]); y += 1
ymax = y
for g, ylab in rows:
    y0 = ypos[g]
    b = best[best.gene == g]
    nn = len(b)
    jit = rng.uniform(-0.28, 0.28, nn)
    axb.scatter(b.pident, y0 + jit, s=3 if nn > 100 else 8, c=GREY, alpha=0.25 if nn > 100 else 0.5, linewidths=0, rasterized=True, zorder=2)
    rr = rec[rec.gene == g]
    nc = rr[rr.o_bact_pid >= 90]; nonc = rr[rr.o_bact_pid < 90]
    axb.scatter(nonc.pident_seed, y0 + rng.uniform(-0.12, 0.12, len(nonc)), s=22, c=VERM, ec="white", lw=0.4, zorder=4)
    axb.scatter(nc.pident_seed, y0 + rng.uniform(-0.12, 0.12, len(nc)), s=30, c=BLUE, marker="x", lw=1.2, zorder=5)
    axb.text(1.02, y0, f"n = {nn:,}", transform=axb.get_yaxis_transform(), va="center", ha="left", fontsize=7)
    print(g, nn, f"non-fungal nearest {len(nonc)}, excluded as likely contamination {len(nc)}")
    for k, gen in enumerate(labels_wanted.get(g, [])):
        s = rr[rr.organism.str.startswith(gen)]
        if s.empty:
            s = b[b.organism.str.startswith(gen)].rename(columns={"pident": "pident_seed"})
        if s.empty:
            print("no label", g, gen); continue
        p = s.loc[s.pident_seed.idxmax()]
        up, dx, ha = dict(Mycena=(True, 1, "left"), Hebeloma=(True, -1, "right"), Monosporozyma=(False, -3, "right")).get(gen, (k % 2 == 0, 0, "center"))
        axb.annotate(gen, (p.pident_seed, y0), xytext=(p.pident_seed + dx, y0 + (0.42 if up else -0.42)), fontsize=6, style="italic", ha=ha,
                     va="bottom" if up else "top", zorder=6, arrowprops=dict(arrowstyle="-", lw=0.4, color="#555", shrinkA=0, shrinkB=2))
axb.set_ylim(ymax - 0.4, -0.7); axb.set_xlim(0, 70)
axb.set_yticks(ticks); axb.set_yticklabels([])
for t, lab in zip(ticks, tl):
    if "(" in lab:
        axb.text(-0.015, t, r"$\it{acs}$ (pKerr14a)", transform=axb.get_yaxis_transform(), ha="right", va="center", fontsize=8)
    else:
        axb.text(-0.015, t, lab, transform=axb.get_yaxis_transform(), ha="right", va="center", fontsize=8, style="italic")
axb.set_xlabel("Identity of best fungal hit to seed protein (%)")
axb.grid(axis="x", color="#e6e6e6", lw=0.5); axb.set_axisbelow(True)
for s in ("top", "right"): axb.spines[s].set_visible(False)
axb.tick_params(axis="y", length=0)
h = [Line2D([], [], marker="o", ls="", ms=3, color=GREY, alpha=0.6, label="forward search: best hit per assembly"),
     Line2D([], [], marker="o", ls="", ms=5, color=VERM, mec="white", label="reverse search: closest relatives non-fungal (bacterial or plant)"),
     Line2D([], [], marker="x", ls="", ms=5, color=BLUE, mew=1.2, label="excluded as likely contamination (≥90% identical to a protein of a non-rhizobial soil bacterium)")]
fig.legend(handles=h, loc="lower left", bbox_to_anchor=(0.02, 0.0), fontsize=6.5, frameon=False, ncol=1, handletextpad=0.3, labelspacing=0.3)
for ext in ("png", "svg", "pdf"):
    fig.savefig(OUT + "Fig2_screen." + ext, dpi=int(os.environ.get("FIG_DPI", "300")))
print("written Fig2_screen.{png,svg,pdf}")
