"""Figure 1 — workflow of the screen (scheme only; 2026-10-12aa).
The dot plot that used to be panel (b) is now Figure 2 (make_fig2.py). Previous combined version: make_fig1.py.bak_20261012aa.
Output: figures/Fig1_workflow.{png,svg,pdf}
"""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import textwrap
OUT = f"{GMO_ROOT}/paper/figures/"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "font.size": 7.5, "mathtext.fontset": "custom", "mathtext.rm": "sans", "mathtext.it": "sans:italic", "mathtext.bf": "sans:bold",
                     "pdf.fonttype": 42, "svg.fonttype": "none"})
HI = "#fde5d3"
fig = plt.figure(figsize=(6.7, 4.6))
ax = fig.add_axes([0.01, 0.01, 0.98, 0.98]); ax.set_xlim(0, 100); ax.set_ylim(0, 60); ax.axis("off")

def box(x, y, w, h, t, wrap, fc="#f2f2f2", fs=7, bold=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.8", fc=fc, ec="#444", lw=0.7))
    txt = "\n".join(l for p in t.split("\n") for l in (textwrap.wrap(p, wrap, break_on_hyphens=False) or [""]))
    ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=fs, linespacing=1.15, fontweight="bold" if bold else "normal")
def arr(x0, y0, x1, y1):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="-|>", lw=0.8, color="#333", shrinkA=0, shrinkB=0))
def ln(xs, ys):
    ax.plot(xs, ys, color="#333", lw=0.8, solid_capstyle="butt")

# seeds
box(15, 53, 70, 6, "11 T-DNA opine synthesis genes (56 seed proteins from Ti/Ri plasmids and plant cT-DNA)", 95)
# left column: protein-level screen
LX, LW, LH = 1, 46, 6.2
ys = [44.5, 36.5, 28.5, 20.5]
box(LX, ys[0], LW, LH, "Protein-level screen: 6,123 annotated fungal proteomes", 52)
box(LX, ys[1], LW, LH, "1,277 proteins checked by reverse search against NCBI nr", 52)
box(LX, ys[2], LW, LH, "116 proteins with non-fungal closest relatives", 52)
box(LX, ys[3], LW, LH, "50 candidates of rhizobial origin:\n44 T-DNA candidates ($\\it{cus}$, $\\it{chs}$, $\\it{sus}$, $\\it{nos}$) +\n6 non-T-DNA candidates ($\\it{acs}$-like GDPD, $\\it{odh}$)", 95, fc=HI)
for y in ys[:3]:
    arr(LX + LW / 2, y, LX + LW / 2, y - 1.8)
arr(30, 53, 30, ys[0] + LH)
# right column: genome-level screen
RX, RW = 52, 47
ax.add_patch(FancyBboxPatch((RX, 20.5), RW, 30.2, boxstyle="round,pad=0,rounding_size=0.8", fc="white", ec="#444", lw=0.7))
ax.text(RX + RW / 2, 49.6, "Genome-level screen", ha="center", va="top", fontsize=7.5, fontweight="bold")
box(RX + 2, 35.5, RW - 4, 10, "2,215 genome assemblies (517 Boletales, 1,386 other orders, 324 $\\it{Fusarium}$ and relatives; 12 in both)", 52)
box(RX + 2, 22.5, RW - 4, 10, "Searched for the protein families of the candidates, all T-DNA genes and plasmid $\\it{vir}$/$\\it{tra}$/$\\it{rep}$ genes (HMM profiles, blastx, miniprot)", 52)
arr(70, 53, 70, 50.7)
arr(LX + LW, 23.6, RX, 23.6)
# bottom row
arr(12, ys[3], 12, 11.5)
ln([RX + RW / 2, RX + RW / 2], [20.5, 14.5]); ln([RX + RW / 2, 18], [14.5, 14.5]); arr(18, 14.5, 18, 11.5)
box(1, 1.5, 22, 10, "Gene models and catalogue: 2,864 loci", 24)
box(25, 1.5, 24, 10, "Tests against contamination: fungal flanks, independent assemblies, read depth, no plasmid genes", 26)
box(51, 1.5, 22, 10, "Lineage assignment by phylogeny (ML trees, AU tests)", 24)
box(75, 1.5, 24, 10, "Seven lineages of bacterial origin: $\\it{sus}$-like, $\\it{cus}$-like, $\\it{chs}$-like, acc-type GDPD, pKerr14a-type, $\\it{nos}$-subfamily ODH, $\\it{odh}$ island", 30, fc=HI, fs=6.0)
arr(23, 6.5, 25, 6.5); arr(49, 6.5, 51, 6.5); arr(73, 6.5, 75, 6.5)
for ext in ("png", "svg", "pdf"):
    fig.savefig(OUT + "Fig1_workflow." + ext, dpi=int(os.environ.get("FIG_DPI", "300")))
print("written Fig1_workflow.{png,svg,pdf}")
