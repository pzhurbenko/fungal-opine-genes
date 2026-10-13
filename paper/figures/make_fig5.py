import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import re, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
# Figure S11 (formerly Fig. 5) - selection on the lineages of bacterial origin (dN/dS).
# Sources: analysis/dnds/{SUS,CUS}.* and analysis/dnds/extra/<set>/<set>.* (journal 2026-10-12, tags dnds_extra, yn00_gap_caveat).
# acc-type GDPD (pseudogenes only) and nos-ODH of Rhypophila (no set computed) are absent -> not shown.
D = f"{GMO_ROOT}/analysis/dnds/"
X = D + "extra/"
OUT = f"{GMO_ROOT}/paper/supplementary/figures/"
VERM, BLUE, TEAL, GREEN, ORANGE, PURPLE, BROWN = "#D55E00", "#0072B2", "#56B4E9", "#009E73", "#E69F00", "#CC79A7", "#8C6D31"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "font.size": 7.5, "pdf.fonttype": 42, "svg.fonttype": "none", "mathtext.default": "regular"})
PAIR = re.compile(r'^\s*\d+\s+\d+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+([\d.]+)\s+[-\d.]+\s+\+-\s+[\d.]+\s+([-\d.]+)')

def pairs(path, ds_cap=None):
    """Pairwise yn00 estimates; same filter as parse_paml.py (dS > 0.01, omega < 5). Optional dS ceiling for saturation."""
    w, s = [], []
    for line in open(path):
        m = PAIR.match(line)
        if m and float(m.group(2)) > 0.01 and float(m.group(1)) < 5:
            w.append(float(m.group(1))); s.append(float(m.group(2)))
    w, s = np.array(w), np.array(s)
    total = len(w)
    if ds_cap is not None:
        w, s = w[s <= ds_cap], s[s <= ds_cap]
    return w, s, total

def omega_m0(path):
    return float(re.search(r'omega \(dN/dS\) =\s*([\d.]+)', open(path).read()).group(1))

HIST = [
    ("$\\it{sus}$-like\n(Boletales)", VERM, D + "SUS.yn00.out", D + "SUS.M0.out", 31, None),
    ("$\\it{cus}$-like\n(all orders)", BLUE, D + "CUS.yn00.out", D + "CUS.M0.out", 27, None),
    ("$\\it{cus}$-like\n(outside Boletales)", TEAL, X + "cus_nonBol/cus_nonBol.yn00.out", X + "cus_nonBol/cus_nonBol.M0.out", 20, 2.0),
    ("$\\it{chs}$-like\n($\\it{Fusarium}$)", GREEN, X + "chs/chs.yn00.out", X + "chs/chs.M0.out", 6, None),
]
DOTS = [
    ("pKerr14a-type\n$\\it{acs}$-like", ORANGE, X + "acs_pKerr14a/acs_pKerr14a.yn00.out", 6,
     [(omega_m0(X + "acs_pKerr14a/acs_pKerr14a.M0.out"), "99 codons"),
      (omega_m0(X + "acs_pKerr14a/loose/acs_pKerr14a.M0.out"), "219 codons")], ""),
    ("$\\it{nos}$-subfamily ODH\n(Boletales)", PURPLE, X + "nosODH_Bol/nosODH_Bol.yn00.out", 4,
     [(omega_m0(X + "nosODH_Bol/nosODH_Bol.M0.out"), None)], ""),
    ("odh island\n($\\it{Fusarium}$)", BROWN, X + "odh/odh.yn00.out", 5,
     [(omega_m0(X + "odh/odh.M0.out"), None)], "unreliable (dS $\\approx$ 0.02)"),
]

fig = plt.figure(figsize=(18 / 2.54, 11.0 / 2.54))
gs = fig.add_gridspec(2, 12, hspace=0.62, wspace=1.1, left=0.065, right=0.985, top=0.9, bottom=0.09)
TOP = [slice(0, 3), slice(3, 6), slice(6, 9), slice(9, 12)]
BOT = [slice(0, 4), slice(4, 8), slice(8, 12)]
bins = np.arange(0, 1.45, 0.04)
letters = iter("abcdefg")

for col, (lab, colr, ypath, mpath, ngenes, cap) in enumerate(HIST):
    ax = fig.add_subplot(gs[0, TOP[col]])
    w, s, total = pairs(ypath, cap)
    med, w0 = float(np.median(w)), omega_m0(mpath)
    ax.hist(w, bins=bins, color=colr, alpha=0.75, ec="white", lw=0.4)
    ymax = ax.get_ylim()[1] * 1.55; ax.set_ylim(0, ymax)
    ax.axvline(1, color="k", ls="--", lw=0.8, ymax=0.60)
    ax.axvline(med, color=colr, lw=1.2, ymax=0.60)
    ax.axvline(w0, color=colr, lw=1.2, ls=":", ymax=0.60)
    note = f"{ngenes} genes, {len(w)} pairs\nmedian = {med:.3f}\nM0 $\\omega$ = {w0:.2f}"
    if cap is not None:
        note += f"\n{total - len(w)}/{total} pairs with dS > {cap:.0f}\nexcluded (saturation)"
    ax.text(0.98, 0.99, note, ha="right", va="top", fontsize=6.4, color=colr, transform=ax.transAxes, linespacing=1.3)
    ax.set_title(lab, fontsize=7.8, color=colr, loc="left", pad=4, linespacing=1.25)
    ax.set_xlabel("Pairwise dN/dS (yn00)", fontsize=7)
    ax.set_xlim(0, 1.45); ax.set_xticks([0, 0.5, 1.0]); ax.tick_params(labelsize=7)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    ax.text(-0.17, 1.14, next(letters), transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")
    if col == 0: ax.set_ylabel("Number of gene pairs", fontsize=7)

rng = np.random.default_rng(7)
for col, (lab, colr, ypath, ngenes, m0s, flag) in enumerate(DOTS):
    ax = fig.add_subplot(gs[1, BOT[col]])
    w, s, total = pairs(ypath)
    med = float(np.median(w))
    ax.scatter(w, 0.35 + rng.uniform(-0.14, 0.14, len(w)), s=16, color=colr, alpha=0.85, ec="white", lw=0.4, zorder=3)
    ax.axvline(1, color="k", ls="--", lw=0.8, ymax=0.52)
    ax.axvline(med, color=colr, lw=1.2, zorder=2, ymax=0.52)
    styles = [":", (0, (1, 3))]
    txt = [f"{ngenes} genes, {len(w)} pairs", f"median = {med:.3f}"]
    for i, (w0, tag) in enumerate(m0s):
        ax.axvline(w0, color=colr, lw=1.1, ls=styles[i % 2], zorder=2, ymax=0.52)
        txt.append(f"M0 $\\omega$ = {w0:.2f}" + (f" ({tag})" if tag else ""))
    if flag: txt.append(flag)
    ax.text(0.98, 0.99, "\n".join(txt), ha="right", va="top", fontsize=6.4, color=colr, transform=ax.transAxes, linespacing=1.3)
    ax.set_title(lab, fontsize=7.8, color=colr, loc="left", pad=4, linespacing=1.25)
    ax.set_xlabel("Pairwise dN/dS (yn00)", fontsize=7)
    ax.set_xlim(0, 1.6); ax.set_xticks([0, 0.5, 1.0, 1.5]); ax.tick_params(labelsize=7)
    ax.set_ylim(0, 1.45); ax.set_yticks([])
    for sp in ("top", "right", "left"): ax.spines[sp].set_visible(False)
    ax.text(-0.17, 1.14, next(letters), transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")
    if col == 0: ax.set_ylabel("Gene pairs", fontsize=7)

for e in ("png", "svg", "pdf"):
    fig.savefig(OUT + "FigureS11_dnds." + e, dpi=300)
print("written", OUT + "FigureS11_dnds.png")
