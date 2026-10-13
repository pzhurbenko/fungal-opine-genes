#!/usr/bin/env python3
"""Figure 5: distribution of the five lineages across genera of Boletaceae (grouped by subfamily) and other Boletales.
Input: Fig2_genus_counts.tsv (written by analysis/boletales_mash_517/make_fig2_genus.py).
Subfamily layout: Fig5_subfamilies_Tremble2024.tsv (Tremble et al. 2024, Mycologia 116:392-408; subfamily names as in Wu et al. 2014).
Run (PlantCad2): cd paper/figures && MPLBACKEND=Agg python make_fig5_genera.py
Outputs: Fig5_boletales_genera.{png,svg,pdf}"""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import os
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mp

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 6.5, 'axes.linewidth': 0.5,
                     'mathtext.fontset': 'dejavusans', 'svg.fonttype': 'none', 'pdf.fonttype': 42})
FS = 6.5
OUT = f'{GMO_ROOT}/paper/figures'
LIN = [('sus-like', 'sus-like', '#d62728'), ('cus-like', 'cus-like', '#1f77b4'),
       ('acc-type_GDPD', 'acc-type\nGDPD', '#2ca02c'),
       ('pKerr14a-type_acs-like', 'pKerr14a-type\nacs-like', '#9467bd'),
       ('nos-subfamily_ODH', 'nos-subfamily\nODH', '#ff7f0e')]
# genus -> subfamily from Fig5_subfamilies_Tremble2024.tsv (single source: Tremble et al. 2024, Fig. 2 of the preprint;
# genera absent from that tree, a single tentatively identified sample and the polyphyletic Gastroboletus -> 'not assigned')
SUBFAM_ORDER = ['Boletoideae', 'Suillelloideae', 'Leccinoideae', 'Xerocomoideae', 'Zangioideae',
                'Austroboletoideae', 'Chalciporoideae', 'Phylloboletelloideae']
_sf = {r['genus']: r['subfamily'] for r in csv.DictReader(open(f'{OUT}/Fig5_subfamilies_Tremble2024.tsv'), delimiter='\t')}
SUB = [(n, ' '.join(g for g, sf in _sf.items() if sf == n)) for n in SUBFAM_ORDER]
NA_GEN = [g for g, sf in _sf.items() if sf == 'not assigned']
UNRES = 'Boletaceae, genus unresolved'
NA_LABEL = 'Boletaceae, subfamily not assigned'
OTHER = 'other Boletales'

data = {}
for r in csv.DictReader(open(f'{OUT}/Fig2_genus_counts.tsv'), delimiter='\t'):
    data[r['row']] = r
def rec(r):
    n = int(r['n_clusters'])
    k = [(int(r[f'{l}_k_any']), int(r[f'{l}_k_intact'])) for l, _, _ in LIN]
    return n, k
known = {g for _, s in SUB for g in s.split()} | set(NA_GEN) | {UNRES}
other_rows = [g for g, r in data.items() if r['family'] != 'Boletaceae']
extra = [g for g, r in data.items() if r['family'] == 'Boletaceae' and g not in known]
print('Boletaceae rows not in the list (put in not assigned):', extra)
for g in known:
    if g not in data: print('WARNING: listed genus absent from tsv:', g)
def order(gs):
    gs = [g for g in gs if g in data]
    return sorted(gs, key=lambda g: (-int(data[g]['n_clusters']), g))
def has_gene(g):
    return any(int(data[g][f'{l}_k_any']) for l, _, _ in LIN)
OTHER_ROWS = {}   # collapsed row label -> (n genera, n clusters)
def split(gs):
    gs = order(gs)
    keep = [g for g in gs if has_gene(g)]
    rest = [g for g in gs if not has_gene(g)]
    if rest:
        lab = f'other genera ({len(rest)}): none\u200b' + '\u200b' * len(OTHER_ROWS)   # unique key, same display
        n = sum(int(data[g]['n_clusters']) for g in rest)
        OTHER_ROWS[lab] = (len(rest), n)
        data[lab] = {'n_clusters': str(n), **{f'{l}_k_any': '0' for l, _, _ in LIN}, **{f'{l}_k_intact': '0' for l, _, _ in LIN}}
        keep.append(lab)
    return keep
groups = []
for n, s in SUB:
    gs = split(s.split())
    # make the collapsed label unique per subfamily (dict keys) while keeping the display text
    groups.append((n, gs))
na = split(NA_GEN + extra) + [UNRES]
groups.append((NA_LABEL, na))
on, ok = 0, [[0, 0] for _ in LIN]
for g in other_rows:
    n, k = rec(data[g]); on += n
    for i, (a, b) in enumerate(k): ok[i][0] += a; ok[i][1] += b
print('other Boletales pooled from', len(other_rows), 'rows;', on, 'clusters')
R = {}
for n, gs in groups:
    for g in gs: R[g] = rec(data[g])
R[OTHER] = (on, [tuple(x) for x in ok])
KMAX = max(a for n, k in R.values() for a, b in k)

# layout (inches)
W = 170 / 25.4
PITCH, HDR, GAP = 0.125, 0.125, 0.07
TOP, BOT = 0.5, 0.93
items = []   # (kind, label, y_centre)
y = TOP
for n, gs in groups:
    y += GAP
    items.append(('hdr', n, y + HDR / 2)); y += HDR
    for g in gs:
        items.append(('row', g, y + PITCH / 2)); y += PITCH
y += GAP
items.append(('hdr', OTHER, y + PITCH / 2)); y += PITCH  # placeholder, replaced below
items.pop()
items.append(('row', OTHER, y + PITCH / 2)); y += PITCH
H = y + BOT
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.axis('off')
X_BR, X_HDR, X_GEN, X_N = 0.12, 0.2, 0.3, 2.45
X0, STEP = 3.0, 0.77
RMAX = 0.068
# column headers aligned to the top
for j, (lin, short, col) in enumerate(LIN):
    ax.text(X0 + j * STEP, 0.12, short, ha='center', va='top', fontsize=FS, color=col, weight='bold', linespacing=1.0)
ax.text(X_N, 0.12, 'Species', ha='right', va='top', fontsize=FS, weight='bold')
ax.plot([0.05, W - 0.05], [TOP - 0.02] * 2, color='#333333', lw=0.5)
stripe = 0
for kind, lab, yc in items:
    if kind == 'hdr':
        ax.text(X_HDR, yc, lab, ha='left', va='center', fontsize=FS, weight='bold', color='#333333')
        continue
    n, k = R[lab]
    if stripe % 2 == 0:
        ax.add_patch(mp.Rectangle((0.03, yc - PITCH / 2), W - 0.06, PITCH, fc='#f0f0f0', ec='none', zorder=0))
    stripe += 1
    if lab == OTHER:
        ax.text(X_HDR, yc, lab, ha='left', va='center', fontsize=FS, weight='bold')
    elif lab == UNRES:
        ax.text(X_GEN, yc, lab, ha='left', va='center', fontsize=FS, color='#555555')
    elif lab in OTHER_ROWS:
        ax.text(X_GEN, yc, lab.replace('\u200b', ''), ha='left', va='center', fontsize=FS, color='#555555')
    else:
        ax.text(X_GEN, yc, lab, ha='left', va='center', fontsize=FS, style='italic')
    ax.text(X_N, yc, str(n), ha='right', va='center', fontsize=FS)
    for j, ((a, b), (_, _, col)) in enumerate(zip(k, LIN)):
        if not a: continue
        r = RMAX * (a / KMAX) ** 0.5
        cx = X0 + j * STEP
        ax.add_patch(mp.Circle((cx, yc), r, fc='white', ec=col, lw=0.6, zorder=3))
        if b == a:
            ax.add_patch(mp.Circle((cx, yc), r, fc=col, ec=col, lw=0.6, zorder=4))
        elif b:
            # y axis is inverted: angles measured in data space, so mirror to start at 12 o'clock going clockwise
            ax.add_patch(mp.Wedge((cx, yc), r, -90 - 360 * b / a, -90, fc=col, ec='none', zorder=4))
# subfamily brackets
i = 0
for n, gs in groups:
    ys = [yc for kind, lab, yc in items if lab in gs and kind == 'row'] if gs else []
    hy = [yc for kind, lab, yc in items if kind == 'hdr' and lab == n][0]
    y1, y2 = hy - HDR / 2 + 0.02, ys[-1] + PITCH / 2 - 0.02
    ax.plot([X_BR, X_BR], [y1, y2], color='#333333', lw=0.9, solid_capstyle='butt')
# legend
ly = H - BOT + 0.15
ax.text(0.1, ly, 'Circle area: number of species with the gene', fontsize=FS, va='center')
for i, k in enumerate((1, 3, KMAX)):
    xx = 2.75 + i * 0.5
    rk = RMAX * (k / KMAX) ** 0.5
    ax.add_patch(mp.Circle((xx, ly), rk, fc='#666666', ec='#666666', lw=0.6))
    ax.text(xx + rk + 0.04, ly, str(k), fontsize=FS, va='center')
ly2 = ly + 0.22
ax.text(0.1, ly2, 'Fill: share of these species with an intact gene', fontsize=FS, va='center')
for i, (fr, lab) in enumerate(((1, 'all intact'), (0.5, 'some intact'), (0, 'pseudogenes / fragments only'))):
    xx = 2.75 + i * 0.85
    ax.add_patch(mp.Circle((xx, ly2), 0.055, fc='white', ec='#666666', lw=0.6))
    if fr: ax.add_patch(mp.Wedge((xx, ly2), 0.055, -90 - 360 * fr, -90, fc='#666666', ec='none'))
    ax.text(xx + 0.09, ly2, lab, fontsize=FS, va='center')
ax.text(0.1, ly2 + 0.22, 'Genera are grouped by subfamily; species = species-level clusters; genera without any gene are pooled per subfamily',
        fontsize=FS - 0.5, va='center', color='#333333')
ax.text(0.1, ly2 + 0.40, '(other genera: number of genera; species = their species); other Boletales = all other families and unresolved.',
        fontsize=FS - 0.5, va='center', color='#333333')
for ext, kw in (('png', dict(dpi=int(os.environ.get("FIG_DPI", "300")))), ('svg', {}), ('pdf', {})):
    fig.savefig(f'{OUT}/Fig5_boletales_genera.{ext}', **kw)
print('rows', stripe, 'KMAX', KMAX, 'fig mm', 170, round(H * 25.4))
