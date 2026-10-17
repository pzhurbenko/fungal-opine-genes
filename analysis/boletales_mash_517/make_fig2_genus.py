#!/usr/bin/env python3
"""Figure 2: distribution of the five acquired lineages across genera of Boletales.
Left: Mash similarity tree of 358 species-level clusters (517 assemblies) collapsed to genera of
Boletaceae and to families of the other Boletales (cladogram). Right: per-row presence of the five lineages
(circle area ~ number of clusters with the gene; filled fraction ~ share of those with an intact gene).
Run (env paper_build): cd analysis/boletales_mash_517 && MPLBACKEND=Agg python make_fig2_genus.py
Inputs: Boletales_517_mashtree_rooted.nwk, clusters.tsv, genomes.tsv, paper/tables/TableS1_all_lineages.tsv.
Outputs: paper/figures/Fig6_boletales_genera.{png,svg,pdf}, paper/figures/Fig2_genus_counts.tsv.
Row genus = genus of the accepted cluster name (not the NCBI name of the individual assembly)."""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import csv
from collections import Counter, defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mp
from ete3 import Tree

G = f'{GMO_ROOT}'
D = f'{G}/analysis/boletales_mash_517'
OUT = f'{G}/paper/figures'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 6.5, 'axes.linewidth': 0.5,
                     'mathtext.fontset': 'dejavusans', 'svg.fonttype': 'none', 'pdf.fonttype': 42})
FS = 6.5
# lineage colours as in make_fig4.py (COL['sus'], 'cus', 'acc', 'kerr', 'nos')
LIN = [('sus-like (T-DNA type)', 'sus-like', '#d62728'),
       ('cus-like (T-DNA type)', 'cus-like', '#1f77b4'),
       ('acc-type GDPD', 'acc-type\nGDPD', '#2ca02c'),
       ('pKerr14a-type acs-like', 'pKerr14a-\ntype\nacs-like', '#9467bd'),
       ('nos-subfamily opine dehydrogenase (Boletales)', 'nos-\nsubfamily\nODH', '#ff7f0e')]
# clusters whose cluster_name is 'extra' (assemblies added by hand; see make_figures_517.py `extra`)
EXTRA = {'C015': 'Boletus', 'C149': 'Boletus', 'C217': 'Retiboletus', 'C251': 'Boletus'}
# family of the non-Boletaceae genera: NCBI Taxonomy (efetch on the 517 organism tax_ids, 2026-10-12);
# Hydnomerulius is absent from the NCBI names of the assemblies (cluster name only) -> Paxillaceae (NCBI)
FAM = {'Scleroderma': 'Sclerodermataceae', 'Paxillus': 'Paxillaceae', 'Hydnomerulius': 'Paxillaceae',
       'Suillus': 'Suillaceae', 'Rhizopogon': 'Rhizopogonaceae', 'Leucogyrophana': 'Rhizopogonaceae',
       'Pisolithus': 'Pisolithaceae', 'Phlebopus': 'Boletinellaceae', 'Astraeus': 'Astraeaceae',
       'Melanogaster': 'Melanogastraceae', 'Serpula': 'Serpulaceae', 'Gyrodon': 'Gyrodontaceae',
       'Gyroporus': 'Gyroporaceae', 'Hygrophoropsis': 'Hygrophoropsidaceae', 'Coniophora': 'Coniophoraceae',
       'Chroogomphus': 'Gomphidiaceae'}
UNRES_BOL = 'Boletaceae, genus unresolved'
UNRES_ORD = 'Boletales, family unresolved'

rd = lambda f: list(csv.DictReader(open(f), delimiter='\t'))

# ---- clusters -> row group
clu_of, cgenus, cgroup, cfam = {}, {}, {}, {}
for r in rd(f'{D}/clusters.tsv'):
    c = r['cluster_id']
    for g in r['genomes'].split(';'):
        clu_of[g] = c
    if r['cluster_name'] == 'extra':
        gen = EXTRA[c]
    else:
        w = [x for x in r['cluster_name'].split() if x not in ('cf.', 'aff.')]
        gen = w[0]
    cgenus[c] = gen
    if gen == 'Boletaceae':
        cgroup[c], cfam[c] = UNRES_BOL, 'Boletaceae'
    elif gen == 'Boletales':
        cgroup[c], cfam[c] = UNRES_ORD, 'unresolved'
    elif gen in FAM:
        cgroup[c], cfam[c] = FAM[gen], FAM[gen]
    else:
        cgroup[c], cfam[c] = gen, 'Boletaceae'
assert len(cgroup) == 358, len(cgroup)

# ---- carriers (Table S1, placed loci of the five lineages, Boletales assemblies)
car = defaultdict(lambda: defaultdict(set))   # lineage -> cluster -> statuses
miss = []
for r in rd(f'{G}/paper/tables/TableS1_all_lineages.tsv'):
    if r['lineage'] not in [l[0] for l in LIN] or r['order'] != 'Boletales':
        continue
    if r['gca'] not in clu_of:
        miss.append(r['gca']); continue
    car[r['lineage']][clu_of[r['gca']]].add(r['status'])
if miss:
    print('WARNING Boletales carriers outside the 517 set:', sorted(set(miss)))

# ---- tree: one leaf per cluster, then one leaf per row group
t = Tree(f'{D}/Boletales_517_mashtree_rooted.nwk', format=1)
leaves = t.get_leaves()
assert len(leaves) == 517 and all(l.name in clu_of for l in leaves)
def largest_pure(tree, key, val):
    """leaf name from the largest clade whose leaves all have feature key == val"""
    best = None
    for n in tree.traverse('postorder'):
        if all(getattr(x, key) == val for x in n.get_leaves()) and (best is None or len(n) > len(best)):
            best = n
    return best.get_leaves()[0].name
# step 1: 517 assemblies -> 358 species clusters (representative = leaf of the largest pure clade)
for l in leaves:
    l.add_feature('clu', clu_of[l.name])
cmem = defaultdict(list)
for l in leaves:
    cmem[l.clu].append(l)
cl_mono = sum(len(ls) == 1 or set(t.get_common_ancestor(ls).get_leaves()) == set(ls) for ls in cmem.values())
crep = {c: largest_pure(t, 'clu', c) for c in cgroup}
t.prune(list(crep.values()), preserve_branch_length=True)
for l in t.get_leaves():
    l.name = l.clu
    l.add_feature('grp', cgroup[l.clu])
assert len(t) == 358
print('species clusters monophyletic on the 517 tree:', cl_mono, 'of 358')
# step 2: 358 clusters -> rows (genus of Boletaceae / family of other Boletales)
gleaves = defaultdict(list)
for l in t.get_leaves():
    gleaves[l.grp].append(l)
mono = {grp: len(ls) == 1 or set(t.get_common_ancestor(ls).get_leaves()) == set(ls) for grp, ls in gleaves.items()}
rep = {grp: largest_pure(t, 'grp', grp) for grp in gleaves}
t.prune(list(rep.values()), preserve_branch_length=False)
name2grp = {v: k for k, v in rep.items()}
for l in t.get_leaves():
    l.name = name2grp[l.name]
t.ladderize(direction=1)

# ---- per-row counts
ncl = Counter(cgroup.values())
nasm = Counter(cgroup[clu_of[g]] for g in clu_of)
gen_in = defaultdict(set)
for c, grp in cgroup.items():
    gen_in[grp].add(cgenus[c])
rows = [l.name for l in t.get_leaves()]
cnt = {}
for grp in rows:
    for lin, _, _ in LIN:
        cs = [c for c in car[lin] if cgroup[c] == grp]
        cnt[grp, lin] = (len(cs), sum('intact' in car[lin][c] for c in cs))

with open(f'{OUT}/Fig2_genus_counts.tsv', 'w') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    hdr = ['row', 'family', 'genera_in_row', 'n_clusters', 'n_assemblies', 'monophyletic_on_mash_tree']
    for s in ('sus-like', 'cus-like', 'acc-type_GDPD', 'pKerr14a-type_acs-like', 'nos-subfamily_ODH'):
        hdr += [f'{s}_k_any', f'{s}_k_intact']
    w.writerow(hdr)
    for grp in rows:
        fam = 'Boletaceae' if grp in gen_in and all(cfam[c] == 'Boletaceae' for c in cgroup if cgroup[c] == grp) else \
            ('unresolved' if grp == UNRES_ORD else grp)
        rec = [grp, fam, ','.join(sorted(gen_in[grp])), ncl[grp], nasm[grp],
               'NA' if grp in (UNRES_BOL, UNRES_ORD) else ('yes' if mono[grp] else 'no')]
        for lin, _, _ in LIN:
            rec += list(cnt[grp, lin])
        w.writerow(rec)

# ---- drawing
N = len(rows)
W_IN, H_IN = 170 / 25.4, 228 / 25.4
TOP_IN, BOT_IN = 0.42, 0.88               # header (column titles) and legend/footnote
fig = plt.figure(figsize=(W_IN, H_IN))
ax = fig.add_axes([0.005, BOT_IN / H_IN, 0.99, 1 - (TOP_IN + BOT_IN) / H_IN])
ax_h_in = H_IN - TOP_IN - BOT_IN
unit = ax_h_in / N                         # inches per row (1 data unit in y)
XW = (0.99 * W_IN) / unit                  # axes width in row units -> equal aspect
ax.set_xlim(0, XW); ax.set_ylim(N - 0.5, -0.5); ax.axis('off')

# cladogram: leaves aligned, internal node x by height in edges
TREE_W = 1.25 / unit
for n in t.traverse('postorder'):
    n.add_feature('h', 0 if n.is_leaf() else 1 + max(c.h for c in n.children))
H = t.h
yi = {g: i for i, g in enumerate(rows)}
def xpos(n): return TREE_W * (1 - n.h / H)
def ypos(n):
    if n.is_leaf(): return yi[n.name]
    ys = [ypos(c) for c in n.children]
    n.add_feature('y', (min(ys) + max(ys)) / 2)
    return n.y
ypos(t)
lw = 0.45
for n in t.traverse():
    if n.is_leaf():
        continue
    cy = [c.y if not c.is_leaf() else yi[c.name] for c in n.children]
    ax.plot([xpos(n)] * 2, [min(cy), max(cy)], color='#333333', lw=lw, solid_capstyle='butt')
    for c, y in zip(n.children, cy):
        ax.plot([xpos(n), xpos(c)], [y, y], color='#333333', lw=lw, solid_capstyle='butt')
ax.plot([xpos(t) - 0.6, xpos(t)], [t.y, t.y], color='#333333', lw=lw)

# row labels
LAB_X = TREE_W + 0.4
COL0 = LAB_X + 2.12 / unit                 # first lineage column centre
CW = 0.63 / unit                           # column spacing
for g in rows:
    y = yi[g]
    unres = g in (UNRES_BOL, UNRES_ORD)
    fam_row = g not in gen_in or not all(cfam[c] == 'Boletaceae' for c in cgroup if cgroup[c] == g)
    star = '' if (unres or mono[g]) else '*'
    if unres:
        txt, style, wt = g, 'normal', 'normal'
    elif fam_row:
        txt, style, wt = g, 'normal', 'bold'
    else:
        txt, style, wt = g, 'italic', 'normal'
    ax.text(LAB_X, y, txt, fontsize=FS, va='center', ha='left', style=style, weight=wt,
            color='#555555' if unres else 'black')
    # count suffix (upright), placed after the name using the rendered extent
    ax.annotate(f' ({ncl[g]}/{nasm[g]}){star}', xy=(1, 0.5), xycoords=ax.texts[-1], fontsize=FS,
                va='center', ha='left', color='#555555' if unres else 'black')
    if y % 2 == 0:
        ax.add_patch(mp.Rectangle((COL0 - CW * 0.5, y - 0.5), CW * 5, 1, fc='#f3f3f3', ec='none', zorder=0))

# circles
KMAX = max(k for (g, l), (k, ki) in cnt.items())
RMAX = 0.62                                 # radius (row units) for KMAX
for j, (lin, short, col) in enumerate(LIN):
    x = COL0 + j * CW
    ax.text(x, -0.8, short, rotation=0, ha='center', va='bottom', fontsize=FS, color=col, weight='bold',
            linespacing=1.0)
    for g in rows:
        k, ki = cnt[g, lin]
        if not k:
            continue
        y = yi[g]
        r = RMAX * (k / KMAX) ** 0.5
        cx = x - 0.13 / unit
        ax.add_patch(mp.Circle((cx, y), r, fc='white', ec=col, lw=0.6, zorder=3))
        if ki == k:
            ax.add_patch(mp.Circle((cx, y), r, fc=col, ec=col, lw=0.6, zorder=4))
        elif ki:
            ax.add_patch(mp.Wedge((cx, y), r, 90, 90 + 360 * ki / k, fc=col, ec='none', zorder=4))
        ax.text(cx + RMAX + 0.25, y, f'{k}/{ncl[g]}', fontsize=FS - 0.5, va='center',
                ha='left', color='#333333', zorder=5)

# legend / footnotes (figure coordinates, bottom)
fl = fig.add_axes([0.005, 0.0, 0.99, (BOT_IN - 0.02) / H_IN])
lh = BOT_IN - 0.02
fl.set_xlim(0, 0.99 * W_IN / unit); fl.set_ylim(lh / unit, 0); fl.axis('off')
yy = 0.8
fl.text(0.3, yy, 'Circle area: number of species clusters with the gene (any status)', fontsize=FS, va='center')
xs = 0.3 + 4.0 / unit
for i, k in enumerate(sorted({1, 2, KMAX})):
    xx = xs + i * 0.35 / unit
    fl.add_patch(mp.Circle((xx, yy), RMAX * (k / KMAX) ** 0.5, fc='#666666', ec='#666666', lw=0.6))
    fl.text(xx + 0.1 / unit, yy, str(k), fontsize=FS, va='center')
yy2 = yy + 1.9
fl.text(0.3, yy2, 'Fill: share of these clusters with an intact gene', fontsize=FS, va='center')
xs2 = 0.3 + 3.0 / unit
for i, (fr, lab) in enumerate(((1, 'all intact'), (0.5, 'some intact'), (0, 'pseudogenes / fragments only'))):
    xx = xs2 + i * 0.85 / unit
    rr = RMAX * 0.75
    fl.add_patch(mp.Circle((xx, yy2), rr, fc='white', ec='#666666', lw=0.6))
    if fr:
        fl.add_patch(mp.Wedge((xx, yy2), rr, 90, 90 + 360 * fr, fc='#666666', ec='none'))
    fl.text(xx + 0.1 / unit, yy2, lab, fontsize=FS, va='center')
fl.text(0.3, yy2 + 1.9, 'k/n: species clusters with the gene / species clusters in the row.', fontsize=FS, va='center')
fl.text(0.3, yy2 + 3.4, 'Rows: genus of Boletaceae (italic) or family (bold) (species clusters / assemblies).  '
        '* genus (or family) not monophyletic on the similarity tree.', fontsize=FS, va='center')
fl.text(0.3, yy2 + 4.9, 'Cladogram from the Mash genome-similarity tree (branch lengths not shown); '
        'a similarity tree, not a phylogeny.', fontsize=FS, va='center')

for ext, kw in (('png', dict(dpi=300)), ('svg', {}), ('pdf', {})):
    fig.savefig(f'{OUT}/Fig6_boletales_genera.{ext}', **kw)
print('rows', N, 'fig mm', 170, round(H_IN * 25.4), 'row height pt', round(unit * 72, 2), 'KMAX', KMAX)
print('non-monophyletic:', [g for g in rows if not mono[g] and g not in (UNRES_BOL, UNRES_ORD)])
