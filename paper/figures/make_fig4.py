"""Figure 4: opine islands in Boletaceae genomes (panel a, common kb scale), the
strain-specific insertion on the B. edulis contig JBBPRY010001778 (panel b, same scale) and the
odh-pepP-orf3 island of Fusarium poae Fp133 vs. its donor locus in Agrobacterium larrymoorei (panel c, own scale).
Run: MPLBACKEND=Agg python make_fig4.py
Data: analysis/island_order/island_gene_order.tsv, analysis/islands/islands_genes.tsv (gene roles),
analysis/islands/island_JAHRGV010000066_genes.tsv (JAHRGV010000066 genes),
analysis/fusarium_odh_island/{scan/loci.tsv, fig4/fp133_window_genes.tsv, donor/NZ_JAVIZC010000002.1.gb} (panel c),\npaper/tables/TableS1_loci_catalog.tsv (final lineage), loci_pipeline/round2/loci/ACS_JBBPRY010001778_L1/augustus.gff3 (exons)."""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import os
import csv, re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mp

G = f'{GMO_ROOT}'
A = f'{G}/analysis'
OUT = f'{G}/paper/figures'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 6.5, 'axes.linewidth': 0.5,
                     'mathtext.fontset': 'dejavusans', 'hatch.linewidth': 0.5,
                     'svg.fonttype': 'none', 'pdf.fonttype': 42})
FS = 6.5
COL = {'sus': '#d62728', 'cus': '#1f77b4', 'acc': '#2ca02c', 'kerr': '#9467bd', 'mis': '#17becf',
       'unpl': '#9a9a3a', 'nos': '#ff7f0e', 'tn': '#000000', 'pepp': '#8c564b', 'g3': '#e377c2',
       'fung': '#d9d9d9', 'none': '#ffffff'}
LIN = {'other fungal GDPD (position unresolved)': 'unpl', 'acc-type GDPD': 'acc', 'cus-like (T-DNA type)': 'cus',
       'mis lineage': 'mis', 'mis lineage (not placed)': 'mis', 'sus-like (T-DNA type)': 'sus',
       'sus-like (T-DNA type) (not placed)': 'sus', 'pKerr14a-type acs-like': 'kerr'}

def tsv(p):
    return list(csv.DictReader(open(p), delimiter='\t'))

lin = {r['locus_id']: r['lineage_final'] for r in tsv(f'{G}/paper/tables/TableS1_loci_catalog.tsv')}
role = {}
for r in tsv(f'{A}/islands/islands_genes.tsv'):
    role[(r['contig'], r['gene'])] = (r['role'], r['top_title'])
for r in tsv(f'{A}/islands/island_JAHRGV010000066_genes.tsv'):
    role.setdefault(('JAHRGV010000066', r['gene']), (r['role'], r['top_title']))
gene_rows = {}
for r in tsv(f'{A}/island_order/island_gene_order.tsv'):
    gene_rows.setdefault(r['island_contig'], []).append(r)
# 066: neighbours g18/g19/g27 are in the islands/island_JAHRGV010000066_genes.tsv only (outside the island table)
for r in tsv(f'{A}/islands/island_JAHRGV010000066_genes.tsv'):
    if r['gene'] in ('g27',):
        gene_rows['JAHRGV010000066'].append({'island_contig': 'JAHRGV010000066', 'gene': r['gene'], 'start': r['start'],
            'end': r['end'], 'strand': r['strand'], 'exons': r['exons'], 'kind': 'neighbor', 'product': '', 'plasmid_gene': ''})

def build(contig, w0, w1):
    """list of dict(start,end,strand,cat,intact,hatch,key)"""
    rows = [r for r in gene_rows[contig] if int(r['end']) >= w0 and int(r['start']) <= w1]
    reg = [r for r in rows if r['kind'].startswith('registry')]
    out = []
    for r in reg:
        lid = r['gene'][2:]
        st = r['kind'].split(':')[-1]
        if r['kind'].startswith('registry:') and r['gene'][2:] not in lin: cat = 'none'
        else: cat = LIN[lin[lid]]
        out.append(dict(start=int(r['start']), end=int(r['end']), strand=r['strand'] or '+', cat=cat,
                        intact=(st == 'intact'), hatch=False, key=lid, dbg=lin.get(lid, '')))
    for r in rows:
        if r['kind'].startswith('registry'): continue
        s, e = int(r['start']), int(r['end'])
        ov = max((min(e, x['end']) - max(s, x['start']) for x in out if x['key'] and x['key'].count('_') >= 2), default=0)
        if ov > 0.3 * (e - s) and ov > 0: continue   # same gene covered by the registry locus
        rl, title = role.get((contig, r['gene']), (None, ''))
        t = (title or '').lower() + ' ' + (r.get('product') or '').lower()
        hatch = False
        if 'dehydrogenase' in t and ('nopaline' in t or 'octopine' in t): cat = 'nos'
        elif 'transposase' in t or 'integrase' in t: cat = 'tn'
        elif contig == 'JAHRGV010000066' and r['gene'] in ('g20', 'g26', 'g27'): cat, hatch = 'fung', True
        elif rl in ('Fungi', 'anchor', 'bacterial_like'): cat = 'fung'
        else: cat = 'none'
        out.append(dict(start=s, end=e, strand=r['strand'] or '+', cat=cat, intact=True, hatch=hatch, key=r['gene'], dbg=''))
    return sorted(out, key=lambda x: x['start'])

def pale(c, f=0.55):
    c = np.array(matplotlib.colors.to_rgb(c)); return tuple(c + (1 - c) * f)

def arrow(ax, g, x0, h=0.17, head=0.45, y0=0.0):
    s, e = g['start'] - x0, g['end'] - x0
    s, e = s / 1000, e / 1000
    L = e - s; hl = min(head, 0.55 * L)
    if g['strand'] == '-': pts = [(e, y0 - h), (s + hl, y0 - h), (s, y0), (s + hl, y0 + h), (e, y0 + h)]
    else: pts = [(s, y0 - h), (e - hl, y0 - h), (e, y0), (e - hl, y0 + h), (s, y0 + h)]
    c = COL[g['cat']]
    if g['cat'] == 'none': fc, ec, ls = 'white', '#8c8c8c', '-'
    elif g['cat'] == 'fung': fc, ec, ls = c, '#8c8c8c', '-'
    elif g['intact']: fc, ec, ls = c, c, '-'
    else: fc, ec, ls = pale(c), c, (0, (2, 1.2))
    if g['hatch']:
        ax.add_patch(mp.Polygon(pts, closed=True, fc='white', ec='#555555', hatch='//////', lw=0.6, zorder=3))
    else:
        ax.add_patch(mp.Polygon(pts, closed=True, fc=fc, ec=ec, lw=0.7, ls=ls, zorder=3))

# ---------------- layout (cm) ----------------
CM = 1 / 2.54
FW, FH = 17.6, 22.0
fig = plt.figure(figsize=(FW * CM, FH * CM), dpi=int(os.environ.get("FIG_DPI", "300")))
LM, W60 = 3.2, 12.4            # left margin, cm for 60 kb
KB = W60 / 60.0                 # cm per kb  (common to all tracks)
def axcm(x, ytop, w, h, **kw):
    return fig.add_axes([x / FW, 1 - (ytop + h) / FH, w / FW, h / FH], **kw)

tracks = [  # contig, window, species (genus, species), gca, labels
 ('JAHRGV010000066', (1144500, 1204500), ('Butyriboletus', 'roseoflavus'), 'GCA_019364595.1'),
 ('JAHRGV010000067', (1426000, 1480000), ('Butyriboletus', 'roseoflavus'), 'GCA_019364595.1'),
 ('JAHRGV010000044', (1144500, 1194500), ('Butyriboletus', 'roseoflavus'), 'GCA_019364595.1'),
 ('JBBPJK010000045', (0, 36000), ('Boletus', 'variipes'), 'GCA_038088715.1'),
 ('JBBPWX010000076', (1300, 38400), ('Tylopilus', 'sp. UT-M0000335'), 'GCA_041722695.1'),
]
# label spec: key -> (text, x offset kb of text centre relative to gene centre, side +1 above/-1 below, level 0/1)
LAB = {
 'JAHRGV010000066': {'g20': ('repeat', 0, 1, 0), 'ACS_JAHRGV010000066_L1': ('GDPD, unplaced\n(intact)', -1.5, -1, 0),
                     'g23': ('IS481\ntransposase', 0, 1, 1), 'CUS_JAHRGV010000066_L1': ('$\\it{cus}$-like\n(pseudogene)', 2.5, -1, 0),
                     'g26': ('repeat', 0, 1, 0), 'g27': ('repeat', 0, 1, 0)},
 'JAHRGV010000067': {'ACS_JAHRGV010000067_L1': ('acc-type GDPD\n(pseudogene)', -2.5, 1, 0),
                     'CUS_JAHRGV010000067_L1': ('$\\it{mis}$-lineage\ngene', 1.0, -1, 0), 'g30': ('nopaline-dehydrogenase-like', 6.0, 1, 0),
                     'SUS_JAHRGV010000067_L1': ('$\\it{sus}$-like\nfragment', 0, -1, 0)},
 'JAHRGV010000044': {'CUS_JAHRGV010000044_L1': ('$\\it{cus}$-like\n(pseudogene)', 0, 1, 0),
                     'CUS_JAHRGV010000044_L2': ('$\\it{cus}$-like\n(intact)', 0, 1, 0)},
 'JBBPJK010000045': {'CUS_JBBPJK010000045_L1': ('$\\it{mis}$-lineage\n(pseudogene)', -1.5, 1, 0),
                     'ACS_JBBPJK010000045_L1': ('acc-type GDPD\n(pseudogene)', 0, 1, 0)},
 'JBBPWX010000076': {'ACS_JBBPWX010000076_S2L1': ('acc-type GDPD\n(pseudogene)', -1.0, 1, 0),
                     'CUS_JBBPWX010000076_S2L1': ('$\\it{mis}$-lineage\n(fragment)', 3.0, -1, 0)},
}
PITCH = 1.8
y = 0.45
fig.text(0.1 / FW, 1 - 0.1 / FH, '(a)', fontsize=9, fontweight='bold', va='top')
axes_a = []
for i, (ctg, (w0, w1), (gn, sp), gca) in enumerate(tracks):
    wk = (w1 - w0) / 1000
    ax = axcm(LM, y, wk * KB, PITCH - 0.05)
    ax.set_xlim(0, wk); ax.set_ylim(-1, 1); ax.axis('off')
    ax.plot([0, wk], [0, 0], color='#bbbbbb', lw=0.5, zorder=1)
    gs = build(ctg, w0, w1)
    for g in gs:
        arrow(ax, g, w0, h=0.17)
        lab = LAB[ctg].get(g['key'])
        if lab:
            txt, dx, side, lv = lab
            xc = (g['start'] + g['end']) / 2000 - w0 / 1000
            ty = side * (0.52 + 0.0 * lv) if lv == 0 else side * 0.52
            ax.plot([xc, xc + dx * 0.5], [side * 0.19, side * (0.40)], color='#666666', lw=0.4, zorder=2)
            ax.text(xc + dx, side * 0.42, txt, ha='center', va='bottom' if side > 0 else 'top', fontsize=FS, linespacing=0.95,
                    color='black', zorder=4)
    # left label
    fig.text((LM - 0.2) / FW, 1 - (y + (PITCH - 0.05) / 2) / FH,
             (f'$\\it{{{gn}}}$\n$\\it{{{sp}}}$' if not sp.startswith('sp.') else f'$\\it{{{gn}}}$\nsp. UT-M0000335')
             + f'\n{ctg}\n{w0/1000:,.1f}\u2013{w1/1000:,.1f} kb',
             fontsize=FS, ha='right', va='center', linespacing=1.05)
    fig.text((LM + wk * KB + 0.15) / FW, 1 - (y + (PITCH - 0.05) / 2 - 0.0) / FH, f'{wk:.1f} kb', fontsize=FS, ha='left', va='center', fontweight='bold')
    axes_a.append(ax)
    y += PITCH
# ruler (common scale)
axr = axcm(LM, y - 0.05, 60 * KB, 0.5); axr.patch.set_visible(False)
axr.set_xlim(0, 60); axr.set_ylim(0, 1); axr.set_yticks([])
for s_ in ('top', 'right', 'left'): axr.spines[s_].set_visible(False)
axr.set_xticks(range(0, 61, 10)); axr.set_xticks(range(0, 61, 5), minor=True)
axr.tick_params(labelsize=FS, width=0.5, length=2.5, pad=1); axr.tick_params(which='minor', width=0.4, length=1.5)
axr.set_xlabel('Distance from the left end of each shown region, kb (common scale, all tracks)', fontsize=FS, labelpad=1)
y_b = y - 0.05 + 0.5 + 0.8
fig.text(0.1 / FW, 1 - (y_b - 0.55) / FH, '(b)', fontsize=9, fontweight='bold', va='top')
y_b += 0.15

# ---------------- panel b ----------------
wb = 7.415                      # contig length, kb
KBB = 0.85                      # cm per kb in panel b (own, larger scale)
MH = 1.8
y_m = y_b
axm = axcm(LM, y_m, wb * KBB, MH)
axm.set_xlim(0, wb); axm.set_ylim(0.7, 3.0); axm.axis('off')
# gene row (y 1.9-2.8), segment row
gff = [l.split('\t') for l in open(f'{G}/analysis/loci_pipeline/round2/loci/ACS_JBBPRY010001778_L1/augustus.gff3') if '\tCDS\t' in l]
ex = sorted((int(r[3]), int(r[4])) for r in gff if int(r[3]) > 2000)
gs_, ge_ = 2690, 6606
axm.plot([gs_ / 1000, ge_ / 1000], [2.35, 2.35], color=COL['kerr'], lw=0.6, zorder=2)
for s, e in ex: axm.add_patch(mp.Rectangle((s / 1000, 2.0), (e - s + 1) / 1000, 0.7, fc=COL['kerr'], ec=COL['kerr'], lw=0.4, zorder=3))
axm.add_patch(mp.Rectangle((0.051, 2.0), 0.938, 0.7, fc='white', ec='#555555', hatch='//////', lw=0.6, zorder=3))
axm.add_patch(mp.Rectangle((7.202, 2.0), 0.213, 0.7, fc='white', ec='#555555', hatch='//////', lw=0.6, zorder=3))
axm.text(0.52, 2.8, 'repeat', ha='center', va='bottom', fontsize=FS)
axm.text(wb + 0.4, 2.35, 'pKerr14a-type $\\it{acs}$-like, intact (5 exons)', ha='left', va='center', fontsize=FS, clip_on=False)
axm.plot([6.65, wb + 0.3], [2.35, 2.35], color='#666666', lw=0.4, clip_on=False)
# segment row
axm.add_patch(mp.Rectangle((0, 0.95), 1.313, 0.5, fc='#bfe3b4', ec='none'))
axm.add_patch(mp.Rectangle((7.202, 0.95), 0.213, 0.5, fc='#bfe3b4', ec='none'))
axm.add_patch(mp.Rectangle((1.313, 0.95), 5.889, 0.5, fc='#f3a58a', ec='none'))
axm.text(wb + 0.4, 1.5, 'segment 1,313–7,202 bp (5.9 kb) absent from\nfive other assemblies of the species', fontsize=FS, va='center', clip_on=False, color='#a8472a')
axm.text(wb + 0.4, 0.83, 'flanks (green) present in other assemblies', fontsize=FS, va='center', clip_on=False, color='#3d7a30')
axm.set_clip_on(False)
axm.axvline(1.313, color='k', ls=':', lw=0.8, ymax=1)
# x axis for panel b
axx = axcm(LM, y_m + MH + 0.0, wb * KBB, 0.5)
axx.set_xlim(0, wb); axx.set_ylim(0, 1); axx.set_yticks([])
for s_ in ('top', 'right', 'left'): axx.spines[s_].set_visible(False)
axx.set_xticks(range(0, 8, 2)); axx.set_xticks(range(0, 8), minor=True); axx.tick_params(labelsize=FS, width=0.5, length=2.5, pad=1)
axx.set_xlabel('Position on contig, kb', fontsize=FS, labelpad=1)
fig.text((LM - 0.2) / FW, 1 - (y_m + MH * 0.55) / FH, '$\\it{Boletus\\ edulis}$\nJBBPRY010001778\n7.4 kb', fontsize=FS, ha='right', va='center', linespacing=1.05)

# ---------------- panel c: F. poae Fp133 odh island vs. donor locus ----------------
FA = f'{A}/fusarium_odh_island'
y_c = y_m + MH + 1.25
fig.text(0.1 / FW, 1 - y_c / FH, '(c)', fontsize=9, fontweight='bold', va='top')
WC = 26.4; OFS = 13.0                  # window, kb; cassette starts OFS kb from the left window edge
KBC = 12.4 / WC                        # cm per kb in panel c (common to the three tracks)
HC = 5.0
axc = axcm(LM, y_c + 0.25, WC * KBC, HC)
axc.set_xlim(0, WC); axc.set_ylim(HC, 0); axc.axis('off')
YT = [1.45, 4.0]                       # baselines (cm from the top of the axes): donor, chr2
HH = 0.17                              # arrow half-height, cm
def view(ctg_or_none, s, e, st, w0, flip):
    """genomic bp -> view bp (bp from the left window edge); flip = reverse-complemented view"""
    if flip: s, e = w1c_[0] - e, w1c_[0] - s; st = '+' if st == '-' else '-'
    else: s, e = s - w0, e - w0
    return s, e, st
def gdict(s, e, st, cat, intact=True):
    return dict(start=s, end=e, strand=st, cat=cat, intact=intact, hatch=False)
def row(spec, y0):
    for (s, e, st, cat, intact) in spec: arrow(axc, gdict(s * 1000, e * 1000, st, cat, intact), 0, h=HH, head=0.25, y0=y0)
    axc.plot([0, WC], [y0, y0], color='#bbbbbb', lw=0.5, zorder=1)
# donor locus (GenBank NZ_JAVIZC010000002.1)
from Bio import SeqIO
rec = SeqIO.read(f'{FA}/donor/NZ_JAVIZC010000002.1.gb', 'genbank')
d0 = 123871 - int(OFS * 1000); d1 = d0 + int(WC * 1000)
CAT_D = {'WP_309770756.1': 'g3', 'WP_309770757.1': 'pepp', 'WP_309770758.1': 'pepp', 'WP_309770759.1': 'nos', 'WP_309770762.1': 'nos'}
dspec = []; dpos = {}
for f in rec.features:
    if f.type != 'CDS': continue
    s, e = int(f.location.start) + 1, int(f.location.end)
    if s < d0 or e > d1: continue
    q = f.qualifiers; pid = q.get('protein_id', [''])[0]; prod = q.get('product', [''])[0].lower()
    st = '+' if f.location.strand == 1 else '-'
    if pid in CAT_D: cat = CAT_D[pid]
    elif 'transposase' in prod: cat = 'tn'
    else: cat = 'none'
    dspec.append(((s - d0) / 1000, (e - d0) / 1000, st, cat, 'pseudo' not in q))
    dpos[pid] = ((s - d0) / 1000, (e - d0) / 1000)
row([(a, b, st, c, ok) for a, b, st, c, ok in dspec], YT[0])
# Fp133 copies (cassette coordinates: scan/loci.tsv, GCA_046056315.1)
LOC = {r['locus_id']: r for r in tsv(f'{FA}/scan/loci.tsv') if r['gca'] == 'GCA_046056315.1'}
fpw = tsv(f'{FA}/fig4/fp133_window_genes.tsv')
CASS = {'FPOAC2_07382', 'FPOAC2_07383', 'FPOAC2_10362', 'FPOAC2_10363', 'FPOAC2_10364'}
def fp_track(ctg, ids, flip, y0):
    cs = {k: LOC[f'GCA_046056315.1_{v}'] for k, v in ids.items()}      # gene name -> locus row
    g3s = int(cs['g3']['start']) if not flip else int(cs['g3']['end'])
    if flip: w0_ = None; ref = g3s + int(OFS * 1000); w1c_[0] = ref; w0 = ref - int(WC * 1000)
    else: w0 = g3s - int(OFS * 1000); w1 = w0 + int(WC * 1000)
    lo, hi = (w0, w0 + int(WC * 1000))
    spec = []; pos = {}
    for r in fpw:
        if r['contig'] != ctg or r['locus_tag'] in CASS: continue
        s, e = int(r['start']), int(r['end'])
        if s < lo or e > hi: continue
        a, b, st = view(ctg, s, e, r['strand'], w0, flip); spec.append((a / 1000, b / 1000, st, 'fung', True))
    for k, r in cs.items():
        a, b, st = view(ctg, int(r['start']), int(r['end']), r['strand'], w0, flip)
        cat = {'odh': 'nos', 'pepP': 'pepp', 'g3': 'g3'}[k]
        spec.append((a / 1000, b / 1000, st, cat, r['status'] == 'intact')); pos[k] = (a / 1000, b / 1000, r['status'])
    row(spec, y0)
    return pos, (lo, hi)
w1c_ = [0]
pos2, win2 = fp_track('CM100255.1', {'odh': 'odh1', 'pepP': 'pepP1', 'g3': 'g31'}, False, YT[1])
# homology ribbons (donor - chr2 - chr3), cassette genes only
def ribbon(xa, ya, xb, yb, cat):
    c = pale(COL[cat], 0.6)
    axc.add_patch(mp.Polygon([(xa[0], ya + HH), (xa[1], ya + HH), (xb[1], yb - HH), (xb[0], yb - HH)], closed=True, fc=c, ec='none', alpha=0.75, zorder=0))
dg = {'g3': (dpos['WP_309770756.1'][0], dpos['WP_309770756.1'][1], 'g3'), 'pepP': (dpos['WP_309770757.1'][0], dpos['WP_309770758.1'][1], 'pepp'), 'odh': (dpos['WP_309770759.1'][0], dpos['WP_309770759.1'][1], 'nos')}
for k in ('g3', 'pepP', 'odh'):
    ribbon(dg[k][:2], YT[0], pos2[k][:2], YT[1], dg[k][2])
# gene labels
def glab(x, y, t, va='bottom', **kw): axc.text(x, y, t, ha=kw.pop('ha', 'center'), va=va, fontsize=FS, linespacing=0.95, **kw)
def lead(x, y1, y2): axc.plot([x, x], [y1, y2], color='#666666', lw=0.4, zorder=2)
cx = lambda k: (dg[k][0] + dg[k][1]) / 2
glab(cx('g3') - 1.7, YT[0] - 0.62, '$\\it{orf3}$-like\nWP_309770756'); lead(cx('g3'), YT[0] - 0.2, YT[0] - 0.6)
glab(cx('pepP'), YT[0] - 1.05, '$\\it{pepP}$-like\nWP_309770757/8'); lead(cx('pepP'), YT[0] - 0.2, YT[0] - 1.03)
glab(cx('odh') + 1.1, YT[0] - 0.62, '$\\it{odh}$\nWP_309770759'); lead(cx('odh'), YT[0] - 0.2, YT[0] - 0.6)
for k in ('g3', 'pepP', 'odh'):
    a_, b_, stt = pos2[k]; glab((a_ + b_) / 2, YT[1] + 0.2, {'g3': 'orf3'}.get(k, k), va='top', fontstyle='italic')
glab((dpos['WP_309770762.1'][0] + dpos['WP_309770762.1'][1]) / 2, YT[0] - 0.25, '$\\it{odh}$-like\nWP_309770762')
tn = [(a, b) for a, b, st, c, ok in dspec if c == 'tn']
left_is = [t for t in tn if t[1] < dg['g3'][0]]; right_is = [t for t in tn if t[0] > dg['odh'][1]]
for (a, b) in left_is: glab((a + b) / 2, YT[0] - 0.25, 'transposase')
xs0, xs1 = min(a for a, b in right_is), max(b for a, b in right_is)
axc.plot([xs0, xs1], [YT[0] + 0.3, YT[0] + 0.3], color='#444444', lw=0.6)
glab(xs0, YT[0] + 0.36, 'IS3, IS110, IS66-family transposases\n(4 of 5 are pseudogenes)', va='top', ha='left')
# scale bar and track labels
axc.plot([0.3, 5.3], [HC - 0.25, HC - 0.25], color='k', lw=0.8); axc.text(2.8, HC - 0.32, '5 kb', ha='center', va='bottom', fontsize=FS)
def tlab(y0, txt): fig.text((LM - 0.2) / FW, 1 - (y_c + 0.25 + y0) / FH, txt, fontsize=FS, ha='right', va='center', linespacing=1.05)
tlab(YT[0], '$\\it{Agrobacterium}$\n$\\it{larrymoorei}$ plasmid\nNZ_JAVIZC010000002.1\n%.1f–%.1f kb' % (d0 / 1000, d1 / 1000))
tlab(YT[1], '$\\it{Fusarium\\ poae}$ Fp133\nchr. 2 (CM100255.1)\n%.1f–%.1f kb\n%.1f kb to chr. end →' % (win2[0] / 1000, win2[1] / 1000, int(LOC['GCA_046056315.1_odh1']['dist_right']) / 1000))
print('panel c', win2, pos2, dpos.get('WP_309770756.1'))
y_c_end = y_c + 0.25 + HC

# ---------------- legend (bottom) ----------------
leg = [('$\\it{sus}$-like', 'sus', 'sus'), ('$\\it{cus}$-like', 'cus', 'cus'), ('acc-type GDPD', 'acc', 'acc'),
       ('pKerr14a-type $\\it{acs}$-like', 'kerr', 'kerr'), ('$\\it{mis}$ lineage', 'mis', 'mis'),
       ('GDPD, unplaced', 'unpl', 'unpl'), ('opine dehydrogenase (nos/ocs-like)', 'nos', 'nos'),
       ('transposase / integrase', 'tn', 'tn')]
def sw(ax, x, yy, cat, hatch=False, dashed=False):
    c = COL.get(cat, 'white')
    if hatch: ax.add_patch(mp.Rectangle((x, yy - 0.13), 0.5, 0.26, fc='white', ec='#555555', hatch='//////', lw=0.6))
    elif cat == 'fung': ax.add_patch(mp.Rectangle((x, yy - 0.13), 0.5, 0.26, fc=c, ec='#8c8c8c', lw=0.6))
    elif cat == 'none': ax.add_patch(mp.Rectangle((x, yy - 0.13), 0.5, 0.26, fc='white', ec='#8c8c8c', lw=0.6))
    elif dashed: ax.add_patch(mp.Rectangle((x, yy - 0.13), 0.5, 0.26, fc=pale(c), ec=c, lw=0.7, ls=(0, (2, 1.2))))
    else: ax.add_patch(mp.Rectangle((x, yy - 0.13), 0.5, 0.26, fc=c, ec=c, lw=0.7))
LX = LM + wb * KBB + 5.8     # legend left edge, cm (right of panel-b annotations)
axl = axcm(LX - 0.0, y_b + 0.0, FW - LX - 0.1, 0.1)
axl.axis('off')
fig_leg_y0 = y_c_end + 0.2
axl = axcm(0.3, fig_leg_y0, FW - 0.6, FH - fig_leg_y0 - 0.1)
axl.set_xlim(0, FW - 0.6); axl.set_ylim(0, FH - fig_leg_y0 - 0.1); axl.axis('off')
H = FH - fig_leg_y0 - 0.1
cols = [(0.0, leg[:3] + [('fungal gene (anchor / neighbour)', 'fung', 'fung')]),
        (4.6, leg[3:6] + [('no homologs / not classified', 'none', 'none')]),
        (9.0, leg[6:8]),
        (13.4, None)]
yy0 = H - 0.35
for x0, items in cols[:3]:
    for j, (t, cat, _) in enumerate(items):
        yy = yy0 - j * 0.38
        sw(axl, x0, yy, cat); axl.text(x0 + 0.65, yy, t, fontsize=FS, va='center')
# style key
x0 = 9.0; yy = yy0 - 2 * 0.38
sw(axl, x0, yy, 'cus'); sw(axl, x0 + 0.62, yy, 'cus', dashed=True); axl.text(x0 + 1.4, yy, 'intact / pseudogene or fragment', fontsize=FS, va='center')
sw(axl, 14.0, yy0, 'x', hatch=True); axl.text(14.65, yy0, 'repeat element', fontsize=FS, va='center')
axl.text(14.0, yy0 - 0.38, 'arrow = strand', fontsize=FS, va='top', color='#444444')
sw(axl, x0, yy0 - 3 * 0.38, 'pepp'); axl.text(x0 + 0.65, yy0 - 3 * 0.38, '$\\it{pepP}$ (Xaa-Pro aminopeptidase)', fontsize=FS, va='center')
sw(axl, x0, yy0 - 4 * 0.38, 'g3'); axl.text(x0 + 0.65, yy0 - 4 * 0.38, '$\\it{orf3}$ (hypothetical)', fontsize=FS, va='center')
fig.savefig(f'{OUT}/Fig4_islands.png', dpi=int(os.environ.get("FIG_DPI", "300")))
fig.savefig(f'{OUT}/Fig4_islands.svg')
fig.savefig(f'{OUT}/Fig4_islands.pdf')
print('ok', y, y_b, y_m, fig_leg_y0)
