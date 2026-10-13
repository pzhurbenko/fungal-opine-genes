"""Supplementary Figure S13: the odh-pepP-orf3 island of Fusarium poae Fp133 on both chromosome copies (chr. 2 and chr. 3)
vs. the homologous locus of the Agrobacterium larrymoorei plasmid. Formerly the full panel (c) of Figure 4.
Run: MPLBACKEND=Agg python make_figS13_fp133.py"""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import csv, re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mp

G = f'{GMO_ROOT}'
A = f'{G}/analysis'
OUT = f'{G}/paper/supplementary/figures'
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

CM = 1 / 2.54
LM = 3.2
FW, FH = 17.0, 10.4
fig = plt.figure(figsize=(FW * CM, FH * CM), dpi=300)
def axcm(x, ytop, w, h, **kw):
    return fig.add_axes([x / FW, 1 - (ytop + h) / FH, w / FW, h / FH], **kw)
y_c = 0.0
FA = f'{A}/fusarium_odh_island'
WC = 26.4; OFS = 13.0                  # window, kb; cassette starts OFS kb from the left window edge
KBC = 12.4 / WC                        # cm per kb in panel c (common to the three tracks)
HC = 8.4
axc = axcm(LM, y_c + 0.25, WC * KBC, HC)
axc.set_xlim(0, WC); axc.set_ylim(HC, 0); axc.axis('off')
YT = [1.45, 4.35, 6.95]                # baselines (cm from the top of the axes): donor, chr2, chr3
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
pos3, win3 = fp_track('CM100256.1', {'odh': 'odh2', 'pepP': 'pepP2', 'g3': 'g32'}, True, YT[2])
# homology ribbons (donor - chr2 - chr3), cassette genes only
def ribbon(xa, ya, xb, yb, cat):
    c = pale(COL[cat], 0.6)
    axc.add_patch(mp.Polygon([(xa[0], ya + HH), (xa[1], ya + HH), (xb[1], yb - HH), (xb[0], yb - HH)], closed=True, fc=c, ec='none', alpha=0.75, zorder=0))
dg = {'g3': (dpos['WP_309770756.1'][0], dpos['WP_309770756.1'][1], 'g3'), 'pepP': (dpos['WP_309770757.1'][0], dpos['WP_309770758.1'][1], 'pepp'), 'odh': (dpos['WP_309770759.1'][0], dpos['WP_309770759.1'][1], 'nos')}
for k in ('g3', 'pepP', 'odh'):
    ribbon(dg[k][:2], YT[0], pos2[k][:2], YT[1], dg[k][2])
    ribbon(pos2[k][:2], YT[1], pos3[k][:2], YT[2], dg[k][2])
# gene labels
def glab(x, y, t, va='bottom', **kw): axc.text(x, y, t, ha=kw.pop('ha', 'center'), va=va, fontsize=FS, linespacing=0.95, **kw)
def lead(x, y1, y2): axc.plot([x, x], [y1, y2], color='#666666', lw=0.4, zorder=2)
cx = lambda k: (dg[k][0] + dg[k][1]) / 2
glab(cx('g3') - 1.7, YT[0] - 0.62, '$\\it{orf3}$-like\nWP_309770756'); lead(cx('g3'), YT[0] - 0.2, YT[0] - 0.6)
glab(cx('pepP'), YT[0] - 1.05, '$\\it{pepP}$-like\nWP_309770757/8'); lead(cx('pepP'), YT[0] - 0.2, YT[0] - 1.03)
glab(cx('odh') + 1.1, YT[0] - 0.62, '$\\it{odh}$\nWP_309770759'); lead(cx('odh'), YT[0] - 0.2, YT[0] - 0.6)
for k in ('g3', 'pepP', 'odh'):
    a_, b_, stt = pos3[k]; glab((a_ + b_) / 2, YT[2] + 0.2, {'g3': 'orf3'}.get(k, k), va='top', fontstyle='italic')
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
tlab(YT[2], '$\\it{Fusarium\\ poae}$ Fp133\nchr. 3 (CM100256.1)\n%.1f–%.1f kb\n(reverse complement)\n← %.1f kb to chr. end' % (win3[0] / 1000, win3[1] / 1000, int(LOC['GCA_046056315.1_g32']['dist_right']) / 1000))
print('panel c', win2, win3, pos2, pos3, dpos.get('WP_309770756.1'))
# legend
leg = [('$\\it{pepP}$ (Xaa-Pro aminopeptidase)', 'pepp'), ('$\\it{orf3}$ (hypothetical)', 'g3'), ('opine dehydrogenase (nos/ocs-like)', 'nos'),
       ('transposase / integrase', 'tn'), ('fungal neighbour gene', 'fung'), ('no homologs / not classified', 'none')]
yl = y_c + 0.25 + HC + 0.2
axl = axcm(0.3, yl, FW - 0.6, FH - yl - 0.1); axl.set_xlim(0, FW - 0.6); axl.set_ylim(0, FH - yl - 0.1); axl.axis('off')
Hh = FH - yl - 0.1
for j, (t, cat) in enumerate(leg):
    x0 = 0.0 if j < 3 else 6.0; yy = Hh - 0.35 - (j % 3) * 0.38
    c = COL[cat]
    axl.add_patch(mp.Rectangle((x0, yy - 0.13), 0.5, 0.26, fc=('white' if cat == 'none' else c), ec=('#8c8c8c' if cat in ('none', 'fung') else c), lw=0.7))
    axl.text(x0 + 0.65, yy, t, fontsize=FS, va='center')
axl.text(12.0, Hh - 0.35, 'dashed outline = pseudogene\narrow = strand', fontsize=FS, va='center', color='#444444')
# outline = pseudogene; arrow = strand', fontsize=FS, va='center', color='#444444')
fig.savefig(f'{OUT}/FigureS13_fp133_island_two_copies.png', dpi=300)
fig.savefig(f'{OUT}/FigureS13_fp133_island_two_copies.svg')
fig.savefig(f'{OUT}/FigureS13_fp133_island_two_copies.pdf')
print('ok')
