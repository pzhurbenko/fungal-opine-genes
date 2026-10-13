"""Supplementary Figure S10: all opine islands of Boletaceae = contigs carrying >=2 loci of bacterial type
(Table S12), one track per contig, common kb scale (as Figure 4a).
Run: MPLBACKEND=Agg python make_figS10_islands.py
Data: paper/supplementary/TableS12_Gene_content_of_the_islands.tsv (genes, focal loci), TableS1b_Complete_catalogue_of_loci.tsv
(final lineage/status/species by locus_id), analysis/islands/islands_genes.tsv (transposase / dehydrogenase annotation of neighbours).
Repeats are NOT drawn: repeat annotation exists only for the five contigs of Figure 4a (hand-curated), not for all contigs."""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import csv, re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mp

G = f'{GMO_ROOT}'
OUT = f'{G}/paper/supplementary/figures'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 6.5, 'axes.linewidth': 0.5,
                     'mathtext.fontset': 'dejavusans', 'svg.fonttype': 'none', 'pdf.fonttype': 42})
FS = 6.5
COL = {'other': '#bbbbbb', 'sus': '#d62728', 'cus': '#1f77b4', 'acc': '#2ca02c', 'kerr': '#9467bd', 'mis': '#17becf',
       'unpl': '#9a9a3a', 'nos': '#ff7f0e', 'tn': '#000000', 'fung': '#d9d9d9'}
def LINCAT(l):
    if l.startswith('sus-like'): return 'sus'
    if l.startswith('cus-like'): return 'cus'
    if l.startswith('acc-type'): return 'acc'
    if l.startswith('pKerr14a'): return 'kerr'
    if l.startswith(('mis lineage', 'mis group')): return 'mis'
    if l.startswith('GDPD of unresolved'): return 'unpl'
    return 'other'
def tsv(p): return list(csv.DictReader(open(p), delimiter='\t'))
def pale(c, f=0.55):
    c = np.array(matplotlib.colors.to_rgb(c)); return tuple(c + (1 - c) * f)

S12 = tsv(f'{G}/paper/supplementary/TableS12_Gene_content_of_the_islands.tsv')
CAT = {r['locus_id']: r for r in tsv(f'{G}/paper/supplementary/TableS1b_Complete_catalogue_of_loci.tsv')}
ROLE = {(r['contig'], r['gene']): r['top_title'] for r in tsv(f'{G}/analysis/islands/islands_genes.tsv')}
asm_org = {}
for r in CAT.values():
    if r['assembly'] and len(r['organism']) > len(asm_org.get(r['assembly'], '')): asm_org[r['assembly']] = r['organism']

by = {}
for r in S12: by.setdefault(r['contig'], []).append(r)
# nos-subfamily opine dehydrogenase loci (blastx, not in Table S12): focal loci from Table S1b (locus_id ODH_<contig>.1)
ISL = {}
for r in tsv(f'{G}/analysis/islands/islands_genes.tsv'): ISL.setdefault(r['contig'], {}).setdefault(r['gene'], r)
for lid, c in CAT.items():
    if not c['lineage'].startswith('nos-subfamily opine dehydrogenase (Boletales)'): continue
    ctg = c['contig'].rsplit('.', 1)[0]
    if ctg not in by:   # contig absent from Table S12: other loci of S1b + neighbour genes of the islands table
        by[ctg] = []
        for o in CAT.values():
            if o['contig'] == ctg and o['locus_id'] != lid and o['lineage'] not in ('not a homolog', 'native fungal GDPD') and o['start']:
                by[ctg].append(dict(contig=ctg, assembly=o['assembly'], gene='T_' + o['locus_id'], start=o['start'], end=str(int(float(o['end']))), strand=o['strand'], gene_class='focal locus'))
        for gname, g in ISL.get(ctg, {}).items():
            by[ctg].append(dict(contig=ctg, assembly=c['assembly'], gene=gname, start=g['start'], end=g['end'], strand=g['strand'], gene_class='neighboring gene'))
    by[ctg].append(dict(contig=ctg, assembly=c['assembly'], gene='T_' + lid, start=c['start'], end=str(int(float(c['end']))), strand=c['strand'], gene_class='focal locus'))
PAD = 12000
islands = []; skipped = []
for ctg, rows in by.items():
    foc = [r for r in rows if r['gene_class'] == 'focal locus']
    # distinct coordinates only (a locus assigned to two families at identical coordinates is one gene)
    dist = {(int(r['start']), int(r['end'])) for r in foc}
    if len(dist) < 2: skipped.append((ctg, len(foc), len(dist))); continue
    lo = max(min(int(r['start']) for r in rows), min(s for s, e in dist) - PAD)
    hi = min(max(int(r['end']) for r in rows), max(e for s, e in dist) + PAD)
    genes = []
    for r in foc:
        lid = r['gene'][2:]; c = CAT[lid]
        cat = 'nos' if lid.startswith('ODH_') else LINCAT(c['lineage'])
        genes.append(dict(start=int(r['start']), end=int(r['end']), strand=r['strand'], cat=cat, intact=c['status'] == 'intact', lid=lid, lin=c['lineage']))
    for r in rows:
        if r['gene_class'] == 'focal locus': continue
        s, e = int(r['start']), int(r['end'])
        if e < lo or s > hi: continue
        if any(min(e, g['end']) - max(s, g['start']) > 0.3 * (e - s) for g in genes): continue
        t = ROLE.get((ctg, r['gene']), '').lower()
        if 'transposase' in t or 'integrase' in t: cat = 'tn'
        elif 'dehydrogenase' in t and ('nopaline' in t or 'octopine' in t): cat = 'nos'
        else: cat = 'fung'
        genes.append(dict(start=s, end=e, strand=r['strand'], cat=cat, intact=True, lid='', lin=''))
    genes.sort(key=lambda g: g['cat'] == 'nos')   # nos-ODH drawn last (on top of any overlapping gene)
    asm = rows[0]['assembly']
    islands.append(dict(contig=ctg, asm=asm, org=asm_org[asm], lo=lo, hi=hi, genes=genes,
                        lins=sorted({g['cat'] for g in genes if g['lid']})))
islands.sort(key=lambda d: (d['org'].replace(' NGJ', ''), d['asm'], d['contig']))
print('skipped', skipped)
for d in islands: print(d['org'], d['asm'], d['contig'], d['lo'], d['hi'], (d['hi'] - d['lo']) / 1000, d['lins'], len(d['genes']))

def arrow(ax, g, x0, h=0.2, head=0.35, y0=0.0):
    s, e = (g['start'] - x0) / 1000, (g['end'] - x0) / 1000
    L = e - s; hl = min(head, 0.55 * L)
    if g['strand'] == '-': pts = [(e, y0 - h), (s + hl, y0 - h), (s, y0), (s + hl, y0 + h), (e, y0 + h)]
    else: pts = [(s, y0 - h), (e - hl, y0 - h), (e, y0), (e - hl, y0 + h), (s, y0 + h)]
    c = COL.get(g['cat'], '#888888')
    if g['cat'] == 'fung': fc, ec, ls = c, '#8c8c8c', '-'
    elif g['cat'] == 'other': fc, ec, ls = '#bbbbbb', '#555555', (0, (2, 1.2))
    elif g['intact'] or g['cat'] == 'tn': fc, ec, ls = c, c, '-'
    else: fc, ec, ls = pale(c), c, (0, (2, 1.2))
    ax.add_patch(mp.Polygon(pts, closed=True, fc=fc, ec=ec, lw=0.7, ls=ls, zorder=3))

CM = 1 / 2.54
FW = 17.6; LM = 4.2; WMAX = max((d['hi'] - d['lo']) for d in islands) / 1000
KB = 11.0 / WMAX                 # cm per kb, common
PITCH = 1.15; TOP = 0.5
FH = TOP + PITCH * len(islands) + 1.6 + 2.0
fig = plt.figure(figsize=(FW * CM, FH * CM), dpi=300)
def axcm(x, ytop, w, h): return fig.add_axes([x / FW, 1 - (ytop + h) / FH, w / FW, h / FH])
y = TOP
for d in islands:
    wk = (d['hi'] - d['lo']) / 1000
    ax = axcm(LM, y, wk * KB, PITCH - 0.1); ax.set_xlim(0, wk); ax.set_ylim(-1, 1); ax.axis('off')
    ax.plot([0, wk], [0, 0], color='#bbbbbb', lw=0.5, zorder=1)
    for g in d['genes']: arrow(ax, g, d['lo'], h=0.4, head=0.35)
    gn, sp = (d['org'].replace(' NGJ', '').split(' ', 1) + [''])[:2]
    spl = f'$\\it{{{gn}}}$ ' + (f'$\\it{{{sp}}}$' if sp and not sp.startswith('sp.') else sp)
    fig.text((LM - 0.2) / FW, 1 - (y + (PITCH - 0.1) / 2) / FH,
             f"{spl}\n{d['asm']}\n{d['contig']}\n{d['lo']/1000:,.1f}–{d['hi']/1000:,.1f} kb",
             fontsize=FS - 0.5, ha='right', va='center', linespacing=1.05)
    fig.text((LM + wk * KB + 0.15) / FW, 1 - (y + (PITCH - 0.1) / 2) / FH, f'{wk:.1f} kb', fontsize=FS, ha='left', va='center', fontweight='bold')
    y += PITCH
axs = axcm(LM, y + 0.1, 5 * KB, 0.1); axs.axis('off'); axs.set_xlim(0, 5); axs.set_ylim(0, 1)
axs.plot([0, 5], [0.5, 0.5], color='k', lw=0.9); axs.text(2.5, 0.0, '5 kb', ha='center', va='top', fontsize=FS)
# legend
yl = y + 0.9
axl = axcm(0.3, yl, FW - 0.6, FH - yl - 0.1); axl.set_xlim(0, FW - 0.6); axl.set_ylim(0, FH - yl - 0.1); axl.axis('off')
H = FH - yl - 0.1
items = [('$\\it{sus}$-like', 'sus'), ('$\\it{cus}$-like', 'cus'), ('acc-type GDPD', 'acc'),         ('GDPD, position unresolved', 'unpl'), ('$\\it{mis}$ group (other bacterial type)', 'mis'),
         ('opine dehydrogenase (nos-like)', 'nos'), ('other bacterial-type ORF', 'other'), ('transposase / integrase', 'tn'), ('fungal neighbour gene', 'fung')]
def sw(x, yy, cat, dashed=False):
    c = COL[cat]
    if cat == 'fung': axl.add_patch(mp.Rectangle((x, yy - 0.13), 0.5, 0.26, fc=c, ec='#8c8c8c', lw=0.6))
    elif dashed or cat == 'other': axl.add_patch(mp.Rectangle((x, yy - 0.13), 0.5, 0.26, fc=pale(c), ec=c, lw=0.7, ls=(0, (2, 1.2))))
    else: axl.add_patch(mp.Rectangle((x, yy - 0.13), 0.5, 0.26, fc=c, ec=c, lw=0.7))
for j, (t, cat) in enumerate(items):
    x0 = 0.0 + 5.9 * (j // 4 if j < 8 else 2); yy = H - 0.3 - 0.4 * (j % 4 if j < 8 else 4)
    if j == 8: x0, yy = 0.0, H - 0.3 - 0.4 * 4
    sw(x0, yy, cat); axl.text(x0 + 0.65, yy, t, fontsize=FS, va='center')
x0 = 11.4; yy = H - 0.3
sw(x0, yy, 'cus'); sw(x0 + 0.62, yy, 'cus', dashed=True); axl.text(x0 + 1.4, yy, 'intact / pseudogene,\nfragment or remnant', fontsize=FS, va='center')
axl.text(x0, yy - 0.75, 'arrow = strand; repeats not annotated', fontsize=FS, va='center', color='#444444')
for e in ('png', 'svg', 'pdf'): fig.savefig(f'{OUT}/FigureS10_islands.{e}', dpi=300)
print('ok', FW, FH)
