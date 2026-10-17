#!/usr/bin/env python3
"""Graphical abstract v2, minimal (2026-10-15): cropped loci, no scale bars/legend, numbers, patchy distribution cladogram: plasmid ring (T-DNA vs other opine loci) -> fungal genomes (two real loci coloured by source) -> fate.
Detailed variant kept as make_graphical_abstract_detailed.py.
Run: python make_graphical_abstract_v2.py   (conda env paper_build)
Data: paper/tables/TableS1_loci_catalog.tsv, analysis/islands/islands_genes.tsv, analysis/islands/island_JAHRGV010000066_genes.tsv,
analysis/island_order/island_gene_order.tsv, paper/figures/{Fig2_genus_counts,Fig5_subfamilies_Tremble2024}.tsv, analysis/fusarium_odh_island/{scan/loci.tsv, fig4/fp133_window_genes.tsv},
Hard-coded: gene schematics of the plasmid track, 'fungal sus-like (56 loci ...)' label, table values, footer."""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import csv, collections, textwrap, sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mp
from matplotlib import font_manager as fm
from matplotlib.font_manager import FontProperties
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

_FD = '/usr/share/fonts/urw-base35'
_FILES = [f'{_FD}/NimbusSans-{f}.otf' for f in ('Regular', 'Bold', 'Italic', 'BoldItalic')]
if all(os.path.exists(p) for p in _FILES):
    for p in _FILES: fm.fontManager.addfont(p)
    FONT = 'Nimbus Sans'
else:
    FONT = 'DejaVu Sans'
    print(f'NOTE: {_FD}/NimbusSans-*.otf not found; falling back to DejaVu Sans (text widths differ slightly)')
plt.rcParams.update({'font.family': FONT, 'svg.fonttype': 'none', 'hatch.linewidth': 0.8})

G = f'{GMO_ROOT}'; A = f'{G}/analysis'
OUT = f'{G}/paper/figures/graphical_abstract'; os.makedirs(OUT, exist_ok=True)
COL = {'sus': '#d62728', 'cus': '#1f77b4', 'acc': '#2ca02c', 'kerr': '#9467bd', 'mis': '#17becf', 'unpl': '#9a9a3a',
       'nos': '#ff7f0e', 'tn': '#000000', 'pepp': '#8c564b', 'g3': '#e377c2', 'fung': '#d9d9d9', 'none': '#ffffff',
       'chs': '#E69F00', 'acs': '#6E6E6E'}
LIN = {'other fungal GDPD (position unresolved)': 'unpl', 'acc-type GDPD': 'acc', 'cus-like (T-DNA type)': 'cus',
       'mis lineage': 'mis', 'mis lineage (not placed)': 'mis', 'sus-like (T-DNA type)': 'sus',
       'sus-like (T-DNA type) (not placed)': 'sus', 'pKerr14a-type acs-like': 'kerr'}
GREY = '#6e6e6e'; TXTC = '#222222'; LRED = '#F4B6AE'; MG = '#A0A0A0'; LG = '#D9D9D9'

W, H = 2200, 1100
fig = plt.figure(figsize=(22, 11), dpi=100, facecolor='white')
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.axis('off')
fig.canvas.draw(); R = fig.canvas.get_renderer()
PX = 100 / 72.0

def tw(s, fs, weight='normal', style='normal'):
    fp = FontProperties(family=FONT, size=fs, weight=weight, style=style)
    w, h, d = R.get_text_width_height_descent(s, fp, False)
    return w

def T(x, y, s, fs=15, color='k', weight='normal', style='normal', ha='left', va='center'):
    ax.text(x, y, s, fontsize=fs, color=color, fontweight=weight, fontstyle=style, ha=ha, va=va)
    return x + tw(s, fs, weight, style)

def runs(x, y, parts, fs=15, color='k', weight='normal', ha='left'):
    """parts: list of (text, 'n'|'i'|'b'), optional per-run colour as 3rd item; returns end x"""
    ws = [tw(p[0], fs, weight, 'italic' if p[1] == 'i' else 'normal') for p in parts]
    if ha == 'center': x = x - sum(ws) / 2
    cx = x
    for p, w in zip(parts, ws):
        T(cx, y, p[0], fs, p[2] if len(p) > 2 else color, weight, 'italic' if p[1] == 'i' else 'normal')
        cx += w
    return cx

def pale(c, f=0.55):
    c = np.array(matplotlib.colors.to_rgb(c)); return tuple(c + (1 - c) * f)

def garrow(x0, x1, y, h, strand, cat, intact=True, hatch=False, head=16):
    L = abs(x1 - x0); hl = min(head, 0.55 * L); s, e = min(x0, x1), max(x0, x1)
    if strand == '-': pts = [(e, y - h), (s + hl, y - h), (s, y), (s + hl, y + h), (e, y + h)]
    else: pts = [(s, y - h), (e - hl, y - h), (e, y), (e - hl, y + h), (s, y + h)]
    c = COL[cat]
    if hatch: kw = dict(fc='white', ec='#555555', hatch='//////', lw=1.0)
    elif cat == 'none': kw = dict(fc='white', ec='#8c8c8c', lw=1.0)
    elif cat == 'fung': kw = dict(fc=c, ec='#8c8c8c', lw=1.0)
    elif intact: kw = dict(fc=c, ec=c, lw=1.0)
    else: kw = dict(fc=pale(c), ec=c, lw=1.2, ls=(0, (2, 1.2)))
    ax.add_patch(mp.Polygon(pts, closed=True, zorder=3, **kw))

# ---------------- data (copied from make_fig4.py) ----------------
def tsv(p): return list(csv.DictReader(open(p), delimiter='\t'))
lin = {r['locus_id']: r['lineage_final'] for r in tsv(f'{G}/paper/tables/TableS1_loci_catalog.tsv')}
role = {}
for r in tsv(f'{A}/islands/islands_genes.tsv'): role[(r['contig'], r['gene'])] = (r['role'], r['top_title'])
for r in tsv(f'{A}/islands/island_JAHRGV010000066_genes.tsv'): role.setdefault(('JAHRGV010000066', r['gene']), (r['role'], r['top_title']))
gene_rows = {}
for r in tsv(f'{A}/island_order/island_gene_order.tsv'): gene_rows.setdefault(r['island_contig'], []).append(r)
for r in tsv(f'{A}/islands/island_JAHRGV010000066_genes.tsv'):
    if r['gene'] in ('g27',):
        gene_rows['JAHRGV010000066'].append({'island_contig': 'JAHRGV010000066', 'gene': r['gene'], 'start': r['start'], 'end': r['end'],
            'strand': r['strand'], 'exons': r['exons'], 'kind': 'neighbor', 'product': '', 'plasmid_gene': ''})

def build(contig, w0, w1):
    rows = [r for r in gene_rows[contig] if int(r['end']) >= w0 and int(r['start']) <= w1]
    out = []
    for r in [r for r in rows if r['kind'].startswith('registry')]:
        lid = r['gene'][2:]; st = r['kind'].split(':')[-1]
        cat = 'none' if (r['kind'].startswith('registry:') and lid not in lin) else LIN[lin[lid]]
        out.append(dict(start=int(r['start']), end=int(r['end']), strand=r['strand'] or '+', cat=cat, intact=(st == 'intact'), hatch=False, key=lid))
    for r in rows:
        if r['kind'].startswith('registry'): continue
        s, e = int(r['start']), int(r['end'])
        ov = max((min(e, x['end']) - max(s, x['start']) for x in out if x['key'] and x['key'].count('_') >= 2), default=0)
        if ov > 0.3 * (e - s) and ov > 0: continue
        rl, title = role.get((contig, r['gene']), (None, ''))
        t = (title or '').lower() + ' ' + (r.get('product') or '').lower(); hatch = False
        if 'dehydrogenase' in t and ('nopaline' in t or 'octopine' in t): cat = 'nos'
        elif 'transposase' in t or 'integrase' in t: cat = 'tn'
        elif contig == 'JAHRGV010000066' and r['gene'] in ('g20', 'g26', 'g27'): cat, hatch = 'fung', True
        elif rl in ('Fungi', 'anchor', 'bacterial_like'): cat = 'fung'
        else: cat = 'none'
        out.append(dict(start=s, end=e, strand=r['strand'] or '+', cat=cat, intact=True, hatch=hatch, key=r['gene']))
    return sorted(out, key=lambda x: x['start'])


# ---------------- source colours ----------------
RED, TEAL, LTEAL = '#D55E00', '#1B7F8C', '#5FA8B3'
SRC = {'sus': RED, 'cus': RED, 'mis': RED, 'chs': RED, 'acc': TEAL, 'nos': TEAL, 'kerr': TEAL, 'pepp': LTEAL, 'g3': LTEAL,
       'fung': '#d9d9d9', 'none': '#ffffff', 'tn': '#000000', 'unpl': '#ffffff'}
def sarrow(x0, x1, y, h, strand, cat, intact=True, head=18):
    L = abs(x1 - x0); hl = min(head, 0.55 * L); s, e = min(x0, x1), max(x0, x1)
    if strand == '-': pts = [(e, y - h), (s + hl, y - h), (s, y), (s + hl, y + h), (e, y + h)]
    else: pts = [(s, y - h), (e - hl, y - h), (e, y), (e - hl, y + h), (s, y + h)]
    c = SRC[cat]
    if cat in ('none', 'unpl'): kw = dict(fc='white', ec='#9a9a9a', lw=1.0)
    elif cat == 'fung': kw = dict(fc=c, ec='#9a9a9a', lw=1.0)
    elif intact: kw = dict(fc=c, ec=c, lw=1.0)
    else: kw = dict(fc='white', ec=c, lw=1.6, ls=(0, (2, 1.2)))
    ax.add_patch(mp.Polygon(pts, closed=True, zorder=3, **kw))


import csv as _csv, collections


# ---------------- three text styles ----------------
C1, C3 = '#1a1a1a', '#5a5a5a'
def S1(x, y, s, **k): return T(x, y, s, 24, C1, 'bold', **k)                 # block titles
def S2(x, y, s, **k): return T(x, y, s, 18, C1, 'normal', **k)               # text
def S3(x, y, s, **k): return T(x, y, s, 17, C3, 'normal', 'italic', **k)     # captions and labels
def R1(x, y, parts, **k): return runs(x, y, parts, fs=24, color=C1, weight='bold', **k)
def R2(x, y, parts, **k): return runs(x, y, parts, fs=18, color=C1, **k)

# ================= LEFT: plasmid ring =================
cx, cy, r = 300, 480, 190
R1(60, 62, [('Agrobacterium', 'i'), (' Ti/Ri plasmid', 'n')])
ax.add_patch(mp.Circle((cx, cy), r, fc='none', ec='#C4C4C4', lw=14, zorder=1))
ax.add_patch(mp.Wedge((cx, cy), r + 12, 215, 325, width=24, fc=RED, ec='none', zorder=2))
for ang in (215, 325):
    a = np.radians(ang); ax.plot([cx + (r - 24) * np.cos(a), cx + (r + 26) * np.cos(a)], [cy + (r - 24) * np.sin(a), cy + (r + 26) * np.sin(a)], color=C1, lw=3, zorder=3)
ax.add_patch(mp.Wedge((cx, cy), r + 12, 95, 150, width=24, fc=TEAL, ec='none', zorder=2))
S2(cx, cy - r - 48, 'T-DNA', ha='center')
runs(cx, cy - r + 62, [('sus', 'i'), (',  ', 'n'), ('cus', 'i'), (',  ', 'n'), ('chs', 'i')], fs=17, color=C3, ha='center')
S3(cx, cy - r + 98, 'opine synthesis genes', ha='center')
S2(cx - 60, cy + r + 50, 'other opine loci', ha='center')
runs(cx - 60, cy + r + 84, [('acc', 'i'), (',  ', 'n'), ('nos', 'i'), ('-like,  ', 'n'), ('odh', 'i')], fs=17, color=C3, ha='center')
S3(cx - 60, cy + r + 116, 'opine uptake and catabolism', ha='center')

# ================= arrow =================
AX0, AX1 = cx + r + 42, 745
ax.add_patch(mp.FancyArrowPatch((AX0, cy), (AX1, cy), arrowstyle='-|>', mutation_scale=38, color=C1, lw=4.5, shrinkA=0, shrinkB=0))
S1((AX0 + AX1) / 2, cy - 72, 'Horizontal', ha='center'); S1((AX0 + AX1) / 2, cy - 38, 'gene transfer', ha='center')

# ================= RIGHT =================
XL = 790; XT1 = XL + 880; XR = 2140
S1(XL, 62, 'Fungal genomes')
def numline(x, y, parts):
    for s in parts: x = S2(x, y, s)
numline(XL, 140, ['Screened:  6,123 annotated proteomes   ·   2,215 genome assemblies   ·   41 orders'])
numline(XL, 178, ['Found:  7 lineages of bacterial origin   ·   254 loci   ·   121 assemblies   ·   78 species   ·   8 orders'])

RULES = (222, 560)
x_ = S2(XL, 262, 'Examples of acquired loci'); runs(x_ + 26, 262, [('Butyriboletus roseoflavus', 'i'), (' (top),  ', 'n'), ('Fusarium poae', 'i'), (' (bottom)', 'n')], fs=17, color=C3)
HA = 16
def solidish(x0, x1, y, h, strand, cat, intact): sarrow(x0, x1, y, h, strand, 'fung' if cat in ('none', 'unpl') else cat, intact, head=14)
ctg = 'JAHRGV010000067'; w0, w1 = 1428500, 1474500
YA = 340; PXA = (XT1 - XL) / ((w1 - w0) / 1000)
ax.plot([XL, XT1], [YA, YA], color='#bbbbbb', lw=1.2, zorder=1)
GA = build(ctg, 1426000, 1480000)
for g in GA:
    a = XL + (g['start'] - w0) / 1000 * PXA; b = XL + (g['end'] - w0) / 1000 * PXA
    if b < XL or a > XT1: continue
    solidish(max(a, XL), min(b, XT1), YA, HA, g['strand'], g['cat'], g['intact'])
GC = lambda k: next((XL + ((g['start'] + g['end']) / 2 - w0) / 1000 * PXA) for g in GA if g['key'] == k)
S3(GC('ACS_JAHRGV010000067_L1') - 8, YA - 30, 'acc-like', ha='center')
S3(GC('CUS_JAHRGV010000067_L1') + 12, YA + 34, 'mis-like', ha='center')
S3(GC('g30') + 22, YA - 30, 'nos-like', ha='center')
S3(GC('SUS_JAHRGV010000067_L1'), YA + 34, 'sus-like', ha='center')
FA = f'{A}/fusarium_odh_island'
LOC = {r_['locus_id']: r_ for r_ in tsv(f'{FA}/scan/loci.tsv') if r_['gca'] == 'GCA_046056315.1'}
fpw = tsv(f'{FA}/fig4/fp133_window_genes.tsv')
CASS = {'FPOAC2_07382', 'FPOAC2_07383', 'FPOAC2_10362', 'FPOAC2_10363', 'FPOAC2_10364'}
cs = {k: LOC[f'GCA_046056315.1_{v}'] for k, v in {'odh': 'odh1', 'pepP': 'pepP1', 'g3': 'g31'}.items()}
lo = int(cs['g3']['start']) - 6500; hi = int(cs['odh']['end']) + 6500
YB = 465; PXB = (XT1 - XL) / ((hi - lo) / 1000)
ax.plot([XL, XT1], [YB, YB], color='#bbbbbb', lw=1.2, zorder=1)
for r_ in fpw:
    if r_['contig'] != 'CM100255.1' or r_['locus_tag'] in CASS: continue
    s, e = int(r_['start']), int(r_['end'])
    if e < lo or s > hi: continue
    solidish(max(XL + (s - lo) / 1000 * PXB, XL), min(XL + (e - lo) / 1000 * PXB, XT1), YB, HA, r_['strand'], 'fung', True)
for k, cat in (('odh', 'nos'), ('pepP', 'pepp'), ('g3', 'g3')):
    r_ = cs[k]; a = XL + (int(r_['start']) - lo) / 1000 * PXB; b = XL + (int(r_['end']) - lo) / 1000 * PXB
    solidish(a, b, YB, HA, r_['strand'], cat, r_['status'] == 'intact')
    S3((a + b) / 2 + {'g3': -30, 'pepP': 0, 'odh': 28}[k], YB + 34, {'g3': 'orf3', 'pepP': 'pepP-like', 'odh': 'odh-like'}[k], ha='center')
# key
def ksw(x, y, fc, ec, dashed=False): ax.add_patch(mp.Rectangle((x, y - 9), 28, 18, fc=fc, ec=ec, lw=1.6, ls=(0, (2, 1.2)) if dashed else '-'))
kx = XT1 + 50
ksw(kx, YA, RED, RED); S2(kx + 40, YA, 'intact gene')
ksw(kx, YA + 36, 'white', RED, True); S2(kx + 40, YA + 36, 'disrupted or partial gene')
ksw(kx, YA + 72, '#d9d9d9', '#9a9a9a'); S2(kx + 40, YA + 72, 'fungal gene')

# ================= map =================
G5 = {r_['row']: r_ for r_ in _csv.DictReader(open(f'{G}/paper/figures/Fig2_genus_counts.tsv'), delimiter='\t')}
SF = {r_['genus']: r_['subfamily'] for r_ in _csv.DictReader(open(f'{G}/paper/figures/Fig5_subfamilies_Tremble2024.tsv'), delimiter='\t')}
GEN = [('Boletoideae', ['Boletus', 'Tylopilus']), ('Suillelloideae', ['Pulveroboletus', 'Butyriboletus']),
       ('Leccinoideae', ['Leccinum']), ('Xerocomoideae', ['Boletellus', 'Phylloporus'])]
for sf, gs in GEN:
    for g_ in gs: assert SF.get(g_) == sf, (g_, SF.get(g_))
LINS = [('sus-like_k_any', 'sus-like', RED), ('cus-like_k_any', 'cus-like', RED), ('acc-type_GDPD_k_any', 'acc-like GDPD', TEAL)]
MY0 = 660; DY = 34; XSF, XGN, XD0, DXC, XN = XL, XL + 215, XL + 480, 185, XL + 990
x_ = S2(XL, 602, 'Patchy distribution in Boletaceae'); S3(x_ + 14, 602, '(circle: species carrying the gene)')
S3(XSF, MY0, 'subfamily'); S3(XN, MY0, 'species screened')
for j, (_, lab, col) in enumerate(LINS): S3(XD0 + j * DXC, MY0, lab, ha='center')
y = MY0 + 36
for sf, gs in GEN:
    y0 = y
    for g_ in gs:
        row_ = G5[g_]; n = int(row_['n_clusters'])
        S3(XGN, y, g_, va='center')
        for j, (key, lab, col) in enumerate(LINS):
            k_ = int(row_[key]); xc = XD0 + j * DXC
            if k_:
                ax.add_patch(mp.Circle((xc, y), 11, fc=col, ec=col, zorder=3)); S2(xc + 20, y, str(k_), va='center')
            else: ax.add_patch(mp.Circle((xc, y), 3.5, fc='white', ec='#b0b0b0', lw=1.2, zorder=3))
        S2(XN + 60, y, f'{n}', va='center', ha='center')
        y += DY
    ax.plot([XSF + 195, XSF + 195], [y0 - 10, y - DY + 10], color='#c8c8c8', lw=1.5)
    S2(XSF, (y0 + y - DY) / 2, sf, va='center')
    y += 10

# ---- equal margins: measure the content box, shift the view so margins are equal, rules end at the content edge
fig.canvas.draw(); REN2 = fig.canvas.get_renderer()
xs0, xs1, ys0, ys1 = [], [], [], []
for art in list(ax.texts) + list(ax.patches) + list(ax.lines):
    bb = art.get_window_extent(REN2)
    if bb.width == 0 and bb.height == 0: continue
    xs0.append(bb.x0); xs1.append(bb.x1); ys0.append(bb.y0); ys1.append(bb.y1)
bx0, bx1 = min(xs0), max(xs1); by_top, by_bot = H - max(ys1), H - min(ys0)   # data coords (y inverted)
for yr in RULES: ax.plot([XL, bx1], [yr, yr], color='#b8b8b8', lw=1.6, zorder=0)
mx = (W - (bx1 - bx0)) / 2; my = (H - (by_bot - by_top)) / 2
ax.set_xlim(bx0 - mx, bx1 + mx); ax.set_ylim(by_bot + my, by_top - my)
print(f'content box x {bx0:.0f}-{bx1:.0f}, y {by_top:.0f}-{by_bot:.0f}; margins x {mx:.0f}, y {my:.0f}')
fig.savefig(f'{OUT}/GA_v2_minimal.svg', dpi=100, facecolor='white')
fig.savefig(f'{OUT}/GA_v2_minimal_full.png', dpi=100, facecolor='white')
# crop PNG/TIFF to the content plus a uniform 60-px margin (SVG keeps the full canvas)
from PIL import Image
M = 60
im = Image.open(f'{OUT}/GA_v2_minimal_full.png')
box = (int(round(mx - M)), int(round(my - M)), int(round(W - mx + M)), int(round(H - my + M)))
im2 = im.crop(box).convert('RGB')   # no alpha channel, white background
im2.save(f'{OUT}/GA_v2_minimal.png', dpi=(300, 300)); im2.save(f'{OUT}/GA_v2_minimal.tiff', compression='tiff_lzw', dpi=(300, 300))
os.remove(f'{OUT}/GA_v2_minimal_full.png')
print('written v2, cropped to', im2.size)
