"""Supplementary Figure S12: strain-specific insertion on the Boletus edulis contig JBBPRY010001778
(read depth, 5.9-kb segment absent from the other assemblies, gene model, repeat|insertion junction).
Formerly panel (b) of Figure 4. Run: MPLBACKEND=Agg python make_figS12_bedulis.py
Data: analysis/report_A/data/bedulis_contig_depth.tsv, analysis/loci_pipeline/round2/loci/ACS_JBBPRY010001778_L1/augustus.gff3"""
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
OUT = f'{G}/paper/figures'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 6.5, 'axes.linewidth': 0.5,
                     'mathtext.fontset': 'dejavusans', 'hatch.linewidth': 0.5,
                     'svg.fonttype': 'none', 'pdf.fonttype': 42})
FS = 6.5
OUT = f'{G}/paper/supplementary/figures'
COL = {'kerr': '#9467bd'}
CM = 1 / 2.54
FW, FH = 17.0, 7.6
fig = plt.figure(figsize=(FW * CM, FH * CM), dpi=300)
LM = 3.0
def axcm(x, ytop, w, h, **kw):
    return fig.add_axes([x / FW, 1 - (ytop + h) / FH, w / FW, h / FH], **kw)
y_b = 0.6
dep = [int(l.split('\t')[2]) for l in open(f'{A}/report_A/data/bedulis_contig_depth.tsv')]
n = len(dep); xk = np.arange(1, n + 1) / 1000
w = 200; sm = np.convolve(dep, np.ones(w) / w, mode='same')
wb = n / 1000
KBB = 0.85                      # cm per kb in panel b (own, larger scale)
DH, MH = 2.6, 2.4
axd = axcm(LM, y_b, wb * KBB, DH)
axd.fill_between(xk, sm, color='#bcd9ea', lw=0); axd.plot(xk, sm, color='#3a6ea5', lw=0.6)
axd.axhline(91.3, ls='--', color='#555555', lw=0.6)
axd.set_xlim(0, wb); axd.set_ylim(0, 140); axd.set_xticks([])
axd.set_yticks([0, 50, 100]); axd.tick_params(labelsize=FS, width=0.5, length=2, pad=1)
axd.set_ylabel('Read depth, ×', fontsize=FS, labelpad=1)
for s_ in ('top', 'right'): axd.spines[s_].set_visible(False)
axd.text(wb + 0.4, 91.3, 'genome median 91×', fontsize=FS, color='#555555', va='center', clip_on=False)
axd.axvline(1.313, color='k', ls=':', lw=0.8)
axd.annotate('', xy=(1.313, 130), xytext=(2.4, 130), arrowprops=dict(arrowstyle='-|>', lw=0.6, color='k', mutation_scale=5), annotation_clip=False)
axd.text(wb + 0.4, 130, 'repeat|insertion junction at 1,313 bp is spanned\nby 19 reads (9 read pairs): assembly is not chimeric',
         fontsize=FS, va='center', clip_on=False)
axd.set_title('')
y_m = y_b + DH + 0.45
axm = axcm(LM, y_m, wb * KBB, MH)
axm.set_xlim(0, wb); axm.set_ylim(0, 3.0); axm.axis('off')
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
axm.text(wb + 0.4, 1.2, 'segment 1,313–7,202 bp (5.9 kb) absent from\nfive other assemblies of the species', fontsize=FS, va='center', clip_on=False, color='#a8472a')
axm.text(wb + 0.4, 0.2, 'flanks (green) present in other assemblies', fontsize=FS, va='center', clip_on=False, color='#3d7a30')
axm.set_clip_on(False)
axm.axvline(1.313, color='k', ls=':', lw=0.8, ymax=1)
# x axis for panel b
axx = axcm(LM, y_m + MH + 0.0, wb * KBB, 0.5)
axx.set_xlim(0, wb); axx.set_ylim(0, 1); axx.set_yticks([])
for s_ in ('top', 'right', 'left'): axx.spines[s_].set_visible(False)
axx.set_xticks(range(0, 8, 2)); axx.set_xticks(range(0, 8), minor=True); axx.tick_params(labelsize=FS, width=0.5, length=2.5, pad=1)
axx.set_xlabel('Position on contig, kb', fontsize=FS, labelpad=1)
fig.text((LM - 0.2) / FW, 1 - (y_m + MH * 0.55) / FH, '$\\it{Boletus\\ edulis}$\nJBBPRY010001778\n7.4 kb', fontsize=FS, ha='right', va='center', linespacing=1.05)

fig.savefig(f'{OUT}/FigureS12_bedulis_insertion.png', dpi=300)
fig.savefig(f'{OUT}/FigureS12_bedulis_insertion.svg')
fig.savefig(f'{OUT}/FigureS12_bedulis_insertion.pdf')
print('ok')
