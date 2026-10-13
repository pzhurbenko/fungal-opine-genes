#!/usr/bin/env python3
"""Supplementary Figure: AU topology tests of all lineages, drawn from Table S11
(supplementary/TableS11_Approximately_unbiased_AU_topology_tests.tsv). Split off Fig. 3 on 2026-10-12 (decision 6.3b).
Run: python make_figS_AU.py
Blocks = runs of rows starting with tree_number 0 (one block per AU run); rows with tree_number 0 are the ML reference
and are listed first (p = 1 by construction unless the reference tree was not the best)."""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import os, csv, re, textwrap, glob, math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
TSV = glob.glob(os.path.join(PAPER, 'supplementary', 'TableS11_*.tsv'))[0]
OUT = os.path.join(PAPER, 'supplementary', 'figures', 'FigureS8_AU_tests')

plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
                     'font.size': 6.0, 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
FS = 6.0; W = 17.5; CM = 1 / 2.54
RED, BLUE, GREEN, PURP = '#D7301F', '#0072B2', '#009E73', '#7B3294'
ORANGE, PINK, INDIGO = '#E69F00', '#CC79A7', '#332288'
REJ, OK = '#222222', '#E6AB02'
# colour of a block by its family and by the lineage named in the constraints
def block_colour(fam, hyps):
    f = fam.lower(); h = ' '.join(x[0] for x in hyps).lower()
    if f.startswith('ocd'): return RED
    if 'cucumopine' in f: return BLUE
    if f.startswith('gdpd'):
        return GREEN if h.count('acc-type') > h.count('pkerr14a-type') else PURP
    if 'chrysopine' in f: return ORANGE
    if 'nos subfamily' in f: return PINK
    if 'odh island' in f: return INDIGO
    return '#444444'

def block_title(fam, hyps):
    f = fam; h = ' '.join(x[0] for x in hyps).lower()
    if f.startswith('OCD'): return 'sus-like (OCD-fold family)'
    if f.startswith('cucumopine'):
        if 'outside boletales' in h: return 'cus-like outside Boletales (cucumopine synthase family; separate AU run)'
        return 'cus-like (cucumopine synthase family)'
    if f.startswith('GDPD'): return 'GDPD family (acc-type GDPD and pKerr14a-type acs-like)'
    if f.startswith('chrysopine'): return 'chs-like (chrysopine synthase; M24 aminopeptidase family)'
    if 'nos subfamily' in f: return 'nos-subfamily opine dehydrogenase (Boletales; Rhypophila homolog, not assigned)'
    if 'odh island' in f: return 'odh island (Fusarium; opine dehydrogenase family)'
    return f

def read_blocks():
    rows = list(csv.DictReader(open(TSV), delimiter='\t'))
    blocks = []
    for r in rows:
        n = int(r['tree_number'])
        if n == 0:
            blocks.append([r['family'], []])
        blocks[-1][1].append((r['tested_topology_constraint'], float(r['deltaL']), float(r['p_AU']), n, r['verdict']))
    return blocks

WRAP = 78
def main():
    blocks = read_blocks()
    RH = 0.27  # row height (cm) per text line
    rows = []  # (y, kind, ...)
    y = 0.45
    for fam, hyps in blocks:
        y += 0.42; rows.append((y, 'g', block_title(fam, hyps), block_colour(fam, hyps)))
        for lab, dl, p, n, verdict in hyps:
            lines = textwrap.wrap(f'{n}. {lab}', WRAP)
            y += RH * len(lines) + 0.06
            rows.append((y, 'h', lines, dl, p, n, verdict))
        y += 0.12
    H = y + 1.5
    fig = plt.figure(figsize=(W * CM, H * CM))
    ax0 = fig.add_axes([0, 0, 1, 1]); ax0.set_xlim(0, W); ax0.set_ylim(0, H); ax0.axis('off')
    # p-value axis (log scale) in cm coordinates: x from 10.6 to 13.3
    axl, axw = 10.6, 2.7
    ax = fig.add_axes([axl / W, 1.0 / H, axw / W, (H - 1.4) / H])
    ax.set_xscale('log'); ax.set_xlim(2.0, 1e-83); ax.set_ylim(H - 1.0, 0.4)
    for s in ('top', 'right', 'left'): ax.spines[s].set_visible(False)
    ax.patch.set_visible(False); ax.set_yticks([]); ax.tick_params(labelsize=FS, length=2.5, width=0.6); ax.spines['bottom'].set_linewidth(0.6)
    ax.set_xticks([1, 1e-20, 1e-40, 1e-60, 1e-80]); ax.set_xticks([], minor=True)
    ax.set_xlabel('p-AU (log scale)', fontsize=FS, labelpad=1.5)
    ax.axvline(0.05, color='#333333', lw=0.8, ls='--', zorder=1)
    ax0.text(axl + axw * (-math.log10(0.05) + math.log10(2.0)) / (83 + math.log10(2.0)), H - 0.08, 'p = 0.05', fontsize=FS, ha='center', va='top')
    def Y(yc): return H - yc  # cm from top -> axes-cm from bottom
    ax0.text(0.05, H - 0.1, 'AU tests of the monophyly constraints (Table S11)', fontsize=FS + 1.5, fontweight='bold', va='top')
    ax0.text(13.75, H - 0.1, 'ΔlnL', fontsize=FS, fontweight='bold', va='top')
    ax0.text(15.5, H - 0.1, 'p-AU', fontsize=FS, fontweight='bold', va='top')
    for r in rows:
        if r[1] == 'g':
            yc, _, title, col = r
            ax0.text(0.3, Y(yc), title, fontsize=FS + 0.5, fontweight='bold', color=col, va='center')
            continue
        yc, _, lines, dl, p, n, verdict = r
        ymid = yc - RH * (len(lines) - 1) / 2
        for i, ln in enumerate(lines):
            ax0.text(0.6, Y(yc - RH * (len(lines) - 1 - i)), ln, fontsize=FS, va='center')
        ref = (n == 0)
        rej = p < 0.05
        c = '#888888' if ref else (REJ if rej else OK)
        ax.plot([1, max(p, 1e-83)], [ymid, ymid], color='#BBBBBB', lw=0.6, zorder=2)
        ax.plot([max(p, 1e-83)], [ymid], 'o' if not ref else 's', ms=4.0, color=c, mec='k' if (not rej or ref) else c, mew=0.4, zorder=3, clip_on=False)
        ps = f'{p:.3f}' if p >= 0.01 else f'{p:.1e}'
        ps = re.sub(r'e-0?(\d+)', lambda m: '×10$^{-' + m.group(1) + '}$', ps)
        ax0.text(13.75, Y(ymid), f'{dl:g}', fontsize=FS, va='center')
        ax0.text(15.5, Y(ymid), ps, fontsize=FS, va='center')
    ax0.legend(handles=[Line2D([], [], marker='s', ls='', color='#888888', mec='k', mew=0.4, ms=4.0, label='unconstrained ML tree (reference)'),
                        Line2D([], [], marker='o', ls='', color=OK, mec='k', mew=0.4, ms=4.0, label='constraint not rejected (p ≥ 0.05)'),
                        Line2D([], [], marker='o', ls='', color=REJ, ms=4.0, label='constraint rejected (p < 0.05)')],
               loc='lower left', bbox_to_anchor=(0.01, 0.002), ncol=1, columnspacing=1.0, fontsize=FS, frameon=False, handletextpad=0.3, labelspacing=0.2)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    for ext, kw in (('png', dict(dpi=300)), ('svg', {}), ('pdf', {})):
        fig.savefig(f'{OUT}.{ext}', **kw)
    print('blocks', len(blocks), 'rows', sum(len(b[1]) for b in blocks), 'size cm', W, round(H, 2))

if __name__ == '__main__':
    main()
