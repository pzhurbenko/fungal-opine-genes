#!/usr/bin/env python3
import re, glob, statistics
from scipy.stats import chi2
def lnls(f):
    txt = open(f).read()
    models = re.findall(r'Model (\d+):[^\n]*\n(?:.*\n)*?lnL\(ntime:\s*\d+\s+np:\s*(\d+)\):\s*([-\d.]+)', txt)
    return {int(m): (int(np), float(l)) for m, np, l in models}
out = open('dnds_summary.tsv', 'w'); out.write('family\tn_seqs\tcodon_sites\tM0_omega\tM7_lnL\tM8_lnL\tLRT_M7vsM8\tp_M7vsM8\tM8a_lnL\tLRT_M8avsM8\tp_M8avsM8\tM8_positive_sites_BEB>0.95\tyn00_median_dNdS\tyn00_q90_dNdS\tyn00_pairs\n')
for fam in ['SUS', 'CUS']:
    try: m0 = open(f'{fam}.M0.out').read()
    except FileNotFoundError: continue
    w0 = re.search(r'omega \(dN/dS\) =\s*([\d.]+)', m0); n, sites = re.search(r'^\s*(\d+)\s+(\d+)', open(f'{fam}.codon.phy').readline()).groups()
    row = {'family': fam, 'n_seqs': n, 'codon_sites': int(sites)//3, 'M0_omega': w0.group(1) if w0 else ''}
    try:
        L = lnls(f'{fam}.M7M8.out'); l7, l8 = L[7][1], L[8][1]; lrt = 2*(l8-l7); row.update(M7_lnL=l7, M8_lnL=l8, LRT_M7vsM8=f'{lrt:.2f}', p_M7vsM8=f'{chi2.sf(lrt,2):.3g}')
        beb = re.search(r'Bayes Empirical Bayes.*?\n(.*?)\n\n', open(f'{fam}.M7M8.out').read(), re.S)
        row['M8_positive_sites_BEB>0.95'] = (len(re.findall(r'\*', beb.group(1))) if beb else '') if chi2.sf(lrt, 2) < 0.05 else 0
    except Exception as e: row.update(M7_lnL='', M8_lnL='', LRT_M7vsM8='', p_M7vsM8='', **{'M8_positive_sites_BEB>0.95': ''})
    try:
        l8a = float(re.search(r'lnL\(ntime:\s*\d+\s+np:\s*\d+\):\s*([-\d.]+)', open(f'{fam}.M8a.out').read()).group(1)); lrt2 = 2*(row['M8_lnL']-l8a); row.update(M8a_lnL=l8a, LRT_M8avsM8=f'{lrt2:.2f}', p_M8avsM8=f'{chi2.sf(lrt2,1):.3g}')
    except Exception: row.update(M8a_lnL='', LRT_M8avsM8='', p_M8avsM8='')
    vals = []
    try:
        for l in open(f'{fam}.yn00.out'):
            m = re.match(r'^\s*\d+\s+\d+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+([\d.]+)\s+[-\d.]+\s+\+-\s+[\d.]+\s+([-\d.]+)', l)
            if m and float(m.group(2)) > 0.01 and float(m.group(1)) < 5: vals.append(float(m.group(1)))
    except FileNotFoundError: pass
    row.update(yn00_median_dNdS=f'{statistics.median(vals):.3f}' if vals else '', yn00_q90_dNdS=f'{sorted(vals)[int(0.9*len(vals))]:.3f}' if vals else '', yn00_pairs=len(vals))
    out.write('\t'.join(str(row.get(k, '')) for k in ['family','n_seqs','codon_sites','M0_omega','M7_lnL','M8_lnL','LRT_M7vsM8','p_M7vsM8','M8a_lnL','LRT_M8avsM8','p_M8avsM8','M8_positive_sites_BEB>0.95','yn00_median_dNdS','yn00_q90_dNdS','yn00_pairs']) + '\n')
    print(row)
out.close()
