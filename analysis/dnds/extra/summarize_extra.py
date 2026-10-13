import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
#!/usr/bin/env python3
# summary.tsv for the extra dN/dS sets; yn00 filter as in ../parse_paml.py (pairs with dS>0.01 and omega<5)
import re, statistics, os
from scipy.stats import chi2
X = f'{GMO_ROOT}/analysis/dnds/extra'
SETS = ['chs', 'odh', 'pepP', 'g3', 'nosODH_Bol', 'cus_nonBol', 'cus_nonBol_sp', 'acs_pKerr14a']
def sites(p):
    try: return int(open(p).readline().split()[1]) // 3
    except Exception: return ''
def m0(p):
    try: m = re.search(r'omega \(dN/dS\) =\s*([\d.]+)', open(p).read()); return m.group(1) if m else ''
    except FileNotFoundError: return ''
def lnl(p, model=None):
    txt = open(p).read()
    if model is None: return float(re.search(r'lnL\(ntime:\s*\d+\s+np:\s*\d+\):\s*([-\d.]+)', txt).group(1))
    d = {int(m): float(l) for m, l in re.findall(r'Model (\d+):[^\n]*\n(?:.*\n)*?lnL\(ntime:\s*\d+\s+np:\s*\d+\):\s*([-\d.]+)', txt)}
    return d[model]
def yn(p):
    allp, keep = [], []
    if p.endswith('.tsv'):
        import csv
        for r in csv.DictReader(open(p), delimiter='\t'):
            if not r['omega']: continue
            w, dn, ds = float(r['omega']), float(r['dN']), float(r['dS']); allp.append((w, dn, ds))
            if ds > 0.01 and w < 5: keep.append((w, dn, ds))
    try:
        for l in ([] if p.endswith('.tsv') else open(p)):
            m = re.match(r'^\s*\d+\s+\d+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+([\d.]+)\s+([-\d.]+)\s+\+-\s+[\d.]+\s+([-\d.]+)', l)
            if m:
                w, dn, ds = map(float, m.groups()); allp.append((w, dn, ds))
                if ds > 0.01 and w < 5: keep.append((w, dn, ds))
    except FileNotFoundError: pass
    f = lambda v: f'{statistics.median(v):.3f}' if v else ''
    r = lambda v: f'{min(v):.3f}-{max(v):.3f}' if v else ''
    W = [x[0] for x in keep]; S = [x[2] for x in keep]
    return dict(pairs_estimable=len(allp), pairs_used=len(keep), dNdS_median=f(W), dNdS_range=r(W), dS_median=f(S), dS_range=r(S),
                pairs_dS_gt2=sum(1 for x in keep if x[2] > 2), dNdS_median_dS_le2=f([x[0] for x in keep if x[2] <= 2]))
cols = ['set', 'n_CDS', 'codons_strict', 'codons_le20gap', 'M0_omega_strict', 'M0_omega_le20gap', 'LRT_M7vsM8', 'p_M7vsM8', 'p_M8avsM8', 'BEB_sites_P>0.95']
ycols = ['pairs_estimable', 'pairs_used', 'dNdS_median', 'dNdS_range', 'dS_median', 'dS_range', 'pairs_dS_gt2', 'dNdS_median_dS_le2']
out = open(f'{X}/summary.tsv', 'w')
out.write('\t'.join(cols + [f'yn00_strict_{c}' for c in ycols] + [f'yn00_le20gap_{c}' for c in ycols]) + '\n')
for s in SETS:
    d = f'{X}/{s}'; row = {'set': s, 'n_CDS': sum(1 for l in open(f'{d}/{s}.cds.fna') if l[0] == '>'),
         'codons_strict': sites(f'{d}/{s}.codon.phy'), 'codons_le20gap': sites(f'{d}/loose/{s}.codon.phy'),
         'M0_omega_strict': m0(f'{d}/{s}.M0.out'), 'M0_omega_le20gap': m0(f'{d}/loose/{s}.M0.out')}
    try:
        l7, l8 = lnl(f'{d}/{s}.M7M8.out', 7), lnl(f'{d}/{s}.M7M8.out', 8); t = 2 * (l8 - l7)
        row.update(LRT_M7vsM8=f'{t:.2f}', p_M7vsM8=f'{chi2.sf(max(t,0), 2):.3g}')
        t2 = 2 * (l8 - lnl(f'{d}/{s}.M8a.out')); row['p_M8avsM8'] = f'{chi2.sf(max(t2,0), 1):.3g}'
        txt = open(f'{d}/{s}.M7M8.out').read(); beb = txt[txt.rfind('Bayes Empirical Bayes'):]
        beb = beb[:beb.find('The grid')] if 'The grid' in beb else beb
        row['BEB_sites_P>0.95'] = len(re.findall(r'^\s+\d+ \S\s+[\d.]+\*', beb, re.M))
    except Exception: pass
    a = yn(f'{d}/{s}.yn00.out'); pw = f'{d}/loose/{s}.yn00_pairwise.tsv'
    # yn00 deletes every gapped column, so on the <=20%-gap alignment it reproduces the strict result; use pairwise deletion where computed
    b = yn(pw) if os.path.exists(pw) else yn(f'{d}/loose/{s}.yn00.out')
    out.write('\t'.join([str(row.get(c, '')) for c in cols] + [str(a[c]) for c in ycols] + [str(b[c]) for c in ycols]) + '\n')
out.close(); print(open(f'{X}/summary.tsv').read())
