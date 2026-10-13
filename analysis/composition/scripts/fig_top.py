import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import csv, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT = f'{GMO_ROOT}/analysis/composition'
r = list(csv.DictReader(open(f'{OUT}/per_locus.tsv'), delimiter='\t'))
col = {'sus-like (T-DNA type)': '#d62728', 'cus-like (T-DNA type)': '#1f77b4', 'acc-type GDPD': '#2ca02c', 'pKerr14a-type acs-like': '#ff7f0e'}
lab = {'sus-like (T-DNA type)': 'SUS-like (T-DNA type)', 'cus-like (T-DNA type)': 'CUS-like (T-DNA type)', 'acc-type GDPD': 'ACS-type GDPD (acc-type)', 'pKerr14a-type acs-like': 'pKerr14a-type ACS-like'}
f = lambda v: float(v) if v != '' else np.nan
fig, ax = plt.subplots(1, 2, figsize=(11, 5), sharex=True, sharey=True)
for a, mode, title in [(ax[0], 'own', 'Host-genome background (>=30 genes)'), (ax[1], 'pool', 'Pooled Boletales background (Boletales loci)')]:
    for l in col:
        xs = [(f(x['bg_med_gc3']) if mode == 'own' else f(x['bolpool_med_gc3'])) for x in r if x['lineage'] == l and (x['bg_small_flag'] == 'ok' if mode == 'own' else x['order'] == 'Boletales')]
        ys = [f(x['gc3']) for x in r if x['lineage'] == l and (x['bg_small_flag'] == 'ok' if mode == 'own' else x['order'] == 'Boletales')]
        a.scatter(xs, ys, c=col[l], label=f'{lab[l]} (n={len(xs)})', s=36, edgecolor='k', linewidth=0.4, alpha=0.85)
    a.plot([0.2, 0.9], [0.2, 0.9], 'k--', lw=0.8); a.set_title(title, fontsize=10)
    a.set_xlabel('Median GC3 of background genes'); a.legend(fontsize=7, loc='upper left')
ax[0].set_ylabel('GC3 of locus CDS')
plt.tight_layout()
for e in ('png', 'svg'): plt.savefig(f'{OUT}/fig_gc3_locus_vs_bg.{e}', dpi=200)
# top deviations: primary background = own genome if >=30 genes else pooled Boletales
rows = []
for x in r:
    if x['bg_small_flag'] == 'ok': bgm, z, p, cp, d = 'own_genome', f(x['gc3_z']), f(x['gc3_pctile']), f(x['codon_dist_pctile']), f(x['delta_gc3'])
    elif x['order'] == 'Boletales': bgm, z, p, cp, d = 'boletales_pool', f(x['bolpool_gc3_z']), f(x['bolpool_gc3_pctile']), f(x['bolpool_codon_dist_pctile']), f(x['bolpool_delta_gc3'])
    else: continue
    rows.append((abs(z), x['locus_id'], x['lineage'][:3], x['organism'], x['order'], bgm, x['n_bg_genes'] if bgm == 'own_genome' else x['bolpool_n_genes'], f(x['gc3']), d, z, p, cp))
rows.sort(reverse=True)
with open(f'{OUT}/top15_deviation.tsv', 'w') as o:
    o.write('locus_id\tlineage\torganism\torder\tbackground\tn_bg\tgc3\tdelta_gc3\tgc3_z\tgc3_pctile\tcodon_dist_pctile\n')
    for t in rows[:15]: o.write('\t'.join([t[1], t[2], t[3], t[4], t[5], str(t[6])] + ['%.3f' % v for v in t[7:]]) + '\n')
print(len(rows))
