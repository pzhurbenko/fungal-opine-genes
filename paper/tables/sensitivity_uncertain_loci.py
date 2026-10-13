#!/usr/bin/env python3
"""Sensitivity of the Table 2 counts and the Figure 5 inputs to excluding 'uncertain' loci (integration_status.tsv).
Does NOT change production files. Figure 5 inputs are recomputed by running a copy of analysis/boletales_mash_517/make_fig2_genus.py
with the uncertain loci removed from TableS1_all_lineages.tsv and all outputs redirected to a temporary directory.
Run (env PlantCad2): python integration_status.py && python sensitivity_uncertain_loci.py
Outputs: sensitivity_uncertain_loci.tsv, intact_but_uncertain_loci.tsv (this directory).
Columns: all_loci | confirmed_supported_only."""
import csv, os, re, subprocess, sys, tempfile, glob
from collections import defaultdict
P = os.path.dirname(os.path.abspath(__file__)); G = P + '/../..'
rd = lambda f: list(csv.DictReader(open(f), delimiter='\t'))
st = {r['locus_id']: r for r in rd(P + '/integration_status.tsv')}
s1a = {r['locus_id']: r for r in rd(P + '/../supplementary/TableS1a_Loci_of_bacterial_origin.tsv')}
def keep_a(l): return st[l]['status'] != 'uncertain'

sp = lambda s: re.sub(r' \(cluster [^)]*\)', '', s)
rows = []
lins = sorted({r['lineage'] for r in s1a.values()}) + ['TOTAL']
for lin in lins:
    for metric in ('loci', 'intact_loci', 'species', 'assemblies'):
        v = []
        for keep in (lambda l: True, keep_a):
            ls = [r for l, r in s1a.items() if (lin == 'TOTAL' or r['lineage'] == lin) and keep(l)]
            if metric == 'loci': v.append(len(ls))
            elif metric == 'intact_loci': v.append(sum(r['status'] == 'intact' for r in ls))
            elif metric == 'species': v.append(len({sp(r['species']) for r in ls}))
            else: v.append(len({r['assembly'] for r in ls}))
        rows.append(['Table 2: ' + lin, metric] + v)
# Figure 5 inputs
src = open(G + '/analysis/boletales_mash_517/make_fig2_genus.py').read()
tmp = tempfile.mkdtemp(prefix='fig5sens_')
res = {}
for tag, keep in (('all', lambda l: True), ('cs', keep_a)):
    d = tmp + '/' + tag; os.makedirs(d)
    s1 = rd(G + '/paper/tables/TableS1_all_lineages.tsv')
    with open(d + '/S1.tsv', 'w') as f:
        w = csv.DictWriter(f, list(s1[0]), delimiter='\t', lineterminator='\n'); w.writeheader()
        w.writerows(r for r in s1 if r['locus_id'] not in st or keep(r['locus_id']))
    code = src.replace("OUT = f'{G}/paper/figures'", "OUT = %r" % d).replace("{G}/paper/tables/TableS1_all_lineages.tsv", d + '/S1.tsv')
    assert d in code
    open(d + '/run.py', 'w').write(code)
    subprocess.run([sys.executable, d + '/run.py'], check=True, cwd=G + '/analysis/boletales_mash_517', env=dict(os.environ, MPLBACKEND='Agg'), stdout=subprocess.DEVNULL)
    res[tag] = rd(d + '/Fig2_genus_counts.tsv')
LN = ['sus-like', 'cus-like', 'acc-type_GDPD', 'pKerr14a-type_acs-like', 'nos-subfamily_ODH']
for ln in LN + ['ALL5']:
    for m, f in (('Fig5: genus rows with the lineage', lambda r, l: int(int(r[l + '_k_any']) > 0)),
                 ('Fig5: genus rows with an intact copy', lambda r, l: int(int(r[l + '_k_intact']) > 0)),
                 ('Fig5: sum over rows of species clusters with the lineage', lambda r, l: int(r[l + '_k_any'])),
                 ('Fig5: sum over rows of clusters with an intact copy', lambda r, l: int(r[l + '_k_intact']))):
        v = [sum(f(r, l) for r in res[t] for l in (LN if ln == 'ALL5' else [ln])) for t in ('all', 'cs')]
        rows.append(['Fig5: ' + ln, m.split(': ')[1], *v])
with open(P + '/sensitivity_uncertain_loci.tsv', 'w') as f:
    f.write('lineage\tmetric\tall_loci\tconfirmed_supported_only\n')
    for r in rows: f.write('\t'.join(map(str, r)) + '\n')
# intact but uncertain + dN/dS membership
dn = defaultdict(list)
for fa in glob.glob(G + '/analysis/dnds/*.cds.fna') + glob.glob(G + '/analysis/dnds/extra/*/*.cds.fna'):
    for l in open(fa):
        if l.startswith('>'): dn[l[1:].split()[0]].append(os.path.basename(fa).replace('.cds.fna', ''))
with open(P + '/intact_but_uncertain_loci.tsv', 'w') as f:
    f.write('locus_id\tlineage\tspecies\tassembly\tcontig_length_bp\trules_matched\tin_dnds_sets\n')
    n = 0
    for l, r in st.items():
        if r['status'] == 'uncertain' and s1a[l]['status'] == 'intact':
            n += 1; f.write('\t'.join([l, r['lineage'], r['species'], r['assembly'], s1a[l]['contig_length_bp'], r['rules_matched'], ';'.join(dn.get(l, [])) or 'no']) + '\n')
print(n, 'intact-but-uncertain loci;', sum(1 for l, r in st.items() if r['status'] == 'uncertain' and s1a[l]['status'] == 'intact' and l in dn), 'in dN/dS sets')
print(open(P + '/sensitivity_uncertain_loci.tsv').read())
