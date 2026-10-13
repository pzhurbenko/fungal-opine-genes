#!/usr/bin/env python
"""Join flank genes with best DIAMOND hits (plasmid catalog 7105 prot / 42 plasmids; 389 representatives).
backbone_ge60 = 1 if best catalog hit is vir/tra/rep_partition/mobile with identity >= 60 %."""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import csv
D = f'{GMO_ROOT}/analysis/fusarium_odh_island/rhypophila_flank'
def best(p):
    b = {}
    for r in csv.reader(open(p), delimiter='\t'):
        if r[0] not in b or float(r[11]) > float(b[r[0]][11]): b[r[0]] = r
    return b
bc, bq = best(f'{D}/flank50_vs_plasmid.tsv'), best(f'{D}/flank50_vs_queries389.tsv')
BB = {'vir', 'tra', 'rep_partition', 'mobile'}
rows = list(csv.DictReader(open(f'{D}/flank_genes.tsv'), delimiter='\t'))
cols = list(rows[0]) + ['cat_best_hit', 'cat_category', 'cat_plasmid', 'cat_organism', 'cat_product', 'cat_pident', 'cat_aln_len', 'cat_evalue',
                        'rep389_best_hit', 'rep389_pident', 'backbone_ge60']
with open(f'{D}/flank50_hits_table.tsv', 'w') as f:
    f.write('\t'.join(cols) + '\n')
    for r in rows:
        h = bc.get(r['protein_id']); q = bq.get(r['protein_id'])
        if h:
            a = h[1].split('|')
            r.update(cat_best_hit=a[0], cat_category=a[1], cat_plasmid=a[2], cat_organism=a[3], cat_product=a[5], cat_pident=h[2], cat_aln_len=h[3], cat_evalue=h[10])
        if q: r.update(rep389_best_hit=q[1], rep389_pident=q[2])
        r['backbone_ge60'] = int(bool(h) and h[1].split('|')[1] in BB and float(h[2]) >= 60)
        f.write('\t'.join(str(r.get(c, '')) for c in cols) + '\n')
