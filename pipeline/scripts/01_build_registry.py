#!/usr/bin/env python3
"""Stage 01 — registry of loci from all_blast_stats.tsv.

One locus = one family x one contig x one cluster of HSPs closer than registry.merge_gap.
Protein-record subjects (blastp hits to annotated proteins, e.g. KAG..., XP_...) go to a
separate table: they already have gene models in NCBI.

Outputs (out/registry/):
  loci.tsv              one row per locus
  contigs.tsv           unique WGS contigs with the loci they carry
  protein_subjects.tsv  subjects that are protein accessions (no annotation needed)
"""
import os, re, sys
from collections import defaultdict
from common import load_config, read_tsv, write_tsv, IS_WGS_CONTIG, log

cfg = load_config()
P = cfg['_paths']; os.makedirs(P['registry'], exist_ok=True)
rows = read_tsv(P['blast_stats'])
log('blast_stats rows:', len(rows))

fam_norm = lambda f: f.upper()
by_key = defaultdict(list); prot = []
for r in rows:
    acc = r['Subj_Accession'].strip(); fam = fam_norm(r['T-DNA_homolog'])
    if IS_WGS_CONTIG.match(acc):
        by_key[(fam, acc.split('.')[0])].append(r)
    else:
        prot.append({'fam': fam, 'protein_acc': acc, 'organism': r['Species'], 'genus': r['Genera'],
                     'query': r['Query_name'], 'subj_start': r['Subj_start'], 'subj_end': r['Subj_end'],
                     'identity': r['Identity'], 'evalue': r['E_value'], 'query_cover': r['Query_cover'],
                     'was_deleted': r['Was_deleted'], 'reason': r['Reason'], 'tree_clade': r['Tree_clade'],
                     'division': r['Division'], 'class': r['Class'], 'order': r['Order']})

gap = cfg['registry']['merge_gap']
loci = []; contig_fams = defaultdict(set); rep_contigs = set(r['Rep_seq_mmseqs'].split('.')[0] for r in rows)
for (fam, acc), hs in sorted(by_key.items()):
    hs = sorted(hs, key=lambda r: int(float(r['Subj_start'])))
    clusters = []
    for r in hs:
        s, e = int(float(r['Subj_start'])), int(float(r['Subj_end']))
        s, e = min(s, e), max(s, e)
        if clusters and s - clusters[-1]['end'] <= gap:
            c = clusters[-1]; c['end'] = max(c['end'], e); c['rows'].append(r)
        else:
            clusters.append({'start': s, 'end': e, 'rows': [r]})
    for i, c in enumerate(clusters, 1):
        rr = c['rows']; r0 = rr[0]
        ident = max(float(x['Identity'] or 0) for x in rr)
        ev = min(float(x['E_value'] or 1) for x in rr)
        deleted = [x['Was_deleted'] for x in rr]
        loci.append({
            'locus_id': f'{fam}_{acc}_L{i}', 'fam': fam, 'contig': acc, 'genome': r0['WGS_Accession'],
            'organism': r0['Species'], 'genus': r0['Genera'], 'division': r0['Division'], 'class': r0['Class'],
            'order': r0['Order'], 'family_tax': r0['Family'],
            'hsp_start': c['start'], 'hsp_end': c['end'], 'n_hsp': len(rr),
            'best_identity': f'{ident:.2f}', 'best_evalue': f'{ev:.2e}',
            'queries': ';'.join(sorted(set(x['Query_name'] for x in rr))),
            'strategies': ';'.join(sorted(set(x['Blast_strategy'] for x in rr))),
            'was_deleted': 'all' if all(d == 'yes' for d in deleted) else ('some' if 'yes' in deleted else 'no'),
            'reasons': ';'.join(sorted(set(x['Reason'] for x in rr if x['Reason']))),
            'tree_clade': ';'.join(sorted(set(x['Tree_clade'] for x in rr if x['Tree_clade']))),
            'is_rep_contig': 'yes' if acc in rep_contigs else 'no',
        })
        contig_fams[acc].add(fam)

for l in loci:
    l['other_fams_on_contig'] = ','.join(sorted(contig_fams[l['contig']] - {l['fam']}))
write_tsv(os.path.join(P['registry'], 'loci.tsv'), loci)
contigs = [{'contig': a, 'fams': ','.join(sorted(f)), 'n_loci': sum(1 for l in loci if l['contig'] == a),
            'organism': next(l['organism'] for l in loci if l['contig'] == a),
            'hsp_min': min(l['hsp_start'] for l in loci if l['contig'] == a),
            'hsp_max': max(l['hsp_end'] for l in loci if l['contig'] == a)} for a, f in sorted(contig_fams.items())]
write_tsv(os.path.join(P['registry'], 'contigs.tsv'), contigs)
write_tsv(os.path.join(P['registry'], 'protein_subjects.tsv'), prot)
from collections import Counter
log('loci:', len(loci), dict(Counter(l['fam'] for l in loci)))
log('contigs:', len(contigs), ' multi-family contigs:', sum(1 for c in contigs if ',' in c['fams']))
log('protein subjects:', len(prot))
