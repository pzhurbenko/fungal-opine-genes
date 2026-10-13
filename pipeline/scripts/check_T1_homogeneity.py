#!/usr/bin/env python3
"""QC — are the T1 reference sets homogeneous (all true T-DNA genes) or do they contain chromosomal paralogs?

For every T1 protein of a family:
  id_to_canonical   normalised identity (shared residues / shorter length) to the canonical T-DNA protein of the family
  max_id_other_T1   best identity to any other T1 protein (an outlier has a low value)
  mean_id_T1        mean identity to the other T1 proteins
  id_to_paralog     identity to a known CHROMOSOMAL paralog of the same fold (if one is configured)
  pfam              Pfam architecture from analysis/domains/sensitive/<FAM>_domain_arch.tsv (if present)
Canonical / paralog accessions are fetched from NCBI protein. Output: guides/<FAM>_T1_check.tsv (sorted by id_to_canonical).
"""
import os, re, sys, statistics
from common import load_config, read_fasta, write_fasta, write_tsv, read_tsv, eutils, log
from Bio import Align
from Bio.Align import substitution_matrices

cfg = load_config(); P = cfg['_paths']
# canonical T-DNA proteins (most frequent T2 best hits / literature) and chromosomal paralogs of the same fold
CANON = {'ACS': 'CUX65485.1',        # putative agrocinopine synthase, pTi A. tumefaciens Kerr 14
         'CUS': 'BAB16139.1',        # riorf20 = cucumopine synthase, pRi R. rhizogenes
         'SUS': 'WP_174080856.1',    # pRi T-DNA sus-like (0.98 to Agro_WP_269710210 ref)
         'ORF358': 'WP_269710219.1'} # pRi orf358-like
PARALOG = {'ACS': 'glycerophosphodiester phosphodiesterase[Title] AND Agrobacterium fabrum[Organism] AND refseq[filter]',
           'SUS': 'WP_114836055.1',   # chromosomal ornithine cyclodeaminase A. tumefaciens
           'CUS': None, 'ORF358': None}
al = Align.PairwiseAligner(); al.substitution_matrix = substitution_matrices.load('BLOSUM62')
al.open_gap_score = -10; al.extend_gap_score = -0.5; al.mode = 'local'
def pid(a, b):
    x = al.align(a, b)[0]; i = sum(1 for p, q in zip(x[0], x[1]) if p == q and p != '-'); return i / min(len(a), len(b))
def fetch(acc_or_query):
    if '[' in acc_or_query:
        ids = eutils(cfg, 'esearch.fcgi', {'db': 'protein', 'term': acc_or_query, 'retmax': 1, 'retmode': 'json'}, as_json=True)['esearchresult']['idlist']
        if not ids: return None
        acc_or_query = ids[0]
    txt = eutils(cfg, 'efetch.fcgi', {'db': 'protein', 'id': acc_or_query, 'rettype': 'fasta', 'retmode': 'text'})
    h, s = txt.split('\n', 1); return h[1:], s.replace('\n', '')

for fam in cfg['families']:
    t1 = [(h, s.replace('*', '')) for h, s in read_fasta(os.path.join(P['curated_dir'], f'{fam}.fasta')) if h.startswith(('Agro_', 'Plant_'))]
    canon = fetch(CANON[fam]); par = fetch(PARALOG[fam]) if PARALOG.get(fam) else None
    pf = os.path.join(cfg['project_dir'], 'analysis', 'domains', 'sensitive', f'{fam}_domain_arch.tsv')
    arch = {r['clean_id']: r['architecture'] for r in read_tsv(pf)} if os.path.exists(pf) else {}
    rows = []
    for h, s in t1:
        others = [pid(s, s2) for h2, s2 in t1 if h2 != h]
        rows.append({'fam': fam, 'ref': h[:90], 'len': len(s), 'id_to_canonical': f'{pid(s, canon[1]):.2f}',
                     'max_id_other_T1': f'{max(others):.2f}', 'mean_id_T1': f'{statistics.mean(others):.2f}',
                     'id_to_paralog': f'{pid(s, par[1]):.2f}' if par else '', 'pfam': arch.get(h, '')})
    rows.sort(key=lambda r: float(r['id_to_canonical']))
    write_tsv(os.path.join(P['guides'], f'{fam}_T1_check.tsv'), rows)
    log(f'{fam}: canonical={canon[0][:60]} | paralog={par[0][:60] if par else "-"}')
    for r in rows: print(f"{fam:6} {r['ref'][:70]:70} len={r['len']:4} canon={r['id_to_canonical']} maxT1={r['max_id_other_T1']} meanT1={r['mean_id_T1']} par={r['id_to_paralog'] or '-':5} {r['pfam'][:35]}")
