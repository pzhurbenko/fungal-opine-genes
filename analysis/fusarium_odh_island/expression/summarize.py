#!/usr/bin/env python3
# Collect results_<run>_<target>.json -> summary.tsv (one row per run x island gene), same columns as chs_fusarium/expression/summary.tsv
import json, glob, os, statistics
GENES = {'FPOAC2_10362': 'odh (CM100256 copy)', 'FPOAC2_10363': 'pepP (CM100256 copy)', 'FPOAC2_10364': 'g3 (CM100256 copy)',
         'FPOAC2_07383': 'pepP (CM100255 copy, partial)', 'FPOAC2_07382': 'g3 (CM100255 copy, fragment)'}
H = ['run', 'gene', 'locus_tag', 'mapped_pct', 'n_reads_used', 'reads_on_target', 'target_TPM', 'target_reads_per_kb', 'pct_all_genes_strictly_below', 'pct_all_genes_leq',
     'genes_with_zero_reads_pct', 'nb_TPM_min', 'nb_TPM_median', 'nb_TPM_max', 'nb_zero_genes', 'nb_n_genes', 'median_TPM_genes_with_reads',
     'target_TPM_over_median_expressed', 'actin_TPM', 'EF1a_TPM', 'GAPDH_TPM', 'tubulinB_TPM', 'lib_sense_frac', 'target_breadth', 'target_bins_max_share']
rows = []
for f in sorted(glob.glob('results_*_FPOAC2_*.json')):
    d = json.load(open(f)); T = d['target']; run = d['run']; lt = T['lt']
    c = {x['label']: x['tpm'] for x in d['controls']}
    rows.append([run, GENES.get(lt, lt), lt, d['mapped_pct'], d['n_primary_reads'], T['n'], T['tpm'], round(1000 * T['n'] / T['len'], 2), T['pct_lt'], T['pct_le'],
                 round(100 * (d['n_genes'] - d['genes_with_reads']) / d['n_genes'], 1), d['nb_tpm_min'], d['nb_tpm_median'], d['nb_tpm_max'], d['nb_zero'], len(d['neighbours']),
                 d['median_tpm_expressed'], round(T['tpm'] / d['median_tpm_expressed'], 3) if d['median_tpm_expressed'] else '', c.get('actin', ''), c.get('EF1a', ''), c.get('GAPDH', ''),
                 c.get('tubulin beta', ''), d['strand_sense_frac_genes_ge200'], d['target_breadth'], max(d['target_bin_share'])])
rows.sort(key=lambda r: (list(GENES).index(r[2]) if r[2] in GENES else 9, r[0]))
with open('summary.tsv', 'w') as fh:
    fh.write('\t'.join(H) + '\n')
    for r in rows: fh.write('\t'.join(map(str, r)) + '\n')
# region counts (incl. unannotated odh copy on CM100255)
reg = {}
for l in open('region_counts.tsv'):
    run, R, m1, al = l.rstrip('\n').split('\t'); reg.setdefault(R, {})[run] = (int(m1.split('=')[1]), int(al.split('=')[1]))
print('runs:', len({r[0] for r in rows}))
for lt, nm in GENES.items():
    v = [r for r in rows if r[2] == lt]
    if v: print(f"{nm:32s} reads {min(r[5] for r in v)}-{max(r[5] for r in v)}  TPM {min(r[6] for r in v)}-{max(r[6] for r in v)}  median TPM {statistics.median(r[6] for r in v):.2f}  pct {min(r[8] for r in v)}-{max(r[8] for r in v)}  breadth {min(r[23] for r in v)}-{max(r[23] for r in v)}")
for R, d in reg.items():
    print(f"region {R}: MAPQ>=1 reads {min(v[0] for v in d.values())}-{max(v[0] for v in d.values())}, all primary {min(v[1] for v in d.values())}-{max(v[1] for v in d.values())} over {len(d)} runs")
v = [r for r in rows if r[2] == 'FPOAC2_10362']
print('controls (TPM, odh rows): actin', min(r[18] for r in v), '-', max(r[18] for r in v), '; EF1a', min(r[19] for r in v), '-', max(r[19] for r in v), '; median expressed', min(r[16] for r in v), '-', max(r[16] for r in v))
