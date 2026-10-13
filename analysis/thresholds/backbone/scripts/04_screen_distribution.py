#!/usr/bin/env python3
"""Distribution of the backbone signal in the real step-2 screen (2,033 assemblies), read from the existing
analysis/refset/screen/step2/step2_loci.tsv (read-only). For every contig: max identity of a locus whose best subject
is a backbone protein (vir/tra/rep_partition/mobile; best_pident = max HSP identity of the best chain, the value the
60 % rule is applied to) and the number of such loci. Outputs:
  screen_backbone_contigs.tsv  (one row per contig with >= 1 backbone-best locus)
  screen_backbone_hist.tsv     (5 % bins: n contigs, n genomes, median contig length, n contigs < 20 kb,
                                n contigs with >= 2 backbone loci, n contigs carrying a family-lineage locus)
  screen_backbone_scan.tsv     (threshold 30..100: n contigs / genomes flagged)
"""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import csv, statistics
from collections import defaultdict

B = f'{GMO_ROOT}/analysis/thresholds/backbone'
L = f'{GMO_ROOT}/analysis/refset/screen/step2/step2_loci.tsv'
BB = {'vir', 'tra', 'rep_partition', 'mobile'}
FAMC = {'ACS', 'CUS', 'SUS', 'ORF358'}
c = {}
fam_on = set()
for r in csv.DictReader(open(L), delimiter='\t'):
    k = (r['gca'], r['contig'])
    if r['best_cat'] in FAMC: fam_on.add(k)
    if r['best_cat'] not in BB or not r['best_pident']: continue
    p = float(r['best_pident'])
    d = c.setdefault(k, dict(gca=r['gca'], organism=r['organism'], tax_class=r['tax_class'], contig=r['contig'],
                             contig_len=int(r['contig_len'] or 0), max_pid=0.0, n_loci=0, pids=[], best=''))
    d['n_loci'] += 1; d['pids'].append(p)
    if p > d['max_pid']: d['max_pid'] = p; d['best'] = f"{r['best_cat']}:{r['best_lineage']}"
with open(f'{B}/screen_backbone_contigs.tsv', 'w') as fo:
    fo.write('gca\torganism\ttax_class\tcontig\tcontig_len\tmax_backbone_pid\tn_backbone_loci\tn_backbone_loci_ge60\tbest_label\tfamily_locus_on_contig\n')
    for k, d in sorted(c.items(), key=lambda x: -x[1]['max_pid']):
        fo.write(f"{d['gca']}\t{d['organism']}\t{d['tax_class']}\t{d['contig']}\t{d['contig_len']}\t{d['max_pid']:.1f}\t"
                 f"{d['n_loci']}\t{sum(p >= 60 for p in d['pids'])}\t{d['best']}\t{int(k in fam_on)}\n")
bins = defaultdict(list)
for k, d in c.items(): bins[min(int(d['max_pid'] // 5) * 5, 95)].append((k, d))
with open(f'{B}/screen_backbone_hist.tsv', 'w') as fo:
    fo.write('pid_bin\tn_contigs\tn_genomes\tmedian_contig_len\tn_contigs_lt20kb\tn_contigs_ge2_bb_loci\tn_with_family_locus\n')
    for b in sorted(bins):
        xs = bins[b]
        fo.write(f"{b}-{b + 5}\t{len(xs)}\t{len({k[0] for k, _ in xs})}\t{statistics.median(d['contig_len'] for _, d in xs):.0f}\t"
                 f"{sum(d['contig_len'] < 20000 for _, d in xs)}\t{sum(d['n_loci'] >= 2 for _, d in xs)}\t{sum(k in fam_on for k, _ in xs)}\n")
with open(f'{B}/screen_backbone_scan.tsv', 'w') as fo:
    fo.write('threshold\tn_contigs_flagged\tn_genomes_flagged\tn_flagged_contigs_with_family_locus\n')
    for t in range(30, 101, 2):
        ks = [k for k, d in c.items() if d['max_pid'] >= t]
        fo.write(f"{t}\t{len(ks)}\t{len({k[0] for k in ks})}\t{sum(k in fam_on for k in ks)}\n")
print('contigs with backbone-best loci', len(c))
