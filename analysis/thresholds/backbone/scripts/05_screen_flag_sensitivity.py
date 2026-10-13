#!/usr/bin/env python3
"""How the step-2 contamination flags of family loci (best_cat ACS/CUS/SUS/ORF358) would change if CONTAM_PID
(aggregate_step2.py line 18; applied at lines 129, 172-173) were 40..80 instead of 60. Read-only on step2_loci.tsv.
  contig  : the locus' contig carries >= 1 backbone-best locus with pid >= t
  genome  : the genome has >= 3 such contigs
  short   : the family locus' own best pid >= t and contig < 20 kb (same constant reused)
Outputs: screen_flag_by_threshold.tsv, screen_flagged_family_loci_60.tsv
"""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import csv
from collections import defaultdict

B = f'{GMO_ROOT}/analysis/thresholds/backbone'
L = f'{GMO_ROOT}/analysis/refset/screen/step2/step2_loci.tsv'
BB = {'vir', 'tra', 'rep_partition', 'mobile'}; FAM = {'ACS', 'CUS', 'SUS', 'ORF358'}
bbmax = defaultdict(float); fam = []
for r in csv.DictReader(open(L), delimiter='\t'):
    if not r['best_pident']: continue
    p = float(r['best_pident']); k = (r['gca'], r['contig'])
    if r['best_cat'] in BB: bbmax[k] = max(bbmax[k], p)
    elif r['best_cat'] in FAM: fam.append((k, p, int(r['contig_len'] or 0), r['best_cat'], r['organism'], r['contam_flag'], r['start'], r['end']))
out = open(f'{B}/screen_flag_by_threshold.tsv', 'w')
out.write('threshold\tn_contigs_flagged\tn_genomes_ge3_contigs\tn_family_loci\tfam_contig_flag\tfam_genome_flag\tfam_short_contig_flag\tfam_any_flag\n')
for t in range(40, 85, 5):
    gcount = defaultdict(int)
    for (g, c), p in bbmax.items():
        if p >= t: gcount[g] += 1
    nc = sum(1 for p in bbmax.values() if p >= t); gflag = {g for g, n in gcount.items() if n >= 3}
    cf = gf = sf = anyf = 0
    rows60 = []
    for k, p, clen, cat, org, cflag, s, e in fam:
        a = bbmax.get(k, 0) >= t; b = (not a) and k[0] in gflag; c = (not a and not b) and p >= t and clen < 20000
        cf += a; gf += b; sf += c; anyf += a or b or c
        if t == 60 and (a or b or c): rows60.append((k[0], org, k[1], clen, s, e, cat, p, 'contig' if a else 'genome' if b else 'short_contig', cflag))
    out.write(f'{t}\t{nc}\t{len(gflag)}\t{len(fam)}\t{cf}\t{gf}\t{sf}\t{anyf}\n')
    if t == 60:
        with open(f'{B}/screen_flagged_family_loci_60.tsv', 'w') as fo:
            fo.write('gca\torganism\tcontig\tcontig_len\tstart\tend\tcat\tbest_pid\tflag_recomputed\tflag_in_step2_loci\n')
            for r in rows60: fo.write('\t'.join(map(str, r)) + '\n')
out.close()
