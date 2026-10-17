#!/usr/bin/env python3
"""Human-readable summary of step-2 screening (reads step2_*.tsv). Prints to stdout; redirect to step2_summary.txt."""
import csv
from collections import Counter, defaultdict
loci = list(csv.DictReader(open('step2_loci.tsv'), delimiter='\t'))
fc = list(csv.DictReader(open('step2_family_calls.tsv'), delimiter='\t'))
ct = list(csv.DictReader(open('step2_contamination.tsv'), delimiter='\t'))
gs = list(csv.DictReader(open('step2_genome_summary.tsv'), delimiter='\t'))
print(f"Assemblies processed: {len(gs)}; loci (hit regions): {len(loci)}")
# ---- family calls: unique loci, confident only
uniq = {}
for r in fc:
    k = (r['gca'], r['contig'], r['start'], r['family'])
    if r['call'] == 'ambiguous': continue
    if k not in uniq or float(r['score']) > float(uniq[k]['score']): uniq[k] = r
print("\n== Our families: confident calls (subgroup > OUTGROUP, gap>=20, score>=30) ==")
print(f"{'family':8s} {'loci':>8s} {'genomes':>8s} {'new loci':>10s} {'new genomes':>14s} {'of them contam.':>15s}")
by_fam = defaultdict(list)
for r in uniq.values(): by_fam[r['family']].append(r)
for fam in ['ACS', 'CUS', 'SUS', 'ORF358']:
    rs = by_fam[fam]; new = [r for r in rs if r['is_new'] == 'yes']
    known_g = {r['gca'] for r in rs if r['is_new'] == 'no'}; new_g = {r['gca'] for r in new} - known_g
    print(f"{fam:8s} {len(rs):8d} {len({r['gca'] for r in rs}):8d} {len(new):10d} {len(new_g):14d} {sum(1 for r in new if r['contam_flag']!='no'):15d}")
    sub = Counter(r['call'] for r in rs); subn = Counter(r['call'] for r in new)
    print("         subgroups: " + ", ".join(f"{k} {v} (new {subn[k]})" for k, v in sub.most_common()))
    cls = Counter(r['tax_class'] for r in new)
    print("         new by class: " + ", ".join(f"{k or '?'} {v}" for k, v in cls.most_common(8)))
# ---- new loci lists for rare families
for fam, n in [('SUS', 40), ('ACS', 40)]:
    new = sorted([r for r in by_fam[fam] if r['is_new'] == 'yes'], key=lambda r: -float(r['score']))
    print(f"\n== New loci {fam} (top {n} by score) ==")
    for r in new[:n]:
        print(f"  {r['organism'][:30]:30s} {r['tax_class'][:16]:16s} {r['contig']:18s} {r['start']:>9s}-{r['end']:<9s} {r['call']:10s} score {r['score']:>4s} gap {r['gap']:>4s} dmnd {r['best_subject_cat']:13s} {r['best_pident']:>5s}% contam={r['contam_flag']}")
# ---- multi-family genomes
gf = defaultdict(set)
for r in uniq.values(): gf[r['gca']].add(r['family'])
multi = {g: f for g, f in gf.items() if len(f) >= 2}
print(f"\n== Genomes with >=2 families (confident calls): {len(multi)} ==")
org = {r['gca']: r['organism'] for r in fc}
newmulti = [g for g in multi if any(r['is_new'] == 'yes' for r in uniq.values() if r['gca'] == g)]
for g in sorted(multi, key=lambda g: -len(multi[g]))[:40]:
    print(f"  {org[g][:34]:34s} {g:16s} {','.join(sorted(multi[g])):22s} {'(has new loci)' if g in newmulti else ''}")
# ---- T-DNA lineages and opine outliers outside our families, excluding contaminated contigs/genomes
print("\n== T-DNA lineages outside the 4 families: hits >=45% id, cov>=0.5, no contamination ==")
tl = defaultdict(list)
for r in loci:
    if r['best_cat'] == 'TDNA' and float(r['best_pident'] or 0) >= 45 and float(r['subj_cov']) >= 0.5 and r['contam_flag'] == 'no':
        tl[r['best_lineage'][:30]].append(r)
for k, v in sorted(tl.items(), key=lambda x: -len(x[1])):
    ex = sorted(v, key=lambda r: -float(r['best_pident']))[:4]
    print(f"  {k:30s} {len(v):4d}  genomes {len({r['gca'] for r in v}):3d}  e.g. " + "; ".join(f"{r['organism'][:22]} {r['best_pident']}% {r['contig_len']}bp" for r in ex))
print("\n== Opine catabolism/transport (opine_related) >=45%, cov>=0.5, no contamination ==")
ol = defaultdict(list)
for r in loci:
    if r['best_cat'] == 'opine_related' and float(r['best_pident'] or 0) >= 45 and float(r['subj_cov']) >= 0.5 and r['contam_flag'] == 'no':
        ol[r['best_lineage'][:30]].append(r)
for k, v in sorted(ol.items(), key=lambda x: -len(x[1]))[:15]:
    ex = sorted(v, key=lambda r: -float(r['best_pident']))[:3]
    print(f"  {k:30s} {len(v):4d}  genomes {len({r['gca'] for r in v}):3d}  max {max(float(r['best_pident']) for r in v):.0f}%  e.g. " + "; ".join(f"{r['organism'][:22]} {r['best_pident']}%" for r in ex))
# ---- contamination
print("\n== Contamination: assemblies with contigs carrying plasmid-backbone genes (vir/tra/rep/mobile) >=60% ==")
g = defaultdict(list)
for r in ct: g[(r['gca'], r['organism'])].append((int(r['contig_len'] or 0), int(r['n_backbone_hits_ge60']), int(r['n_hi_id_hits'])))
for (gca, o), v in sorted(g.items(), key=lambda x: -sum(c[0] for c in x[1])):
    print(f"  {o[:32]:32s} {gca:16s} contigs {len(v):3d}  sum {sum(c[0] for c in v)/1e3:7.0f} kb  max {max(c[0] for c in v)/1e3:6.0f} kb  backbone hits {sum(c[1] for c in v):3d}  all hits>=60% {sum(c[2] for c in v):3d}")
