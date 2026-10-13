#!/usr/bin/env python3
"""Per-locus integration status (confirmed / supported / uncertain) for the 254 placed loci of Table S1a.
Inputs (final-numbered supplementary TSVs, paper/supplementary/): TableS1a, TableS8 (read depth), TableS9 (independent assemblies),
TableS12 (genes on the island contigs). Output: integration_status.tsv, integration_status_summary.tsv (this directory).
Run (env PlantCad2): python integration_status.py
Rules (a locus gets the first matching class):
 confirmed  = any of
   F2  fungal genes on both flanks: S1a genomic_context == 'fungal context' (assumed to mean fungal neighbours on both sides),
       or S12 lists annotated genes on the same contig on both sides of the locus
   R   reproduced with flanks in an independent assembly (S9: locus covered >=0.5 and both flanks covered >=0.5)
   D2  read depth ratio (S8 locus_depth_ratio) 0.5-2 AND junction reads (S8) > 0 on both flanks
 Lineage-specific evidence (2026-10-13 refinement; S13a/S13b have no neighbour-gene columns, only contig geometry and assembly level):
   CHR  locus on a chromosome-level sequence (S13a / S13b column assembly_level == 'Chromosome') at >= 1 kb from both contig ends, not SHORT/EDGE: treated as F2 (fungal flanks)
   R2   the same locus (chs-like; odh/pepP/orf3 genes of the island; same species) found in >= 2 assemblies of one species (S13a/S13b: species, assembly), not SHORT/EDGE
   (nos-ODH: S10b/S12 - annotated island genes on both sides = F2 via S12; otherwise LONG below)
 uncertain (first override; applies if not confirmed) =
   SHORT  contig shorter than locus length + 2 kb (1 kb on either side)
   EDGE   locus touches a contig end (start <= 1 or end >= contig length; S8 'edge' flanks)
 supported = any of
   LONG contig >= 100 kb, locus not at a contig end, genomic_context not bacterial, and depth ratio not outside 0.5-2 (a fungal-sized contig is itself evidence)
   F1  fungal gene on one flank only: S1a genomic_context == 'mixed', or S12 genes on one side only
   D1  depth ratio 0.5-2 without junction reads on both sides
 uncertain (remaining) =
   DEPTH  depth ratio < 0.5 or > 2 without other support (no F/R/D rule met)
   NONB   no neighbour information and no other support
"""
import csv, os, re
from collections import Counter, defaultdict
P = os.path.dirname(os.path.abspath(__file__)); S = P + '/../supplementary/'
rd = lambda f: list(csv.DictReader(open(S + f), delimiter='\t'))
s1a = rd('TableS1a_Loci_of_bacterial_origin.tsv')
s8 = {r['locus_id']: r for r in rd('TableS8_Read_depth_at_the_loci.tsv')}
s9 = rd('TableS9_Reproducibility_of_loci_in_independent_assemblies.tsv')
s12 = rd('TableS12_Gene_content_of_the_islands.tsv')
lvl = {}
for fn in ('TableS13a_chs_like_loci_in_Fusarium_and_relatives.tsv', 'TableS13b_Loci_of_the_odh_like_island_in_Fusarium_and_relati.tsv'):
    for r in rd(fn): lvl[r['locus_id']] = (r['assembly_level'], r['species'], r['assembly'], r.get('island_gene', 'chs'))
asm_by = defaultdict(set)   # (gene key, species) -> assemblies
gk = lambda lid, g: 'chs' if g == 'chs' else g.split()[0].replace('orf3', 'g3')
for lid, (al, sp, asm, g) in lvl.items(): asm_by[(gk(lid, g), sp)].add(asm)
fl = lambda x: float(x) if x not in ('', 'nan', 'NA', None) else None

repro = set()
for r in s9:
    lc, l, rr = fl(r['locus_covered_fraction']), fl(r['left_flank_covered_fraction']), fl(r['right_flank_covered_fraction'])
    if lc is not None and l is not None and rr is not None and lc >= .5 and l >= .5 and rr >= .5:
        repro.update(x for x in re.split(r'[;,\s]+', r['loci']) if x)

# S12: sides with genes, per S1a locus (same assembly+contig, overlapping coordinates of the focal row)
byc = defaultdict(list)
for r in s12: byc[(r['assembly'], r['contig'])].append(r)
s12sides = {}
for o in s1a:
    rows = byc.get((o['assembly'], o['contig']))
    if not rows: continue
    a, b = int(o['start']), int(o['end'])
    foc = [r for r in rows if r['gene_class'] == 'focal locus' and int(r['start']) <= b and int(r['end']) >= a]
    if not foc: continue
    nb = [r for r in rows if r['gene_class'] != 'focal locus']
    s12sides[o['locus_id']] = (any(int(r['end']) < a for r in nb), any(int(r['start']) > b for r in nb))

out = []
for o in s1a:
    lid = o['locus_id']; a, b, L = int(o['start']), int(o['end']), int(o['contig_length_bp'])
    d = s8.get(lid); ratio = fl(d['locus_depth_ratio']) if d else None
    jl = jr = None
    if d:
        jl = d['reads_spanning_left_flank']; jr = d['reads_spanning_right_flank']
        jl = None if jl == 'edge' else fl(jl); jr = None if jr == 'edge' else fl(jr)
    ctx = o['genomic_context']
    left = right = None
    if lid in s12sides: left, right = s12sides[lid]
    elif ctx == 'fungal context': left = right = True
    elif ctx == 'mixed': left, right = True, False   # one flank only (which side is not recorded)
    rules = []; conf = []; sup = []; unc = []
    if left and right: conf.append('F2')
    if lid in repro: conf.append('R')
    inr = ratio is not None and .5 <= ratio <= 2
    if inr and jl and jr and jl > 0 and jr > 0: conf.append('D2')
    short = L < (b - a + 1) + 2000
    edge = a <= 1 or b >= L or (d is not None and 'edge' in (d['reads_spanning_left_flank'], d['reads_spanning_right_flank']))
    if lid in lvl and not short and not edge:
        al, sp, asm, g = lvl[lid]
        if al == 'Chromosome' and a > 1000 and L - b > 1000: conf.append('CHR')
        if len(asm_by[(gk(lid, g), sp)]) >= 2: conf.append('R2')
    outdepth = ratio is not None and not inr
    if L >= 100000 and not edge and not short and ctx != 'bacterial context' and not outdepth: sup.append('LONG')
    if (left or right) and not (left and right): sup.append('F1')
    if inr and 'D2' not in conf: sup.append('D1')
    if short: unc.append('SHORT')
    if edge: unc.append('EDGE')
    if ratio is not None and not inr and not conf and not sup: unc.append('DEPTH')
    if left is None and ctx in ('no neighbors', 'short contig', '') and not conf and not sup and ratio is None: unc.append('NONB')
    if conf: st, rules = 'confirmed', conf
    elif short or edge: st, rules = 'uncertain', unc + sup
    elif sup: st, rules = 'supported', sup
    else:
        st = 'uncertain'; rules = unc or ['NO_SUPPORT']
    out.append(dict(locus_id=lid, lineage=o['lineage'], species=o['species'], assembly=o['assembly'], status=st,
                    rules_matched=';'.join(rules), locus_state=o['status'], contig_length_bp=L))

cols = ['locus_id', 'lineage', 'species', 'assembly', 'status', 'rules_matched']
with open(P + '/integration_status.tsv', 'w') as f:
    w = csv.DictWriter(f, cols, delimiter='\t', lineterminator='\n', extrasaction='ignore'); w.writeheader(); w.writerows(out)
cnt = Counter(); 
for r in out:
    cnt[(r['lineage'], r['status'], 'all')] += 1
    if r['locus_state'] == 'intact': cnt[(r['lineage'], r['status'], 'intact')] += 1
lins = sorted({r['lineage'] for r in out}) + ['TOTAL']
with open(P + '/integration_status_summary.tsv', 'w') as f:
    f.write('lineage\tn_loci\tconfirmed\tsupported\tuncertain\tintact_loci\tintact_confirmed\tintact_supported\tintact_uncertain\n')
    for l in lins:
        g = lambda st, k: sum(v for (a, s, kk), v in cnt.items() if (l == 'TOTAL' or a == l) and s == st and kk == k)
        n = sum(g(s, 'all') for s in ('confirmed', 'supported', 'uncertain')); ni = sum(g(s, 'intact') for s in ('confirmed', 'supported', 'uncertain'))
        f.write('\t'.join(map(str, [l, n] + [g(s, 'all') for s in ('confirmed', 'supported', 'uncertain')] + [ni] + [g(s, 'intact') for s in ('confirmed', 'supported', 'uncertain')])) + '\n')
print(open(P + '/integration_status_summary.tsv').read())
print(Counter(r['rules_matched'] for r in out).most_common(20))
