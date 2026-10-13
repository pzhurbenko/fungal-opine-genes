#!/usr/bin/env python3
"""Per-locus integration status (confirmed / supported / uncertain) for the 254 placed loci of Table S1a.
Inputs (final-numbered supplementary TSVs, paper/supplementary/): TableS1a, TableS8 (read depth), TableS9 (independent assemblies),
TableS12 (genes on the island contigs). Output: integration_status.tsv, integration_status_summary.tsv (this directory).
Run (env PlantCad2): python integration_status.py
Rules (a locus gets the first matching class). 2026-10-13ak, review v14 point 1: F2/R2/CHR redefined (see below).
 confirmed  = any of
   F2  fungal-classified gene on BOTH sides of the locus, computed per side from per-neighbour classifications
       (analysis/loci_pipeline/round2/neighbors_class.tsv majority_kingdom == Fungi, top_org genus different from the carrier species' genus (self-only hits do not count) and n_hits >= 1; analysis/islands/islands_genes.tsv role == Fungi, which
       overrides neighbors_class when both cover the same gene, so self_only / bacterial_like / TARGET genes do not count;
       analysis/fusarium_odh_island/reverse_search_summary.tsv for the flanking genes of the F. poae Fp133 odh island, positions from
       fig4/fp133_window_genes.tsv). Contig names are matched without the '.1' version suffix. Genes overlapping the locus or any other
       registered locus of the contig (S1a, S12 focal rows, islands_genes TARGET/registry rows) are ignored. S1a genomic_context is NOT
       evidence of two flanks any more. chs-like loci have no per-neighbour classification and get no F2.
   R   reproduced with flanks in an independent assembly (S9: locus covered >=0.5 and both flanks covered >=0.5)
   D2  read depth ratio (S8 locus_depth_ratio) 0.5-2 AND junction reads (S8) > 0 on both flanks
   CHR (any assembly level, incl. 'Complete Genome'; the criterion is the sequence role) locus on an assembled-molecule sequence with Assigned-Molecule-Location/Type == Chromosome in the NCBI assembly report
       (analysis/data_checks_20261013/assembly_reports/; S13a/S13b sequences), >= 1 kb from both contig ends, not SHORT/EDGE.
       Unlocalized / unplaced scaffolds and sequences without a downloaded report are not verified and get no CHR (they may get LONG).
 uncertain (first override; applies if not confirmed) =
   SHORT  contig shorter than locus length + 2 kb (1 kb on either side)
   EDGE   locus touches a contig end (start <= 1 or end >= contig length; S8 'edge' flanks)
 supported = any of
   R2   a gene of the same lineage found in >=2 assemblies of the species (S13a/S13b: species, assembly), not SHORT/EDGE;
        the identity of the insertion (same flanks) is NOT verified
   LONG contig >= 100 kb, locus not at a contig end, genomic_context not bacterial, and depth ratio not outside 0.5-2 (a fungal-sized contig is itself evidence)
   F1  fungal-classified gene on one side only (per-side data); if no per-side data exist, S1a genomic_context 'fungal context' or 'mixed'
       counts as F1 (at least one fungal neighbour)
   D1  depth ratio 0.5-2 without junction reads on both sides
 uncertain (remaining) =
   DEPTH  depth ratio < 0.5 or > 2 without other support (no F/R/D rule met)
   NONB   no neighbour information and no other support
"""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import csv, os, re, glob
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

# 2026-10-13ak, review v14 point 1: per-side fungal neighbours from per-neighbour classification tables (not from S1a genomic_context)
G = P + '/../../analysis/'
nrm = lambda c: re.sub(r'\.\d+$', '', c)
zones = defaultdict(list)                      # normalised contig -> [(start, end)] of registered loci (ignored as neighbours)
for o in s1a: zones[nrm(o['contig'])].append((int(o['start']), int(o['end'])))
for r in s12:
    if r['gene_class'] == 'focal locus': zones[nrm(r['contig'])].append((int(r['start']), int(r['end'])))
genes = defaultdict(dict)                      # normalised contig -> {(start, end): set of classes 'F'/'N'/'B'}
def addg(contig, s, e, cls): genes[nrm(contig)].setdefault((s, e), set()).add(cls)
carrier = {nrm(o['contig']): re.sub(r' \(cluster [^)]*\)', '', o['species']).split()[0].lower() for o in s1a if o['species'].split()}
for r in csv.DictReader(open(G + 'loci_pipeline/round2/neighbors_class.tsv'), delimiter='\t'):
    if r['kind'] != 'neighbor': continue
    q = r['query'].split('|'); m = re.match(r'^(\d+)-(\d+)[+-]$', q[2]) if len(q) > 2 else None
    if not m: continue
    contig = re.sub(r'_(S\d)?[LA]\d+$', '', r['locus_id'].split('_', 1)[1])
    # 2026-10-13ak follow-up: a neighbour is fungal only if its top hit is from another genus than the carrier species (self-only hits do not count) and n_hits >= 1
    cls = {'Fungi': 'F', 'Bacteria': 'B'}.get(r['majority_kingdom'], 'N')
    if cls == 'F' and (int(r['n_hits'] or 0) < 1 or (carrier.get(nrm(contig)) and r['top_org'].split()[:1] and r['top_org'].split()[0].lower() == carrier[nrm(contig)])): cls = 'N'
    addg(contig, int(m.group(1)), int(m.group(2)), cls)
igo = defaultdict(dict)
for r in csv.DictReader(open(G + 'islands/islands_genes.tsv'), delimiter='\t'):
    s_, e_ = int(r['start']), int(r['end'])
    if r['role'] == 'TARGET' or r['registry_locus']: zones[nrm(r['contig'])].append((s_, e_)); continue
    igo[nrm(r['contig'])][(s_, e_)] = {'Fungi': 'F', 'bacterial_like': 'B'}.get(r['role'], 'N')
for c, d in igo.items():
    for k, v in d.items(): genes[c][k] = {v}   # islands_genes (self_only aware) overrides neighbors_class for the same gene
# Fp133 left flank (KAL3595991) has a 'synthetic construct' best bacterial hit (87.9% over 58 aa) but 47/50 top hits are fungal: kept as fungal
# odh island of F. poae Fp133 (= GCA_046056315.1, CM100256.1): classified flanking genes
win = {r['protein']: r for r in csv.DictReader(open(G + 'fusarium_odh_island/fig4/fp133_window_genes.tsv'), delimiter='\t')}
for r in csv.DictReader(open(G + 'fusarium_odh_island/reverse_search_summary.tsv'), delimiter='\t'):
    if 'flanking' in r['gene'] and r['query'] in win:
        w = win[r['query']]
        addg(w['contig'], int(w['start']), int(w['end']), 'F' if int(r['top50_non_Fusarium_fungal']) > int(r['top50_non_Fusarium_other']) else 'N')
def sides(o):
    c = nrm(o['contig']); a, b = int(o['start']), int(o['end'])
    gl = [(k, v) for k, v in genes.get(c, {}).items() if not any(k[0] <= ze and k[1] >= zs for zs, ze in zones[c])]
    if not gl: return None
    fung = [k for k, v in gl if 'F' in v and 'B' not in v]
    return any(k[1] < a for k in fung), any(k[0] > b for k in fung)

# CHR: assembled-molecule Chromosome sequences from the NCBI assembly reports (six chromosome-level Fusarium assemblies)
chrom = set()
for fn in glob.glob(G + 'data_checks_20261013/assembly_reports/*_assembly_report.txt'):
    for l in open(fn):
        f = l.rstrip('\n').split('\t')
        if not l.startswith('#') and len(f) > 4 and f[1] == 'assembled-molecule' and f[3] == 'Chromosome': chrom.add(nrm(f[4]))

out = []
for o in s1a:
    lid = o['locus_id']; a, b, L = int(o['start']), int(o['end']), int(o['contig_length_bp'])
    d = s8.get(lid); ratio = fl(d['locus_depth_ratio']) if d else None
    jl = jr = None
    if d:
        jl = d['reads_spanning_left_flank']; jr = d['reads_spanning_right_flank']
        jl = None if jl == 'edge' else fl(jl); jr = None if jr == 'edge' else fl(jr)
    ctx = o['genomic_context']
    sd = sides(o)
    if sd is not None: left, right = sd
    elif ctx in ('fungal context', 'mixed'): left, right = True, False   # no per-side data: at least one fungal neighbour = F1 only
    else: left = right = None
    rules = []; conf = []; sup = []; unc = []
    if left and right: conf.append('F2')
    if lid in repro: conf.append('R')
    inr = ratio is not None and .5 <= ratio <= 2
    if inr and jl and jr and jl > 0 and jr > 0: conf.append('D2')
    short = L < (b - a + 1) + 2000
    edge = a <= 1 or b >= L or (d is not None and 'edge' in (d['reads_spanning_left_flank'], d['reads_spanning_right_flank']))
    if lid in lvl and not short and not edge:
        al, sp, asm, g = lvl[lid]
        if nrm(o['contig']) in chrom and a > 1000 and L - b > 1000: conf.append('CHR')
        if len(asm_by[(gk(lid, g), sp)]) >= 2: sup.append('R2')
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
