#!/usr/bin/env python3
"""Step 1-2: collect contigs for the backbone-threshold sensitivity test.

(a) bacterial: the 11 proteins with fate == 'probable contamination' (paper/tables/TableS3b_fate.tsv);
    the contig carrying each protein is found in the NCBI GFF3 of its assembly (downloaded with `datasets`)
    and extracted whole. Assembly-wide GC and N50 and the number of annotated genes on the contig are recorded.
(b) fungal: intact loci of paper/tables/TableS1_all_lineages.tsv, stratified by lineage: up to K loci per lineage
    drawn at random (seed=1) among loci whose whole contig was already fetched by the loci pipeline
    (analysis/loci_pipeline/contigs/<contig>.fa); one locus per contig.
Outputs (work/): contigs.fna (ids '<sample>__<contig>'), manifest.tsv
Usage: python -I 01_collect_contigs.py   (paths are hard-coded below)
"""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import csv, glob, gzip, os, random, re
from collections import defaultdict

GMO = f'{GMO_ROOT}'
OUT = f'{GMO}/analysis/thresholds/backbone/work'
DL = f'{GMO}/analysis/thresholds/backbone/downloads/contam9/ncbi_dataset/data'
S3B = f'{GMO}/paper/tables/TableS3b_fate.tsv'
S1 = f'{GMO}/paper/tables/TableS1_all_lineages.tsv'
CONTIGS = f'{GMO}/analysis/loci_pipeline/contigs'
K_PER_LINEAGE = 5
SEED = 1


def read_fasta(path):
    name, seq = None, []
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as fh:
        for l in fh:
            if l.startswith('>'):
                if name: yield name, ''.join(seq)
                name, seq = l[1:].split()[0], []
            else:
                seq.append(l.strip())
    if name: yield name, ''.join(seq)


def gc(s):
    s = s.upper(); g = s.count('G') + s.count('C'); at = s.count('A') + s.count('T')
    return 100.0 * g / (g + at) if g + at else float('nan')


def n50(lens):
    lens = sorted(lens, reverse=True); half = sum(lens) / 2; c = 0
    for x in lens:
        c += x
        if c >= half: return x


os.makedirs(OUT, exist_ok=True)
man = []; seqs = {}

# ---------- (a) bacterial ----------
rows = [r for r in csv.DictReader(open(S3B), delimiter='\t') if r['fate'] == 'probable contamination']
by_gca = defaultdict(list)
for r in rows: by_gca[r['gca']].append(r)
for gca, rl in by_gca.items():
    d = f'{DL}/{gca}'
    fna = glob.glob(f'{d}/*_genomic.fna')[0]
    prot2contig = {}; genes_on = defaultdict(set)
    for l in open(f'{d}/genomic.gff'):
        if l.startswith('#'): continue
        f = l.rstrip('\n').split('\t')
        if len(f) < 9: continue
        if f[2] == 'gene': genes_on[f[0]].add(f[8].split(';')[0])
        if f[2] == 'CDS':
            m = re.search(r'protein_id=([^;]+)', f[8])
            if m: prot2contig[m.group(1)] = f[0]
    want = {prot2contig[r['protein']] for r in rl}
    lens, gc_num, gc_den = [], 0, 0
    for name, s in read_fasta(fna):
        lens.append(len(s)); su = s.upper()
        gc_num += su.count('G') + su.count('C'); gc_den += su.count('G') + su.count('C') + su.count('A') + su.count('T')
        if name in want: seqs[f'bacterial__{name}'] = s
    asm_gc = 100.0 * gc_num / gc_den; asm_n50 = n50(lens)
    for r in rl:
        c = prot2contig[r['protein']]
        man.append(dict(sample='bacterial', contig=c, source=f"{r['protein']}|{r['gene']}|{r['o_bact_org']} {r['o_bact_pid']}%",
                        lineage='probable contamination', gca=gca, organism=r['organism'],
                        annot_genes_on_contig=len(genes_on[c]), assembly_gc=f'{asm_gc:.2f}', assembly_n50=asm_n50))

# ---------- (b) fungal ----------
rows = [r for r in csv.DictReader(open(S1), delimiter='\t')
        if r['status'] == 'intact' and r['lineage'] != 'not a homolog' and os.path.exists(f"{CONTIGS}/{r['contig']}.fa")]
by_lin = defaultdict(list)
for r in sorted(rows, key=lambda r: r['locus_id']): by_lin[r['lineage']].append(r)
rng = random.Random(SEED); used = set()
for lin in sorted(by_lin):
    cand = by_lin[lin][:]; rng.shuffle(cand); k = 0
    for r in cand:
        if k >= K_PER_LINEAGE: break
        if r['contig'] in used: continue
        used.add(r['contig']); k += 1
        s = next(read_fasta(f"{CONTIGS}/{r['contig']}.fa"))[1]
        seqs[f"fungal__{r['contig']}"] = s
        man.append(dict(sample='fungal', contig=r['contig'], source=r['locus_id'], lineage=lin, gca=r['gca'],
                        organism=r['organism'], annot_genes_on_contig='', assembly_gc='', assembly_n50=''))

with open(f'{OUT}/contigs.fna', 'w') as fo:
    for k, s in seqs.items():
        fo.write(f'>{k}\n')
        for i in range(0, len(s), 80): fo.write(s[i:i + 80] + '\n')
for m in man:
    s = seqs[f"{m['sample']}__{m['contig']}"]; m['length'] = len(s); m['GC'] = f'{gc(s):.2f}'
cols = ['sample', 'contig', 'source', 'lineage', 'gca', 'organism', 'length', 'GC', 'assembly_gc', 'assembly_n50', 'annot_genes_on_contig']
with open(f'{OUT}/manifest.tsv', 'w') as fo:
    w = csv.DictWriter(fo, cols, delimiter='\t', extrasaction='ignore'); w.writeheader(); w.writerows(man)
print('bacterial rows', sum(m['sample'] == 'bacterial' for m in man), 'fungal rows', sum(m['sample'] == 'fungal' for m in man),
      'unique contigs', len(seqs))
