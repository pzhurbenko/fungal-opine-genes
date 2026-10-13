#!/usr/bin/env python3
# usage: count_expr.py RUN SPECIES(fol|van) BAM GFF TARGET_LOCUS CTRL_TSV(locus<TAB>label)
# per-gene read counts (exon overlap, MAPQ>=1, primary), TPM-like, percentile, neighbours +-10, controls, coverage evenness, strand
import sys, re, bisect, json, statistics
import numpy as np, pysam
run, sp, bam, gff, target, ctrl = sys.argv[1:7]
genes = {}; rna2g = {}
for l in open(gff):
    if l[0] == '#': continue
    f = l.rstrip('\n').split('\t')
    if len(f) < 9: continue
    a = dict(x.split('=', 1) for x in f[8].split(';') if '=' in x)
    if f[2] == 'gene':
        genes[a['ID']] = dict(id=a['ID'], lt=a.get('locus_tag', a.get('Name', a['ID'])), chr=f[0], s=int(f[3]), e=int(f[4]), st=f[6], ex=[], mrna=False, bt=a.get('gene_biotype', ''))
    elif f[2] == 'mRNA':
        rna2g[a['ID']] = a['Parent']; genes[a['Parent']]['mrna'] = True
    elif f[2] == 'exon':
        g = rna2g.get(a['Parent'])
        if g: genes[g]['ex'].append((int(f[3]), int(f[4])))
G = [g for g in genes.values() if g['mrna'] and g['ex']]
for g in G:  # union of exons
    ex = sorted(set(g['ex'])); m = []
    for s, e in ex:
        if m and s <= m[-1][1] + 1: m[-1][1] = max(m[-1][1], e)
        else: m.append([s, e])
    g['ex'] = [tuple(x) for x in m]; g['len'] = sum(e - s + 1 for s, e in m)
bf = pysam.AlignmentFile(bam, 'rb')
tot_prim = 0; tot_mapped = 0
for st in bf.get_index_statistics(): pass
n_all = int(pysam.view('-c', '-F', '0x900', bam).strip()); n_map = int(pysam.view('-c', '-F', '0x904', bam).strip())
for g in G:
    c = sg = ag = 0
    for r in bf.fetch(g['chr'], g['s'] - 1, g['e']):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_qcfail or r.mapping_quality < 1: continue
        hit = False
        for bs, be in r.get_blocks():  # 0-based half-open
            for s, e in g['ex']:
                if bs < e and be > s - 1: hit = True; break
            if hit: break
        if not hit: continue
        c += 1
        sense = (r.is_reverse == (g['st'] == '-'))
        if r.is_read2: sense = not sense
        if sense: sg += 1
        else: ag += 1
    g['n'] = c; g['sense'] = sg; g['anti'] = ag
tot = sum(g['n'] for g in G)
for g in G:
    g['rpk'] = 1000 * g['n'] / g['len']
rs = sum(g['rpk'] for g in G)
for g in G: g['tpm'] = 1e6 * g['rpk'] / rs if rs else 0.0
tp = sorted(g['tpm'] for g in G)
def pct(x): return (100 * bisect.bisect_left(tp, x) / len(tp), 100 * bisect.bisect_right(tp, x) / len(tp))
expr = [g['tpm'] for g in G if g['n'] > 0]
# library strandedness from genes with >=200 reads
big = [g for g in G if g['n'] >= 200]
ss = sum(g['sense'] for g in big); aa = sum(g['anti'] for g in big)
strand_frac = ss / (ss + aa) if ss + aa else float('nan')
T = next(g for g in G if g['lt'] == target)
bychr = sorted((g for g in G if g['chr'] == T['chr']), key=lambda g: g['s'])
i = bychr.index(T)
left = bychr[max(0, i - 10):i]; right = bychr[i + 1:i + 11]
overl = [g['lt'] for g in genes.values() if g['chr'] == T['chr'] and g['id'] != T['id'] and g['s'] <= T['e'] and g['e'] >= T['s']]
# coverage evenness over exons of target (10 bins of concatenated exon positions)
pos = [p for s, e in T['ex'] for p in range(s - 1, e)]
idx = {p: k for k, p in enumerate(pos)}
cov = np.zeros(len(pos), dtype=np.int32)
for r in bf.fetch(T['chr'], T['s'] - 1, T['e']):
    if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < 1: continue
    for bs, be in r.get_blocks():
        for p in range(max(bs, T['s'] - 1), min(be, T['e'])):
            k = idx.get(p)
            if k is not None: cov[k] += 1
bins = np.array_split(cov, 10)
binmean = [float(b.mean()) for b in bins]
cs = sum(binmean)
bin_share = [round(b / cs, 3) if cs else 0 for b in binmean]
def row(g): return dict(lt=g['lt'], chr=g['chr'], s=g['s'], e=g['e'], st=g['st'], len=g['len'], n=g['n'], tpm=round(g['tpm'], 3), sense=g['sense'], anti=g['anti'], pct_lt=round(pct(g['tpm'])[0], 1), pct_le=round(pct(g['tpm'])[1], 1))
ctrls = []
for l in open(ctrl):
    if l.startswith('#') or not l.strip(): continue
    lt, lab = l.rstrip('\n').split('\t')[:2]
    g = next((x for x in G if x['lt'] == lt), None)
    if g: d = row(g); d['label'] = lab; ctrls.append(d)
nb = left + right
out = dict(run=run, sp=sp, n_primary_reads=n_all, n_mapped_primary=n_map, mapped_pct=round(100 * n_map / n_all, 2), n_genes=len(G), reads_on_genes=tot,
           frac_reads_in_genes=round(tot / n_map, 3) if n_map else 0, genes_with_reads=len(expr), median_tpm_expressed=round(statistics.median(expr), 2) if expr else 0,
           median_tpm_all=round(statistics.median(tp), 3), strand_sense_frac_genes_ge200=round(strand_frac, 3), target=row(T), target_overlapping_genes=overl,
           target_breadth=round(float((cov > 0).mean()), 3), target_mean_depth=round(float(cov.mean()), 2), target_bin_share=bin_share,
           neighbours=[row(g) for g in nb], nb_tpm_min=round(min(g['tpm'] for g in nb), 2), nb_tpm_max=round(max(g['tpm'] for g in nb), 2),
           nb_tpm_median=round(statistics.median(g['tpm'] for g in nb), 2), nb_zero=sum(1 for g in nb if g['n'] == 0), controls=ctrls)
json.dump(out, open(f'results_{run}_{target}.json', 'w'), indent=1)
with open(f'gene_counts_{run}_{target}.tsv', 'w') as fh:
    fh.write('locus_tag\tchr\tstart\tend\tstrand\texon_len\treads\tsense\tantisense\ttpm\n')
    for g in G: fh.write(f"{g['lt']}\t{g['chr']}\t{g['s']}\t{g['e']}\t{g['st']}\t{g['len']}\t{g['n']}\t{g['sense']}\t{g['anti']}\t{g['tpm']:.4f}\n")
print(run, 'mapped%', out['mapped_pct'], 'target n', T['n'], 'tpm', round(T['tpm'], 3), 'pct', pct(T['tpm']))
