#!/usr/bin/env python
"""Composition (GC, GC3, RSCU) of T-DNA-type fungal loci vs host-genome background genes.
Outputs to analysis/composition/: per_locus.tsv, summary.tsv, background_genes.tsv.gz,
reference_bacterial.tsv, cds_loci.fna, failed.tsv, bg_sources.tsv, fig_gc3_locus_vs_bg.{png,svg}
"""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import csv, gzip, os, re, sys, json, itertools
from collections import defaultdict, Counter
import numpy as np

ROOT = f'{GMO_ROOT}'
OUT = f'{ROOT}/analysis/composition'
R2 = f'{ROOT}/analysis/loci_pipeline/round2'
LINEAGES = ['sus-like (T-DNA type)', 'cus-like (T-DNA type)', 'acc-type GDPD', 'pKerr14a-type acs-like']
MIN_BG = 30
MIN_CDS = 300          # nt, background genes (without stop)
POOL_CAP = 300         # genes per genome in pooled Boletales background
TE_RE = re.compile(r'transpos|retrotranspos|reverse transcriptase|gag-pol|gag[- ]polyprotein|\bpol polyprotein|integrase|'
                   r'copia|gypsy|\bLTR\b|mutator|helitron|tc1|hAT family|piggybac|endonuclease.*reverse|ty3|ty1|'
                   r'retroelement|retrovirus|\bRNase H\b', re.I)
rng = np.random.default_rng(12345)

# ---------------- genetic code
B = 'TCAG'
AA = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODONS = [a + b + c for a in B for b in B for c in B]
CODE = dict(zip(CODONS, AA))
SENSE = [c for c in CODONS if CODE[c] != '*']
CIDX = {c: i for i, c in enumerate(CODONS)}
FAM = defaultdict(list)
for c in SENSE:
    FAM[CODE[c]].append(CIDX[c])
RS_CODONS = [i for aa, l in FAM.items() if len(l) > 1 for i in l]   # exclude Met, Trp
RS_FAMS = [(aa, l) for aa, l in FAM.items() if len(l) > 1]
RS_POS = {c: j for j, c in enumerate(RS_CODONS)}
COMP = str.maketrans('ACGTacgtN', 'TGCAtgcaN')


def rc(s): return s.translate(COMP)[::-1]


def read_fasta(p):
    h, buf = None, []
    op = gzip.open if p.endswith('.gz') else open
    with op(p, 'rt') as f:
        for l in f:
            if l[0] == '>':
                if h is not None: yield h, ''.join(buf)
                h, buf = l[1:].strip(), []
            else:
                buf.append(l.strip())
    if h is not None: yield h, ''.join(buf)


def check_cds(s):
    """return (cds_without_stop, problem or None)"""
    s = s.upper()
    if len(s) % 3: return None, 'len_not_mult3'
    cod = [s[i:i + 3] for i in range(0, len(s), 3)]
    if any(c not in CODE for c in cod): return None, 'ambiguous_base'
    if CODE[cod[-1]] == '*': cod = cod[:-1]
    if any(CODE[c] == '*' for c in cod): return None, 'internal_stop'
    if not cod: return None, 'empty'
    return ''.join(cod), None


def counts(cds):
    v = np.zeros(64, dtype=np.int32)
    for i in range(0, len(cds), 3): v[CIDX[cds[i:i + 3]]] += 1
    return v


def gc_gc3(cds):
    n = len(cds); gc = sum(cds.count(x) for x in 'GC') / n
    t = cds[2::3]; gc3 = sum(t.count(x) for x in 'GC') / len(t)
    return gc, gc3


# ---------------- RSCU / distance
def rscu_matrix(cnt):
    """cnt: (n,64) -> (rscu (n,59), mask (n,59)); rscu undefined (mask 0) for amino acids absent from the gene."""
    cnt = np.atleast_2d(cnt).astype(float)
    n = cnt.shape[0]
    R = np.zeros((n, len(RS_CODONS))); M = np.zeros_like(R)
    for aa, l in RS_FAMS:
        tot = cnt[:, l].sum(1)
        ok = tot > 0
        for i in l:
            j = RS_POS[i]
            R[ok, j] = cnt[ok, i] / (tot[ok] / len(l))
            M[ok, j] = 1
    return R, M


def wpearson(X, M, y):
    """masked Pearson r of each row of X (n,k) vs y (n,k or k) using mask M (n,k)."""
    if y.ndim == 1: y = np.broadcast_to(y, X.shape)
    w = M.sum(1)
    mx = (X * M).sum(1) / w; my = (y * M).sum(1) / w
    dx = (X - mx[:, None]) * M; dy = (y - my[:, None]) * M
    den = np.sqrt((dx ** 2).sum(1) * (dy ** 2).sum(1))
    with np.errstate(invalid='ignore', divide='ignore'):
        return (dx * dy).sum(1) / den


def bg_rscu_vec(total_cnt):
    R, _ = rscu_matrix(total_cnt)
    return R[0]


def bg_distances(cnt):
    """1 - r of each bg gene vs leave-one-out bg RSCU."""
    tot = cnt.sum(0).astype(float)
    R, M = rscu_matrix(cnt)
    # leave-one-out RSCU for each gene (chunked)
    d = np.empty(cnt.shape[0])
    for a in range(0, cnt.shape[0], 2000):
        sl = slice(a, a + 2000)
        loo = tot[None, :] - cnt[sl]
        Rl, _ = rscu_matrix(loo)
        d[sl] = 1 - wpearson(R[sl], M[sl], Rl)
    return d


# ---------------- inputs
cat = list(csv.DictReader(open(f'{ROOT}/paper/tables/TableS1_loci_catalog.tsv'), delimiter='\t'))
sel = [x for x in cat if x['lineage_final'] in LINEAGES and x['status'] == 'intact']
lt = {x['locus_id']: x for x in csv.DictReader(open(f'{R2}/loci_table.tsv'), delimiter='\t')}

# registry loci (all, any status) by gca/contig for overlap exclusion
reg = defaultdict(list)    # (gca, contig) -> [(s,e)]
for x in cat:
    if not x['start'] or not x['end']: continue
    reg[(x['gca'], x['contig'])].append((int(x['start']), int(x['end'])))
for x in lt.values():
    if not x['contig_start'] or not x['contig_end']: continue
    reg[(x['gca'], x['contig'])].append((int(x['contig_start']), int(x['contig_end'])))


def overlaps_reg(gca, contig, s, e):
    contig = re.sub(r'\.\d+$', '', contig)
    return any(s <= b and a <= e for a, b in reg.get((gca, contig), []))


# ---------------- locus CDS
models = {}
for fam in ['SUS', 'CUS', 'ACS']:
    for h, s in read_fasta(f'{R2}/proteins/{fam}_models.fna'):
        models[h.split('|')[0]] = s
cds_loc, failed, cds_src = {}, [], {}
EXTRA = {'GCA_058411305.1': 'GCA_058411305.1_ASM5841130v1_genomic.fna', 'GCA_964275055.1': 'GCA_964275055.1_MOSE0_genomic.fna'}
extra_dir = f'{ROOT}/data/fungal_genomes/extra2/x/ncbi_dataset/data'
for x in sel:
    lid = x['locus_id']
    if lid.endswith('_A1') and x['gca'] in EXTRA:
        g = x['gca']; seqs = dict(read_fasta(f'{extra_dir}/{g}/{EXTRA[g]}'))
        key = [k for k in seqs if k.split()[0].startswith(x['contig'])][0]
        s = seqs[key][int(x['start']) - 1:int(x['end'])]
        s = rc(s) if x['strand'] == '-' else s
        src = 'extra2_genome_fasta(single exon, catalog coords)'
    elif lid in models:
        s, src = models[lid], 'proteins/<FAM>_models.fna'
    else:
        p = f'{R2}/loci/{lid}/target_cds.fna'
        if os.path.exists(p):
            s = next(read_fasta(p))[1]; src = 'loci/<id>/target_cds.fna'
        else:
            failed.append((lid, 'no_cds_found')); continue
    cds, prob = check_cds(s)
    if prob: failed.append((lid, prob)); continue
    cds_loc[lid] = cds; cds_src[lid] = src

with open(f'{OUT}/cds_loci.fna', 'w') as f:
    for k, v in cds_loc.items(): f.write(f'>{k}\n{v}\n')

# ---------------- background
genomes = sorted({x['gca'] for x in sel})
order_of = {x['gca']: x['order'] for x in sel}
org_of = {x['gca']: x['organism'] for x in sel}
bg = {}          # gca -> list of (gene_id, cds)
bg_src = {}
dl = f'{OUT}/cds_dl'
EXTRA_ANN = set(EXTRA)
by_gca_loci = defaultdict(list)
for k, x in lt.items(): by_gca_loci[x['gca']].append(k)

te_nb = set()
for r in csv.DictReader(open(f'{R2}/neighbors_class.tsv'), delimiter='\t'):
    if r['kind'] == 'neighbor' and r['is_TE'] == 'yes':
        te_nb.add((r['locus_id'], r['query'].split('|')[1]))


def bg_from_datasets(g):
    out, nte, nov, nbad = {}, 0, 0, 0
    for h, s in read_fasta(f'{dl}/{g}/ncbi_dataset/data/{g}/cds_from_genomic.fna'):
        m = re.match(r'lcl\|(\S+?)_cds_\S+', h)
        contig = m.group(1)
        loc = re.search(r'\[location=([^\]]+)\]', h).group(1)
        nums = [int(n) for n in re.findall(r'\d+', loc)]
        prod = (re.search(r'\[protein=([^\]]*)\]', h) or [0, ''])[1]
        lt_ = (re.search(r'\[locus_tag=([^\]]*)\]', h) or [0, h.split()[0]])[1]
        if TE_RE.search(prod): nte += 1; continue
        if overlaps_reg(g, contig, min(nums), max(nums)): nov += 1; continue
        if 'join' not in loc and False: pass
        cds, prob = check_cds(s)
        if prob or len(cds) < MIN_CDS: nbad += 1; continue
        if lt_ not in out or len(cds) > len(out[lt_]): out[lt_] = cds     # longest isoform per locus_tag
    return list(out.items()), f'NCBI datasets cds_from_genomic.fna; excluded TE-product={nte}, overlap registry loci={nov}, bad/short(<{MIN_CDS}nt)={nbad}'


def bg_from_gff_extra(g):
    gff = f'{extra_dir}/{g}/genomic.gff'
    seqs = {k.split()[0]: v for k, v in read_fasta(f'{extra_dir}/{g}/{EXTRA[g]}')}
    tr = defaultdict(list); info = {}
    for l in open(gff):
        if l[0] == '#': continue
        p = l.rstrip('\n').split('\t')
        if len(p) < 9 or p[2] != 'CDS' or 'pseudo=true' in p[8]: continue
        a = dict(kv.split('=', 1) for kv in p[8].split(';') if '=' in kv)
        k = a.get('Parent', a.get('ID'))
        tr[k].append((p[0], int(p[3]), int(p[4]), p[6], int(p[7]) if p[7] != '.' else 0))
        info[k] = a
    out, nte, nov, nbad = [], 0, 0, 0
    for k, fs in tr.items():
        prod = info[k].get('product', '')
        if TE_RE.search(prod): nte += 1; continue
        c = fs[0][0]; st = fs[0][3]
        s0, e0 = min(f[1] for f in fs), max(f[2] for f in fs)
        if overlaps_reg(g, c, s0, e0): nov += 1; continue
        fs = sorted(fs, key=lambda f: f[1], reverse=(st == '-'))
        seq = ''.join(rc(seqs[c][f[1] - 1:f[2]]) if st == '-' else seqs[c][f[1] - 1:f[2]] for f in fs)
        seq = seq[fs[0][4]:]
        cds, prob = check_cds(seq)
        if prob or len(cds) < MIN_CDS: nbad += 1; continue
        out.append((k, cds))
    return out, f'GFF CDS + genome FASTA (extra2); excluded TE-product={nte}, overlap registry loci={nov}, bad/short={nbad}'


def bg_from_augustus(g):
    out, seen, nte, nov, nbad, nloci = [], set(), 0, 0, 0, 0
    for lid in by_gca_loci[g]:
        d = f'{R2}/loci/{lid}'
        if not os.path.exists(f'{d}/augustus.gff3'): continue
        nloci += 1
        ctxs = dict(read_fasta(f'{d}/context.fa'))
        hdr = open(f'{d}/context.fa').readline()
        ls = int(re.search(r'local_start=(\d+)', hdr).group(1))
        seq = next(iter(ctxs.values())).upper()
        contig = lt[lid]['contig']
        genes = defaultdict(list); flags = defaultdict(set)
        for l in open(f'{d}/augustus.gff3'):
            if l[0] == '#': continue
            p = l.rstrip('\n').split('\t')
            if len(p) < 9: continue
            m = re.search(r'Parent=([^;]+)', p[8])
            if p[2] == 'CDS': genes[m.group(1)].append((int(p[3]), int(p[4]), p[6], int(p[7])))
            elif p[2] in ('start_codon', 'stop_codon'): flags[m.group(1)].add(p[2])
        for t, fs in genes.items():
            if flags[t] != {'start_codon', 'stop_codon'}: continue
            gid = t.split('.')[0]
            st = fs[0][2]
            s0, e0 = min(f[0] for f in fs) + ls - 1, max(f[1] for f in fs) + ls - 1
            key = (contig, s0, e0, st)
            if key in seen: continue
            seen.add(key)
            if (lid, gid) in te_nb: nte += 1; continue
            if overlaps_reg(g, contig, s0, e0): nov += 1; continue
            fs = sorted(fs, key=lambda f: f[0], reverse=(st == '-'))
            sq = ''.join(rc(seq[f[0] - 1:f[1]]) if st == '-' else seq[f[0] - 1:f[1]] for f in fs)
            sq = sq[fs[0][3]:]
            cds, prob = check_cds(sq)
            if prob or len(cds) < MIN_CDS: nbad += 1; continue
            out.append((f'{contig}:{s0}-{e0}{st}', cds))
    return out, f'AUGUSTUS neighbor genes of {nloci} locus dirs (complete genes: start+stop codon); excluded TE-flagged(diamond)={nte}, overlap registry loci={nov}, bad/short={nbad}'


for g in genomes:
    if g in EXTRA_ANN:
        bg[g], bg_src[g] = bg_from_gff_extra(g); bg_src[g] = ('annotation_GFF_extra2', bg_src[g])
    elif os.path.exists(f'{dl}/{g}/ncbi_dataset/data/{g}/cds_from_genomic.fna'):
        bg[g], s = bg_from_datasets(g); bg_src[g] = ('ncbi_annotation_CDS', s)
    else:
        bg[g], s = bg_from_augustus(g); bg_src[g] = ('augustus_neighbors', s)
    print(g, org_of[g], bg_src[g][0], len(bg[g]), flush=True)

with open(f'{OUT}/bg_sources.tsv', 'w') as f:
    f.write('gca\torganism\torder\tbg_source\tn_bg_genes\tnotes\n')
    for g in genomes: f.write(f'{g}\t{org_of[g]}\t{order_of[g]}\t{bg_src[g][0]}\t{len(bg[g])}\t{bg_src[g][1]}\n')

# ---------------- background stats
def prep(genes):
    cnt = np.array([counts(c) for _, c in genes])
    gcs = np.array([gc_gc3(c) for _, c in genes])
    L = np.array([len(c) // 3 for _, c in genes])
    return dict(cnt=cnt, gc=gcs[:, 0], gc3=gcs[:, 1], L=L, tot=cnt.sum(0).astype(float))


def pooled(gs):
    sub = []
    for g in gs:
        gl = bg[g]
        if len(gl) > POOL_CAP:
            idx = rng.choice(len(gl), POOL_CAP, replace=False); gl = [gl[i] for i in sorted(idx)]
        sub += gl
    return sub


bol = [g for g in genomes if order_of[g] == 'Boletales']
pool_genes = pooled(bol)
S = {g: prep(bg[g]) for g in genomes if len(bg[g]) >= 5}
S_pool = prep(pool_genes)
for st in list(S.values()) + [S_pool]:
    st['d'] = bg_distances(st['cnt'])
    st['rs'] = bg_rscu_vec(st['tot'][None, :])


def eval_locus(cds, st):
    gc, gc3 = gc_gc3(cds)
    cnt = counts(cds)[None, :]
    R, M = rscu_matrix(cnt)
    d = float(1 - wpearson(R, M, st['rs'])[0])
    L = len(cds) // 3
    dl_ = np.abs(st['L'] - L) / L
    idx = np.where(dl_ <= 0.25)[0]
    if len(idx) < 100: idx = np.argsort(dl_)[:min(100, len(dl_))]
    ref = st['d'][idx]
    pct_codon = 100 * (np.sum(ref < d) + 0.5 * np.sum(ref == d)) / len(ref)
    n = len(st['gc3'])
    pct3 = 100 * (np.sum(st['gc3'] < gc3) + 0.5 * np.sum(st['gc3'] == gc3)) / n
    pctgc = 100 * (np.sum(st['gc'] < gc) + 0.5 * np.sum(st['gc'] == gc)) / n
    sd = st['gc3'].std(ddof=1) if n > 1 else np.nan
    return dict(med_gc=np.median(st['gc']), q1_gc=np.percentile(st['gc'], 25), q3_gc=np.percentile(st['gc'], 75),
                med_gc3=np.median(st['gc3']), q1_gc3=np.percentile(st['gc3'], 25), q3_gc3=np.percentile(st['gc3'], 75),
                gc3_pct=pct3, gc_pct=pctgc, gc3_z=(gc3 - st['gc3'].mean()) / sd, codon_dist=d, codon_dist_pct=pct_codon,
                n_len_matched=len(idx))


# ---------------- per locus
fmt = lambda v: ('%.4f' % v) if isinstance(v, (float, np.floating)) else str(v)
rows = []
cols = ['locus_id', 'family', 'lineage', 'organism', 'order', 'gca', 'cds_len_nt', 'cds_source', 'gc', 'gc3',
        'bg_source', 'n_bg_genes', 'bg_small_flag', 'bg_med_gc', 'bg_q1_gc', 'bg_q3_gc', 'bg_med_gc3', 'bg_q1_gc3', 'bg_q3_gc3',
        'delta_gc3', 'gc3_pctile', 'gc3_z', 'gc_pctile', 'codon_dist_1minusr_rscu', 'codon_dist_pctile', 'n_len_matched_bg',
        'bolpool_n_genes', 'bolpool_med_gc3', 'bolpool_delta_gc3', 'bolpool_gc3_pctile', 'bolpool_gc3_z',
        'bolpool_codon_dist_pctile']
for x in sel:
    lid = x['locus_id']; g = x['gca']
    if lid not in cds_loc: continue
    cds = cds_loc[lid]; gc, gc3 = gc_gc3(cds)
    r = dict(locus_id=lid, family=x['family'], lineage=x['lineage_final'], organism=x['organism'], order=x['order'], gca=g,
             cds_len_nt=len(cds) + 3, cds_source=cds_src[lid], gc=gc, gc3=gc3, bg_source=bg_src[g][0], n_bg_genes=len(bg[g]))
    if g in S and len(bg[g]) >= 5:
        e = eval_locus(cds, S[g])
        r.update(bg_small_flag='bg_small(<30)' if len(bg[g]) < MIN_BG else 'ok', bg_med_gc=e['med_gc'], bg_q1_gc=e['q1_gc'],
                 bg_q3_gc=e['q3_gc'], bg_med_gc3=e['med_gc3'], bg_q1_gc3=e['q1_gc3'], bg_q3_gc3=e['q3_gc3'],
                 delta_gc3=gc3 - e['med_gc3'], gc3_pctile=e['gc3_pct'], gc3_z=e['gc3_z'], gc_pctile=e['gc_pct'],
                 codon_dist_1minusr_rscu=e['codon_dist'], codon_dist_pctile=e['codon_dist_pct'], n_len_matched_bg=e['n_len_matched'])
    else:
        r['bg_small_flag'] = 'bg_unavailable(<5)'
    if x['order'] == 'Boletales':
        e = eval_locus(cds, S_pool)
        r.update(bolpool_n_genes=len(pool_genes), bolpool_med_gc3=e['med_gc3'], bolpool_delta_gc3=gc3 - e['med_gc3'],
                 bolpool_gc3_pctile=e['gc3_pct'], bolpool_gc3_z=e['gc3_z'], bolpool_codon_dist_pctile=e['codon_dist_pct'])
    rows.append(r)
with open(f'{OUT}/per_locus.tsv', 'w') as f:
    f.write('\t'.join(cols) + '\n')
    for r in rows: f.write('\t'.join(fmt(r.get(c, '')) for c in cols) + '\n')
with open(f'{OUT}/failed.tsv', 'w') as f:
    f.write('locus_id\treason\n')
    for a, b in failed: f.write(f'{a}\t{b}\n')
with gzip.open(f'{OUT}/background_genes.tsv.gz', 'wt') as f:
    f.write('gca\tgene\tcds_len_nt\tgc\tgc3\n')
    for g in genomes:
        for gid, c in bg[g]:
            a, b = gc_gc3(c); f.write(f'{g}\t{gid}\t{len(c)}\t{a:.4f}\t{b:.4f}\n')

# ---------------- summary
def summ(label, rs, mode):
    rs = [r for r in rs]
    if mode == 'own':
        rs = [r for r in rs if r['bg_small_flag'] == 'ok']
        k = dict(d='delta_gc3', p='gc3_pctile', c='codon_dist_pctile')
    else:
        rs = [r for r in rs if r['order'] == 'Boletales']
        k = dict(d='bolpool_delta_gc3', p='bolpool_gc3_pctile', c='bolpool_codon_dist_pctile')
    n = len(rs)
    if n == 0: return [label, mode, 0] + [''] * 7
    d = np.array([r[k['d']] for r in rs]); p = np.array([r[k['p']] for r in rs]); c = np.array([r[k['c']] for r in rs])
    return [label, mode, n, '%.4f' % np.median(d), int((p > 95).sum()), int((p < 5).sum()), '%.1f' % np.median(p),
            int((c > 95).sum()), int((c < 5).sum()), '%.1f' % np.median(c)]


groups = [('all', lambda r: True)]
groups += [(l, (lambda l: lambda r: r['lineage'] == l)(l)) for l in LINEAGES]
groups += [('Boletales', lambda r: r['order'] == 'Boletales'), ('non-Boletales', lambda r: r['order'] != 'Boletales')]
groups += [(f'{l} | Boletales', (lambda l: lambda r: r['lineage'] == l and r['order'] == 'Boletales')(l)) for l in LINEAGES]
groups += [(f'{l} | non-Boletales', (lambda l: lambda r: r['lineage'] == l and r['order'] != 'Boletales')(l)) for l in LINEAGES]
with open(f'{OUT}/summary.tsv', 'w') as f:
    f.write('group\tbackground\tn_loci\tmedian_delta_gc3\tn_gc3_above_p95\tn_gc3_below_p5\tmedian_gc3_pctile\t'
            'n_codon_dist_above_p95\tn_codon_dist_below_p5\tmedian_codon_dist_pctile\n')
    for lab, fn in groups:
        rs = [r for r in rows if fn(r)]
        for mode in ('own', 'bolpool'):
            f.write('\t'.join(map(str, summ(lab, rs, mode))) + '\n')
print('done; loci', len(rows), 'failed', len(failed))

# ---------------- bacterial references
from Bio import SeqIO
want = {'AAZ50401.1': 'sus pTiBo542', 'AAK90959.2': 'acs pTiC58'}
pl = f'{ROOT}/analysis/refset/plasmids'
with open(f'{OUT}/reference_bacterial.tsv', 'w') as f:
    f.write('ref_id\tdescription\tsource_file\tcds_len_nt\tgc\tgc3\tnote\n')
    for fn in sorted(os.listdir(pl)):
        if not fn.endswith('.gb'): continue
        for rec in SeqIO.parse(f'{pl}/{fn}', 'genbank'):
            for ft in rec.features:
                if ft.type != 'CDS': continue
                pid = ft.qualifiers.get('protein_id', [''])[0]
                prod = ft.qualifiers.get('product', [''])[0]
                tag = None
                if pid in want: tag = want[pid]
                elif fn == 'NC_002575.1.gb' and 'cucumopine synthase-related' in prod: tag = 'mis-like (pRi1724) ' + ft.qualifiers.get('locus_tag', [''])[0]
                if tag is None: continue
                s = str(ft.extract(rec.seq))
                cds, prob = check_cds(s)
                if prob: f.write(f'{pid}\t{tag}\t{fn}\t{len(s)}\t\t\t{prob}\n'); continue
                a, b = gc_gc3(cds)
                f.write(f'{pid}\t{tag}\t{fn}\t{len(s)}\t{a:.4f}\t{b:.4f}\t{prod}\n')
                print(pid, tag, fn, len(s), a, b)
