#!/usr/bin/env python3
"""Empty-site analysis. Usage: step6.py <Q_gca> <locus_id> <outprefix> <T_gca> [<T_gca> ...]   (run in annot env, threads=16)"""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import sys, os, subprocess, collections, csv, re
S = f'{GMO_ROOT}/analysis/chs_fusarium/stage2'
P = f'{GMO_ROOT}/data/fungal_proteomes/fungi_annotated/ncbi_dataset/data'
G = f'{GMO_ROOT}/data/fungal_genomes/fusarium/fus/ncbi_dataset/data'
MOB = re.compile(r'transpos|integrase|retro|helicase|reverse transcriptase', re.I)
Q, LOC, OUT, TS = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
W = f'{S}/step6/{OUT}'; os.makedirs(W, exist_ok=True)
def run(c): subprocess.run(c, shell=True, check=True, executable='/bin/bash')
def faa(path):
    d = {}; h = None
    for l in open(path):
        if l[0] == '>': h = l[1:].split()[0]; d[h] = []
        else: d[h].append(l.strip())
    return {k: ''.join(v) for k, v in d.items()}
def genes(gca):
    g = {}; m = {}; prot = {}
    for l in open(f'{G}/{gca}/genomic.gff'):
        if l[0] == '#': continue
        f = l.rstrip('\n').split('\t')
        if len(f) < 9: continue
        a = dict(kv.split('=', 1) for kv in f[8].split(';') if '=' in kv)
        if f[2] in ('gene', 'pseudogene'): g[a['ID']] = dict(contig=f[0], start=int(f[3]), end=int(f[4]), strand=f[6], name=a.get('Name', a['ID']), prod='', ex=0, prot='')
        elif f[2] in ('mRNA', 'transcript'): m[a['ID']] = a.get('Parent')
        elif f[2] == 'CDS':
            gid = m.get(a.get('Parent'), a.get('Parent'))
            if gid in g:
                G_ = g[gid]; G_['prod'] = a.get('product', G_['prod']); G_['ex'] += 1
                if not G_['prot']: G_['prot'] = a.get('Name', ''); prot[G_['prot']] = gid
    bc = collections.defaultdict(list)
    for gid, x in g.items(): x['id'] = gid; bc[x['contig']].append(x)
    for c in bc: bc[c].sort(key=lambda x: x['start'])
    return g, bc, prot
# flanks of the query locus
fl = [r for r in csv.DictReader(open(f'{S}/neighbors.tsv'), delimiter='\t') if r['locus_id'] == LOC and -5 <= int(r['rel']) <= 5]
qp = faa(f'{P}/{Q}/protein.faa')
with open(f'{W}/flank.faa', 'w') as f:
    for r in fl: f.write(f">{r['protein']}\n{qp[r['protein']]}\n")
run(f'diamond makedb --in {P}/{Q}/protein.faa -d {W}/Q --quiet')
rows = []; between = []
for T in TS:
    tp = f'{P}/{T}/protein.faa'; run(f'diamond makedb --in {tp} -d {W}/{T} --quiet')
    run(f'diamond blastp -q {W}/flank.faa -d {W}/{T} -o {W}/{T}.f.tsv --very-sensitive -e 1e-5 -k 3 --threads 16 --quiet --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore')
    best = {}
    for l in open(f'{W}/{T}.f.tsv'):
        q, s, pid, L, ql, sl, ev, bs = l.split()
        if q not in best: best[q] = (s, float(pid), int(L) / int(ql), float(bs))
    tg = faa(tp); open(f'{W}/{T}.hits.faa', 'w').write(''.join(f'>{v[0]}\n{tg[v[0]]}\n' for v in best.values()))
    run(f'diamond blastp -q {W}/{T}.hits.faa -d {W}/Q -o {W}/{T}.r.tsv --very-sensitive -e 1e-5 -k 1 --threads 16 --quiet --outfmt 6 qseqid sseqid pident')
    rev = {l.split()[0]: l.split()[1] for l in open(f'{W}/{T}.r.tsv')}
    gg, bc, pr = genes(T)
    pos = {}
    for r in fl:
        q = r['protein']
        if q not in best: rows.append([T, r['rel'], r['gene'], q, 'no_hit'] + [''] * 7); continue
        s, pid, cov, bs = best[q]; mut = rev.get(s) == q; g = gg.get(pr.get(s))
        if g: pos[int(r['rel'])] = (g['contig'], g['start'], g['end'], mut)
        rows.append([T, r['rel'], r['gene'], q, s, round(pid, 1), round(cov, 2), 'mutual' if mut else 'not_mutual', g['contig'] if g else '', g['start'] if g else '', g['end'] if g else '', g['strand'] if g else ''])
    # interval between innermost mutual hits left (rel<0) and right (rel>0) on the same contig
    L_ = [(abs(k), v) for k, v in pos.items() if k < 0 and v[3]]; R_ = [(k, v) for k, v in pos.items() if k > 0 and v[3]]
    L_.sort(key=lambda x: x[0]); R_.sort(key=lambda x: x[0])
    msg = f'{T}: left mutual hits={len(L_)}, right mutual hits={len(R_)}'
    ctgs = collections.Counter(v[0] for _, v in L_ + R_); msg += f'; contigs={dict(ctgs)}'
    if L_ and R_ and L_[0][1][0] == R_[0][1][0]:
        c = L_[0][1][0]; a = L_[0][1]; b = R_[0][1]
        lo, hi = sorted([(a[1], a[2]), (b[1], b[2])])[0][1], sorted([(a[1], a[2]), (b[1], b[2])])[1][0]
        inn = [x for x in bc[c] if x['start'] > lo and x['end'] < hi]
        msg += f'; innermost left(rel -{L_[0][0]}) and right(rel {R_[0][0]}) flank hits on {c}: gap {hi - lo} bp, {len(inn)} annotated genes between'
        for x in inn: between.append([T, c, x['name'], x['start'], x['end'], x['prod'], x['ex'], int(bool(MOB.search(x['prod'])))])
    print(msg)
with open(f'{S}/step6_{OUT}_flanks.tsv', 'w') as f:
    f.write('target\trel\tquery_gene\tquery_protein\ttarget_protein\tpident\tqcov\tmutual\tcontig\tstart\tend\tstrand\n')
    for r in rows: f.write('\t'.join(map(str, r)) + '\n')
with open(f'{S}/step6_{OUT}_between.tsv', 'w') as f:
    f.write('target\tcontig\tgene\tstart\tend\tproduct\texons\tmobile\n')
    for r in between: f.write('\t'.join(map(str, r)) + '\n')
