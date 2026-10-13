#!/usr/bin/env python
"""+-50 kb flanks around nos-like ODH genes of Rhypophila (3 assemblies; NCBI annotation).
Outputs: flank_genes.tsv (all annotated genes in window), flank50.faa (their proteins), flank50_windows.fna (window DNA)."""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import re, sys
from collections import defaultdict
D = f'{GMO_ROOT}/analysis/fusarium_odh_island/rhypophila_flank'
DL = f'{D}/dl/ncbi_dataset/data'
LOCI = [('GCA_054392455.1', 'KAM7187585.1', 'Rhypophila sp. PSN 637'),
        ('GCA_054392445.1', 'KAM7213959.1', 'Rhypophila decipiens'),
        ('GCA_033297075.1', 'KAK4208163.1', 'Rhypophila decipiens')]
FL = 50000


def read_fasta(p):
    h, b = None, []
    for l in open(p):
        if l[0] == '>':
            if h: yield h, ''.join(b)
            h, b = l[1:].split()[0], []
        else: b.append(l.strip())
    if h: yield h, ''.join(b)


import glob
out_t = open(f'{D}/flank_genes.tsv', 'w')
out_t.write('gca\torganism\todh_protein\tcontig\tcontig_len\twin_start\twin_end\tprotein_id\tlocus_tag\tstart\tend\tstrand\tn_exons\tdist_to_odh_bp\tis_odh\tproduct\n')
out_f = open(f'{D}/flank50.faa', 'w'); out_w = open(f'{D}/flank50_windows.fna', 'w')
for g, odh, org in LOCI:
    prot = dict(read_fasta(f'{DL}/{g}/protein.faa'))
    cds = defaultdict(list); info = {}; clen = {}
    for l in open(f'{DL}/{g}/genomic.gff'):
        if l.startswith('##sequence-region'):
            p = l.split(); clen[p[1]] = int(p[3]); continue
        if l[0] == '#': continue
        p = l.rstrip('\n').split('\t')
        if len(p) < 9: continue
        if p[2] == 'region' and p[3] == '1': clen.setdefault(p[0], int(p[4]))
        if p[2] != 'CDS': continue
        a = dict(kv.split('=', 1) for kv in p[8].split(';') if '=' in kv)
        pid = a.get('protein_id')
        if not pid: continue
        cds[pid].append((p[0], int(p[3]), int(p[4]), p[6]))
        info[pid] = (a.get('locus_tag', ''), a.get('product', '').replace('%2C', ','))
    c0 = cds[odh][0][0]; s0 = min(f[1] for f in cds[odh]); e0 = max(f[2] for f in cds[odh])
    ws, we = max(1, s0 - FL), min(clen[c0], e0 + FL)
    genome = glob.glob(f'{DL}/{g}/*_genomic.fna')[0]
    for h, s in read_fasta(genome):
        if h == c0: out_w.write(f'>{g}__{c0}:{ws}-{we}\n{s[ws - 1:we]}\n'); break
    genes = []
    for pid, fs in cds.items():
        if fs[0][0] != c0: continue
        s, e = min(f[1] for f in fs), max(f[2] for f in fs)
        if e < ws or s > we: continue
        d = 0 if pid == odh else (s - e0 if s > e0 else (s0 - e if e < s0 else 0))
        genes.append((s, e, fs[0][3], len(fs), pid, d))
    for s, e, st, nx, pid, d in sorted(genes):
        lt, prod = info[pid]
        out_t.write(f'{g}\t{org}\t{odh}\t{c0}\t{clen[c0]}\t{ws}\t{we}\t{pid}\t{lt}\t{s}\t{e}\t{st}\t{nx}\t{d}\t{int(pid == odh)}\t{prod}\n')
        out_f.write(f'>{pid} {g}|{c0}:{s}-{e}{st}|{prod}\n{prot[pid]}\n')
    print(g, c0, clen[c0], ws, we, len(genes), 'genes')
