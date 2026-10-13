#!/usr/bin/env python3
# nos-ODH Boletales: CDS from miniprot models (query = own-locus protein from nos_check/boletales_odh.faa, identity 1.0) on locus region +-2 kb
import re
from Bio import SeqIO
from Bio.Seq import Seq
reg = {r.id: str(r.seq) for r in SeqIO.parse('regions.fa', 'fasta')}
best = {}; cds = {}
for l in open('mp.gff'):
    f = l.rstrip('\n').split('\t')
    if len(f) < 9: continue
    a = dict(x.split('=', 1) for x in f[8].split(';') if '=' in x)
    if f[2] == 'mRNA' and float(a['Identity']) == 1.0 and f[0] not in best: best[f[0]] = a['ID']
    if f[2] == 'CDS': cds.setdefault(a['Parent'], []).append((int(f[3]), int(f[4]), f[6]))
fc = open('../nosODH_Bol.cds.fna', 'w'); fp = open('../nosODH_Bol.pep.faa', 'w')
for loc, mid in best.items():
    seg = sorted(cds[mid]); st = seg[0][2]
    s = ''.join(reg[loc][a-1:b] for a, b, _ in seg)
    if st == '-': s = str(Seq(s).reverse_complement())
    p = str(Seq(s).translate())
    if p.endswith('*'): p = p[:-1]; s = s[:-3]
    print(loc, len(s), len(s) % 3, p.count('*'), p[:10])
    if '*' not in p and len(s) % 3 == 0: fc.write(f'>{loc}\n{s}\n'); fp.write(f'>{loc}\n{p}\n')
