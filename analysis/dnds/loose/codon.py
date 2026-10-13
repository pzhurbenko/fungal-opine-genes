#!/usr/bin/env python3
# back-translate protein alignment to codon alignment (PHYLIP for PAML + FASTA)
import sys
from Bio import SeqIO
fam = sys.argv[1]
cds = {r.id: str(r.seq).upper() for r in SeqIO.parse(f'{fam}.cds.fna', 'fasta')}
aln = list(SeqIO.parse(f'{fam}.pep.aln', 'fasta'))
rows = []
for r in aln:
    s = cds[r.id]; out = ''; k = 0
    for a in str(r.seq):
        if a == '-': out += '---'
        else: out += s[3*k:3*k+3]; k += 1
    rows.append((r.id, out))
L = len(rows[0][1])
# drop codon columns with any gap (conservative for PAML)
keep = [i for i in range(0, L, 3) if sum(1 for x in rows if x[1][i:i+3] == "---" or "N" in x[1][i:i+3]) <= 0.2*len(rows)]
with open(f'{fam}.codon.phy', 'w') as f:
    f.write(f' {len(rows)} {3*len(keep)}\n')
    for n, s in rows: f.write(n + '  ' + ''.join(s[i:i+3] for i in keep) + '\n')
with open(f'{fam}.codon.fna', 'w') as f:
    for n, s in rows: f.write(f'>{n}\n' + ''.join(s[i:i+3] for i in keep) + '\n')
print(fam, len(rows), 'seqs', 3*len(keep), 'codon-aligned bp kept of', L)
