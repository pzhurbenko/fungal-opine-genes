import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
#!/usr/bin/env python3
# B3: codon alignments for fungal intact genes of SUS T-DNA clade and CUS cus clade (Boletales), from round-2 CDS.
import csv, subprocess
from Bio import SeqIO
from Bio.Seq import Seq
A = f'{GMO_ROOT}/analysis'
car = list(csv.DictReader(open(f'{A}/funnel/carriers_tree_based.tsv'), delimiter='\t'))
sets = {'SUS': [r['locus_id'] for r in car if r['lineage'].startswith('SUS') and r['status'] == 'intact'],
        'CUS': [r['locus_id'] for r in car if r['lineage'].startswith('CUS') and r['status'] == 'intact']}
for fam, ids in sets.items():
    cds = {r.id.split('|')[0]: r for r in SeqIO.parse(f'{A}/loci_pipeline/round2/proteins/{fam}_models.fna', 'fasta')}
    prot = {r.id.split('|')[0]: r for r in SeqIO.parse(f'{A}/loci_pipeline/round2/proteins/{fam}_models.faa', 'fasta')}
    ok = []
    for i in ids:
        if i not in cds: continue
        s = str(cds[i].seq).upper()
        if len(s) % 3: continue
        p = str(Seq(s).translate()).rstrip('*')
        q = str(prot[i].seq).rstrip('*')
        if '*' in p or not (p.startswith(q[:len(p)]) or q.startswith(p)) or abs(len(p)-len(q)) > 2: continue
        ok.append((i, s, p))
    print(fam, 'intact', len(ids), 'usable CDS', len(ok))
    with open(f'{fam}.cds.fna', 'w') as f:
        for i, s, p in ok: f.write(f'>{i}\n{s}\n')
    with open(f'{fam}.pep.faa', 'w') as f:
        for i, s, p in ok: f.write(f'>{i}\n{p}\n')
