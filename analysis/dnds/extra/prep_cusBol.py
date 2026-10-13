import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
#!/usr/bin/env python3
# cus_Bol: 14 intact cus-like genes of Boletales (TableS1a), CDS/protein from ../CUS.{cds.fna,pep.faa} (same source as the all-cus run, prep.py).
# cus_Bol_sp: one copy per species (cluster tag stripped), longest CDS, ties by locus_id.
import csv, re, os
from Bio import SeqIO
S=f'{GMO_ROOT}/paper/supplementary/TableS1a_Loci_of_bacterial_origin.tsv'
R=[r for r in csv.DictReader(open(S),delimiter='\t') if r['lineage']=='cus-like' and r['status']=='intact' and r['order']=='Boletales']
cds={r.id:str(r.seq) for r in SeqIO.parse('../CUS.cds.fna','fasta')}; pep={r.id:str(r.seq) for r in SeqIO.parse('../CUS.pep.faa','fasta')}
ids=[r['locus_id'] for r in R]; assert len(ids)==14 and all(i in cds for i in ids)
sp={r['locus_id']:re.sub(r' \(cluster \w+\)','',r['species']) for r in R}
best={}
for i in sorted(ids):
    if sp[i] not in best or len(cds[i])>len(cds[best[sp[i]]]): best[sp[i]]=i
for name,sel in (('cus_Bol',ids),('cus_Bol_sp',sorted(best.values()))):
    os.makedirs(name,exist_ok=True)
    with open(f'{name}/{name}.cds.fna','w') as fc, open(f'{name}/{name}.pep.faa','w') as fp:
        for i in sel: fc.write(f'>{i}\n{cds[i]}\n'); fp.write(f'>{i}\n{pep[i]}\n')
    print(name,len(sel),[sp[i] for i in sel])
