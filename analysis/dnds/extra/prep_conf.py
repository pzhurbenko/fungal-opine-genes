#!/usr/bin/env python3
# *_conf sets (2026-10-13): parent sets minus loci with integration status 'uncertain' in paper/tables/integration_status.tsv.
# CDS/protein copied from the parent set files (same source). cus_Bol_sp_conf: one copy per species (cluster tag stripped), longest CDS, ties by locus_id (rule of prep_cusBol.py).
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import csv, re, os
from Bio import SeqIO
st={r['locus_id']:r['status'] for r in csv.DictReader(open(f'{GMO_ROOT}/paper/tables/integration_status.tsv'),delimiter='\t')}
spc={r['locus_id']:re.sub(r' \(cluster \w+\)','',r['species']) for r in csv.DictReader(open(f'{GMO_ROOT}/paper/tables/integration_status.tsv'),delimiter='\t')}
def load(p): return {r.id:str(r.seq) for r in SeqIO.parse(p,'fasta')}
def keep(parent):
    cds=load(f'{parent}/{parent}.cds.fna'); pep=load(f'{parent}/{parent}.pep.faa')
    ids=[i for i in cds if st[i]!='uncertain']; return ids,cds,pep
def write(name,ids,cds,pep):
    os.makedirs(name,exist_ok=True)
    with open(f'{name}/{name}.cds.fna','w') as fc, open(f'{name}/{name}.pep.faa','w') as fp:
        for i in ids: fc.write(f'>{i}\n{cds[i]}\n'); fp.write(f'>{i}\n{pep[i]}\n')
    print(name,len(ids),ids)
for parent in ('cus_Bol','acs_pKerr14a','nosODH_Bol'):
    ids,cds,pep=keep(parent)
    name={'cus_Bol':'cus_Bol_conf','acs_pKerr14a':'acs_pKerr14a_conf','nosODH_Bol':'nosODH_Bol_conf'}[parent]
    write(name,ids,cds,pep)
    if parent=='cus_Bol':
        best={}
        for i in sorted(ids):
            if spc[i] not in best or len(cds[i])>len(cds[best[spc[i]]]): best[spc[i]]=i
        write('cus_Bol_sp_conf',sorted(best.values()),cds,pep)
