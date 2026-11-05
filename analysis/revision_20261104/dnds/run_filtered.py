#!/usr/bin/env python3
"""Integration-filter sensitivity: original CDS, same alignment/model, two M0 starts.

Run with PlantCad2 Python on an LSF worker. Original files are never modified.
"""
import csv, hashlib, json, re, subprocess, sys, os
from pathlib import Path
from Bio import SeqIO, Phylo

ROOT = Path(__file__).resolve().parents[3]
ENV = Path(os.environ.get('GMO_CONDA_ENVS', str(Path.home()/'miniforge3/envs')))
HERE = Path(__file__).resolve().parent
status = {r['locus_id']: r for r in csv.DictReader(
    open(ROOT/'paper/tables/integration_status.tsv'), delimiter='\t')}
name = sys.argv[1]
parent = {'SUS_conf': ROOT/'analysis/dnds/SUS',
          'cus_nonBol_conf': ROOT/'analysis/dnds/extra/cus_nonBol/cus_nonBol'}[name]
out = HERE/name
out.mkdir(exist_ok=True)
cds = {r.id: r for r in SeqIO.parse(str(parent)+'.cds.fna', 'fasta')}
pep = {r.id: r for r in SeqIO.parse(str(parent)+'.pep.faa', 'fasta')}
keep = [i for i in cds if status[i]['status'] != 'uncertain']
assert len(keep) == {'SUS_conf':15, 'cus_nonBol_conf':19}[name]
qc = []
for i in cds:
    dna = str(cds[i].seq).upper()
    aa = str(pep[i].seq)
    assert len(dna)%3 == 0 and set(dna) <= set('ACGTN'), i
    translated = str(cds[i].seq.translate())
    assert translated.rstrip('*') == aa.rstrip('*'), i
    assert '*' not in translated.rstrip('*'), i
    qc.append(dict(locus_id=i,status=status[i]['status'],included=i in keep,
                   cds_bp=len(dna),ambiguous_bases=dna.count('N'), sha256=hashlib.sha256(dna.encode()).hexdigest()))
(out/'input_audit.json').write_text(json.dumps(qc,indent=2)+'\n')
SeqIO.write([cds[i] for i in keep],out/(name+'.cds.fna'),'fasta')
SeqIO.write([pep[i] for i in keep],out/(name+'.pep.faa'),'fasta')
def run(cmd,wd=out,log=None):
    with open(wd/(log or (Path(str(cmd[0])).name+'.log')),'w') as f:
        subprocess.run([str(c) for c in cmd],cwd=wd,stdout=f,stderr=subprocess.STDOUT,check=True)
with open(out/(name+'.pep.aln'),'w') as f, open(out/'mafft.log','w') as err:
    subprocess.run([str(ENV/'annot/bin/mafft'),'--localpair','--maxiterate','1000',
                    '--thread','2',str(out/(name+'.pep.faa'))],stdout=f,stderr=err,check=True)
rows=[]
for r in SeqIO.parse(out/(name+'.pep.aln'),'fasta'):
    k=0; row=[]; s=str(cds[r.id].seq).upper()
    for a in str(r.seq):
        if a=='-':row.append('---')
        else:row.append(s[3*k:3*k+3]);k+=1
    assert k*3 in (len(s),len(s)-3),r.id
    rows.append((r.id,row))
results=[]
for mode,frac in [('strict',0),('loose',.2)]:
    wd=out/mode;wd.mkdir(exist_ok=True)
    columns=[j for j in range(len(rows[0][1])) if
             sum(row[j]=='---' or 'N' in row[j] for _,row in rows)<=frac*len(rows)
             and all(row[j] not in ('TAA','TGA','TAG') for _,row in rows)]
    (wd/(name+'.codon.phy')).write_text(f'{len(rows)} {len(columns)*3}\n'+''.join(
        i+'  '+''.join(row[j] for j in columns)+'\n' for i,row in rows))
    (wd/(name+'.codon.fna')).write_text(''.join('>'+i+'\n'+''.join(row[j] for j in columns)+'\n' for i,row in rows))
    (wd/'retained_columns.json').write_text(json.dumps(columns)+'\n')
    if mode=='strict':
        run([ENV/'phylo_env/bin/iqtree','-s',name+'.codon.fna','-st','DNA','-m','GTR+G',
             '-T','2','-seed','20261104','-pre',name+'.tree','-redo'],wd)
        tree=Phylo.read(wd/(name+'.tree.treefile'),'newick')
        for c in tree.find_clades():c.confidence=None
    Phylo.write(tree,wd/(name+'.tree.nwk'),'newick')
    treepath=wd/(name+'.tree.nwk')
    treepath.write_text(f'{len(rows)} 1\n'+treepath.read_text())
    for start in (.4,1.5):
        d=wd/f'M0_start{start}';d.mkdir(exist_ok=True)
        ctl=f'''seqfile = ../{name}.codon.phy
treefile = ../{name}.tree.nwk
outfile = M0.out
noisy = 0
verbose = 0
runmode = 0
seqtype = 1
CodonFreq = 2
model = 0
NSsites = 0
icode = 0
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = {start}
cleandata = 0
'''
        (d/'M0.ctl').write_text(ctl)
        run([ENV/'PlantCad2/bin/codeml','M0.ctl'],d)
        txt=(d/'M0.out').read_text()
        omega=float(re.search(r'omega \(dN/dS\)\s*=\s*(\S+)',txt)[1])
        likelihood=float(re.search(r'lnL\(.*?\):\s*([-\d.]+)',txt)[1])
        results.append(dict(dataset=name,mode=mode,n=len(rows),codons=len(columns),
                            start_omega=start,omega=omega,lnL=likelihood))
        (out/'M0_results.json').write_text(json.dumps(results,indent=2)+'\n')
        print(results[-1],flush=True)
for mode in ('strict','loose'):
    pair=[r for r in results if r['mode']==mode]
    assert abs(pair[0]['lnL']-pair[1]['lnL'])<.01,(mode,pair)
print('COMPLETE',name,flush=True)
