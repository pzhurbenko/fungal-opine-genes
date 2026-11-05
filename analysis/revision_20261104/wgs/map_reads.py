#!/usr/bin/env python3
"""Whole-assembly mapping of checksum-verified matched-strain genomic reads."""
import csv, hashlib, json, subprocess, sys, os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent
ENV=Path(os.environ.get('GMO_CONDA_ENVS', str(Path.home()/'miniforge3/envs')))
run=sys.argv[1]
config={
 'SRR23927955':('GCA_046056315.1','FP133','sr','data/fungal_genomes/fusarium/fus/ncbi_dataset/data/GCA_046056315.1/GCA_046056315.1_ASM4605631v1_genomic.fna'),
 'SRR23927956':('GCA_046056315.1','FP133','map-ont','data/fungal_genomes/fusarium/fus/ncbi_dataset/data/GCA_046056315.1/GCA_046056315.1_ASM4605631v1_genomic.fna'),
 'ERR12742700':('GCA_964275055.1','CLIB189','sr','data/fungal_genomes/extra2/x/ncbi_dataset/data/GCA_964275055.1/GCA_964275055.1_MOSE0_genomic.fna'),
 'SRR305351':('GCA_000149955.2','Fol4287','sr','data/fungal_genomes/fusarium/fus/ncbi_dataset/data/GCA_000149955.2/GCA_000149955.2_ASM14995v2_genomic.fna')}
assembly,strain,preset,rel=config[run]
rows=[]
for f in HERE.glob('*.tsv'):
    rows.extend(r for r in csv.DictReader(open(f),delimiter='\t') if r.get('run_accession')==run)
r=rows[0];paths=[HERE/'reads'/u.split('/')[-1] for u in r['fastq_ftp'].split(';')]
for p,expected in zip(paths,r['fastq_md5'].split(';')):
    assert p.with_name(p.name+'.md5ok').read_text().strip()==expected,p
d=HERE/run;d.mkdir(exist_ok=True)
fasta=ROOT/rel
(d/'provenance.json').write_text(json.dumps(dict(run=run,assembly=assembly,strain=strain,
    reference=str(fasta),reference_sha256=hashlib.sha256(fasta.read_bytes()).hexdigest(),
    metadata=r,mapping_preset=preset),indent=2)+'\n')
samtools=ENV/'annot/bin/samtools';minimap=ENV/'minimap2_env/bin/minimap2'
with open(d/'minimap2.log','w') as err,open(d/'samtools_sort.log','w') as err2:
    mapper=subprocess.Popen([str(minimap),'-ax',preset,'-t','12',str(fasta)]+[str(p) for p in paths],stdout=subprocess.PIPE,stderr=err)
    sorter=subprocess.run([str(samtools),'sort','-@','3','-m','1G','-o',str(d/'reads.bam'),'-'],stdin=mapper.stdout,stderr=err2)
    mapper.stdout.close();assert mapper.wait()==0 and sorter.returncode==0
subprocess.run([str(samtools),'index','-@','4',str(d/'reads.bam')],check=True)
with open(d/'flagstat.txt','w') as f:subprocess.run([str(samtools),'flagstat',str(d/'reads.bam')],stdout=f,check=True)
subprocess.run([str(ENV/'PlantCad2/bin/python'),str(HERE/'summarize_reads.py'),run],check=True)
print('COMPLETE',run,flush=True)
