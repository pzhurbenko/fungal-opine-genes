#!/usr/bin/env python3
"""Download a specified ENA run; require exact byte count and MD5."""
import csv, hashlib, subprocess, sys
from pathlib import Path
root=Path(__file__).resolve().parent
run=sys.argv[1]
rows=[]
for f in root.glob('*.tsv'):
    rows.extend(r for r in csv.DictReader(open(f),delimiter='\t') if r.get('run_accession')==run)
assert rows,run
r=rows[0]
assert r['library_strategy']=='WGS',r
d=root/'reads';d.mkdir(exist_ok=True)
for url,n,md5 in zip(r['fastq_ftp'].split(';'),r['fastq_bytes'].split(';'),r['fastq_md5'].split(';')):
    p=d/url.split('/')[-1]
    if not (p.exists() and p.stat().st_size==int(n)):
        success=False
        for scheme in ('ftp://','https://'):
            result=subprocess.run(['curl','-fsSL','--retry','5','--retry-all-errors','--connect-timeout','30','--speed-time','120',
                        '--speed-limit','1000','-C','-','-o',str(p),scheme+url])
            if result.returncode==0:
                success=True;break
        assert success,('Download failed',url)
    assert p.stat().st_size==int(n),(p,p.stat().st_size,n)
    h=hashlib.md5()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    assert h.hexdigest()==md5,(p,h.hexdigest(),md5)
    (d/(p.name+'.md5ok')).write_text(md5+'\n')
print('COMPLETE',run,flush=True)
