#!/usr/bin/env python3
"""ENA byte-range download with eight streams and a final archive MD5 check."""
import csv,concurrent.futures,hashlib,os,time,urllib.request,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
run=sys.argv[1];row=None
for file in ROOT.glob('*.tsv'):
    for r in csv.DictReader(open(file),delimiter='\t'):
        if r.get('run_accession')==run:row=r
assert row and row['library_strategy']=='WGS'
for url,n,md5 in zip(row['fastq_ftp'].split(';'),row['fastq_bytes'].split(';'),row['fastq_md5'].split(';')):
    n=int(n);target=ROOT/'reads'/url.split('/')[-1]
    parts=target.with_name(target.name+'.parts');parts.mkdir(parents=True,exist_ok=True)
    if target.with_name(target.name+'.md5ok').exists():continue
    size=(n+7)//8
    # Preserve any already downloaded prefix as part of the first chunk.
    if target.exists():
        with open(target,'rb') as src:
            for i in range(8):
                available=min(size,max(0,target.stat().st_size-i*size))
                if not available:break
                part=parts/f'{i:02}'
                if not part.exists():
                    src.seek(i*size)
                    with open(part,'wb') as dst:dst.write(src.read(available))
    def one(i):
        start=i*size;end=min(n,start+size)-1;dest=parts/f'{i:02}'
        for attempt in range(15):
            done=dest.stat().st_size if dest.exists() else 0
            if done==end-start+1:return dest
            assert done<end-start+1,(dest,done)
            try:
                req=urllib.request.Request('https://'+url,headers={'Range':f'bytes={start+done}-{end}'})
                with urllib.request.urlopen(req,timeout=60) as response:
                    assert response.status==206 and response.headers['Content-Range'].startswith(f'bytes {start+done}-{end}/')
                    with open(dest,'ab') as f:
                        while True:
                            b=response.read(1024*1024)
                            if not b:break
                            f.write(b)
                print('chunk',i,'bytes',dest.stat().st_size,flush=True)
            except Exception as ex:
                print('retry',i,attempt,str(ex),flush=True);time.sleep(min(10,attempt+1))
        raise RuntimeError(('failed range',i))
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:files=list(ex.map(one,range(8)))
    tmp=target.with_name(target.name+'.complete');h=hashlib.md5()
    with open(tmp,'wb') as f:
        for p in files:
            with open(p,'rb') as src:
                for b in iter(lambda:src.read(4*1024*1024),b''):f.write(b);h.update(b)
    assert tmp.stat().st_size==n and h.hexdigest()==md5,(tmp,h.hexdigest())
    os.replace(tmp,target);target.with_name(target.name+'.md5ok').write_text(md5+'\n')
print('COMPLETE',run,flush=True)
