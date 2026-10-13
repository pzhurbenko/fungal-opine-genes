import sys,os
from collections import OrderedDict
clean,tsv,wd=sys.argv[1:4]
cdir=os.path.join(wd,'clusters'); os.makedirs(cdir,exist_ok=True)
seqs={};h=None;buf=[]
for line in open(clean):
    line=line.rstrip('\n')
    if line.startswith('>'):
        if h: seqs[h]=''.join(buf)
        h=line[1:];buf=[]
    else: buf.append(line)
if h: seqs[h]=''.join(buf)
clusters=OrderedDict()
for line in open(tsv):
    rep,mem=line.rstrip('\n').split('\t')
    clusters.setdefault(rep,[]).append(mem)
summ=open(os.path.join(wd,'cluster_summary.tsv'),'w')
summ.write('cluster_id\trep\tsize\tmembers\n')
singles=[];multi=0
for i,(rep,mems) in enumerate(clusters.items(),1):
    cid=f'cl{i:03d}'
    summ.write(f'{cid}\t{rep}\t{len(mems)}\t{",".join(mems)}\n')
    if len(mems)>=2:
        multi+=1
        with open(os.path.join(cdir,f'{cid}.fasta'),'w') as o:
            for mm in mems: o.write(f'>{mm}\n{seqs[mm]}\n')
    else: singles.append(rep)
summ.close()
open(os.path.join(wd,'singletons.txt'),'w').write('\n'.join(singles)+('\n' if singles else ''))
print(f'  clusters: {len(clusters)} total | {multi} multi(>=2) | {len(singles)} singletons')
