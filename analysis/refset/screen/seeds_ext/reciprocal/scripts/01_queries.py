import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..', '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import pandas as pd, numpy as np, os
B=f'{GMO_ROOT}/'
d=pd.read_csv(B+'analysis/refset/screen/seeds_ext/seeds_ext_best.tsv',sep='\t')
d=d[d.dup_dropped==0].copy()
RARE=['sus','acs','acs_pKerr14a','cus','mis','chs']
CTRL=['mas1','ags','ocs','nos','vis']
rows=[]
r=d[d.gene.isin(RARE)]
rare_keys=set(zip(r.gca,r.protein))
for (g,p),x in r.groupby(['gca','protein']):
    x=x.sort_values('pident',ascending=False)
    rows.append(dict(protein=p,gca=g,organism=x.organism.iloc[0],gene=','.join(sorted(set(x.gene))),pident=x.pident.iloc[0],cov_seed=x.cov_seed.iloc[0],group='rare'))
for gene in CTRL:
    x=d[(d.gene==gene)&~d.set_index(['gca','protein']).index.isin(rare_keys)].drop_duplicates(['gca','protein'])
    x=x.sort_values(['pident','bits','protein'],ascending=[False,False,True])
    top=x.head(50); rest=x.iloc[50:]
    rnd=rest.sample(n=100,random_state=1)
    for grp,t in (('control_top',top),('control_random',rnd)):
        for _,y in t.iterrows():
            rows.append(dict(protein=y.protein,gca=y.gca,organism=y.organism,gene=gene,pident=y.pident,cov_seed=y.cov_seed,group=grp))
for _,y in d[d.gene=='mas2'].iterrows():
    rows.append(dict(protein=y.protein,gca=y.gca,organism=y.organism,gene='mas2',pident=y.pident,cov_seed=y.cov_seed,group='control_top'))
q=pd.DataFrame(rows)
# unique query id
q['qid']=q.protein
dup=q.protein.duplicated(keep=False)
q.loc[dup,'qid']=q.gca+'|'+q.protein
print('rare',(q.group=='rare').sum(),'dup protein ids',dup.sum())
seqs={}
for g,x in q.groupby('gca'):
    want=set(x.protein); cur=None; buf=[]
    got={}
    for l in open(f'{B}data/fungal_proteomes/fungi_annotated/ncbi_dataset/data/{g}/protein.faa'):
        if l[0]=='>':
            if cur in want: got[cur]=''.join(buf)
            cur=l[1:].split()[0]; buf=[]
        else: buf.append(l.strip())
    if cur in want: got[cur]=''.join(buf)
    for p in want: seqs[(g,p)]=got[p]
with open('queries.faa','w') as f:
    for _,y in q.iterrows(): f.write(f'>{y.qid}\n{seqs[(y.gca,y.protein)]}\n')
q['qlen']=[len(seqs[(a,b)]) for a,b in zip(q.gca,q.protein)]
q[['qid','protein','gca','organism','gene','pident','cov_seed','group','qlen']].to_csv('queries.tsv',sep='\t',index=False)
print(q.groupby(['group']).size(), len(q), q.qid.nunique())
