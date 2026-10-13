#!/usr/bin/env python3
"""Species clusters (connected components at mash d<0.03) over Boletales assemblies; carrier counts per lineage by cluster/genus."""
import csv, os
from collections import defaultdict, Counter
gen={}
for g,org,genus,f in csv.reader(open('genomes.tsv'),delimiter='\t'): gen[os.path.basename(f)]=(g,org,genus)
G={v[0]:v for v in gen.values()}
adj=defaultdict(set); dist={}
for a,b,d,p,sh in csv.reader(open('dist.tsv'),delimiter='\t'):
    ga,gb=gen[os.path.basename(a)][0],gen[os.path.basename(b)][0]; d=float(d)
    if ga!=gb: dist[(ga,gb)]=d
    if ga!=gb and d<0.03: adj[ga].add(gb); adj[gb].add(ga)
seen=set(); clusters=[]
for g in G:
    if g in seen: continue
    st=[g]; comp=set()
    while st:
        x=st.pop()
        if x in comp: continue
        comp.add(x); st.extend(adj[x]-comp)
    seen|=comp; clusters.append(sorted(comp))
clusters.sort(key=lambda c:(-len(c),c[0]))
# cluster name = most common organism among members that are not 'extra'; corrected label where NCBI names disagree
pres={r['gca']:r for r in csv.DictReader(open('../distribution/presence_step2.tsv'),delimiter='\t')}
lins=[k for k in next(iter(pres.values())).keys() if ':' in k]
def status(g,l):
    v=pres.get(g,{}).get(l,'')
    return v if v not in ('','not_homolog','no_model') else ''
out=open('clusters.tsv','w'); out.write('cluster_id\tcluster_name\tn_genomes\tgenomes\tncbi_names\tname_conflict\t'+'\t'.join(lins)+'\n')
rows=[]
# 517 change: clusters containing genomes of the old 140-set keep their previously accepted name (old cluster name, if not 'extra'),
# so that adding new assemblies does not silently rename carrier clusters (tie-break would otherwise change e.g. B. edulis -> B. chippewaensis)
oldname={}
for r in csv.DictReader(open('../boletales_mash/clusters.tsv'),delimiter='\t'):
    for g in r['genomes'].split(';'): oldname[g]=r['cluster_name']
for i,c in enumerate(clusters,1):
    names=[G[g][1] for g in c if G[g][1]!='extra']
    cn=Counter(n for n in names); 
    # prefer a name that is not 'sp.'/'aff.' and appears in the 2021 reference assemblies
    best=sorted(cn.items(),key=lambda x:(-x[1],x[0]))
    name=best[0][0] if best else G[c[0]][1]
    on=Counter(oldname[g] for g in c if g in oldname and oldname[g]!='extra')
    if on: name=sorted(on.items(),key=lambda x:(-x[1],x[0]))[0][0]
    if 'GCA_019364595.1' in c: name='Butyriboletus roseoflavus'
    conflict='yes' if len({n.split()[0] for n in names})>1 else ''  # first word of NCBI name, as before
    st={l:sorted({status(g,l) for g in c if status(g,l)}) for l in lins}
    out.write('\t'.join([f'C{i:03d}',name,str(len(c)),';'.join(c),';'.join(sorted(set(names))),conflict]+[','.join(st[l]) for l in lins])+'\n')
    rows.append((f'C{i:03d}',name,c,st))
out.close()
# carriers per lineage: clusters with intact / any homolog
summ=open('carriers_by_cluster.tsv','w'); summ.write('lineage\tn_clusters_any\tn_clusters_intact\tn_genera_any\tclusters_intact\n')
for l in lins:
    anyc=[r for r in rows if r[3][l]]; intc=[r for r in rows if 'intact' in r[3][l]]
    genera={r[1].split()[0] for r in anyc}
    summ.write('\t'.join([l,str(len(anyc)),str(len(intc)),str(len(genera)),'; '.join(f"{r[1]}({len(r[2])})" for r in intc)])+'\n')
summ.close()
print('clusters',len(clusters),'multi-genome clusters',sum(1 for c in clusters if len(c)>1))
for r in rows:
    if len(r[2])>1: print(r[0],r[1],len(r[2]),';'.join(sorted({G[g][1] for g in r[2]})))
