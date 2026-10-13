import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import sys,os,glob
from collections import defaultdict,Counter
OUT=f'{GMO_ROOT}/analysis/domains'
CL=f'{GMO_ROOT}/analysis/clustering'
fams=['ACS','CUS','ORF358','SUS']
def load_idmap(fam):
    d={}
    with open(f'{CL}/{fam}/idmap.tsv') as f:
        next(f)
        for line in f:
            a=line.rstrip('\n').split('\t'); d[a[0]]=a[1] if len(a)>1 else ''
    return d
def kind(hdr):
    if hdr.startswith('Plant_'): return 'Plant'
    if hdr.startswith('Agro_'): return 'Agro'
    return 'Fungi'
for fam in fams:
    dt=f'{OUT}/{fam}.domtbl'
    if not os.path.exists(dt): print(f'{fam}: no domtbl'); continue
    idmap=load_idmap(fam)
    # per query: list of (env_from, dom, acc, ievalue, score)
    perq=defaultdict(list)
    for line in open(dt):
        if line.startswith('#'): continue
        f=line.split()
        # domtblout: 0 target(dom name),1 acc,2 tlen,3 query,4 qacc,5 qlen,6 fullEval,7 fullScore,...
        # 11 i-Evalue,12 dom score,... 19 env_from,20 env_to
        dom=f[0]; acc=f[1]; query=f[3]; ievalue=float(f[12]); score=float(f[13])
        env_from=int(f[19]); env_to=int(f[20])
        perq[query].append((env_from,env_to,dom,acc,ievalue,score))
    # all queries (incl those with no hit)
    allq=set(idmap.keys())
    archi=Counter(); domcount=Counter(); nohit=[]
    bykind=defaultdict(Counter)
    rows=[]
    for q in allq:
        hits=sorted(perq.get(q,[]))
        k=kind(idmap[q])
        if not hits:
            nohit.append(q); archi['(no Pfam domain)']+=1; bykind[k]['(none)']+=1
            rows.append((q,k,'(none)','')); continue
        doms=[h[2] for h in hits]
        for d in set(doms): domcount[d]+=1
        arch=' + '.join(doms)
        archi[arch]+=1
        bykind[k][doms[0]]+=1  # primary domain
        rows.append((q,k,arch,';'.join(f'{h[2]}:{h[0]}-{h[1]}(E{h[4]:.0e})' for h in hits)))
    n=len(allq)
    print('='*80); print(f'{fam}: {n} seqs, {len(nohit)} without any Pfam domain')
    print('  top domains (по числу белков):')
    for d,c in domcount.most_common(8): print(f'     {d:22} {c}')
    print('  доменные архитектуры (топ):')
    for a,c in archi.most_common(8): print(f'     {c:4}  {a}')
    print('  первичный домен по типу (Plant/Agro/Fungi):')
    for k in ['Plant','Agro','Fungi']:
        if bykind[k]: print(f'     {k:6}:',dict(bykind[k].most_common(5)))
    # write per-seq table
    with open(f'{OUT}/{fam}_domain_arch.tsv','w') as o:
        o.write('clean_id\tkind\tarchitecture\tdomain_coords\toriginal_header\n')
        for q,k,arch,coords in sorted(rows,key=lambda x:(x[1],x[2])):
            o.write(f'{q}\t{k}\t{arch}\t{coords}\t{idmap[q]}\n')
print('\nper-seq tables: analysis/domains/<FAM>_domain_arch.tsv')
