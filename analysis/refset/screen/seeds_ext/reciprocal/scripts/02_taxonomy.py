"""Resolve organism names (from nr stitle [..] and query genomes) -> kingdom/class/order/genus via NCBI taxdump (current, downloaded 2026-10-06).
After this script, scripts/02b_taxonomy_fix.py is mandatory: it re-resolves non-exact and homonymous names by the taxid of the hit protein (2026-10-13)."""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..', '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import sys, re, json, csv
R=f'{GMO_ROOT}/analysis/refset/screen/seeds_ext/reciprocal/'
hits=R+'diamond_nr.tsv'
BR=re.compile(r'\[([^\[\]]+)\]')
def orgs_from_stitle(s):
    return [o.strip() for o in BR.findall(s)]
def clean(o):
    return re.sub(r'\s+(strain|str\.|isolate|sp\.|subsp\.|var\.|f\. sp\.|pv\.|cf\.).*$','',o).strip()
if __name__=='__main__':
    names=set()
    for l in open(R+'queries.tsv'):
        names.add(l.split('\t')[3])
    for l in open(hits):
        for o in orgs_from_stitle(l.rstrip('\n').split('\t',8)[8]): names.add(o)
    names.discard('organism')
    # candidates per name
    cand={}
    for o in names:
        w=clean(o).split()
        c=[('exact',o.lower()),('cleaned',clean(o).lower())]
        if len(w)>=2: c.append(('binomial',' '.join(w[:2]).lower()))
        if len(w)>=1: c.append(('genus',w[0].lower()))
        cand[o]=c
    want=set(x for c in cand.values() for _,x in c)
    n2t={}
    for l in open(R+'taxdump/names.dmp'):
        f=l.rstrip('\n').split('\t|\t')
        nm=f[1].lower()
        if nm in want and ('scientific' in f[3] or 'synonym' in f[3] or 'equivalent' in f[3] or 'genbank common' not in f[3]):
            t=int(f[0]); n2t.setdefault(nm,[]).append((0 if 'scientific' in f[3] else 1,t))
    par={};rank={}
    for l in open(R+'taxdump/nodes.dmp'):
        f=l.split('\t|\t'); par[int(f[0])]=int(f[1]); rank[int(f[0])]=f[2]
    sci={}
    merged={int(l.split('\t|\t')[0]):int(l.split('\t|\t')[1].split('\t')[0]) for l in open(R+'taxdump/merged.dmp')}
    K={4751:'Fungi',2:'Bacteria',33090:'Viridiplantae',2157:'Archaea',33208:'Metazoa',10239:'Viruses',2759:'Eukaryota_other'}
    def lin(t):
        out=[];seen=0
        while t in par and t!=par[t] and seen<100:
            out.append(t); t=par[t]; seen+=1
        out.append(1); return out
    # scientific names for lineage ranks
    need=set()
    res={}
    for o,c in cand.items():
        t=None;lev=None
        for lv,x in c:
            if x in n2t:
                t=sorted(n2t[x])[0][1]; lev=lv; break
        if t is None: res[o]=dict(kingdom='NA',cls='',order='',genus='',level='none'); continue
        t=merged.get(t,t); L=lin(t)
        k='NA'
        for a in L:
            if a in K and a!=2759: k=K[a]; break
        if k=='NA' and 2759 in L: k='Eukaryota_other'
        d={}
        for a in L:
            if rank.get(a) in('genus','order','class') and rank[a] not in d: d[rank[a]]=a; need.add(a)
        res[o]=dict(kingdom=k,cls=d.get('class'),order=d.get('order'),genus=d.get('genus'),level=lev)
    for l in open(R+'taxdump/names.dmp'):
        f=l.rstrip('\n').split('\t|\t')
        if 'scientific' in f[3] and int(f[0]) in need: sci[int(f[0])]=f[1]
    for o,d in res.items():
        for k in('cls','order','genus'): d[k]=sci.get(d[k],'') if d[k] else ''
        # fallback genus from name
        if not d['genus']: d['genus']=clean(o).split()[0] if o else ''
    json.dump(res,open(R+'organism_taxonomy.json','w'))
    from collections import Counter
    print(len(res),Counter(d['level'] for d in res.values()),Counter(d['kingdom'] for d in res.values()))
