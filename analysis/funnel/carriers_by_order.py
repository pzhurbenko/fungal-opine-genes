"""Carriers of rhizobial-origin candidates by fungal order. Inputs: paper/supplementary/TableS3_*.tsv, analysis/refset/screen/step3/genome_meta.tsv, taxdump.
Date 2026-10-13. Writes carriers_by_order.tsv here.
Run: python3 -I $GMO_ROOT/analysis/funnel/carriers_by_order.py
"""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import csv,re,collections,sys
G=f'{GMO_ROOT}/'
T=G+'analysis/refset/screen/seeds_ext/reciprocal/taxdump/'
par={};rank={}
for l in open(T+'nodes.dmp'):
    f=l.split('\t|\t'); par[int(f[0])]=int(f[1]); rank[int(f[0])]=f[2]
name={}
for l in open(T+'names.dmp'):
    f=l.split('\t|\t')
    if f[3].startswith('scientific name'): name[int(f[0])]=f[1]
merged={int(l.split('\t|\t')[0]):int(l.split('\t|\t')[1].split('\t')[0]) for l in open(T+'merged.dmp')}
def lin(t):
    t=merged.get(t,t); d={}
    while t in par and t!=par[t]:
        d.setdefault(rank[t],name[t]); t=par[t]
    return d
num=lambda a:re.sub(r'\.\d+$','',re.sub(r'^GC[AF]_','',a))
ver=lambda a:int(a.rsplit('.',1)[1])
m=list(csv.DictReader(open(G+'analysis/refset/screen/step3/genome_meta.tsv'),delimiter='\t'))
grp={}
for r in m: grp.setdefault(num(r['gca']),[]).append(r)
keep={n:sorted(v,key=lambda r:(not r['gca'].startswith('GCA_'),-ver(r['gca']),r['gca']))[0] for n,v in grp.items()}
info={}; unk=collections.Counter()
for n,r in keep.items():
    d=lin(int(r['taxid'])); o=d.get('order')
    ph=d.get('phylum','NA')
    if not o:
        o=f"unassigned ({d.get('class') or d.get('subphylum') or ph})"; unk[o]+=1
    info[n]=(o,ph,r['organism'].split()[0].strip('[]'))
print(len(keep),'unassigned',dict(unk),file=sys.stderr)
s3=list(csv.DictReader(open(G+'paper/supplementary/TableS3_Proteins_with_non_fungal_closest_relatives_and_the.tsv'),delimiter='\t'))
cand=[r for r in s3 if r['fate'] in('T-DNA candidate','non-T-DNA (rhizobial) candidate') and r['duplicate_assembly'].strip().lower()!='yes']
print(len(s3),len(cand),collections.Counter(r['fate'] for r in cand),file=sys.stderr)
nprot=collections.Counter(v[0] for v in info.values()); ph={v[0]:v[1] for v in info.values()}
C=collections.defaultdict(lambda:{'a':set(),'t':set(),'n':set(),'g':set()})
miss=0
for r in cand:
    n=num(r['assembly'])
    if n not in info: miss+=1; continue
    o,_,ge=info[n]; c=C[o]; c['a'].add(n); c['g'].add(ge)
    c['t' if r['fate']=='T-DNA candidate' else 'n'].add(n)
print('missing',miss,'carriers',len({n for c in C.values() for n in c['a']}),file=sys.stderr)
rows=[]
for o,c in C.items():
    rows.append((o,ph[o],nprot[o],len(c['a']),round(100*len(c['a'])/nprot[o],2),len(c['t']),len(c['n']),'; '.join(sorted(c['g']))))
rows.sort(key=lambda r:(-r[4],-r[3]))
with open(G+'analysis/funnel/carriers_by_order.tsv','w') as f:
    w=csv.writer(f,delimiter='\t',lineterminator='\n')
    w.writerow('order phylum n_proteomes n_with_candidates share_pct n_TDNA_candidate_proteomes n_nonTDNA_candidate_proteomes genera_with_candidates'.split())
    w.writerows(rows)
for r in rows: print('\t'.join(map(str,r)))
B=nprot['Boletales'];Bc=len(C['Boletales']['a']);tot=sum(nprot.values());tc=sum(len(c['a']) for c in C.values())
print('Boletales',B,Bc,round(100*Bc/B,2)); print('others',tot-B,tc-Bc,round(100*(tc-Bc)/(tot-B),2)); print('all',tot,tc)
