import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import sys,re,collections,os; sys.path.insert(0,'/tmp/au')
from parse import *
T=f'{GMO_ROOT}/analysis/'
def key(l):
    if '__' in l:
        k=l.split('__')[0]; return re.sub(r'_(intact|pseudogene_like|fragment|pseudogene)$','',k)
    if l.startswith('NOS_'): return '_'.join(l.split('_')[:2])
    return l.split('_')[0]
def clades(t):
    out=set()
    def w(n):
        if isinstance(n,str): return frozenset([n])
        s=frozenset().union(*[w(c) for c in n]); out.add(s); return s
    w(t); return out
def mono(S,all_,cl):
    S=frozenset(S); return len(S)<=1 or S in cl or (all_-S) in cl
def comp(path):
    al,g=groups(path)
    return al,g
if __name__=="__main__":
    ml=sys.argv[1]; t=parse(open(ml).read()); al=frozenset(leaves(t)); cl=clades(t)
    byc=collections.defaultdict(set)
    for l in al: byc[key(l)].add(l)
    for k,v in sorted(byc.items()): print(k,len(v),'mono' if mono(v,al,cl) else 'NOT mono')
    for p in sys.argv[2:]:
        a,g=groups(p); print(os.path.basename(p),len(a),len(g))
        for gg in g:
            c=collections.Counter(key(l) for l in gg)
            print('   ',len(gg),dict(c),'ML-mono' if mono(gg,al,cl) else 'ML-not-mono')

CAND=('FUNGI_','CHSlike','ISL','BOL','FUNG')
def report(ml,paths):
    t=parse(open(ml).read()); al=frozenset(leaves(t)); cl=clades(t)
    for p in paths:
        a,g=groups(p)
        print(os.path.basename(p),len(a),'groups',len(g))
        for gg in g:
            c=collections.Counter(key(l) for l in gg)
            part=[l for l in gg if not key(l).startswith(CAND)]
            cand=[l for l in gg if key(l).startswith(CAND)]
            print('   n=%d %s'%(len(gg),dict(c)))
            print('      cand mono in ML: %s ; partner(non-cand, n=%d) mono in ML: %s'%(mono(cand,al,cl) if cand else '-',len(part),mono(part,al,cl) if part else '-'))

def prune_clades(t,remove):
    # clades of tree after removing leaves in `remove`
    out=set()
    def w(n):
        if isinstance(n,str): return frozenset() if n in remove else frozenset([n])
        s=frozenset().union(*[w(c) for c in n]); out.add(s); return s
    allv=w(t); return allv,out
def report2(ml,paths,cand=CAND):
    t=parse(open(ml).read()); al=frozenset(leaves(t)); cl=clades(t)
    cands=frozenset(l for l in al if key(l).startswith(cand))
    pal,pcl=prune_clades(t,cands)
    for p in paths:
        a,g=groups(p)
        for gg in g:
            c=collections.Counter(key(l) for l in gg)
            part=frozenset(l for l in gg if not key(l).startswith(cand))
            cd=frozenset(l for l in gg if key(l).startswith(cand))
            r=dict(n=len(gg),classes=dict(c),cand_n=len(cd),part_n=len(part),
              cand_mono_ML=mono(cd,al,cl) if cd else None, part_mono_ML=mono(part,al,cl) if part else None,
              part_mono_ML_candidates_pruned=mono(part,pal,pcl) if part else None)
            yield os.path.basename(p),r
