import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import os,glob,re
BASE=f'{GMO_ROOT}/analysis/clustering'
fams=['ACS','CUS','ORF358','SUS']
def load_idmap(p):
    m={}
    with open(p) as f:
        next(f)
        for line in f:
            a=line.rstrip('\n').split('\t'); m[a[0]]=a[1] if len(a)>1 else ''
    return m
def parse_mmseqs(p):
    # cluster_summary.tsv: cluster_id rep size members
    mem2cl={}; sizes=[]
    with open(p) as f:
        next(f)
        for line in f:
            cid,rep,size,members=line.rstrip('\n').split('\t')
            sizes.append(int(size))
            for m in members.split(','): mem2cl[m]=cid
    return mem2cl,sizes
def parse_of(fam):
    rd=glob.glob(f'{BASE}/{fam}/of_input/OrthoFinder/Results_*')
    if not rd: return {},[],[]
    rd=rd[0]
    mem2og={}; ogsizes=[]
    txt=os.path.join(rd,'Orthogroups','Orthogroups.txt')
    if os.path.exists(txt):
        for line in open(txt):
            og,_,rest=line.partition(':')
            members=rest.split()
            ogsizes.append(len(members))
            for m in members: mem2og[m]=og
    unassigned=[]
    ua=os.path.join(rd,'Orthogroups','Orthogroups_UnassignedGenes.tsv')
    if os.path.exists(ua):
        lines=open(ua).read().splitlines()
        for line in lines[1:]:
            parts=line.split('\t')
            g=[x for x in parts[1:] if x]
            if g: unassigned.append(g[0]); mem2og.setdefault(g[0],'unassigned')
    return mem2og,ogsizes,unassigned

print(f"{'fam':7} | {'nseq':>4} | mmseqs: {'tot':>3} {'multi':>5} {'singl':>5} | OF: {'OGs':>3} {'in_OG':>5} {'unassigned':>10}")
print('-'*90)
for fam in fams:
    d=f'{BASE}/{fam}'
    idmap=load_idmap(f'{d}/idmap.tsv')
    mem2cl,csizes=parse_mmseqs(f'{d}/cluster_summary.tsv')
    mem2og,ogsizes,unassigned=parse_of(fam)
    nseq=len(idmap)
    mm_tot=len(csizes); mm_multi=sum(1 for s in csizes if s>=2); mm_sing=sum(1 for s in csizes if s==1)
    of_ogs=len(ogsizes); of_inog=sum(ogsizes); of_un=len(unassigned)
    print(f"{fam:7} | {nseq:4} | mmseqs: {mm_tot:3} {mm_multi:5} {mm_sing:5} | OF: {of_ogs:3} {of_inog:5} {of_un:10}")
    # membership table
    out=f'{d}/membership_{fam}.tsv'
    with open(out,'w') as o:
        o.write('clean_id\tmmseqs_cluster\tof_orthogroup\toriginal_header\n')
        for cid in idmap:
            o.write(f'{cid}\t{mem2cl.get(cid,"NA")}\t{mem2og.get(cid,"unassigned")}\t{idmap[cid]}\n')
print('\nmembership tables: <fam>/membership_<fam>.tsv')
print('mmseqs cluster size dist per fam:')
for fam in fams:
    _,cs=parse_mmseqs(f'{BASE}/{fam}/cluster_summary.tsv')
    from collections import Counter
    print(f'  {fam:7}',dict(sorted(Counter(cs).items())))
print('OrthoFinder OG size dist per fam:')
for fam in fams:
    _,og,_=parse_of(fam)
    from collections import Counter
    print(f'  {fam:7}',dict(sorted(Counter(og).items())))
