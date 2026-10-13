import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import csv,os,re,math,glob
from collections import defaultdict,Counter
BASE=f'{GMO_ROOT}'
CL=f'{BASE}/analysis/clustering'
fams=['ACS','CUS','ORF358','SUS']

def load_membership(fam):
    m={}
    for line in open(f'{CL}/{fam}/cluster_summary.tsv'):
        pass
    with open(f'{CL}/{fam}/cluster_summary.tsv') as f:
        next(f)
        for line in f:
            cid,rep,size,members=line.rstrip('\n').split('\t')
            for mm in members.split(','): m[mm]=cid
    return m
def load_idmap(fam):
    d={}
    with open(f'{CL}/{fam}/idmap.tsv') as f:
        next(f)
        for line in f:
            a=line.rstrip('\n').split('\t'); d[a[0]]=a[1] if len(a)>1 else ''
    return d
def load_clade(fam):
    cand=[f'{BASE}/data/tsv/rep_seqs/{fam}_rep_seqs.tsv',
          f'{BASE}/data/tsv/rep_seqs/{fam.lower()}_rep_seqs.tsv']
    p=next(x for x in cand if os.path.exists(x))
    r=list(csv.reader(open(p),delimiter='\t')); h=r[0]
    ia=h.index('Subj_Accession'); it=h.index('Tree_clade')
    acc2cl={}
    for row in r[1:]:
        if len(row)<=max(ia,it): continue
        acc=row[ia].strip(); cl=row[it].strip()
        if acc: acc2cl[acc]=cl
    return acc2cl

# ---- metrics ----
def ari(labels_true,labels_pred):
    # contingency
    from collections import Counter
    ab=Counter(zip(labels_true,labels_pred))
    a=Counter(labels_true); b=Counter(labels_pred); n=len(labels_true)
    comb2=lambda x: x*(x-1)//2
    sum_ij=sum(comb2(v) for v in ab.values())
    sum_a=sum(comb2(v) for v in a.values()); sum_b=sum(comb2(v) for v in b.values())
    exp=sum_a*sum_b/comb2(n) if n>1 else 0
    mx=(sum_a+sum_b)/2
    return (sum_ij-exp)/(mx-exp) if (mx-exp)!=0 else 1.0
def entropy(labels):
    n=len(labels); c=Counter(labels)
    return -sum((v/n)*math.log(v/n) for v in c.values()) if n else 0
def mutual_info(t,p):
    n=len(t); ab=Counter(zip(t,p)); a=Counter(t); b=Counter(p); mi=0
    for (x,y),nij in ab.items():
        mi+=(nij/n)*math.log((nij*n)/(a[x]*b[y]))
    return mi
def vmeasure(t,p):
    Ht=entropy(t); Hp=entropy(p); I=mutual_info(t,p)
    hom=I/Ht if Ht>0 else 1.0      # homogeneity: each cluster -> single clade
    comp=I/Hp if Hp>0 else 1.0     # completeness: each clade -> single cluster
    v=2*hom*comp/(hom+comp) if (hom+comp)>0 else 0
    nmi=I/math.sqrt(Ht*Hp) if Ht>0 and Hp>0 else 1.0
    return hom,comp,v,nmi

print(f"{'fam':7} {'matched':>7} {'unmatch_fungal':>14} {'refs(Plant/Agro)':>16} | {'ARI':>6} {'homog':>6} {'compl':>6} {'V':>6} {'NMI':>6}")
print('-'*100)
allrows={}
for fam in fams:
    memb=load_membership(fam); idmap=load_idmap(fam); acc2cl=load_clade(fam)
    accs=sorted(acc2cl, key=len, reverse=True)  # longer first to avoid partial
    pairs=[]  # (clade, cluster)
    ref=0; unmatched=[]
    detail=[]
    for cid_clean,cluster in memb.items():
        hdr=idmap.get(cid_clean,cid_clean)
        if hdr.startswith('Plant_') or hdr.startswith('Agro_'):
            ref+=1; continue
        # find clade by accession substring
        found=None
        for acc in accs:
            base=re.sub(r'\.\d+$','',acc)
            if acc in hdr or base in hdr:
                found=acc2cl[acc]; break
        if found is None: unmatched.append(hdr[:60])
        else:
            pairs.append((found,cluster)); detail.append((found,cluster,hdr))
    if not pairs:
        print(f"{fam:7} no pairs"); continue
    t=[x[0] for x in pairs]; p=[x[1] for x in pairs]
    a=ari(t,p); hom,comp,v,nmi=vmeasure(t,p)
    print(f"{fam:7} {len(pairs):7} {len(unmatched):14} {ref:16} | {a:6.3f} {hom:6.3f} {comp:6.3f} {v:6.3f} {nmi:6.3f}")
    allrows[fam]=(detail,unmatched)
    # crosstab
    ct=defaultdict(Counter)
    for cl,clu in pairs: ct[cl][clu]+=1
    print(f"   Tree_clade -> mmseqs clusters (n):")
    for cl in sorted(ct):
        items=', '.join(f'{k}:{v}' for k,v in sorted(ct[cl].items(),key=lambda x:-x[1]))
        pure='  [PURE]' if len(ct[cl])==1 else ''
        print(f"     {cl:14} (n={sum(ct[cl].values()):3}) -> {items}{pure}")
    if unmatched:
        print(f"   [!] unmatched fungal reps ({len(unmatched)}): {unmatched[:3]}")
    print()
