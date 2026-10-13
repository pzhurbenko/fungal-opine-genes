#!/usr/bin/env python3
"""Build seeds_ext.faa/.tsv: seeds11 seed + T1 references (loci pipeline guides) / T-DNA cluster members (plasmid catalog).
Dedup: exact sequence duplicates within gene, then 95% identity collapse (diamond all-vs-all, greedy; seed11 seed always kept)."""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import csv,re,subprocess,os,sys,collections
S=f'{GMO_ROOT}/analysis/refset/screen/seeds_ext'; G=f'{GMO_ROOT}/analysis/loci_pipeline/guides'
P=f'{GMO_ROOT}/analysis/refset/plasmids'; S11=f'{GMO_ROOT}/analysis/refset/screen/seeds11'
def fa(path):
    d={};k=None
    for l in open(path):
        l=l.rstrip()
        if l.startswith('>'): k=l[1:]; d[k]=[]
        elif k: d[k].append(l)
    return {k:''.join(v).replace('*','') for k,v in d.items()}
seeds11=fa(S11+'/seeds11.faa'); meta11={r['gene']:r for r in csv.DictReader(open(S11+'/seeds11.tsv'),delimiter='\t')}
s11={k.split('|')[0]:(k,v) for k,v in seeds11.items()}
recs=collections.defaultdict(list)  # gene -> [(id,source,origin,seq,is_seed)]
for g,(k,v) in s11.items(): recs[g].append((meta11[g]['accession'],'bacterial',f'seeds11 seed: {meta11[g]["organism_plasmid"]}',v,True))
# T1 refs from guides
sub={'CUS':{'cus':'cus','mis':'mis'},'SUS':{'sus':'sus'},'ACS':{'acs':'acs','kerr14':'acs_pKerr14a'}}
for fam,m in sub.items():
    hs={r['header']:r['subgroup'] for r in csv.DictReader(open(f'{G}/{fam}_subgroups.tsv'),delimiter='\t') if r['tier']=='T1'}
    for h,seq in fa(f'{G}/{fam}_T1.faa').items():
        sg=next((v for k,v in hs.items() if h==k or h.startswith(k) or k.startswith(h)),None)
        if sg not in m: continue
        src='plant_cTDNA' if h.startswith('T1|Plant_') else 'bacterial' if h.startswith('T1|Agro_') else None
        assert src,h
        acc=re.sub(r'^T1\|(Plant|Agro)_','',h)
        acc=acc.split('_')[0]+('_'+acc.split('_')[1] if acc.startswith(('WP_','MFA','MBA')) else '')
        m2=re.search(r'_([A-Z]{4,6}\d{8,})_(\d+)_',h)
        if m2 and 'shotgun' in h: acc=h.split('Plant_')[1].split('_')[0]+'_'+m2.group(1)+'_'+m2.group(2)
        acc=re.sub(r'[^A-Za-z0-9_.-]','',acc)
        recs[m[sg]].append((acc,src,f'T1 {fam}/{sg}: {h}',seq,False))
# T-DNA clusters
tp=fa(P+'/tdna_proteins_all.faa'); byacc=collections.defaultdict(list)
for k,v in tp.items(): byacc[k.split('|')[0]].append((k,v))
cl={'mas1':'AAF77131.1','mas2':'WP_262512305.1','ags':'WP_080864006.1','nos':'WP_010974826.1','ocs':'AAF77127.1'}
mem=collections.defaultdict(set)
for r in csv.reader(open(P+'/tdna_clu_cluster.tsv'),delimiter='\t'): mem[r[0].split('|')[0]].add(r[1])
for g,rep in cl.items():
    for hd in sorted(mem[rep]):
        a=hd.split('|')[0]; seq=tp.get(hd) or byacc[a][0][1]
        recs[g].append((a,'bacterial',f'tdna cluster {rep}: {hd}',seq,False))
GENES=['cus','mis','sus','acs','acs_pKerr14a','mas1','mas2','ags','nos','ocs','vis','chs']
rows=[];out=[];stat=[]
for g in GENES:
    print(g,len(recs[g]),file=sys.stderr)
    n0=len(recs[g]); seen={}; ex=[]
    for r in recs[g]:  # exact dups: keep first (seed first)
        if r[3] in seen: continue
        seen[r[3]]=1; ex.append(r)
    n1=len(ex)
    # unique ids
    ids=collections.Counter(); uniq=[]
    for r in ex:
        ids[r[0]]+=1; i=r[0] if ids[r[0]]==1 else f'{r[0]}.v{ids[r[0]]}'; uniq.append((i,)+r[1:])
    tmp=f'{S}/tmp_{g}.faa'; open(tmp,'w').write(''.join(f'>{r[0]}\n{r[3]}\n' for r in uniq))
    keep=set()
    if len(uniq)>1:
        subprocess.run(f'diamond makedb --in {tmp} -d {tmp}.db --quiet && diamond blastp -q {tmp} -d {tmp}.db -o {tmp}.out --very-sensitive -e 1e-5 --max-target-seqs 1000 --outfmt 6 qseqid sseqid pident length qlen slen --quiet',shell=True,check=True,executable='/bin/bash')
        hit=collections.defaultdict(set)
        for l in open(tmp+'.out'):
            q,s,pid,al,ql,sl=l.split('\t'); 
            if float(pid)>=95 and int(al)/min(int(ql),int(sl))>=0.9 and int(al)/max(int(ql),int(sl))>=0.9: hit[q].add(s)
        order=sorted(uniq,key=lambda r:(not r[4],-len(r[3])))  # seed first, then longest
        removed=set()
        for r in order:
            if r[0] in removed: continue
            keep.add(r[0])
            for s in hit[r[0]]:
                if s!=r[0] and s not in keep: removed.add(s)
    else: keep={uniq[0][0]}
    for f in (tmp,tmp+'.db.dmnd',tmp+'.out'):
        if os.path.exists(f): os.remove(f)
    fin=[r for r in uniq if r[0] in keep]
    stat.append((g,n0,n1,len(fin)))
    for r in fin:
        sid=r[0]; out.append(f'>{g}|{sid}|{r[1]}\n{r[3]}\n'); rows.append((g,sid,r[1],r[2],len(r[3])))
open(S+'/seeds_ext.faa','w').write(''.join(out))
with open(S+'/seeds_ext.tsv','w') as f:
    f.write('gene\tseed_id\tsource\torigin\tlength\n'); [f.write('\t'.join(map(str,r))+'\n') for r in rows]
for s in stat: print('gene raw=%d after_exact=%d after_95=%d'%s[1:],s[0])
c=collections.Counter((r[0],r[2]) for r in rows)
for g in GENES: print(g,c[(g,'bacterial')],c[(g,'plant_cTDNA')])
print('short<150:',[(r[0],r[1],r[4]) for r in rows if r[4]<150])
