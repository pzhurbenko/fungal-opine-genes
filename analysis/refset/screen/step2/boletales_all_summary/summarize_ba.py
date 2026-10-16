"""Summary of step-2 screen for Boletales_all (376 GCA). Does NOT touch shared step2_*.tsv. Rules replicate aggregate_step2.py."""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import re,os,csv
from collections import defaultdict,Counter
D=f'{GMO_ROOT}/analysis/refset/screen/step2'; OUT=f'{D}/boletales_all_summary'
ACC=f'{GMO_ROOT}/data/fungal_genomes/boletales_all/accessions_20261006.tsv'
CAR=f'{GMO_ROOT}/analysis/funnel/carriers_tree_based.tsv'
GAP=5000; FAMS=['ACS','CUS','SUS']; SUBS=['sus','cus','gdpd_meso','kerr14','acs','mis','cusrel181','rcorf12']
CONTAM_PID=60.0; BACK=('vir','tra','rep_partition','mobile')
def ps(s):
    p=s.split('|')
    if p[0].startswith('FAM_'): return dict(cat=p[0][4:],lin=p[1])
    return dict(cat=p[1],lin=p[4] if len(p)>4 and p[4]!='-' else '')
def scov_of(iv,slen):
    mg=[]
    for a,e in sorted(iv):
        if mg and a<=mg[-1][1]: mg[-1][1]=max(mg[-1][1],e)
        else: mg.append([a,e])
    return sum(e-a+1 for a,e in mg)/slen
meta={}
for r in csv.DictReader(open(ACC),delimiter='\t'): meta[r['accession']]=r
paths={l.split('\t')[0]:l.split('\t')[1] for l in open(f'{D}/boletales_all.tsv')}
def binom(n): 
    w=re.sub(r'[\[\]]','',n).split(); return ' '.join(w[:2]).lower()
known=defaultdict(lambda: defaultdict(list))
for r in csv.DictReader(open(CAR),delimiter='\t'):
    ln=r['lineage']; key='sus' if ln.startswith('SUS') else 'cus' if ln.startswith('CUS') else 'gdpd_meso' if 'gdpd_meso' in ln else 'kerr14' if 'kerr14' in ln else None
    if key: known[key][binom(r['organism_ncbi'])].append((r['organism_ncbi'],r['gca']))
def fasta_lens(fna,need):
    L={};cur=None
    for l in open(fna):
        if l[0]=='>':
            cur=l[1:].split()[0]; L[cur]=0
        else: L[cur]+=len(l.strip())
    return L
calls=[];gen=[];strongs=[];notdone=[]
for g,m in meta.items():
    od=f'{D}/out/{g}'
    if not os.path.exists(f'{od}/done'): notdone.append(g); continue
    hs=defaultdict(list); clen={}
    for l in open(f'{od}/dmnd.tsv'):
        f=l.rstrip('\n').split('\t')
        if len(f)<13: continue
        q=f[0]; clen[q]=int(f[4]); qs,qe=sorted((int(f[6]),int(f[7])))
        hs[(q,'+' if int(f[12])>0 else '-',f[1])].append((qs,qe,int(f[8]),int(f[9]),float(f[2]),float(f[11]),int(f[5])))
    chains=defaultdict(list)
    for (q,st,s),hl in hs.items():
        hl.sort(); cur=None
        for qs,qe,ss,se,pid,bits,slen in hl:
            if cur and qs-cur['end']<=GAP:
                cur['end']=max(cur['end'],qe);cur['bits']+=bits;cur['scov'].append((min(ss,se),max(ss,se)));cur['pid']=max(cur['pid'],pid)
            else:
                if cur: chains[(q,st)].append(cur)
                cur=dict(start=qs,end=qe,subj=s,bits=bits,scov=[(min(ss,se),max(ss,se))],pid=pid,slen=slen)
        chains[(q,st)].append(cur)
    contig_back=Counter();contig_back_strict=Counter();nstrong=0
    chain_by_contig=defaultdict(list)
    for (q,st),cl in chains.items():
        cl.sort(key=lambda c:c['start']); loci=[];cur=None
        for c in cl:
            c['sc']=scov_of(c['scov'],c['slen']); chain_by_contig[q].append((c['start'],c['end'],st,c))
            if cur and c['start']<=cur['end']: cur['end']=max(cur['end'],c['end']);cur['m'].append(c)
            else:
                if cur: loci.append(cur)
                cur=dict(start=cur_s if False else c['start'],end=c['end'],m=[c])
        loci.append(cur)
        for L in loci:
            b=max(L['m'],key=lambda x:x['bits']); p=ps(b['subj'])
            if p['cat'] in BACK and b['pid']>=CONTAM_PID:
                contig_back[q]+=1
                if b['sc']>=0.5: contig_back_strict[q]+=1
            if b['pid']>=45 and b['sc']>=0.5 and p['cat'] in('T1','T2','TDNA','opine_related'):
                nstrong+=1; strongs.append(dict(gca=g,organism=m['organism'],contig=q,contig_len=clen.get(q,''),start=L['start'],end=L['end'],strand=st,subject=b['subj'],cat=p['cat'],lineage=p['lin'],pident=f"{b['pid']:.1f}",subj_cov=f"{b['sc']:.2f}"))
    coord={}
    for l in open(f'{od}/hit_orfs.faa'):
        mm=re.match(r'>(\S+) source=(\S+) coords=(\d+)\.\.(\d+) length=(\d+) frame=(\d+)',l)
        if mm: a,b=int(mm[3]),int(mm[4]); coord[mm[1]]=(mm[2],min(a,b),max(a,b),'+' if int(mm[6])<=3 else '-')
    per=defaultdict(dict)
    for l in open(f'{od}/hmm.tbl'):
        if l[0]=='#': continue
        f=l.split(); per[f[0]][f[2]]=(float(f[5]),float(f[4]))
    cnt=Counter(); amb=0; og_n=0; gc=[]
    for orf,d in per.items():
        for fam in FAMS:
            subs={k:v for k,v in d.items() if k.startswith(fam+'__') and not k.endswith('OUTGROUP')}
            if not subs: continue
            og=d.get(f'{fam}__OUTGROUP',(0,99))[0]
            bs=max(subs,key=lambda k:subs[k][0]); sc=subs[bs][0]; gap=sc-og
            if og>sc: og_n+=1
            elif gap<20 or sc<30: amb+=1
            else: gc.append((fam,bs.split('__')[1],sc,gap,orf))
    need=[coord[o][0] for *_,o in gc if o in coord]
    flen=fasta_lens(paths[g],need) if gc else {}
    for fam,sub,sc,gap,orf in gc:
        cnt[sub]+=1
        c,a,b,st=coord.get(orf,('?',0,0,'?'))
        best=None
        for qs,qe,cst,ch in chain_by_contig.get(c,[]):
            if min(b,qe)-max(a,qs)>=0:
                if best is None or ch['bits']>best[0]['bits']: best=(ch,cst)
        row=dict(gca=g,organism=m['organism'],release_date=m['release_date'],contig=c,contig_len=flen.get(c,clen.get(c,'')),start=a,end=b,strand=st,family=fam,subgroup=sub,score=f'{sc:.0f}',gap=f'{gap:.0f}',orf=orf)
        if best:
            ch,cst=best; p=ps(ch['subj'])
            row.update(dmnd_subject=ch['subj'],dmnd_cat=p['cat'],dmnd_lineage=p['lin'],dmnd_pident=f"{ch['pid']:.1f}",dmnd_scov=f"{ch['sc']:.2f}",dmnd_strand=cst)
        else: row.update(dmnd_subject='',dmnd_cat='',dmnd_lineage='',dmnd_pident='',dmnd_scov='',dmnd_strand='')
        calls.append(row)
    ncall=sum(cnt.values())
    gen.append(dict(gca=g,organism=m['organism'],release_date=m['release_date'],contig_n50=m['contig_n50'],n_calls=ncall,**{s:cnt[s] for s in SUBS},n_ambiguous=amb,n_outgroup=og_n,
        n_strong_diamond_loci=nstrong,backbone_contigs_agg_rule=len(contig_back),backbone_contigs_strict_cov05=len(contig_back_strict),
        contam_flag='yes' if len(contig_back_strict)>0 else 'no'))
# known species comparison
for r in calls:
    s=r['subgroup']
    if s in known:
        k=known[s]; b=binom(r['organism'])
        ex=[o for o,_ in sum(k.values(),[]) if o.lower()==r['organism'].lower()]
        r['known_same_lineage']='exact_name' if ex else 'same_binomial' if b in k else 'new_species_name'
        r['known_names']=';'.join(sorted({o for o,_ in k.get(b,[])}))
    else: r['known_same_lineage']='';r['known_names']=''
def w(name,rows):
    with open(f'{OUT}/{name}','w') as f:
        if not rows: return
        ks=list(rows[0].keys()); f.write('\t'.join(ks)+'\n')
        for r in rows: f.write('\t'.join(str(r.get(k,'')) for k in ks)+'\n')
calls.sort(key=lambda r:(r['subgroup'],-float(r['score'])))
w('calls.tsv',calls); w('genomes.tsv',gen); w('strong_diamond_loci.tsv',strongs)
print('genomes done',len(gen),'not done',notdone)
for s in SUBS:
    print(s,'genomes',len({r['gca'] for r in calls if r['subgroup']==s}),'calls',sum(1 for r in calls if r['subgroup']==s))
print('no confident call',sum(1 for g in gen if g['n_calls']==0),'contam',sum(1 for g in gen if g['contam_flag']=='yes'),'contam agg>=3',sum(1 for g in gen if g['backbone_contigs_agg_rule']>=3),'agg>=1',sum(1 for g in gen if g['backbone_contigs_agg_rule']>=1))
