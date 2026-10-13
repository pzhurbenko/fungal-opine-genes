#!/usr/bin/env python3
"""Aggregate seeds_ext diamond screen (out/<GCA>.tsv, only genomes with out/<GCA>.done).
Outputs: seeds_ext_best.tsv (genome x GENE), seeds_ext_summary.tsv, seeds_ext_candidates.tsv, seeds_ext_hist.tsv.
Best hit per genome x gene = max bitscore over all seeds of the gene; also pident of the best-bitscore hit to bacterial / plant_cTDNA seeds of the gene.
GCA/GCF duplicates (same numeric accession part, version ignored; GCA preferred, then highest version) collapsed (GCA preferred) for summary/candidates (best.tsv has column dup_dropped=0/1 for all).
Candidate rule (as seeds11/step 3): pident>=45, seed coverage>=0.5, pident>=median(best-hit pident of gene)+15.
Contamination flag: step3/dmnd_best.tsv has CAT:vir|CAT:tra|CAT:rep_partition with best_pident>=60."""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import glob,os,csv,re
import numpy as np
D=f'{GMO_ROOT}/analysis/refset/screen'; S=D+'/seeds_ext'
tsv=list(csv.DictReader(open(S+'/seeds_ext.tsv'),delimiter='\t'))
GENES=list(dict.fromkeys(r['gene'] for r in tsv)); SRC={f"{r['gene']}|{r['seed_id']}":r['source'] for r in tsv}
meta={r['gca']:r for r in csv.DictReader(open(D+'/step3/genome_meta.tsv'),delimiter='\t')}
contam=set()
for r in csv.DictReader(open(D+'/step3/dmnd_best.tsv'),delimiter='\t'):
    if r['key'] in ('CAT:vir','CAT:tra','CAT:rep_partition') and float(r['best_pident'])>=60: contam.add(r['gca'])
best={}; genomes=[]
for done in sorted(glob.glob(S+'/out/*.done')):
    g=os.path.basename(done)[:-5]; genomes.append(g)
    for l in open(f'{S}/out/{g}.tsv'):
        q,s,pid,alen,qlen,slen,qs,qe,ss,se,ev,bits=l.rstrip('\n').split('\t'); gene,sid,src=s.split('|')
        bits=float(bits); pid=float(pid); k=(g,gene)
        r=best.setdefault(k,{'gca':g,'gene':gene,'bits':-1,'pident_bact':'','pident_plant':''})
        col='pident_bact' if src=='bacterial' else 'pident_plant'
        if bits>r.get('bits_'+col,-1): r[col]=pid; r['bits_'+col]=bits   # pident of the best-BITSCORE hit within that source
        if bits>r['bits']:
            r.update({'seed_id':sid,'seed_source':src,'protein':q,'pident':pid,'bits':bits,'evalue':ev,
              'cov_seed':round((int(se)-int(ss)+1)/int(slen),3),'cov_query':round((int(qe)-int(qs)+1)/int(qlen),3),'qlen':int(qlen),'contam':int(g in contam)})
num=lambda a:re.sub(r'\.\d+$','',re.sub(r'^GC[AF]_','',a))   # numeric part WITHOUT version
grp={}
for g in genomes: grp.setdefault(num(g),[]).append(g)
ver=lambda a:int(a.rsplit('.',1)[1]) if re.search(r'\.\d+$',a) else 0
keepg={n:sorted(v,key=lambda a:(not a.startswith('GCA_'),-ver(a),a))[0] for n,v in grp.items()}   # prefer GCA, then highest version
keepset=set(keepg.values())
cols=['gca','organism','gene','seed_id','seed_source','protein','pident','bits','evalue','cov_seed','cov_query','qlen','pident_bact','pident_plant','contam','dup_dropped']
def W(name,cols,rows):
    with open(f'{S}/{name}','w') as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',extrasaction='ignore'); w.writeheader(); w.writerows(rows)
allrows=[dict(r,organism=meta.get(r['gca'],{}).get('organism',''),dup_dropped=int(r['gca'] not in keepset)) for r in best.values()]
allrows.sort(key=lambda r:(GENES.index(r['gene']),-r['pident'])); W('seeds_ext_best.tsv',cols,allrows)
rows=[r for r in allrows if not r['dup_dropped']]
summ=[];cand=[];hist=[]; bins=list(range(0,105,5))
for gn in GENES:
    rs=[r for r in rows if r['gene']==gn]; p=np.array([r['pident'] for r in rs])
    med=float(np.median(p)) if len(p) else float('nan')
    c=[r for r in rs if r['pident']>=45 and r['cov_seed']>=0.5 and r['pident']>=med+15]; cand+=c
    summ.append({'gene':gn,'n_seeds':sum(1 for r in tsv if r['gene']==gn),'n_genomes_screened':len(keepset),'n_genomes_hit':len(rs),
      'median_pident':round(med,1),'p95_pident':round(float(np.percentile(p,95)),1) if len(p) else '',
      'p99_pident':round(float(np.percentile(p,99)),1) if len(p) else '','max_pident':float(p.max()) if len(p) else '',
      'n_candidates':len(c),'n_candidates_noncontam':sum(1 for r in c if not r['contam']),'n_contam_genomes_hit':sum(r['contam'] for r in rs)})
    h,_=np.histogram(p,bins=bins+[100.0001]) if len(p) else (np.zeros(len(bins)),0)
    for b,n in zip(bins,h[:len(bins)]): hist.append({'gene':gn,'bin_lo':b,'bin_hi':b+5,'n_genomes':int(n)})
W('seeds_ext_summary.tsv',list(summ[0].keys()),summ)
W('seeds_ext_candidates.tsv',cols,sorted(cand,key=lambda r:(GENES.index(r['gene']),-r['pident'])))
W('seeds_ext_hist.tsv',['gene','bin_lo','bin_hi','n_genomes'],hist)
print('genomes',len(genomes),'after GCA/GCF collapse',len(keepset),'best rows',len(rows),'candidates',len(cand))
