#!/usr/bin/env python3
"""Aggregate step-3 screening over annotated fungal proteomes.
Per genome x HMM: best protein score; per genome x plasmid-protein category/lineage: best diamond hit.
Outputs: genome_meta.tsv, hmm_best.tsv (long), dmnd_best.tsv (long), summary printed."""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import os,json,glob,csv,re,statistics
from collections import defaultdict,Counter
D=f'{GMO_ROOT}/analysis/refset/screen/step3'
meta={}
for l in open(f'{GMO_ROOT}/data/fungal_proteomes/fungi_annotated/ncbi_dataset/data/assembly_data_report.jsonl'):
    d=json.loads(l); acc=d['accession']; org=d.get('organism',{})
    meta[acc]={'organism':org.get('organismName',''),'taxid':org.get('taxId',''),'level':d.get('assemblyInfo',{}).get('assemblyLevel',''),
               'refseq':d.get('assemblyInfo',{}).get('refseqCategory',''),'nprot':d.get('annotationInfo',{}).get('stats',{}).get('geneCounts',{}).get('proteinCoding','')}
# T-DNA lineage name -> short
def hmm_kind(name):
    if '__' in name: return ('family',name)
    return ('tdna',name)
hmm_rows=[]; dmnd_rows=[]
for done in glob.glob(f'{D}/out/*/done'):
    g=os.path.basename(os.path.dirname(done))
    # HMM: per protein collect all hmm hits; family competitive; tdna best
    per=defaultdict(dict)
    for l in open(f'{D}/out/{g}/hmm.tbl'):
        if l.startswith('#'): continue
        f=l.split(); per[f[0]][f[2]]=(float(f[5]),float(f[4]))
    best_by_hmm=defaultdict(lambda:(0,'',99))
    fam_calls=defaultdict(list)
    for prot,d in per.items():
        for h,(sc,ev) in d.items():
            if sc>best_by_hmm[h][0]: best_by_hmm[h]=(sc,prot,ev)
        for fam in ['ACS','CUS','SUS','ORF358']:
            subs={k:v for k,v in d.items() if k.startswith(fam+'__') and not k.endswith('OUTGROUP')}
            if not subs: continue
            og=d.get(f'{fam}__OUTGROUP',(0,99))[0]; bs=max(subs,key=lambda k:subs[k][0]); sc=subs[bs][0]
            call='OUTGROUP' if og>sc else ('ambiguous' if (sc-og<20 or sc<30) else bs.split('__')[1])
            if call not in ('OUTGROUP',): fam_calls[fam].append((call,sc,round(sc-og,1),prot))
    for h,(sc,prot,ev) in best_by_hmm.items():
        hmm_rows.append({'gca':g,'hmm':h,'best_score':sc,'best_evalue':ev,'protein':prot,'n_prot_hits':sum(1 for p in per if h in per[p])})
    for fam,calls in fam_calls.items():
        for call,sc,gap,prot in calls: hmm_rows.append({'gca':g,'hmm':f'{fam}_CALL_{call}','best_score':sc,'best_evalue':gap,'protein':prot,'n_prot_hits':1})
    # diamond: per plasmid category and per gene/product lineage: best pident
    best=defaultdict(lambda:(0,0,'',''))
    for l in open(f'{D}/out/{g}/dmnd.tsv'):
        q,s,pid,alen,qlen,slen,ev,bits=l.rstrip('\n').split('\t')[:8]
        p=s.split('|'); cat=p[1]; lineage=(p[4] if p[4]!='-' else p[5][:30])
        cov=int(alen)/max(int(qlen),int(slen))
        for key in (f'CAT:{cat}',f'LIN:{cat}:{lineage}'):
            if float(bits)>best[key][1]: best[key]=(float(pid),float(bits),q,f'{cov:.2f}')
    for key,(pid,bits,q,cov) in best.items():
        dmnd_rows.append({'gca':g,'key':key,'best_pident':pid,'best_bits':bits,'cov':cov,'protein':q})
with open(f'{D}/genome_meta.tsv','w') as f:
    w=csv.DictWriter(f,fieldnames=['gca','organism','taxid','level','refseq','nprot'],delimiter='\t'); w.writeheader()
    for g,m in meta.items(): w.writerow({'gca':g}|m)
for name,rows in [('hmm_best.tsv',hmm_rows),('dmnd_best.tsv',dmnd_rows)]:
    with open(f'{D}/{name}','w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0].keys()),delimiter='\t'); w.writeheader(); w.writerows(rows)
print('genomes processed:',len({r['gca'] for r in hmm_rows}|{r['gca'] for r in dmnd_rows}),' hmm rows',len(hmm_rows),' dmnd rows',len(dmnd_rows))
