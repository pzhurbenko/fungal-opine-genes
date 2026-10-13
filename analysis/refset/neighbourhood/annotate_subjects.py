#!/usr/bin/env python3
"""Annotate nr-neighbourhood subjects: kingdom/lineage via NCBI Taxonomy (cached), replicon via IPG for bacterial subjects."""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import sys,re,csv,json,os
sys.path.insert(0,f'{GMO_ROOT}/pipeline/scripts'); from common import load_config, eutils, log
cfg=load_config()
rows=list(csv.DictReader(open('subjects.tsv'),delimiter='\t'))
cache_p=f'{GMO_ROOT}/analysis/loci_pipeline/registry/taxonomy_cache.json'
cache=json.load(open(cache_p)) if os.path.exists(cache_p) else {}
lin_cache_p='lineage_cache.json'; lin=json.load(open(lin_cache_p)) if os.path.exists(lin_cache_p) else {}
def lineage(org):
    key=re.sub(r'\s+(strain|str\.|isolate|sp\.|subsp\.).*$','',org).strip()
    if key in lin: return lin[key]
    out={'kingdom':'other','lineage':''}
    try:
        ids=eutils(cfg,'esearch.fcgi',{'db':'taxonomy','term':key,'retmode':'json'},as_json=True)['esearchresult']['idlist']
        if not ids and ' ' in key: ids=eutils(cfg,'esearch.fcgi',{'db':'taxonomy','term':key.split()[0],'retmode':'json'},as_json=True)['esearchresult']['idlist']
        if ids:
            xml=eutils(cfg,'efetch.fcgi',{'db':'taxonomy','id':ids[0],'retmode':'xml'})
            L=re.search(r'<Lineage>(.*?)</Lineage>',xml,re.S); L=L.group(1) if L else ''
            parts=L.split('; ')
            k='Fungi' if 'Fungi' in parts else 'Bacteria' if 'Bacteria' in parts else 'Plants' if 'Viridiplantae' in parts else 'Archaea' if 'Archaea' in parts else 'Metazoa' if 'Metazoa' in parts else ('other_eukaryote' if 'Eukaryota' in parts else 'other')
            out={'kingdom':k,'lineage':L}
    except Exception as e: log('tax fail',org,e)
    lin[key]=out; return out
orgs=sorted({r['organism'] for r in rows if r['organism']}); log('organisms',len(orgs))
for i,o in enumerate(orgs):
    lineage(o)
    if (i+1)%100==0: json.dump(lin,open(lin_cache_p,'w')); log('tax',i+1)
json.dump(lin,open(lin_cache_p,'w'))
for r in rows:
    L=lineage(r['organism']) if r['organism'] else {'kingdom':'?','lineage':''}
    r['kingdom']=L['kingdom']; parts=L['lineage'].split('; ')
    r['class_or_order']=next((p for p in parts if p.endswith(('mycetes','ales'))),'')
    r['is_rhizobiales']='yes' if 'Hyphomicrobiales' in parts or 'Rhizobiales' in parts else 'no'
# replicon via IPG for bacterial WP_/non-WP subjects (Rhizobiales first, then others), cap to keep runtime sane
ipg_p='ipg_cache.json'; ipg=json.load(open(ipg_p)) if os.path.exists(ipg_p) else {}
bact=[r for r in rows if r['kingdom']=='Bacteria']
bact.sort(key=lambda r:(r['is_rhizobiales']!='yes',-float(r['bits'])))
todo=[r['sseqid'] for r in bact if r['sseqid'] not in ipg][:2500]; log('IPG lookups',len(todo))
for i in range(0,len(todo),20):
    chunk=todo[i:i+20]
    try:
        txt=eutils(cfg,'efetch.fcgi',{'db':'protein','id':','.join(chunk),'rettype':'ipg','retmode':'text'})
    except Exception as e: log('ipg fail',e); continue
    cur=None
    for l in txt.strip().split('\n')[1:]:
        f=l.split('\t')
        if len(f)<9: continue
        acc=f[6]  # protein accession in this row
        for c in chunk:
            if c.split('.')[0]==acc.split('.')[0]:
                ipg.setdefault(c,[]).append({'nuc':f[2],'start':f[3],'end':f[4],'org':f[8],'src':f[1]})
    for c in chunk: ipg.setdefault(c,[])
    if (i//20)%20==0: json.dump(ipg,open(ipg_p,'w')); log('ipg',i+20)
json.dump(ipg,open(ipg_p,'w'))
def replicon(recs):
    if not recs: return 'unknown',''
    nucs=[x['nuc'] for x in recs]
    # heuristics from accession prefix: NZ_CP/CP complete -> need title; use nuc accession patterns
    tags=[]
    for n in nucs[:5]:
        if n.startswith(('NC_','NZ_CP','CP','AE','AP','CU','LT','LR','LN','FN','FO','FP','FQ','HE','HG')): tags.append('complete_replicon')
        else: tags.append('wgs_contig')
    return ('complete_replicon' if 'complete_replicon' in tags else 'wgs_contig'), ';'.join(nucs[:3])
for r in rows:
    rep,nuc=replicon(ipg.get(r['sseqid'],[])) if r['kingdom']=='Bacteria' else ('','')
    r['replicon_type']=rep; r['nuc_accs']=nuc
with open('subjects_annotated.tsv','w') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys()),delimiter='\t'); w.writeheader(); w.writerows(rows)
from collections import Counter
log('kingdom x fam:',dict(Counter((r['fam'],r['kingdom']) for r in rows)))
log('DONE')
