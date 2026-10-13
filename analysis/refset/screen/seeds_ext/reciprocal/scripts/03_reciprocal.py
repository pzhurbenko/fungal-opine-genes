import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..', '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import sys, json, re, collections
import numpy as np, pandas as pd
sys.path.insert(0,f'{GMO_ROOT}/analysis/refset/screen/seeds_ext/reciprocal/scripts')
from importlib import import_module
R=f'{GMO_ROOT}/analysis/refset/screen/seeds_ext/reciprocal/'
BR=re.compile(r'\[([^\[\]]+)\]')
tax=json.load(open(R+'organism_taxonomy.json'))
q=pd.read_csv(R+'queries.tsv',sep='\t')
cols=['qseqid','sseqid','pident','length','qlen','slen','evalue','bitscore','stitle']
h=pd.read_csv(R+'diamond_nr.tsv',sep='\t',names=cols,header=None,quoting=3,dtype={'stitle':str})
h['stitle']=h.stitle.fillna('')
# expand: one record per [organism] in stitle (merged identical nr entries)
recs=[]
n_rows_multi=0
for qs,pid,bits,st in zip(h.qseqid,h.pident,h.bitscore,h.stitle):
    os_=list(dict.fromkeys(o.strip() for o in BR.findall(st)))
    if len(os_)>1: n_rows_multi+=1
    if not os_: os_=['']
    for o in os_: recs.append((qs,o,pid,bits))
r=pd.DataFrame(recs,columns=['qid','org','pident','bits'])
t=lambda o,k: tax.get(o,{}).get(k,'') if o else ''
uo=r.org.unique()
mp={o:(tax[o]['kingdom'] if o in tax else 'NA', tax[o]['genus'] if o in tax else '', tax[o]['order'] if o in tax else '') for o in uo}
r['king']=r.org.map(lambda o:mp[o][0]); r['genus']=r.org.map(lambda o:mp[o][1]); r['order']=r.org.map(lambda o:mp[o][2])
r.loc[r.org=='','king']='NA'
st=dict(rows=len(h),records=len(r),rows_multi_org=n_rows_multi,records_no_bracket=int((r.org=='').sum()),
        records_king_NA=int((r.king=='NA').sum()),frac_king_NA=float((r.king=='NA').mean()),
        frac_NA_incl_nobracket=float((r.king=='NA').mean()),
        king_counts=r.king.value_counts().to_dict(),
        fungal_records_no_order=int(((r.king=='Fungi')&(r.order=='')).sum()),
        n_fungal_records=int((r.king=='Fungi').sum()),
        unique_orgs=len(uo),unresolved_orgs=int(sum(1 for o in uo if mp[o][0]=='NA')))
# cross-check with pipeline cache
cache=json.load(open(f'{GMO_ROOT}/analysis/loci_pipeline/registry/taxonomy_cache.json'))
m={'Plants':'Viridiplantae','other_eukaryote':'Eukaryota_other'}
agree=dis=0;ex=[]
for o in uo:
    c=re.sub(r'\s+(strain|str\.|isolate|sp\.).*$','',o).strip()
    if c in cache and mp[o][0]!='NA':
        cv=m.get(cache[c],cache[c])
        if cv==mp[o][0] or (cv in('Metazoa','Eukaryota_other') and mp[o][0] in('Metazoa','Eukaryota_other')): agree+=1
        else: dis+=1; ex.append((o,cache[c],mp[o][0]))
st['cache_crosscheck_agree']=agree; st['cache_crosscheck_disagree']=dis; st['cache_disagree_examples']=ex[:8]
json.dump(st,open(R+'taxonomy_stats.json','w'),indent=1,default=str)
print(json.dumps(st,indent=1,default=str))
r.sort_values(['qid','bits'],ascending=[True,False],kind='stable',inplace=True)
groups={k:v for k,v in r.groupby('qid',sort=False)}
def best(x,kg):
    y=x[x.king==kg]
    if len(y)==0: return ('',np.nan,np.nan)
    z=y.iloc[0]; return (z.org,z.pident,z.bits)
rows=[]
for _,Q in q.iterrows():
    x=groups.get(Q.qid)
    o_org=Q.organism; own=tax.get(o_org,{}); og=own.get('genus') or o_org.split()[0]; oo=own.get('order','')
    row=dict(qid=Q.qid,protein=Q.protein,gca=Q.gca,organism=o_org,gene=Q.gene,group=Q.group,pident_seed=Q.pident,cov_seed=Q.cov_seed,
             own_genus=og,own_class=own.get('cls',''),own_order=oo,n_records=0 if x is None else len(x))
    if x is not None:
        f=x[x.king=='Fungi']
        row['n_fung_genera_all']=f.genus[f.genus!=''].nunique(); row['n_fung_orders_all']=f.order[f.order!=''].nunique()
    else: row['n_fung_genera_all']=row['n_fung_orders_all']=0
    for v in('g','o'):
        if x is None: continue
        if v=='g': y=x[~((x.king=='Fungi')&(x.genus==og))]
        else:
            y=x[~((x.king=='Fungi')&(x.order==oo))] if oo else x[~((x.king=='Fungi')&(x.genus==og))]
        for kg,nm in(('Bacteria','bact'),('Viridiplantae','plant'),('Fungi','fung')):
            a,b,c=best(y,kg); row[f'{v}_{nm}_org']=a; row[f'{v}_{nm}_pid']=b; row[f'{v}_{nm}_bits']=c
        nf=np.nanmax([row[f'{v}_bact_bits'],row[f'{v}_plant_bits']]) if not (np.isnan(row[f'{v}_bact_bits']) and np.isnan(row[f'{v}_plant_bits'])) else 0.0
        fb=0.0 if np.isnan(row[f'{v}_fung_bits']) else row[f'{v}_fung_bits']
        row[f'{v}_h']=nf-fb
        t50=y.head(50).king.value_counts()
        row[f'{v}_top50_bact']=int(t50.get('Bacteria',0)); row[f'{v}_top50_plant']=int(t50.get('Viridiplantae',0)); row[f'{v}_top50_fung']=int(t50.get('Fungi',0))
        row[f'{v}_top50_other']=int(len(y.head(50))-row[f'{v}_top50_bact']-row[f'{v}_top50_plant']-row[f'{v}_top50_fung'])
        bb=0.0 if np.isnan(row[f'{v}_bact_bits']) else row[f'{v}_bact_bits']
        fy=y[(y.king=='Fungi')&(y.bits>=bb)]
        row[f'{v}_nfung_orders_ge_bact']=fy.order[fy.order!=''].nunique()
    rows.append(row)
out=pd.DataFrame(rows)
out.to_csv(R+'reciprocal.tsv',sep='\t',index=False,float_format='%.4g')
print(len(out),'queries; without any hit:',int((out.n_records==0).sum()),'; own_order unknown:',int((out.own_order=='').sum()))
