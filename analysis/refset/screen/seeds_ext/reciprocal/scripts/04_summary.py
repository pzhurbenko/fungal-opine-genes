import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..', '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import pandas as pd, numpy as np
R=f'{GMO_ROOT}/analysis/refset/screen/seeds_ext/reciprocal/'
d=pd.read_csv(R+'reciprocal.tsv',sep='\t',keep_default_na=False,na_values=[''])
d['gene_x']=d.gene.str.split(',')
e=d.explode('gene_x').rename(columns={'gene_x':'g1'})
e['grp']=np.where(e.group=='rare','rare','control')
S=[]
for (g,grp),x in e.groupby(['g1','grp']):
    r=dict(gene=g,group=grp,n=len(x),n_no_hits=int((x.n_records==0).sum()))
    for v,nm in(('g','genus'),('o','order')):
        h=x[f'{v}_h']
        r[f'n_h_gt0_{nm}']=int((h>0).sum()); r[f'median_h_{nm}']=h.median()
        for qn,qv in(('min',0),('q1',.25),('q3',.75),('max',1)): r[f'{qn}_h_{nm}']=h.quantile(qv)
        r[f'median_nfung_orders_ge_bact_{nm}']=x[f'{v}_nfung_orders_ge_bact'].median()
    r['median_n_fung_orders_all']=x.n_fung_orders_all.median(); r['median_n_fung_genera_all']=x.n_fung_genera_all.median()
    S.append(r)
s=pd.DataFrame(S); s.to_csv(R+'summary.tsv',sep='\t',index=False,float_format='%.4g')
pd.set_option('display.width',250); pd.set_option('display.max_columns',50); pd.set_option('display.max_rows',500); pd.set_option('display.max_colwidth',60)
print(s[['gene','group','n','n_h_gt0_genus','n_h_gt0_order','median_h_genus','median_h_order','median_nfung_orders_ge_bact_genus','median_nfung_orders_ge_bact_order']].round(1).to_string(index=False))
print(); print(s[['gene','group','min_h_order','q1_h_order','median_h_order','q3_h_order','max_h_order','median_n_fung_orders_all']].round(1).to_string(index=False))
# also control_top / control_random separately
print()
t=d.assign(gx=d.gene.str.split(',')).explode('gx')
t=t[t.group!='rare'].groupby(['gx','group']).agg(n=('qid','size'),h_gt0_order=('o_h',lambda h:(h>0).sum()),h_gt0_genus=('g_h',lambda h:(h>0).sum()),q1=('o_h',lambda h:h.quantile(.25)),med=('o_h','median'),q3=('o_h',lambda h:h.quantile(.75)))
print(t.round(1).to_string())
p=d[d.o_h>0].copy()
print('\nh>0 (order exclusion):',len(p))
p['best_nonfungal']=np.where(p.o_bact_bits.fillna(0)>=p.o_plant_bits.fillna(0),'B:'+p.o_bact_org.astype(str)+f' ',' P:'+p.o_plant_org.astype(str))
def f(r):
    bb=r.o_bact_bits if not np.isnan(r.o_bact_bits) else -1; pb=r.o_plant_bits if not np.isnan(r.o_plant_bits) else -1
    return ('B:' if bb>=pb else 'P:')+(r.o_bact_org if bb>=pb else r.o_plant_org)+' %.1f%%/%.0f'%((r.o_bact_pid,r.o_bact_bits) if bb>=pb else (r.o_plant_pid,r.o_plant_bits))
p['nonfungal']=p.apply(f,axis=1)
p['fung']=p.apply(lambda r:'%s %.1f%%/%s'%(r.o_fung_org,r.o_fung_pid,'%.0f'%r.o_fung_bits) if r.o_fung_bits==r.o_fung_bits else 'none',axis=1)
cols=['gene','group','organism','protein','nonfungal','fung','o_h','o_top50_bact','o_top50_plant','o_top50_fung','o_nfung_orders_ge_bact']
p[cols].sort_values(['gene','o_h'],ascending=[True,False]).to_csv(R+'h_gt0_order.tsv',sep='\t',index=False,float_format='%.4g')
if len(p)<=150: print(p[cols].sort_values(['gene','o_h'],ascending=[True,False]).round(1).to_string(index=False))
want=['KAG8213181.1','KAG6371153.1','KAG6371150.1','KAG8216553.1','KAJ7289093.1','EEU41297.1','CAL9735043.1','KAO1061018.1']
w=d[d.protein.isin(want)]
w.to_csv(R+'specific_proteins.tsv',sep='\t',index=False,float_format='%.4g')
for _,r in w.iterrows():
    print('\n',r.protein,r.organism,r.gene,r.group,'| own order:',r.own_order,'| n_records',r.n_records,'| fung genera/orders all:',r.n_fung_genera_all,r.n_fung_orders_all)
    for v in('g','o'):
        print(' ',v,'bact',r[f'{v}_bact_org'],r[f'{v}_bact_pid'],r[f'{v}_bact_bits'],'| plant',r[f'{v}_plant_org'],r[f'{v}_plant_pid'],r[f'{v}_plant_bits'],'| fung',r[f'{v}_fung_org'],r[f'{v}_fung_pid'],r[f'{v}_fung_bits'],'| h',r[f'{v}_h'],'| top50 B/P/F/other',r[f'{v}_top50_bact'],r[f'{v}_top50_plant'],r[f'{v}_top50_fung'],r[f'{v}_top50_other'],'| nFungOrd>=bact',r[f'{v}_nfung_orders_ge_bact'])
