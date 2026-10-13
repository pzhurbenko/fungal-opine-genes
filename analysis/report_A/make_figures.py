import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
#!/usr/bin/env python3
import os, csv, glob, re, statistics
os.environ['QT_QPA_PLATFORM']='offscreen'
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt; import matplotlib.patches as mp
import numpy as np, pandas as pd
from ete3 import Tree, TreeStyle, NodeStyle, TextFace, CircleFace, RectFace, faces
A=f'{GMO_ROOT}/analysis'; D='data'; F='fig'
plt.rcParams.update({'font.size':10,'font.family':'DejaVu Sans'})
lt={r['locus_id']:r for r in csv.DictReader(open(f'{A}/loci_pipeline/round2/loci_table.tsv'),delimiter='\t')}
clu={}; cname={}
for r in csv.DictReader(open(f'{D}/clusters.tsv'),delimiter='\t'):
    for g in r['genomes'].split(';'): clu[g]=r['cluster_id']; cname[g]=r['cluster_name']
car=list(csv.DictReader(open(f'{D}/carriers_tree_based.tsv'),delimiter='\t'))
LIN=['SUS T-DNA clade','CUS cus clade','ACS acc-locus (gdpd_meso) clade','ACS kerr14/acc-like ("six") lineage']
LCOL={'SUS T-DNA clade':'#d62728','CUS cus clade':'#1f77b4','ACS acc-locus (gdpd_meso) clade':'#2ca02c','ACS kerr14/acc-like ("six") lineage':'#9467bd'}
LSHORT={'SUS T-DNA clade':'SUS','CUS cus clade':'CUS','ACS acc-locus (gdpd_meso) clade':'ACS-acc','ACS kerr14/acc-like ("six") lineage':'ACS-k14'}
pres={}
for r in car:
    pres.setdefault(r['gca'],{}).setdefault(r['lineage'],set()).add(r['status'])

# ---------- F1: Boletales mash tree with lineage presence ----------
t=Tree(f'{A}/boletales_mash/boletales.mashtree.nwk')
gmeta={}
for g,org,genus,f in csv.reader(open(f'{A}/boletales_mash/genomes.tsv'),delimiter='\t'): gmeta[os.path.basename(f).rsplit('.fna',1)[0]]=(g,org)
extra_names={'GCA_054741165.1':'Boletus edulis (chr, 2026)','GCA_015179015.1':'Boletus edulis BED1','GCA_014925795.1':'Boletus edulis','GCA_022837065.1':'Boletus edulis','GCA_038089135.1':'B. edulis (low-quality asm)','GCA_054948845.1':'Tylopilus felleus (chr, 2026)','GCA_038093275.1':'"Tylopilus felleus" 2024','GCA_038091215.1':'Retiboletus ornatipes 2024','GCA_038093655.1':'Boletus reticuloceps 2024'}
ts=TreeStyle(); ts.show_leaf_name=False; ts.show_scale=False; ts.branch_vertical_margin=1; ts.scale=1200
ts.legend_position=4
for l in LIN: ts.legend.add_face(CircleFace(6,LCOL[l],'circle'),column=0); ts.legend.add_face(TextFace(' '+l+'  (filled = intact, hollow = pseudogene/fragment)',fsize=9),column=1)
for leaf in t.iter_leaves():
    key=leaf.name.rsplit('.fna',1)[0]; g,org=gmeta.get(key,(leaf.name,leaf.name))
    if org=='extra': org=extra_names.get(g,g)
    cn=cname.get(g,''); cn='' if cn=='extra' else cn; lab=f'{org}'+(f'  [{cn}]' if cn and cn.split()[0]!=org.split()[0] else '')+f'  {g}'
    ns=NodeStyle(); ns['size']=0; leaf.set_style(ns)
    leaf.add_face(TextFace(lab,fsize=8,fgcolor='black' if not (cn and cn.split()[0]!=org.split()[0]) else '#b03a2e'),column=0,position='aligned')
    for i,l in enumerate(LIN):
        st=pres.get(g,{}).get(l)
        if st:
            fill = 'intact' in st
            leaf.add_face(CircleFace(5,LCOL[l] if fill else 'white','circle') if fill else RectFace(9,9,LCOL[l],'white'),column=1+i,position='aligned')
        else: leaf.add_face(TextFace('  ',fsize=8),column=1+i,position='aligned')
for n in t.traverse():
    if not n.is_leaf(): ns=NodeStyle(); ns['size']=0; n.set_style(ns)
t.render(f'{F}/F1_boletales_mashtree_presence.png',tree_style=ts,w=1600,units='px',dpi=150); t.render(f'{F}/F1_boletales_mashtree_presence.svg',tree_style=ts,w=1600,units='px'); t.render(f'{F}/F1_boletales_mashtree_presence.pdf',tree_style=ts,w=1600,units='px')
print('F1 done')

# ---------- F2: gene trees ----------
def group_of(n):
    if n.startswith('REF_'):
        sg=n.split('__')[0].split('_',2)[2]
        return 'REF '+sg
    if n.startswith('FUNGI_'):
        p=n.split('__')[0].split('_'); return 'FUNGI '+('cus_actinobacteria' if 'actinobacteria' in n else p[1])
    return n.split('__')[0]
GCOL={'FUNGI sus':'#d62728','FUNGI cus':'#1f77b4','FUNGI mis':'#00838f','FUNGI cus_actinobacteria':'#5a5a5a','FUNGI cusrel181':'#8a7d00',
      'FUNGI gdpd':'#1a9641','FUNGI acs':'#7b3fa0','REF sus':'#b30000','REF cus':'#08519c','REF mis':'#006d6d','REF acs':'#6a1b9a','REF gdpd_meso':'#1b7837',
      'REF kerr14':'#c51b8a','REF rcorf12':'#8a7d00','REF cusrel181':'#8a7d00','REF cus_actinobacteria':'#404040','REF cus_betaproteobacteria':'#404040','REF cus_other_bacteria':'#404040',
      'BACT_acc':'#8c564b','OG':'#000000'}
SIX=['JBBPKL010003317','JBBPKL010018021','JBBPNB010000203','JBBPNB010000771','JBBPOG010012756','JBBPRY010001778']
def render_tree(fam, tf, out, highlights):
    t=Tree(tf,format=1)
    og=[l for l in t.iter_leaves() if l.name.startswith('OG__')]
    if og: t.set_outgroup(og[0])
    ts=TreeStyle(); ts.show_leaf_name=False; ts.show_scale=True; ts.branch_vertical_margin=1; ts.scale=400
    for leaf in t.iter_leaves():
        g=group_of(leaf.name); col=GCOL.get(g,'#444444')
        if g.startswith('FUNGI gdpd'): col=GCOL['FUNGI gdpd']
        if g.startswith('NATIVE'): col='#33691e'
        nm=leaf.name.replace('FUNGI_','').replace('REF_','REF ')[:70]
        if any(s in leaf.name for s in SIX): nm='★ '+nm
        leaf.add_face(TextFace(nm,fsize=9,fgcolor=col),column=0)
        ns=NodeStyle(); ns['size']=0; leaf.set_style(ns)
    for n in t.traverse():
        if n.is_leaf(): continue
        ns=NodeStyle(); ns['size']=0; n.set_style(ns)
        m=re.match(r'^([\d.]+)/(\d+)$',n.name or '')
        if m and int(m.group(2))>=90: n.add_face(TextFace(n.name,fsize=6,fgcolor='#333333'),column=0,position='branch-top')
    for label,sel,color in highlights:
        leaves=[l for l in t.iter_leaves() if sel(l.name)]
        if len(leaves)<2: continue
        mrca=t.get_common_ancestor(leaves); ns=NodeStyle(); ns['bgcolor']=color; ns['size']=0; mrca.set_style(ns)
        mrca.add_face(TextFace(' '+label,fsize=11,fgcolor='black',bold=True),column=0,position='branch-bottom')
    t.render(out,tree_style=ts,w=2400,units='px',dpi=150)
    t.render(out.replace('.png','.svg'),tree_style=ts,w=2400,units='px')
    t.render(out.replace('.png','.pdf'),tree_style=ts,w=2400,units='px')
render_tree('SUS',f'{A}/trees_final/SUS/SUS.mafft.iqtree.treefile',f'{F}/F2_SUS_tree.png',
    [('Клада T-DNA sus: 52 грибных + 24 референса',lambda n:n.startswith('REF_') and '_sus__' in n,'#ffe5e5')])
render_tree('CUS',f'{A}/trees_final/CUS/CUS.mafft.iqtree.treefile',f'{F}/F2_CUS_tree.png',
    [('Клада T-DNA cus: 42 грибных + 21 референс',lambda n:n.startswith('REF_') and '_cus__' in n,'#e5f0ff'),
     ('Актинобактериальная линия (не T-DNA)',lambda n:'cus_actinobacteria' in n and n.startswith('FUNGI'),'#eeeeee'),
     ('mis: грибы + референсы вперемешку',lambda n:n.startswith('FUNGI_mis_'),'#e0f7fa')])
render_tree('ACS',f'{A}/trees_final/ACS/ACS.mafft.iqtree.treefile',f'{F}/F2_ACS_tree.png',
    [('Линия acc-локуса (gdpd_meso): 60 Boletales',lambda n:n.startswith('FUNGI_gdpd_meso'),'#e5ffe5'),
     ('«Шестёрка» + Parasponia',lambda n:any(s in n for s in SIX) or 'PON40589' in n,'#f0e5ff'),
     ('T-DNA acs референсы',lambda n:n.startswith('REF_') and '_acs__' in n and 'PON40589' not in n,'#f7f0ff')])
print('F2 trees done')

# ---------- F2b: AU tests ----------
AU={'SUS (глобальный тест)':[('ML без ограничений',0,1.0),('грибные SUS + T-DNA sus',14.8,0.264),('грибные SUS + альтернативные OCD',176.1,7e-50)],
 'ACS (по линиям)':[('ML без ограничений',0.3,0.579),('шестёрка + T-DNA acs',1723.7,1.4e-60),('шестёрка + kerr14 + Parasponia',17.5,0.415),('шестёрка + acc-линия',62.7,0.132),('шестёрка + нативные GDPD',3267,2.4e-6),('gdpd_meso + acc-линия',0,1.0),('gdpd_meso + T-DNA acs',2121,5.9e-47)],
 'CUS (по линиям, без нативного гена)':[('ML без ограничений',2.6,0.528),('грибные cus + T-DNA cus',16.2,0.315),('грибные cus + T-DNA cus + 4 бакт. реф.',0,1.0),('грибные cus + актино',405.8,1.9e-13),('грибные cus + mis',920.7,3.4e-47),('грибные cus + прочие бактерии',1107,4.9e-4),('все грибные cus монофилетичны',129.4,6.4e-5),('Boletales cus + T-DNA cus (без др. грибов)',519.3,1.3e-3),('Boletales cus + актино',643.7,5e-59),('Boletales cus + mis',1162,3e-77),('Boletales cus + остальные грибные cus',134.8,1.6e-5)]}
fig,axes=plt.subplots(3,1,figsize=(11,13),gridspec_kw={'height_ratios':[3,7,11]})
for ax,(fam,rows) in zip(axes,AU.items()):
    names=[r[0] for r in rows]; dl=[max(r[1],0.1) for r in rows]; p=[r[2] for r in rows]
    cols=['#2ca02c' if x>=0.05 else '#d62728' for x in p]
    y=np.arange(len(rows))[::-1]
    ax.barh(y,dl,color=cols); ax.set_xscale('log'); ax.set_yticks(y); ax.set_yticklabels(names,fontsize=9); ax.set_xlabel('ΔlogL — насколько хуже свободного дерева (лог. шкала)')
    for yi,d,pp in zip(y,dl,p): ax.text(d*1.15,yi,f'p-AU={pp:.2g}',va='center',fontsize=7)
    ax.set_title(fam,fontsize=11,loc='left'); ax.set_xlim(0.08,3e5)
axes[0].legend(handles=[mp.Patch(color='#2ca02c',label='гипотеза не отвергнута (p-AU ≥ 0.05)'),mp.Patch(color='#d62728',label='гипотеза отвергнута (p-AU < 0.05)')],loc='lower right',fontsize=9)
plt.tight_layout(); plt.savefig(f'{F}/F2b_AU_tests.png',dpi=150); plt.savefig(f'{F}/F2b_AU_tests.svg'); plt.close(); print('F2b done')

# ---------- F3: B. edulis contig ----------
dep=[int(l.split('\t')[2]) for l in open(f'{D}/bedulis_contig_depth.tsv')]
x=np.arange(1,len(dep)+1); w=200; sm=np.convolve(dep,np.ones(w)/w,mode='same')
fig,ax=plt.subplots(2,1,figsize=(12,5.5),gridspec_kw={'height_ratios':[3,1.4]},sharex=True)
ax[0].fill_between(x,sm,color='#9ecae1'); ax[0].plot(x,sm,color='#3182bd',lw=1); ax[0].axhline(91.3,ls='--',color='grey'); ax[0].text(7400,93,'медиана генома 91×',ha='right',fontsize=8,color='grey')
ax[0].set_ylabel('глубина покрытия ридами'); ax[0].set_title('Boletus edulis, изолят JBBPRY (GCA_038093215.1), контиг JBBPRY010001778 (7 415 п.н.)',fontsize=10)
ax[0].axvline(1313,color='k',ls=':'); ax[0].text(1400,112,'стык повтор|вставка: 19 ридов и 9 пар\nпересекают границу — сборка не химерная',fontsize=8,va='top')
ax[0].set_ylim(0,125)
ax[1].set_ylim(0,3); ax[1].set_yticks([])
ax[1].add_patch(mp.Rectangle((51,2.1),989-51,0.7,color='#fdae6b')); ax[1].text(520,2.45,'g1: многокопийный повтор (TE)',ha='center',va='center',fontsize=8)
ax[1].add_patch(mp.Rectangle((2690,2.1),6606-2690,0.7,color='#9467bd')); ax[1].text(4650,2.45,'ACS (интактный ген, 5 экзонов, 418 аа)',ha='center',va='center',fontsize=8,color='white')
ax[1].add_patch(mp.Rectangle((7202,2.1),213,0.7,color='#fdae6b'))
ax[1].add_patch(mp.Rectangle((2,1.1),1311,0.7,color='#a1d99b')); ax[1].text(650,1.45,'есть в хромосомной сборке\n(CM139040, 94%)',ha='center',va='center',fontsize=7)
ax[1].add_patch(mp.Rectangle((1313,1.1),5889,0.7,color='#fc9272')); ax[1].text(4250,1.45,'отсутствует во всех 5 других сборках B. edulis (сегмент 5.9 кб)',ha='center',va='center',fontsize=8)
ax[1].add_patch(mp.Rectangle((7202,1.1),213,0.7,color='#a1d99b'))
ax[1].add_patch(mp.Rectangle((0,0.1),7415,0.7,color='#eeeeee')); ax[1].text(3700,0.45,'доля правильно спаренных ридов 0.97 (геном 0.98); MAPQ 57; клиппированные риды только на концах контига',ha='center',va='center',fontsize=7)
ax[1].set_xlabel('позиция на контиге, п.н.'); ax[1].set_xlim(0,7415)
plt.tight_layout(); plt.savefig(f'{F}/F3_Bedulis_contig.png',dpi=150); plt.savefig(f'{F}/F3_Bedulis_contig.svg'); plt.close(); print('F3 done')

# ---------- F4: Butyriboletus island ----------
seq=''.join(l.strip() for l in open(f'{D}/island_JAHRGV010000066_region.fa') if not l.startswith('>')); off=1135000
gw=1000; gc=[(sum(1 for c in seq[i:i+gw] if c in 'GCgc')/max(1,len(seq[i:i+gw]))*100, off+i+gw//2) for i in range(0,len(seq)-gw,500)]
genes=list(csv.DictReader(open(f'{D}/island_JAHRGV010000066_genes.tsv'),delimiter='\t'))
reg={'ACS':(1156808,1157542,'ACS (acs, интактный)'),'CUS':(1163752,1164640,'CUS (cus, псевдоген)')}
fig,ax=plt.subplots(2,1,figsize=(13,6),gridspec_kw={'height_ratios':[1.2,1.8]},sharex=True)
ax[0].plot([g[1] for g in gc],[g[0] for g in gc],color='#3182bd'); ax[0].axhline(44.5,ls='--',color='grey'); ax[0].text(off+500,45,'GC сборки 44.5%',fontsize=8,color='grey')
ax[0].axvspan(1146510,1189199,color='#fff3cd',alpha=0.6); ax[0].set_ylabel('GC, % (окно 1 кб)'); ax[0].set_title('Butyriboletus roseoflavus (GCA_019364595.1), контиг JAHRGV010000066 (2.66 Мб): островок 1 146 510–1 189 199 (42.7 кб)',fontsize=10)
ax[1].axvspan(1146510,1189199,color='#fff3cd',alpha=0.6); ax[1].set_ylim(-0.15,3.5); ax[1].set_yticks([])
ROLE={'bacterial_like':'#e6550d','Fungi':'#31a354','no_hit':'#bdbdbd','anchor':'#31a354','':'#bdbdbd'}
for g in genes:
    s,e=int(g['start']),int(g['end']); lab=g['gene']; role=g['role']; col=ROLE.get(role,'#bdbdbd')
    title=g['top_title'].lower()
    if g['gene'] in ('g20','g26','g27'): col='#fdae6b'; lab+='\nповтор ×3'
    elif 'transposase' in title: col='#e6550d'; lab+='\nIS481-транспозаза'
    elif role=='bacterial_like': lab+='\nбакт.-подобный'
    elif role=='Fungi': lab+='\nгрибной'
    y=2.2 if g['strand']=='+' else 1.5
    k=int(g['gene'][1:]); dy=0.0 if k%2==0 else 0.28
    ax[1].add_patch(mp.Rectangle((s,y),e-s,0.5,color=col)); ax[1].text((s+e)/2,(y+0.62+dy) if y>2 else (y-0.32-dy),lab,ha='center',va='center',fontsize=7)
for k,(s,e,lab) in reg.items():
    ax[1].add_patch(mp.Rectangle((s,0.55),e-s,0.5,color='#9467bd' if k=='ACS' else '#1f77b4')); ax[1].text((s+e)/2,0.25 if k=='ACS' else 0.02,lab,ha='center',va='center',fontsize=7.5)
ax[1].text(1146510,3.15,'граница островка →',fontsize=7); ax[1].text(1189199,3.15,'← граница островка',fontsize=7,ha='right')
ax[1].set_xlabel('координата на контиге, п.н.'); ax[1].set_xlim(off,off+75000)
ax[1].legend(handles=[mp.Patch(color='#9467bd',label='ACS'),mp.Patch(color='#1f77b4',label='CUS'),mp.Patch(color='#e6550d',label='бактериально-подобные гены / транспозаза'),mp.Patch(color='#31a354',label='грибные гены (якоря)'),mp.Patch(color='#fdae6b',label='копии повтора на границах'),mp.Patch(color='#bdbdbd',label='без гомологов')],fontsize=8,loc='upper center',bbox_to_anchor=(0.5,-0.25),ncol=3,frameon=False)
plt.tight_layout(); plt.savefig(f'{F}/F4_Butyriboletus_island.png',dpi=150,bbox_inches='tight'); plt.savefig(f'{F}/F4_Butyriboletus_island.svg',bbox_inches='tight'); plt.close(); print('F4 done')

# ---------- F5: read coverage ----------
cov=list(csv.DictReader(open(f'{D}/locus_coverage.tsv'),delimiter='\t'))
sp=sorted({r['species'] for r in cov}); fig,ax=plt.subplots(figsize=(11,5))
for i,s in enumerate(sp):
    rows=[r for r in cov if r['species']==s and r['status'] in ('intact','pseudogene_like')]
    for r in rows:
        lin=r['fam']+':'+r['subgroup']; col={'SUS:sus':'#d62728','CUS:cus':'#1f77b4','CUS:mis':'#17becf','ACS:gdpd_meso':'#2ca02c','ACS:acs':'#9467bd','ORF358:orf358':'#8c564b','ACS:kerr14':'#e377c2'}.get(lin,'#7f7f7f')
        ax.scatter(i+np.random.uniform(-0.25,0.25),float(r['locus_ratio']),color=col,s=25+min(int(r['contig_len']),60000)/600,alpha=0.8,edgecolor='k',lw=0.3)
ax.axhline(1,color='grey',ls='--'); ax.axhline(0.5,color='grey',ls=':'); ax.set_yscale('log'); ax.set_ylim(0.05,4)
ax.set_xticks(range(len(sp))); ax.set_xticklabels([s.replace('_',' ') for s in sp],rotation=25,ha='right',fontsize=8); ax.set_ylabel('глубина локуса / медиана генома')
ax.set_title('Покрытие локусов ридами SRA (размер точки ~ длина контига; 1 = как весь геном, 0.5 = половина)',fontsize=10)
ax.legend(handles=[mp.Patch(color=c,label=l) for l,c in [('SUS','#d62728'),('CUS cus','#1f77b4'),('CUS mis','#17becf'),('ACS acc-локус','#2ca02c'),('ACS acs/шестёрка','#9467bd'),('ORF358','#8c564b')]],fontsize=8,loc='lower left',ncol=3)
plt.tight_layout(); plt.savefig(f'{F}/F5_read_coverage.png',dpi=150); plt.savefig(f'{F}/F5_read_coverage.svg'); plt.close(); print('F5 done')

# ---------- F6: mash heatmap of 25 assemblies ----------
import seaborn as sns
gm={}; 
for g,s,f in csv.reader(open(f'{D}/indep_genomes.tsv'),delimiter='\t'): gm[os.path.basename(f)]=(g,s)
names=sorted({gm[os.path.basename(a)] for a,b,d,p,sh in csv.reader(open(f'{D}/indep_mash_dist.tsv'),delimiter='\t')})
idx={n:i for i,n in enumerate(names)}; M=np.zeros((len(names),len(names)))
for a,b,d,p,sh in csv.reader(open(f'{D}/indep_mash_dist.tsv'),delimiter='\t'): M[idx[gm[os.path.basename(a)]],idx[gm[os.path.basename(b)]]]=100*(1-float(d))
lab=[f"{s.replace('_',' ')} {g}" for g,s in names]
df=pd.DataFrame(M,index=lab,columns=lab)
cg=sns.clustermap(df,cmap='viridis',vmin=75,vmax=100,figsize=(12,12),annot=False,cbar_kws={'label':'ANI, % (по mash)'})
cg.ax_heatmap.set_xticklabels(cg.ax_heatmap.get_xticklabels(),fontsize=7); cg.ax_heatmap.set_yticklabels(cg.ax_heatmap.get_yticklabels(),fontsize=7)
cg.savefig(f'{F}/F6_mash_ANI_heatmap.png',dpi=130); plt.close('all'); print('F6 done')
