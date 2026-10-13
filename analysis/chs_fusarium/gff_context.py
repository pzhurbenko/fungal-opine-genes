import re,sys,collections
G='gff_dl/ncbi_dataset/data/GCF_000151355.1/genomic.gff'
genes={};mr={};ex=collections.Counter();cds={}
def at(s): return dict(x.split('=',1) for x in s.split(';') if '=' in x)
for l in open(G):
    if l[0]=='#': continue
    p=l.rstrip('\n').split('\t')
    if len(p)<9: continue
    a=at(p[8]); t=p[2]
    if t=='gene': genes[a['ID']]=dict(chr=p[0],s=int(p[3]),e=int(p[4]),st=p[6],name=a.get('Name',''),biotype=a.get('gene_biotype',''),prods=[],prot=[],nex=0)
    elif t in('mRNA','transcript','rna') or t.endswith('RNA'):
        mr[a['ID']]=(a['Parent'],a.get('product',''))
    elif t=='exon': ex[a['Parent']]+=1
    elif t=='CDS' and 'protein_id' in a:
        cds[a['Parent']]=a['protein_id']
for r,(g,prod) in mr.items():
    if g in genes:
        d=genes[g]; d['prods'].append(prod); d['nex']=max(d['nex'],ex[r]); 
        if r in cds: d['prot'].append(cds[r])
def find(pid):
    for k,d in genes.items():
        if pid in d['prot']: return k
bykey={}
for k,d in genes.items(): bykey.setdefault(d['chr'],[]).append(k)
for c in bykey: bykey[c].sort(key=lambda k:genes[k]['s'])
kw=re.compile(r'transposase|integrase|retrotranspos|reverse transcriptase|gag|polyprotein|helitron|MULE|hAT|Tc1|mariner|Gypsy|Copia|transposon|bacteri',re.I)
def show(pid,out):
    k=find(pid); d=genes[k]; lst=bykey[d['chr']]; i=lst.index(k)
    out.write(f'## {pid} {k} {d["chr"]}:{d["s"]}-{d["e"]}({d["st"]}) exons={d["nex"]} prod={d["prods"]} scaffold_genes={len(lst)}\n')
    out.write('rel\tgene\tprotein\tstart\tend\tstrand\texons\tdist_bp\tproduct\n')
    for j in range(max(0,i-10),min(len(lst),i+11)):
        g=genes[lst[j]]; dist=0 if j==i else (g['s']-d['e'] if g['s']>d['e'] else d['s']-g['e'])
        out.write(f'{j-i}\t{g["name"]}\t{",".join(g["prot"])}\t{g["s"]}\t{g["e"]}\t{g["st"]}\t{g["nex"]}\t{dist}\t{"|".join(g["prods"])}\n')
    out.write('# +-50kb window flagged (transposase/bacterial keywords):\n')
    n=0
    for j in lst:
        g=genes[j]
        if j==k: continue
        dist=g['s']-d['e'] if g['s']>d['e'] else d['s']-g['e']
        if dist<=50000 and (g['e']>=d['s']-50000):
            n+=1
            if kw.search(' '.join(g['prods'])+g['name']): out.write(f'FLAG\t{g["name"]}\t{g["s"]}\t{"|".join(g["prods"])}\n')
    out.write(f'# genes in +-50kb window: {n}\n')
    return d
out=open('vanettenii_context.txt','w')
d=show('XP_003047010.1',out)
# acs / GDPD-like and other 'opine_related' proteins
k2=find('XP_003054402.1'); d2=genes[k2]
out.write(f'\n## acs-like (GDPD) XP_003054402.1: {d2["chr"]}:{d2["s"]}-{d2["e"]} same_scaffold={d2["chr"]==d["chr"]} exons={d2["nex"]} prod={d2["prods"]}\n')
lst=bykey[d['chr']]
out.write('\n# all genes within +-50kb:\n')
for j in lst:
    g=genes[j]
    if g['e']>=d['s']-50000 and g['s']<=d['e']+50000:
        out.write(f'{g["name"]}\t{",".join(g["prot"])}\t{g["s"]}-{g["e"]}\t{g["st"]}\tex={g["nex"]}\t{"|".join(g["prods"])}\n')
out.write('# scaffold NW_003315868.1 length 907961 (from ##sequence-region)\n')
out.close()
