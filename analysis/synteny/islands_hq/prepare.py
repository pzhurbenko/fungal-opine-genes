#!/usr/bin/env python3
"""Island synteny vs sister genomes. For each bounded Boletales island: anchors (2 multi-exon fungal genes each side) + island genes
-> miniprot against (a) up to 3 closest genomes LACKING the island lineages and (b) 1 closest genome HAVING them (control)."""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import csv, os, glob, json
from collections import defaultdict
D=f'{GMO_ROOT}'; A=f'{D}/analysis'
lt={r['locus_id']:r for r in csv.DictReader(open(f'{A}/loci_pipeline/round2/loci_table.tsv'),delimiter='\t')}
gen={r[0]:(r[1],r[3]) for r in csv.reader(open(f'{A}/boletales_mash/genomes.tsv'),delimiter='\t')}   # gca -> (org, fna)
pres={r['gca']:r for r in csv.DictReader(open(f'{A}/distribution/presence_step2.tsv'),delimiter='\t')}
# mash distances
base2g={os.path.basename(f):g for g,(o,f) in gen.items()}
dist=defaultdict(dict)
for a,b,d,p,sh in csv.reader(open(f'{A}/boletales_mash/dist.tsv'),delimiter='\t'):
    ga,gb=base2g.get(os.path.basename(a)),base2g.get(os.path.basename(b))
    if ga and gb and ga!=gb: dist[ga][gb]=float(d)
# islands: bounded both sides, Boletales genome, dedup by contig
isl={}
for r in csv.DictReader(open(f'{A}/islands/islands.tsv'),delimiter='\t'):
    L=lt.get(r['locus_id']); 
    if not L or L['gca'] not in gen or r['bounded_left']!='yes' or r['bounded_right']!='yes': continue
    if int(r['island_len'])<3000: continue
    if L['gca'] not in ('GCA_019364595.1','GCA_018397855.1','GCA_003315995.1','GCA_038093215.1','GCA_038091295.1'): continue
    key=(L['gca'],r['contig'])
    e=isl.setdefault(key,{'gca':L['gca'],'org':L['organism'],'contig':r['contig'],'start':int(r['island_start']),'end':int(r['island_end']),'loci':[],'lineages':set(),'gc_island':r['gc_island'],'gc_assembly':r['gc_assembly']})
    e['loci'].append(r['locus_id']); e['lineages'].add(f"{L['fam']}:{L['subgroup']}")
    e['start']=min(e['start'],int(r['island_start'])); e['end']=max(e['end'],int(r['island_end']))
print('islands',len(isl))
def has_lineage(g,lins):
    p=pres.get(g); 
    if not p: return None
    return any(p.get(l,'') not in ('','not_homolog','no_model') for l in lins)
plan=open('plan.tsv','w'); plan.write('island\tsrc_gca\torganism\tcontig\tisland_start\tisland_end\tisland_len\tlineages\ttgt_gca\ttgt_org\tmash_d\ttgt_has_lineage\trole\n')
queries=defaultdict(list)  # tgt_gca -> list of (qid, seq)
meta=open('queries_meta.tsv','w'); meta.write('island\tqid\tkind\tabs_start\tabs_end\tstrand\texons\n')
for key,e in isl.items():
    name=f"{e['gca']}__{e['contig']}"
    # anchors & island genes from the first locus dir with neighbors
    genes=[]
    for lid in e['loci']:
        f=f"{A}/loci_pipeline/round2/loci/{lid}/neighbors.faa"
        if os.path.exists(f):
            cur=None
            for l in open(f):
                if l.startswith('>'):
                    h=l[1:].strip(); parts=h.split('|'); coords=parts[2]; s,t=coords[:-1].split('-'); strand=coords[-1]; ex=int(parts[3].split('=')[1])
                    cur={'id':h,'g':parts[1],'s':int(s),'e':int(t),'strand':strand,'ex':ex,'seq':''}; genes.append(cur)
                else: cur['seq']+=l.strip()
            break
    if not genes: print('no neighbors for',name); continue
    seen=set(); genes=[g for g in genes if not (g['g'] in seen or seen.add(g['g']))]
    left=[g for g in genes if g['e']<e['start'] and g['ex']>=2]; right=[g for g in genes if g['s']>e['end'] and g['ex']>=2]
    # drop repeat-like anchors: proteins >=80% identical (k-mer Jaccard proxy) to another neighbour on the same contig
    def kmers(x,k=5): return {x[i:i+k] for i in range(len(x)-k+1)}
    K={g['g']:kmers(g['seq']) for g in genes}
    def repeatlike(g):
        for h in genes:
            if h['g']==g['g']: continue
            j=len(K[g['g']]&K[h['g']])/max(1,len(K[g['g']]|K[h['g']]))
            if j>=0.5: return True
        return False
    rep={g['g'] for g in genes if repeatlike(g)}
    left=[g for g in left if g['g'] not in rep]; right=[g for g in right if g['g'] not in rep]
    left=sorted(left,key=lambda g:-g['e'])[:4]; right=sorted(right,key=lambda g:g['s'])[:4]
    if rep: print(name,'repeat-like neighbours excluded:',sorted(rep))
    inside=[g for g in genes if g['s']>=e['start'] and g['e']<=e['end']]
    for lid in e['loci']:
        f=f"{A}/loci_pipeline/round2/loci/{lid}/target_protein.faa"
        if os.path.exists(f):
            seq=''.join(l.strip() for l in open(f) if not l.startswith('>'))
            L=lt[lid]; inside.append({'id':lid,'g':'T_'+lid,'s':int(L['contig_start']),'e':int(L['contig_end']),'strand':L['strand'],'ex':int(L['n_exons'] or 1),'seq':seq})
    if len(left)<1 or len(right)<1: print('no anchors',name); continue
    # targets
    lins=sorted(e['lineages'])
    HQ=['GCA_054741165.1','GCA_015179015.1','GCA_054948845.1','GCA_018397855.1','GCA_019364595.1']
    tg=[]
    for g in HQ:
        if g==e['gca']: continue
        d=dist[e['gca']].get(g,-1)
        if d<0.02: continue
        hl=has_lineage(g,lins); tg.append((d,g,'lacking' if hl is False else ('having' if hl else 'unknown')))
    if not tg: print('no targets',name); continue
    qs=[(f"{name}|L{i+1}|{g['g']}",g['seq']) for i,g in enumerate(left)]+[(f"{name}|R{i+1}|{g['g']}",g['seq']) for i,g in enumerate(right)]+[(f"{name}|I|{g['g']}",g['seq']) for g in inside]
    for kind,g in [('L%d'%(i+1),g) for i,g in enumerate(left)]+[('R%d'%(i+1),g) for i,g in enumerate(right)]+[('I',g) for g in inside]:
        meta.write('\t'.join(map(str,[name,f"{name}|{kind}|{g['g']}",kind,g['s'],g['e'],g['strand'],g['ex']]))+'\n')
    for d,g,role in tg:
        plan.write('\t'.join(map(str,[name,e['gca'],e['org'],e['contig'],e['start'],e['end'],e['end']-e['start'],';'.join(lins),g,gen[g][0],f'{d:.4f}',has_lineage(g,lins),role]))+'\n')
        queries[g]+=qs
plan.close(); meta.close()
os.makedirs('q',exist_ok=True); tl=open('targets.tsv','w')
for g,qs in queries.items():
    with open(f'q/{g}.faa','w') as fh:
        for qid,seq in qs: fh.write(f'>{qid}\n{seq}\n')
    tl.write(f'{g}\t{gen[g][1]}\n')
tl.close(); print('target genomes',len(queries))
