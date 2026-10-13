#!/usr/bin/env python3
import csv, glob, os
from collections import defaultdict
qlen={}
for f in glob.glob('q/*.faa'):
    for l in open(f):
        if l.startswith('>'): cur=l[1:].strip(); qlen[cur]=0
        else: qlen[cur]+=len(l.strip())
meta={r['qid']:r for r in csv.DictReader(open('queries_meta.tsv'),delimiter='\t')}
best=defaultdict(dict)  # tgt -> qid -> (ident,cov,contig,s,e)
for gff in glob.glob('gff/*.gff'):
    g=os.path.basename(gff)[:-4]
    for l in open(gff):
        if l.startswith('#') or '\tmRNA\t' not in l: continue
        f=l.rstrip().split('\t'); a=dict(x.split('=',1) for x in f[8].split(';') if '=' in x)
        q,qs,qe=a['Target'].split()[:3]; ident=float(a.get('Identity',0)); cov=(int(qe)-int(qs))/qlen[q]
        if ident<0.4 or cov<0.5: continue
        if q not in best[g] or ident*cov>best[g][q][0]*best[g][q][1]: best[g][q]=(ident,cov,f[0],int(f[3]),int(f[4]))
out=open('island_synteny.tsv','w')
out.write('island\torganism\tisland_len\tlineages\ttgt_gca\ttgt_org\tmash_d\trole\tanchors_found\tanchorL\tanchorR\tsame_contig\ttgt_gap_bp\tsrc_gap_bp\tisland_genes\tisl_hits_between\tisl_hits_elsewhere\tverdict\n')
for r in csv.DictReader(open('plan.tsv'),delimiter='\t'):
    name,g=r['island'],r['tgt_gca']; B=best.get(g,{})
    qs=[q for q in meta if q.startswith(name+'|')]
    L=[q for q in qs if meta[q]['kind'].startswith('L')]; R=[q for q in qs if meta[q]['kind'].startswith('R')]; I=[q for q in qs if meta[q]['kind']=='I']
    hl=[(q,B[q]) for q in L if q in B]; hr=[(q,B[q]) for q in R if q in B]
    def ov(a,b): return a[2]==b[2] and not (a[4]<b[3] or b[4]<a[3])
    repeat_flank=any(ov(x[1],y[1]) for x in hl for y in hr)
    # choose the first L/R pair whose target hits do not overlap (distinct genes)
    pair=None
    for x in hl:
        for y in hr:
            if not ov(x[1],y[1]): pair=(x,y); break
        if pair: break
    aL=f"{pair[0][0].split('|')[1]}:{pair[0][1][2]}:{pair[0][1][3]}-{pair[0][1][4]}({pair[0][1][0]:.2f})" if pair else ('.' if not hl else 'ONLY_REPEAT')
    aR=f"{pair[1][0].split('|')[1]}:{pair[1][1][2]}:{pair[1][1][3]}-{pair[1][1][4]}({pair[1][1][0]:.2f})" if pair else ('.' if not hr else 'ONLY_REPEAT')
    same=bool(pair and pair[0][1][2]==pair[1][1][2]); gap='.'; between=0; elsewhere=0; verdict='anchors_missing'
    if pair:
        srcgap=abs(int(meta[pair[1][0]]['abs_start'])-int(meta[pair[0][0]]['abs_end']))
        if same:
            lo=min(pair[0][1][3],pair[0][1][4],pair[1][1][3],pair[1][1][4]); hi=max(pair[0][1][3],pair[0][1][4],pair[1][1][3],pair[1][1][4]); gap=hi-lo
            for q in I:
                if q in B:
                    if B[q][2]==pair[0][1][2] and lo<=B[q][3]<=hi: between+=1
                    else: elsewhere+=1
            verdict='island_absent_anchors_adjacent' if (gap<0.5*srcgap and between==0) else ('island_present_syntenic' if between>=1 else 'anchors_adjacent_partial')
        else:
            for q in I:
                if q in B: elsewhere+=1
            verdict='anchors_on_different_contigs'
    else:
        srcgap='.'
        if hl and hr: verdict='only_repeat_anchors'
    verdict=verdict+('+repeat_flank' if repeat_flank else '')
    out.write('\t'.join(map(str,[name,r['organism'],r['island_len'],r['lineages'],g,r['tgt_org'],r['mash_d'],r['role'],f'{len(hl)}/{len(L)}+{len(hr)}/{len(R)}',aL,aR,same,gap,srcgap,len(I),between,elsewhere,verdict]))+'\n')
out.close(); print('island_synteny.tsv written')
