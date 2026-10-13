#!/usr/bin/env python3
"""(a) contig-level: for each source contig x target assembly: best target, aligned fraction of contig, identity, locus covered?, flank covered?
   (b) lineage presence: for each target assembly x (fam,subgroup): best miniprot identity/coverage among Boletales models."""
import csv, glob, os, re
from collections import defaultdict
src=list(csv.DictReader(open('source_loci.tsv'),delimiter='\t'))
bycontig=defaultdict(list)
for r in src: bycontig[(r['species'],r['gca'],r['contig'])].append(r)
out=open('contig_reproducibility.tsv','w')
out.write('species\tsrc_gca\tcontig\tcontig_len\tloci\ttgt_gca\tn_tgt_seqs\tbest_tgt\taligned_frac\tidentity\tlocus_covered_frac\tleft_flank_covered\tright_flank_covered\n')
for paf in sorted(glob.glob('paf/*.paf')):
    sp,s,_,t=os.path.basename(paf)[:-4].split('__')
    hits=defaultdict(list)
    for l in open(paf):
        f=l.split('\t'); hits[f[0]].append((int(f[2]),int(f[3]),f[5],int(f[9]),int(f[10])))
    for (sp2,s2,ctg),L in bycontig.items():
        if (sp2,s2)!=(sp,s): continue
        clen=int(L[0]['contig_len']); H=hits.get(ctg,[])
        cov=[0]*(clen+1)
        for qs,qe,tn,m,bl in H:
            for i in range(qs,qe): cov[i]=1
        alfrac=sum(cov)/clen
        tgts=defaultdict(int)
        for qs,qe,tn,m,bl in H: tgts[tn]+=qe-qs
        best=max(tgts,key=tgts.get) if tgts else ''
        ident=sum(m for *_,m,bl in H)/max(1,sum(bl for *_,m,bl in H)) if H else 0
        lo=min(int(r['contig_start']) for r in L); hi=max(int(r['contig_end']) for r in L)
        lc=sum(cov[lo:hi])/max(1,hi-lo)
        lf=sum(cov[max(0,lo-2000):lo])/max(1,min(2000,lo)) if lo>200 else float('nan')
        rf=sum(cov[hi:min(clen,hi+2000)])/max(1,min(2000,clen-hi)) if clen-hi>200 else float('nan')
        out.write('\t'.join(map(str,[sp,s,ctg,clen,';'.join(r['locus_id'] for r in L),t,len(tgts),best,f'{alfrac:.2f}',f'{ident:.3f}',f'{lc:.2f}',f'{lf:.2f}',f'{rf:.2f}']))+'\n')
out.close()
# (b)
qlen={}
for l in open('boletales_models.faa'):
    if l.startswith('>'): cur=l[1:].strip(); qlen[cur]=0
    else: qlen[cur]+=len(l.strip())
gen={r[0]:r[1] for r in csv.reader(open('genomes.tsv'),delimiter='\t')}
pres=open('lineage_presence.tsv','w'); pres.write('tgt_gca\tspecies\tfam\tsubgroup\tn_hits_id50\tbest_identity\tbest_cov\tbest_query\tbest_target\n')
for gff in sorted(glob.glob('gff/*.gff')):
    g=os.path.basename(gff)[:-4]; best={}; n=defaultdict(int)
    for l in open(gff):
        if l.startswith('#') or '\tmRNA\t' not in l: continue
        f=l.rstrip().split('\t'); a=dict(x.split('=',1) for x in f[8].split(';') if '=' in x)
        q,qs,qe=a['Target'].split()[:3]; ident=float(a.get('Identity',0)); cov=(int(qe)-int(qs))/qlen[q]
        fam=q.split('_')[0]; sub=q.split('|')[1]; key=(fam,sub)
        if ident>=0.5 and cov>=0.5: n[key]+=1
        if key not in best or ident*cov>best[key][0]*best[key][1]: best[key]=(ident,cov,q,f'{f[0]}:{f[3]}-{f[4]}')
    for key in sorted(best):
        ident,cov,q,loc=best[key]; pres.write('\t'.join(map(str,[g,gen.get(g,''),key[0],key[1],n[key],f'{ident:.2f}',f'{cov:.2f}',q,loc]))+'\n')
pres.close(); print('summaries written')
