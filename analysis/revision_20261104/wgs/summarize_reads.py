#!/usr/bin/env python3
"""Coverage and physical continuity; coordinates in output are 1-based inclusive.

Short-read gene edges and long-read fungal-anchor spans are reported separately.
No ancestral insertion boundary is inferred from a gene edge.
"""
from pathlib import Path
import csv,json,subprocess,sys,os
import numpy as np
import pysam
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
ENV=Path(os.environ.get('GMO_CONDA_ENVS', str(Path.home()/'miniforge3/envs')))
run=sys.argv[1];wd=HERE/run;meta=json.loads((wd/'provenance.json').read_text())
bam=wd/'reads.bam';samtools=ENV/'annot/bin/samtools';bedtools=ENV/'annot/bin/bedtools'
f=open(HERE/'inputs/S1a_with_legend.tsv');next(f);next(f)
genes=[r for r in csv.DictReader(f,delimiter='\t') if r['assembly']==meta['assembly']]
regions=[dict(region_id=r['locus_id'],kind='gene',contig=r['contig'],start=int(r['start']),end=int(r['end']),
              left_anchor_end=None,right_anchor_start=None) for r in genes]
if meta['strain']=='FP133':
    regions += [dict(region_id='Fp133_chr2_odh_block',kind='three_gene_block',contig='CM100255.1',start=9463190,end=9466883,
                     left_anchor_end=9460447,right_anchor_start=9470427,left_anchor='FPOAC2_07381',right_anchor='FPOAC2_07385'),
                dict(region_id='Fp133_chr3_odh_block',kind='three_gene_block',contig='CM100256.1',start=8420071,end=8423766,
                     left_anchor_end=8413923,right_anchor_start=8426768,left_anchor='FPOAC2_10360',right_anchor='FPOAC2_10366')]
af=pysam.AlignmentFile(bam,'rb');names={r.split('.')[0]:r for r in af.references}
for r in regions:r['contig']=names[r['contig'].split('.')[0]]
def primary(read,q):return not(read.flag & 0xF04) and read.mapping_quality>=q
def overlap(blocks,a,b):return sum(max(0,min(y,b)-max(x,a)) for x,y in blocks)
def median_hist(h):
    total=sum(h.values());cs=0;answers=[]
    for depth,n in sorted(h.items()):
        prev=cs;cs+=n
        for target in [(total-1)//2,total//2]:
            if prev<=target<cs:answers.append(depth)
    assert len(answers)==2
    return sum(answers)/2
output=[];spans=[];genome=[]
for q in [0,20]:
    histfile=wd/f'coverage_hist_MAPQ{q}.tsv'
    with open(histfile,'w') as out,open(wd/f'coverage_MAPQ{q}.log','w') as log:
        filt=subprocess.Popen([str(samtools),'view','-u','-F','3844','-q',str(q),str(bam)],stdout=subprocess.PIPE,stderr=log)
        cov=subprocess.run([str(bedtools),'genomecov','-ibam','stdin','-split'],stdin=filt.stdout,stdout=out,stderr=log)
        filt.stdout.close();assert filt.wait()==0 and cov.returncode==0
    hist={}
    for line in open(histfile):
        fields=line.split()
        if fields[0]=='genome':hist[int(fields[1])]=int(fields[2])
    median=median_hist(hist);total=sum(hist.values())
    genome.append(dict(run=run,min_MAPQ=q,assembly=meta['assembly'],genome_bp=total,median_depth=median,
                       mean_depth=sum(d*n for d,n in hist.items())/total,covered_fraction=1-hist.get(0,0)/total))
    for r in regions:
        a,b=r['start']-1,r['end'];lo=max(0,a-2000);hi=min(af.get_reference_length(r['contig']),b+2000)
        delta=np.zeros(hi-lo+1,dtype=np.int32);left=right=total_reads=proper=0;qualities=[];whole=0
        for read in af.fetch(r['contig'],lo,hi):
            if not primary(read,q):continue
            blocks=read.get_blocks()
            for x,y in blocks:
                x=max(x,lo)-lo;y=min(y,hi)-lo
                if y>x:delta[x]+=1;delta[y]-=1
            if overlap(blocks,a,b):
                total_reads+=1;proper+=read.is_proper_pair;qualities.append(read.mapping_quality)
            # For long reads tolerate small indels while requiring aligned bases
            # on both sides. Record the minimum so the two technologies are explicit.
            minimum=32 if meta['mapping_preset']=='map-ont' else 40
            if overlap(blocks,a-40,a)>=minimum and overlap(blocks,a,a+40)>=minimum:left+=1
            if overlap(blocks,b-40,b)>=minimum and overlap(blocks,b,b+40)>=minimum:right+=1
            if r['left_anchor_end'] is not None:
                x,y=r['left_anchor_end']-100,r['right_anchor_start']-1+100
                if (overlap(blocks,x,x+100)>=80 and overlap(blocks,y-100,y)>=80
                    and overlap(blocks,x,y)>=.8*(y-x)
                    and not any(op in (2,3) and length>500 for op,length in read.cigartuples)):
                    whole+=1
                    spans.append(dict(run=run,region_id=r['region_id'],min_MAPQ=q,read_name=read.query_name,
                         alignment_start=read.reference_start+1,alignment_end=read.reference_end,MAPQ=read.mapping_quality,
                         left_anchor=r['left_anchor'],right_anchor=r['right_anchor'],span_start=x+1,span_end=y,
                         aligned_fraction=overlap(blocks,x,y)/(y-x)))
        depth=np.cumsum(delta[:-1]);middle=depth[a-lo:b-lo];l=depth[:a-lo];right_depth=depth[b-lo:]
        output.append(dict(run=run,assembly=meta['assembly'],strain=meta['strain'],platform=meta['metadata']['instrument_platform'],
            region_id=r['region_id'],region_kind=r['kind'],contig=r['contig'],start=r['start'],end=r['end'],min_MAPQ=q,
            genome_median_depth=median,locus_mean_depth=float(np.mean(middle)),locus_median_depth=float(np.median(middle)),
            depth_ratio=float(np.mean(middle))/median if median else None,covered_fraction=float(np.mean(middle>0)),
            left_flank_median=float(np.median(l)),right_flank_median=float(np.median(right_depth)),
            left_edge_spanning_reads=left,right_edge_spanning_reads=right,min_aligned_bp_each_40bp_side=minimum,
            overlapping_reads=total_reads,median_MAPQ=float(np.median(qualities)) if qualities else None,
            properly_paired_fraction=proper/total_reads if total_reads and meta['mapping_preset']=='sr' else None,
            fungal_anchor_spanning_reads=whole if r['left_anchor_end'] is not None else None,
            left_fungal_anchor=r.get('left_anchor',''),right_fungal_anchor=r.get('right_anchor','')))
def write(name,rows):
    if not rows:return
    with open(wd/name,'w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
write('locus_read_support.tsv',output);write('fungal_anchor_spanning_reads.tsv',spans);write('genome_coverage.tsv',genome)
(wd/'summary.json').write_text(json.dumps(dict(regions=output,genome=genome),indent=2)+'\n')
print(json.dumps(genome),flush=True)
