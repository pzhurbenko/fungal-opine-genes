"""Export CIGAR/placement of the independently counted Fp133 molecules."""
from pathlib import Path
import csv,hashlib,json
import pysam
HERE=Path(__file__).resolve().parent;D=HERE/'SRR23927956'
evidence=list(csv.DictReader(open(D/'fungal_anchor_spanning_reads.tsv'),delimiter='\t'))
chosen=[r for r in evidence if r['min_MAPQ']=='20']
assert len({r['read_name'] for r in chosen})==len(chosen)
bam=pysam.AlignmentFile(D/'reads.bam','rb');rows=[]
for r in chosen:
    contig={'Fp133_chr2_odh_block':'CM100255.1','Fp133_chr3_odh_block':'CM100256.1'}[r['region_id']]
    hit=[x for x in bam.fetch(contig,int(r['span_start'])-1,int(r['span_end'])) if x.query_name==r['read_name'] and not x.flag&0xF04]
    assert len(hit)==1
    x=hit[0]
    rows.append(dict(region=r['region_id'],read=x.query_name,flag=x.flag,reference=x.reference_name,
        start_1based=x.reference_start+1,end_1based=x.reference_end,MAPQ=x.mapping_quality,CIGAR=x.cigarstring,
        query_length=x.query_length,aligned_query_start_0based=x.query_alignment_start,aligned_query_end_0based=x.query_alignment_end,
        NM=x.get_tag('NM'),query_sequence_sha256=hashlib.sha256(x.query_sequence.encode()).hexdigest()))
with open(D/'spanning_alignments.tsv','w') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
print(len(rows),'unique long-read alignments exported')
