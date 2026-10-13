#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../.." && pwd)}"; export GMO_ROOT
# usage: map_rna.sh <SRR> <GCA>   minimap2 splice-aware mapping of RNA-seq to assembly; per-contig and per-locus counts
SRR=$1; GCA=$2; T=32
source $(conda info --base)/etc/profile.d/conda.sh
cd ${GMO_ROOT}/analysis/expression
FNA=$(ls ${GMO_ROOT}/data/fungal_genomes/step2/ncbi_dataset/data/$GCA/*.fna)
conda activate minimap2_env
minimap2 -ax splice -uf -k14 -t $T "$FNA" reads/${SRR}_1.fastq.gz reads/${SRR}_2.fastq.gz 2> bam/$SRR.$GCA.log | $(conda info --base)/envs/annot/bin/samtools sort -@ 8 -m 1500M -o bam/$SRR.$GCA.bam - 2>/dev/null
conda activate annot
samtools index bam/$SRR.$GCA.bam; samtools flagstat bam/$SRR.$GCA.bam > bam/$SRR.$GCA.flagstat
echo "MAP_DONE $SRR $GCA $(date)" >> timing.txt
