#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../.." && pwd)}"; export GMO_ROOT
# usage: map_one.sh <GCA> <SRR>  — minimap2 sr of reads to the WHOLE assembly; coverage summaries
GCA=$1; SRR=$2; T=${T:-32}
source $(conda info --base)/etc/profile.d/conda.sh
cd ${GMO_ROOT}/analysis/read_coverage; mkdir -p bam depth
FNA=$(awk -F'\t' -v g=$GCA '$1==g{print $4}' links.tsv)
R1=$(ls reads/${SRR}_1.fastq.gz); R2=$(ls reads/${SRR}_2.fastq.gz)
for f in $R1 $R2; do gzip -t $f || { echo "CORRUPT $f" >> timing.txt; exit 1; }; done
conda activate minimap2_env
minimap2 -ax sr -t $T "$FNA" $R1 $R2 2> bam/$GCA.minimap2.log | $(conda info --base)/envs/annot/bin/samtools sort -@ 8 -m 1500M -o bam/$GCA.bam - 2> bam/$GCA.sort.log
conda activate annot
samtools index bam/$GCA.bam
samtools flagstat bam/$GCA.bam > bam/$GCA.flagstat
samtools coverage bam/$GCA.bam > depth/$GCA.coverage.tsv
# per-base depth on registry contigs of this assembly
for c in $(awk -F'\t' -v g=$GCA '$12==g && $5!="not_homolog" && $5!="no_model" && $5!=""{print $14}' ../loci_pipeline/round2/loci_table.tsv | sort -u); do
  cid=$(grep -m1 -oP "^>\K${c}\S*" "$FNA"); cid=${cid:-$c}
  samtools depth -a -r $cid bam/$GCA.bam > depth/$GCA.$cid.depth
done
echo "MAP_DONE $GCA $(date)" >> timing.txt
