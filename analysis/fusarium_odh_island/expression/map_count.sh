#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../../.." && pwd)}"; export GMO_ROOT
# usage: map_count.sh RUN -- minimap2 splice:sr (single-end) on F. poae Fp133 (GCA_046056315.1) -> sorted BAM -> count_expr.py per island gene
RUN=$1; T=24
B=${GMO_ROOT}/analysis/fusarium_odh_island/expression; cd $B
D=${GMO_ROOT}/data/fungal_genomes/fusarium/fus/ncbi_dataset/data/GCA_046056315.1; FNA=$D/GCA_046056315.1_ASM4605631v1_genomic.fna
source $(conda info --base)/etc/profile.d/conda.sh
conda activate minimap2_env
minimap2 -ax splice:sr -G 3000 -t $T $FNA reads/${RUN}.fastq.gz 2> bam/$RUN.minimap2.log | $(conda info --base)/envs/annot/bin/samtools sort -@ 8 -m 1G -o bam/$RUN.bam - 2>/dev/null
conda activate annot
samtools index bam/$RUN.bam; samtools flagstat bam/$RUN.bam > bam/$RUN.flagstat
# island copy on CM100256 (odh FPOAC2_10362, pepP 10363, g3 10364) and annotated fragments of the CM100255 copy (pepP 07383, g3 07382)
for TG in FPOAC2_10362 FPOAC2_10363 FPOAC2_10364 FPOAC2_07383 FPOAC2_07382; do
  $(conda info --base)/envs/PlantCad2/bin/python -I count_expr.py $RUN fp133 bam/$RUN.bam $D/genomic.gff $TG ctrl_fp133.tsv
done
# unannotated odh copy on CM100255 (9465822-9466883) and the whole CM100255 island: raw read counts, MAPQ>=1 and all primary
for R in CM100255.1:9465822-9466883 CM100255.1:9463190-9466883 CM100256.1:8420071-8423766; do
  echo -e "$RUN\t$R\tmapq1=$(samtools view -c -F 0x904 -q 1 bam/$RUN.bam $R)\tall=$(samtools view -c -F 0x904 bam/$RUN.bam $R)" >> region_counts.tsv
done
echo "DONE $RUN $(date)"
