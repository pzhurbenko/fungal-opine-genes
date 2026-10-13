#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../../.." && pwd)}"; export GMO_ROOT
# usage: map_count.sh RUN SP(fol|van)  -- minimap2 splice:sr (short-read paired) -> sorted BAM -> count_expr.py
RUN=$1; SP=$2; T=24
B=${GMO_ROOT}/analysis/chs_fusarium/expression; cd $B
source $(conda info --base)/etc/profile.d/conda.sh
if [ $SP = fol ]; then D=genome/fol/ncbi_dataset/data/GCF_000149955.1; FNA=$D/GCF_000149955.1_ASM14995v2_genomic.fna; TG=FOXG_13246
else D=genome/van/ncbi_dataset/data/GCF_000151355.1; FNA=$D/GCF_000151355.1_v2.0_genomic.fna; TG=NECHADRAFT_83559; fi
conda activate minimap2_env
minimap2 -ax splice:sr -G 3000 -t $T $FNA reads/${RUN}_1.fastq.gz reads/${RUN}_2.fastq.gz 2> bam/$RUN.minimap2.log | $(conda info --base)/envs/annot/bin/samtools sort -@ 8 -m 1G -o bam/$RUN.bam - 2>/dev/null
conda activate annot
samtools index bam/$RUN.bam; samtools flagstat bam/$RUN.bam > bam/$RUN.flagstat
$(conda info --base)/envs/PlantCad2/bin/python -I count_expr.py $RUN $SP bam/$RUN.bam $D/genomic.gff $TG genome/ctrl_$SP.tsv
echo "DONE $RUN $(date)"
