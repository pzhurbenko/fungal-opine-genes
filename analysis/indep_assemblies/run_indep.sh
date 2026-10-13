#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../.." && pwd)}"; export GMO_ROOT
# Point 3 (plan A): reproducibility of T-DNA-like loci in INDEPENDENT assemblies of the same species.
# (a) whole registry contigs (gene + flanks) of the source assembly -> minimap2 asm20 vs each other assembly of the species
# (b) all Boletales gene models (4 families) -> miniprot vs every target assembly (presence/absence per lineage)
set -u
source $(conda info --base)/etc/profile.d/conda.sh
D=${GMO_ROOT}/analysis/indep_assemblies; cd $D
T=${T:-16}
conda activate annot; python3 prepare.py || exit 1          # writes species.tsv, contigs/<species>.fa, boletales_models.faa, genomes.tsv
conda activate minimap2_env
mkdir -p paf gff
while IFS=$'\t' read -r sp src tgt fna; do
  [ "$src" = "$tgt" ] && continue
  minimap2 -x asm20 -c -t $T --secondary=no "$fna" contigs/${sp}__${src}.fa > paf/${sp}__${src}__vs__${tgt}.paf 2> paf/${sp}__${src}__vs__${tgt}.log
done < pairs.tsv
conda activate annot
while IFS=$'\t' read -r gca sp fna; do
  miniprot -t $T --gff-only --outs 0.5 -p 0.3 "$fna" boletales_models.faa > gff/${gca}.gff 2> gff/${gca}.log
done < genomes.tsv
python3 summarize.py && echo "INDEP_DONE $(date)" >> timing.txt
