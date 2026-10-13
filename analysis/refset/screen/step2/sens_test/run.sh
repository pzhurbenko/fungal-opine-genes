#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../../../../.." && pwd)}"; export GMO_ROOT
source $(conda info --base)/etc/profile.d/conda.sh; conda activate annot
D=${GMO_ROOT}/analysis/refset/screen/step2
cd $D/sens_test
for G in GCA_019364595.1 GCA_018397855.1; do
  FNA=$(awk -v g=$G -F'\t' '$1==g{print $2}' $D/genomes.tsv)
  for M in very-sensitive ultra-sensitive; do
    T0=$(date +%s)
    diamond blastx -q $FNA -d $D/queries.dmnd -o ${G}.${M}.tsv --threads 2 --$M -F 15 --range-culling --top 10 --evalue 1e-5 --quiet \
      --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qframe
    echo "$G $M hits=$(wc -l < ${G}.${M}.tsv) t=$(( $(date +%s) - T0 ))s" >> timing.txt
  done
done
echo SENS_DONE >> timing.txt
