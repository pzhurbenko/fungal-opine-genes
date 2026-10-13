#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../.." && pwd)}"; export GMO_ROOT
# AU test for one family: constrained ML trees (H_tdna / H_alt) + iqtree -z -au against the unconstrained ML tree.
FAM=$1; T=${2:-32}
source $(conda info --base)/etc/profile.d/conda.sh
D=${GMO_ROOT}/analysis/trees_final/$FAM; cd $D
conda activate annot; python3 ../make_constraints.py $FAM ${FAM}.mafft.trim.aln
conda activate phylo_env
MODEL=$(grep -oP 'Best-fit model according to BIC: \K\S+' ${FAM}.mafft.iqtree.iqtree 2>/dev/null | head -1)
MODEL=${MODEL:-LG+F+G4}; echo "model: $MODEL" >> timing.txt
for H in H_tdna H_alt; do
  [ -s ${FAM}.$H.constraint.nwk ] || continue
  T0=$(date +%s); iqtree -s ${FAM}.mafft.trim.aln -m $MODEL -g ${FAM}.$H.constraint.nwk -T $T -pre ${FAM}.$H -redo -quiet > /dev/null 2>&1
  echo "$FAM iqtree_$H $(( $(date +%s) - T0 ))s" >> timing.txt
done
cat ${FAM}.mafft.iqtree.treefile ${FAM}.H_tdna.treefile ${FAM}.H_alt.treefile > ${FAM}.candidates.nwk 2>/dev/null
T0=$(date +%s); iqtree -s ${FAM}.mafft.trim.aln -m $MODEL -z ${FAM}.candidates.nwk -n 0 -zb 10000 -au -T $T -pre ${FAM}.AU -redo -quiet > /dev/null 2>&1
echo "$FAM au_test $(( $(date +%s) - T0 ))s" >> timing.txt
echo "$FAM FINAL_DONE $(date)" >> timing.txt
