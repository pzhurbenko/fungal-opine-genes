#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../../../.." && pwd)}"; export GMO_ROOT
set -x
source $(conda info --base)/etc/profile.d/conda.sh; conda activate phylo_env
cd ${GMO_ROOT}/analysis/trees_final/SUS/au_sus_native
ALN=../SUS.mafft.trim.aln; MODEL=Q.PFAM+R5
for H in sus_natv_core sus_natv_all; do
  T0=$(date +%s); iqtree -s $ALN -m $MODEL -g SUS.$H.constraint.nwk -T 16 -seed 2 -pre SUS.$H -redo -quiet > /dev/null 2>&1
  echo "$H $(( $(date +%s) - T0 ))s" >> timing.txt
done
cat ../SUS.mafft.iqtree.treefile ../SUS.H_tdna.treefile ../SUS.H_alt.treefile SUS.sus_natv_core.treefile SUS.sus_natv_all.treefile > au_trees.nwk
echo "ML H_tdna H_alt sus_natv_core sus_natv_all" > hyp_order.txt
T0=$(date +%s); iqtree -s $ALN -m $MODEL -z au_trees.nwk -n 0 -zb 10000 -au -T 16 -seed 2 -pre au_sus_native -redo -quiet > /dev/null 2>&1
echo "au $(( $(date +%s) - T0 ))s" >> timing.txt
