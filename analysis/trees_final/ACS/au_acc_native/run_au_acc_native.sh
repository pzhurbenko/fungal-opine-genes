#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../../../.." && pwd)}"; export GMO_ROOT
set -e
source $(conda info --base)/etc/profile.d/conda.sh; conda activate phylo_env
D=${GMO_ROOT}/analysis/trees_final/ACS; cd $D/au_acc_native
MODEL=JTT+R5; T=32
T0=$(date +%s)
iqtree -s ../ACS.mafft.trim.aln -m $MODEL -g ACS.meso_natv.constraint.nwk -T $T -pre ACS.meso_natv -redo -quiet
echo "meso_natv constrain $(( $(date +%s) - T0 ))s" >> timing.txt
cat ../ACS.mafft.iqtree.treefile > au_trees.nwk
for H in six_tdna six_kerr six_acc six_natv meso_acc meso_tdna; do cat ../ACS.$H.treefile >> au_trees.nwk; done
cat ACS.meso_natv.treefile >> au_trees.nwk
echo "order: ML six_tdna six_kerr six_acc six_natv meso_acc meso_tdna meso_natv" > hyp_order.txt
T0=$(date +%s)
iqtree -s ../ACS.mafft.trim.aln -m $MODEL -z au_trees.nwk -n 0 -zb 10000 -au -T $T -pre au_acc_native -redo -quiet
echo "au $(( $(date +%s) - T0 ))s; DONE $(date)" >> timing.txt
