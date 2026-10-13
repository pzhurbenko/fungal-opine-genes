#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../.." && pwd)}"; export GMO_ROOT
# usage: au_lineage_acs.sh constrain <H> | au_lineage_acs.sh au   (run in trees_final/ACS)
source $(conda info --base)/etc/profile.d/conda.sh; conda activate phylo_env
D=${GMO_ROOT}/analysis/trees_final/ACS; cd $D
MODEL=$(grep -oP 'Best-fit model according to BIC: \K\S+' ACS.mafft.iqtree.iqtree | head -1); MODEL=${MODEL:-JTT+R5}
T=${T:-32}
if [ "$1" = constrain ]; then
  H=$2; T0=$(date +%s)
  iqtree -s ACS.mafft.trim.aln -m $MODEL -g ACS.$H.constraint.nwk -T $T -pre ACS.$H -redo -quiet
  echo "ACS lineage_$H $(( $(date +%s) - T0 ))s" >> timing.txt
elif [ "$1" = au ]; then
  HS="six_tdna six_kerr six_acc six_natv meso_acc meso_tdna"
  cat ACS.mafft.iqtree.treefile > ACS.lineage.candidates.nwk
  for H in $HS; do cat ACS.$H.treefile >> ACS.lineage.candidates.nwk; done
  echo "order: ML $HS" > ACS.lineage.AU.order.txt
  T0=$(date +%s)
  iqtree -s ACS.mafft.trim.aln -m $MODEL -z ACS.lineage.candidates.nwk -n 0 -zb 10000 -au -T $T -pre ACS.lineage.AU -redo -quiet
  echo "ACS lineage_au $(( $(date +%s) - T0 ))s; LINEAGE_AU_DONE $(date)" >> timing.txt
fi
