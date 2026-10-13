#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../.." && pwd)}"; export GMO_ROOT
# usage: au_lineage_acs.sh constrain <H> | au_lineage_acs.sh au   (run in trees_final/CUS)
source $(conda info --base)/etc/profile.d/conda.sh; conda activate phylo_env
D=${GMO_ROOT}/analysis/trees_final/CUS; cd $D
MODEL=$(grep -oP 'Best-fit model according to BIC: \K\S+' CUS.mafft.iqtree.iqtree | head -1); MODEL=${MODEL:-JTT+R7}
T=${T:-32}
if [ "$1" = constrain ]; then
  H=$2; T0=$(date +%s)
  iqtree -s CUS.mafft.trim.aln -m $MODEL -g CUS.$H.constraint.nwk -T $T -pre CUS.$H -redo -quiet
  echo "CUS lineage_$H $(( $(date +%s) - T0 ))s" >> timing.txt
elif [ "$1" = au ]; then
  HS="${HS:-cus_tdna cus_tdna4 cus_act cus_mis cus_bact cus_mono bol_tdna bol_act bol_mis bol_othercus}"
  SFX=${SFX:-}; cat CUS.mafft.iqtree.treefile > CUS.lineage$SFX.candidates.nwk
  for H in $HS; do cat CUS.$H.treefile >> CUS.lineage$SFX.candidates.nwk; done
  echo "order: ML $HS" > CUS.lineage$SFX.AU.order.txt
  T0=$(date +%s)
  iqtree -s CUS.mafft.trim.aln -m $MODEL -z CUS.lineage$SFX.candidates.nwk -n 0 -zb 10000 -au -T $T -pre CUS.lineage$SFX.AU -redo -quiet
  echo "CUS lineage${SFX}_au $(( $(date +%s) - T0 ))s; LINEAGE_AU${SFX}_DONE $(date)" >> timing.txt
fi
