#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../../.." && pwd)}"; export GMO_ROOT
# usage: run_au.sh constrain <H>   |   run_au.sh au
# Same scheme as nos_check/run_au.sh: constrained ML (-g), then -z -n 0 -zb 10000 -au.
set -eo pipefail
source $(conda info --base)/etc/profile.d/conda.sh; conda activate phylo_env
D=${GMO_ROOT}/analysis/fusarium_odh_island/au; cd $D
ALN=../tree/odh_tree.trim.aln
MODEL=$(grep -oP 'Best-fit model according to BIC: \K\S+' ../tree/odh_tree.iqtree.iqtree | head -1); MODEL=${MODEL:-Q.PFAM+G4}
T=${T:-32}
HS="${HS:-H1_isl_native H1b_isl_sordario H1c_native_only H2_isl_tdna H2s_isl_seeds H3_isl_larry H4_isl_bact H5_isl_nos}"
if [ "$1" = constrain ]; then
  H=$2; T0=$(date +%s)
  iqtree -s $ALN -m $MODEL -g $H.nwk -T $T -seed 1 -pre $H -redo -quiet
  echo "$H $(( $(date +%s) - T0 ))s" >> timing.txt
elif [ "$1" = au ]; then
  cat ../tree/odh_tree.iqtree.treefile > au_trees.nwk
  echo "ML $HS" | tr ' ' '\n' > au_trees.order.txt
  for H in $HS; do
    grep -q 'Total wall-clock time' $H.iqtree || { echo "ERROR: $H.iqtree not finished"; exit 1; }
    cat $H.treefile >> au_trees.nwk
  done
  T0=$(date +%s)
  iqtree -s $ALN -m $MODEL -z au_trees.nwk -n 0 -zb 10000 -au -T $T -seed 1 -pre au_test -redo -quiet
  echo "au $(( $(date +%s) - T0 ))s; AU_DONE $(date)" >> timing.txt
fi
