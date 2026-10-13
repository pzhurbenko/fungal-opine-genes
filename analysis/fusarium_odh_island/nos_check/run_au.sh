#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../../.." && pwd)}"; export GMO_ROOT
set -eo pipefail
source $(conda info --base)/etc/profile.d/conda.sh; conda activate phylo_env
cd ${GMO_ROOT}/analysis/fusarium_odh_island/nos_check; T=16
for c in c1_bol_in_tdna c2_rh_in_tdna c3_fungal_mono; do
  iqtree -s nos_tree.trim.aln -m Q.PFAM+R4 -g $c.nwk -T $T -seed 1 -pre $c -redo -quiet
done
cat nos_tree.iqtree.treefile c1_bol_in_tdna.treefile c2_rh_in_tdna.treefile c3_fungal_mono.treefile > au_trees.nwk
iqtree -s nos_tree.trim.aln -m Q.PFAM+R4 -z au_trees.nwk -n 0 -zb 10000 -au -T $T -seed 1 -pre au_test -redo -quiet
echo AU_DONE
