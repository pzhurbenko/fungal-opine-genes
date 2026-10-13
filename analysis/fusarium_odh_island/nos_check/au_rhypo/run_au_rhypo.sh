#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../../../.." && pwd)}"; export GMO_ROOT
set -eo pipefail
source $(conda info --base)/etc/profile.d/conda.sh; conda activate phylo_env
cd ${GMO_ROOT}/analysis/fusarium_odh_island/nos_check/au_rhypo; T=32
for c in c6_rh_2sord c7_rh_sord c8_rh_nat; do iqtree -s ../nos_tree.trim.aln -m Q.PFAM+R4 -g $c.nwk -T $T -seed 1 -pre $c -redo -quiet; done
cat ../nos_tree.iqtree.treefile ../c1_bol_in_tdna.treefile ../c2_rh_in_tdna.treefile ../c3_fungal_mono.treefile ../au_native/c4_bol_nat.treefile ../au_native/c5_bol_rh_nat.treefile c6_rh_2sord.treefile c7_rh_sord.treefile c8_rh_nat.treefile > au_trees.nwk
echo "ML c1_bol_in_tdna c2_rh_in_tdna c3_fungal_mono c4_bol_nat c5_bol_rh_nat c6_rh_2sord c7_rh_sord c8_rh_nat" > hyp_order.txt
iqtree -s ../nos_tree.trim.aln -m Q.PFAM+R4 -z au_trees.nwk -n 0 -zb 10000 -au -T $T -seed 1 -pre au_rhypo -redo -quiet
echo AU_DONE
