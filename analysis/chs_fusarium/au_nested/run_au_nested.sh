#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../../.." && pwd)}"; export GMO_ROOT
set -eo pipefail
source $(conda info --base)/etc/profile.d/conda.sh; conda activate phylo_env
cd ${GMO_ROOT}/analysis/chs_fusarium/au_nested; T=16
M=$(grep -oP 'Best-fit model according to BIC: \K\S+' ../chs_tree.iqtree.iqtree)
HS=$(cut -d' ' -f2- hyp_order.txt)
for h in $HS; do iqtree -s ../chs_tree.trim.aln -m $M -g $h.nwk -T $T -seed 1 -pre $h -redo -quiet; done
cat ../chs_tree.iqtree.treefile > au_trees.nwk; for h in $HS; do cat $h.treefile >> au_trees.nwk; done
iqtree -s ../chs_tree.trim.aln -m $M -z au_trees.nwk -n 0 -zb 10000 -au -T $T -seed 1 -pre au_nested -redo -quiet
echo AU_DONE
