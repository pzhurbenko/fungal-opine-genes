#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../../.." && pwd)}"; export GMO_ROOT
set -eo pipefail
source $(conda info --base)/etc/profile.d/conda.sh; conda activate phylo_env
cd ${GMO_ROOT}/analysis/chs_fusarium/au; T=16
M=$(grep -oP 'Best-fit model according to BIC: \K\S+' ../chs_tree.iqtree.iqtree)
for h in chs_tdna chs_bactchs chs_ags chs_nat chs_bact chs_mono; do iqtree -s ../chs_tree.trim.aln -m $M -g $h.nwk -T $T -seed 1 -pre $h -redo -quiet; done
cat ../chs_tree.iqtree.treefile chs_tdna.treefile chs_bactchs.treefile chs_ags.treefile chs_nat.treefile chs_bact.treefile chs_mono.treefile > au_trees.nwk
iqtree -s ../chs_tree.trim.aln -m $M -z au_trees.nwk -n 0 -zb 10000 -au -T $T -seed 1 -pre au_test -redo -quiet
echo AU_DONE
