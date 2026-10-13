#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../../.." && pwd)}"; export GMO_ROOT
set -eo pipefail
source $(conda info --base)/etc/profile.d/conda.sh; conda activate phylo_env
cd ${GMO_ROOT}/analysis/fusarium_odh_island/tree; T=16
mafft --localpair --maxiterate 1000 --thread $T --quiet odh_tree_set.faa > odh_tree.mafft.aln
clipkit odh_tree.mafft.aln -m smart-gap -o odh_tree.trim.aln -l > clipkit.log
iqtree -s odh_tree.trim.aln -m MFP -B 1000 -alrt 1000 -bnni -T $T -seed 1 -pre odh_tree.iqtree -redo -quiet
echo TREE_DONE
