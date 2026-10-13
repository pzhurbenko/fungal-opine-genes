#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../.." && pwd)}"; export GMO_ROOT
set -eo pipefail
source $(conda info --base)/etc/profile.d/conda.sh; conda activate phylo_env
cd ${GMO_ROOT}/analysis/chs_fusarium; T=16
mafft --localpair --maxiterate 1000 --thread $T --quiet chs_tree_set.faa > chs_tree.mafft.aln
clipkit chs_tree.mafft.aln -m smart-gap -o chs_tree.trim.aln -l > clipkit.log
iqtree -s chs_tree.trim.aln -m MFP -B 1000 -alrt 1000 -bnni -T $T -seed 1 -pre chs_tree.iqtree -redo -quiet
