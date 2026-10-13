#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../../.." && pwd)}"; export GMO_ROOT
source $(conda info --base)/etc/profile.d/conda.sh; conda activate phylo_env
cd ${GMO_ROOT}/analysis/fold/foldmason
iqtree -s all_struct_aa.fa -m MFP -B 1000 -alrt 1000 -T 32 -pre struct_aa -redo -quiet
iqtree -s all_struct_3di.fa -m 3DI -B 1000 -alrt 1000 -T 32 -pre struct_3di -redo -quiet
echo "STRUCT_TREE_DONE $(date)" >> timing.txt
