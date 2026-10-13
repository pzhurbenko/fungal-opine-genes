#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../../.." && pwd)}"; export GMO_ROOT
source $(conda info --base)/etc/profile.d/conda.sh; conda activate phylo_env
cd ${GMO_ROOT}/analysis/fold/foldmason; mkdir -p pdbs tmp
for f in ../out/*rank_001*.pdb; do b=$(basename $f); ln -sf $(readlink -f $f) pdbs/${b%%_unrelaxed*}.pdb; done
ls pdbs | wc -l
foldmason easy-msa pdbs all_struct tmp --threads 32 --report-mode 1 > foldmason.log 2>&1
echo "FOLDMASON_DONE $(date)" >> timing.txt
