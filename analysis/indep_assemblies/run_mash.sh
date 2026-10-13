#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../.." && pwd)}"; export GMO_ROOT
source $(conda info --base)/etc/profile.d/conda.sh; conda activate mashtree_env
cd ${GMO_ROOT}/analysis/indep_assemblies; mkdir -p mash
cut -f3 genomes.tsv > mash/list.txt
mash sketch -p 16 -s 20000 -k 21 -o mash/all -l mash/list.txt > mash/sketch.log 2>&1
mash dist -p 16 mash/all.msh mash/all.msh > mash/dist.tsv 2> mash/dist.log
echo "MASH_DONE $(date)" >> timing.txt
