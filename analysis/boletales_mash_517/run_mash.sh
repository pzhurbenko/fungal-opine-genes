#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../.." && pwd)}"; export GMO_ROOT
source $(conda info --base)/etc/profile.d/conda.sh; conda activate mashtree_env
cd ${GMO_ROOT}/analysis/boletales_mash_517
cut -f4 genomes.tsv > list.txt
mash sketch -p 64 -s 20000 -k 21 -o all -l list.txt > sketch.log 2>&1
mash dist -p 64 all.msh all.msh > dist.tsv 2> dist.log
echo "MASH_DONE $(date)" >> timing.txt
