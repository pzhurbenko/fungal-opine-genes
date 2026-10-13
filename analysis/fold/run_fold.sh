#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../.." && pwd)}"; export GMO_ROOT
source $(conda info --base)/etc/profile.d/conda.sh; conda activate colabfold
D=${GMO_ROOT}/analysis/fold; cd $D
echo "=== FOLD START $(date)"
colabfold_batch --num-models 1 --num-recycle 3 --msa-mode mmseqs2_uniref_env --model-type auto \
  input/all.fasta out > colabfold.log 2>&1
echo "=== FOLD DONE $(date): $(ls out/*.pdb 2>/dev/null | wc -l) структур"
