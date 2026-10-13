#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../.." && pwd)}"; export GMO_ROOT
cd ${GMO_ROOT}/analysis/expression
until grep -q ALL_DL_DONE timing.txt 2>/dev/null; do sleep 60; done
while IFS=$'\t' read -r srr sp gca rest; do
  bsub -q medium -n 32 -R "span[hosts=1]" -R "rusage[mem=1500]" -J rna_$srr -e lsf_$srr.err -o lsf_$srr.out "bash $PWD/map_rna.sh $srr $gca" >> chain.log 2>&1
done < links.tsv
