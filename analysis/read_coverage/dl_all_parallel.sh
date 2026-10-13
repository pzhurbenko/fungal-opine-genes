#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../.." && pwd)}"; export GMO_ROOT
cd ${GMO_ROOT}/analysis/read_coverage
while IFS=$'\t' read -r gca sp srr fna ftp byt; do
  ( bash ./download.sh $srr > lsf_dl_$srr.log 2>&1 && bsub -q medium -n 32 -R 'span[hosts=1]' -R 'rusage[mem=1500]' -J map_$gca -e lsf_map_$gca.err -o lsf_map_$gca.out "bash $PWD/map_one.sh $gca $srr" >> lsf_dl_$srr.log 2>&1 ) &
done < links.tsv
wait; echo "ALL_DL_DONE $(date)" >> timing.txt
