#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../.." && pwd)}"; export GMO_ROOT
cd ${GMO_ROOT}/analysis/expression/reads
while IFS=$'\t' read -r srr sp gca ftp byt name title; do
  ( for u in ${ftp//;/ }; do wget -q "ftp://$u" || wget -q "https://$u"; done; for f in ${srr}_*.fastq.gz; do gzip -t $f || echo "CORRUPT $f" >> ../timing.txt; done; echo "DL_DONE $srr $(date)" >> ../timing.txt ) &
done < ../links.tsv
wait; echo "ALL_DL_DONE $(date)" >> ../timing.txt
