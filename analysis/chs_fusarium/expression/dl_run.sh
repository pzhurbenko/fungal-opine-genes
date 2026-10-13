#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../../.." && pwd)}"; export GMO_ROOT
# usage: dl_run.sh RUN  -- first 10M read pairs from ENA (https, URLs from selected_meta.tsv), re-gzipped; gzip -t check
R=$1; B=${GMO_ROOT}/analysis/chs_fusarium/expression; cd $B/reads
U=$(awk -F'\t' -v r=$R '$1==r{print $11}' $B/selected_meta.tsv)
for m in 1 2; do
  ( curl -s -f --retry 3 "https://$(echo $U | cut -d';' -f$m)" | gzip -dc 2>/dev/null | head -n 40000000 | gzip -1 > ${R}_$m.fastq.gz ) &
done; wait
for m in 1 2; do gzip -t ${R}_$m.fastq.gz && echo "OK $R _$m lines=$(zcat ${R}_$m.fastq.gz | wc -l)"; done
