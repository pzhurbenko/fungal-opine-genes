#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../../.." && pwd)}"; export GMO_ROOT
# usage: dl_run.sh RUN -- first 10M single-end reads from ENA (URL from runs.tsv), re-gzipped; gzip -t check
R=$1; B=${GMO_ROOT}/analysis/fusarium_odh_island/expression; cd $B/reads
U=$(awk -F'\t' -v r=$R '$1==r{print $2}' $B/runs.tsv)
curl -s -f --retry 3 "https://$U" | gzip -dc 2>/dev/null | head -n 40000000 | gzip -1 > ${R}.fastq.gz
gzip -t ${R}.fastq.gz && echo "OK $R lines=$(zcat ${R}.fastq.gz | wc -l)"
