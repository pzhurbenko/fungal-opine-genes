#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../.." && pwd)}"; export GMO_ROOT
# usage: download.sh <SRR>  (runs on transfer queue) — fetch ENA fastq.gz files listed in links.tsv
SRR=$1; cd ${GMO_ROOT}/analysis/read_coverage/reads
FTP=$(awk -F'\t' -v s=$SRR '$3==s{print $5}' ../links.tsv)
for u in ${FTP//;/ }; do wget -q -c "ftp://$u" || wget -q -c "https://$u" || exit 1; done
ls -la ${SRR}*; echo "DL_DONE $SRR $(date)" >> ../timing.txt
