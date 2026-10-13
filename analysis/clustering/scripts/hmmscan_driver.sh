#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../../.." && pwd)}"; export GMO_ROOT
set -uo pipefail
source $(conda info --base)/etc/profile.d/conda.sh
conda activate hmmer_pfam
PFAM=${GMO_ROOT}/analysis/pfam/Pfam-A.hmm
CL=${GMO_ROOT}/analysis/clustering
OUT=${GMO_ROOT}/analysis/domains
mkdir -p $OUT
for fam in ACS CUS ORF358 SUS; do
  echo "########## hmmscan $fam ##########"
  hmmscan --cpu 16 --cut_ga \
     --domtblout $OUT/${fam}.domtbl -o $OUT/${fam}.hmmscan.txt \
     $PFAM $CL/$fam/clean.fasta
  echo "  done: $(grep -vc '^#' $OUT/${fam}.domtbl) domain hits"
done
echo "===== HMMSCAN DONE ====="
