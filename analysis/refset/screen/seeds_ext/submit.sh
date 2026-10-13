#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../../../.." && pwd)}"; export GMO_ROOT
# Submit LSF array: 6847 proteomes / PER=50 -> 137 jobs (limit 2000). Run from anywhere.
S=${GMO_ROOT}/analysis/refset/screen/seeds_ext
N=$(wc -l < ${GMO_ROOT}/analysis/refset/screen/step3/proteomes.txt); PER=50
J=$(( (N+PER-1)/PER ))
mkdir -p $S/logs
bsub -q short -J "seedsext[1-$J]" -n 8 -R "span[hosts=1]" -R "rusage[mem=1000]" -o $S/logs/%J_%I.out -e $S/logs/%J_%I.err \
  "PER=$PER THR=8 bash $S/screen_one.sh 0"
