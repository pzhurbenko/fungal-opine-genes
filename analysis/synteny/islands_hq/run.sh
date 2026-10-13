#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../../.." && pwd)}"; export GMO_ROOT
source $(conda info --base)/etc/profile.d/conda.sh; conda activate annot
cd ${GMO_ROOT}/analysis/synteny/islands_hq; mkdir -p gff
while IFS=$'\t' read -r g fna; do miniprot -t 16 --gff-only --outs 0.5 -p 0.3 "$fna" q/$g.faa > gff/$g.gff 2> gff/$g.log; done < targets.tsv
python3 summarize.py && echo "ISL_SYN_DONE $(date)" >> timing.txt
