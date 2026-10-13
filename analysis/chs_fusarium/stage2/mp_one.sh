#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../../.." && pwd)}"; export GMO_ROOT
# miniprot of 3 chs-like seeds on ONE whole genome. Usage: mp_one.sh <list.tsv: GCA fna gff>  (index = LSB_JOBINDEX)
source $(conda info --base)/etc/profile.d/conda.sh; conda activate annot
S=${GMO_ROOT}/analysis/chs_fusarium/stage2
LINE=$(sed -n "${LSB_JOBINDEX}p" $1); GCA=$(echo "$LINE" | cut -f1); FNA=$(echo "$LINE" | cut -f2)
O=$S/mp_out/$GCA; mkdir -p $O; [ -s $O/done ] && exit 0
seqkit fx2tab -n -l -g -i $FNA > $O/contigs.tsv
for sd in chs_A0A2Z2PGA9 EEU41297.1 KNB14377.1; do
  seqkit grep -p $sd $S/seeds3.faa > $O/$sd.faa
  miniprot -t 4 -G 3000 -p 0.3 --outs 0.3 --outc 0.2 -N 20 --gff --trans $FNA $O/$sd.faa > $O/$sd.gff 2> $O/$sd.log || exit 1
done
echo ok > $O/done
