#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../../.." && pwd)}"; export GMO_ROOT
source $(conda info --base)/etc/profile.d/conda.sh
cd ${GMO_ROOT}/analysis/dnds/loose
for FAM in SUS CUS; do
  cp ../$FAM.cds.fna ../$FAM.pep.aln ../$FAM.tree.nwk .; conda activate annot; python3 codon.py $FAM
  conda activate paper_build; python3 ../ctl.py $FAM M0; sed -i 's/cleandata = 0/cleandata = 0/' $FAM.M0.ctl; codeml $FAM.M0.ctl > $FAM.M0.log 2>&1
  printf "seqfile = $FAM.codon.phy\noutfile = $FAM.yn00.out\nverbose = 0\nicode = 0\nweighting = 0\ncommonf3x4 = 0\n" > $FAM.yn00.ctl; yn00 $FAM.yn00.ctl > $FAM.yn00.log 2>&1
  cp 2ML.dN $FAM.2ML.dN; cp 2ML.dS $FAM.2ML.dS
done
echo "LOOSE_DONE $(date)" >> timing.txt
