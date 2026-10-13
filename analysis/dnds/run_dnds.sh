#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../.." && pwd)}"; export GMO_ROOT
source $(conda info --base)/etc/profile.d/conda.sh
cd ${GMO_ROOT}/analysis/dnds
for FAM in SUS CUS; do
  conda activate annot
  mafft --localpair --maxiterate 1000 --quiet --thread 8 $FAM.pep.faa > $FAM.pep.aln
  python3 codon.py $FAM
  conda activate phylo_env
  iqtree -s $FAM.codon.fna -st DNA -m GTR+G -T 8 -pre $FAM.tree -redo -quiet
  python3 - <<PY
from Bio import Phylo
t=Phylo.read('$FAM.tree.treefile','newick')
for c in t.find_clades(): c.confidence=None
Phylo.write(t,'$FAM.tree.nwk','newick')
PY
  conda activate PlantCad2
  for M in M0 M7M8 M8a; do python3 ctl.py $FAM $M; codeml $FAM.$M.ctl > $FAM.$M.log 2>&1; done
  echo "seqfile = $FAM.codon.phy
outfile = $FAM.yn00.out
verbose = 0
icode = 0
weighting = 0
commonf3x4 = 0" > $FAM.yn00.ctl; yn00 $FAM.yn00.ctl > $FAM.yn00.log 2>&1
done
echo "DNDS_DONE $(date)" >> timing.txt
