#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../../.." && pwd)}"; export GMO_ROOT
# [*_conf sets 2026-10-13; identical to run_extra.sh, except N=3: star tree (IQ-TREE needs >=4 taxa), M0 only] Extra dN/dS sets: same procedure as ../run_dnds.sh (strict, no-gap codon columns) and ../loose/run.sh (<=20% gaps).
# n>=4: IQ-TREE GTR+G on strict codon alignment, codeml M0, M7/M8, M8a (strict) and M0 (loose, same tree); yn00 on both. n<4: yn00 only.
source $(conda info --base)/etc/profile.d/conda.sh
X=${GMO_ROOT}/analysis/dnds/extra; D=${GMO_ROOT}/analysis/dnds
yn() { printf "seqfile = $1.codon.phy\noutfile = $1.yn00.out\nverbose = 0\nicode = 0\nweighting = 0\ncommonf3x4 = 0\n" > $1.yn00.ctl; yn00 $1.yn00.ctl > $1.yn00.log 2>&1; }
one() {
  FAM=$1; cd $X/$FAM; N=$(grep -c ">" $FAM.cds.fna)
  conda activate annot
  mafft --localpair --maxiterate 1000 --quiet --thread 2 $FAM.pep.faa > $FAM.pep.aln
  python3 $D/codon.py $FAM
  mkdir -p loose; cp $FAM.cds.fna $FAM.pep.aln loose/; (cd loose; python3 $D/loose/codon.py $FAM)
  if [ $N -ge 4 ]; then
    conda activate phylo_env
    iqtree -s $FAM.codon.fna -st DNA -m GTR+G -T 2 -pre $FAM.tree -redo -quiet
    python3 -c "
from Bio import Phylo
t=Phylo.read('$FAM.tree.treefile','newick')
for c in t.find_clades(): c.confidence=None
Phylo.write(t,'$FAM.tree.nwk','newick')"
    cp $FAM.tree.nwk loose/
  elif [ $N -eq 3 ]; then
    python3 -c "
ids=[l[1:].strip() for l in open('$FAM.cds.fna') if l[0]=='>']
open('$FAM.tree.nwk','w').write('('+','.join(i+':0.1' for i in ids)+');\n')"
    cp $FAM.tree.nwk loose/
  fi
  conda activate PlantCad2
  yn $FAM; (cd loose; yn $FAM)
  if [ $N -ge 3 ]; then
    (cd loose; python3 $D/ctl.py $FAM M0; codeml $FAM.M0.ctl > $FAM.M0.log 2>&1)
    for M in $( [ $N -ge 4 ] && echo M0 M7M8 M8a || echo M0 ); do python3 $D/ctl.py $FAM $M; codeml $FAM.$M.ctl > $FAM.$M.log 2>&1; done
  fi
  echo "$FAM done N=$N $(date)"
}
for F in cus_Bol_conf cus_Bol_sp_conf acs_pKerr14a_conf nosODH_Bol_conf; do one $F > $X/$F.run.log 2>&1 & done
wait; echo "EXTRA_CONF_DONE $(date)" >> $X/timing.txt
