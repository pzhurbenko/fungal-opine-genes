#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../../.." && pwd)}"; export GMO_ROOT
set -uo pipefail
source $(conda info --base)/etc/profile.d/conda.sh
conda activate orthofinder
SP=${GMO_ROOT}/analysis/clustering/scripts
BASE=${GMO_ROOT}
FDIR=$BASE/data/Fasta
OUT=$BASE/analysis/clustering
mkdir -p $OUT
fams=(ACS CUS ORF358 SUS)
files=(ACS_rep_seq_plant_agro.fasta CUS_rep_seq_plant_agro.fasta ORF358_rep_seq_plant_agro.fasta SUS_rep_seq_plant_agro.fasta)
for idx in ${!fams[@]}; do
  fam=${fams[$idx]}; in=$FDIR/${files[$idx]}
  wd=$OUT/$fam; mkdir -p $wd; cd $wd
  echo "########## $fam ##########"
  python $SP/prep.py "$in" clean.fasta idmap.tsv
  rm -rf tmp clu clu_* 
  mmseqs easy-cluster clean.fasta clu tmp -s 7.5 --cluster-steps 5 --cluster-reassign \
     --min-seq-id 0.0 -e 1e-2 -c 0.6 --cov-mode 0 --alignment-mode 3 --threads 16 > mmseqs.log 2>&1
  python $SP/clusterize.py clean.fasta clu_cluster.tsv $wd
  na=0
  for f in $wd/clusters/*.fasta; do [ -e "$f" ] || continue
     mafft --genafpair --maxiterate 1000 --thread 16 "$f" > "${f%.fasta}.aln" 2>/dev/null; na=$((na+1)); done
  echo "  mafft: aligned $na multi-clusters"
  python $SP/split_species.py clean.fasta $wd
  rm -rf $wd/of_input/OrthoFinder
  if orthofinder -f $wd/of_input -S diamond_ultra_sens -og -t 16 -a 4 > orthofinder.log 2>&1; then
     og=$(ls $wd/of_input/OrthoFinder/Results_*/Orthogroups/Orthogroups.txt 2>/dev/null | head -1)
     stats=$(ls $wd/of_input/OrthoFinder/Results_*/Comparative_Genomics_Statistics/Statistics_Overall.tsv 2>/dev/null | head -1)
     echo "  orthofinder OK"
     [ -n "$og" ] && echo "  #orthogroups: $(wc -l < $og)"
  else echo "  orthofinder FAILED (see $wd/orthofinder.log)"; fi
done
echo "===== ALL DONE ====="
