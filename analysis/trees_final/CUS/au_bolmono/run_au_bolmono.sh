#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../../../.." && pwd)}"; export GMO_ROOT
set -eo pipefail
source $(conda info --base)/etc/profile.d/conda.sh; conda activate phylo_env
cd ${GMO_ROOT}/analysis/trees_final/CUS/au_bolmono; T=32
MODEL=$(grep -oP 'Best-fit model according to BIC: \K\S+' ../CUS.mafft.iqtree.iqtree | head -1); MODEL=${MODEL:-JTT+R7}
T0=$(date +%s)
iqtree -s ../CUS.mafft.trim.aln -m $MODEL -g CUS.bol_mono_x.constraint.nwk -T $T -pre CUS.bol_mono_x -redo -quiet
echo "bol_mono_x $(( $(date +%s) - T0 ))s" >> timing.txt
cat ../CUS.mafft.iqtree.treefile ../CUS.bol_tdna_x.treefile ../CUS.bol_act_x.treefile ../CUS.bol_mis_x.treefile ../CUS.bol_othercus_x.treefile CUS.bol_mono_x.treefile > CUS.lineage_bm.candidates.nwk
iqtree -s ../CUS.mafft.trim.aln -m $MODEL -z CUS.lineage_bm.candidates.nwk -n 0 -zb 10000 -au -T $T -pre CUS.lineage_bm.AU -redo -quiet
echo "AU_DONE $(date)" >> timing.txt
