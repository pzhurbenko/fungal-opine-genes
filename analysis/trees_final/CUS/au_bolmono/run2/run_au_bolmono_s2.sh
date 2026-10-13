#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../../../../.." && pwd)}"; export GMO_ROOT
set -eo pipefail
source $(conda info --base)/etc/profile.d/conda.sh; conda activate phylo_env
cd ${GMO_ROOT}/analysis/trees_final/CUS/au_bolmono/run2; T=16
MODEL=JTT+R7; T0=$(date +%s)
iqtree -s ../../CUS.mafft.trim.aln -m $MODEL -g CUS.bol_mono_x.constraint.nwk -T $T -seed 2 -ninit 20 -pre CUS.bol_mono_x_s2 -redo -quiet
echo "bol_mono_x_s2 $(( $(date +%s) - T0 ))s" >> timing.txt
cat ../../CUS.mafft.iqtree.treefile ../../CUS.bol_tdna_x.treefile ../../CUS.bol_act_x.treefile ../../CUS.bol_mis_x.treefile ../../CUS.bol_othercus_x.treefile CUS.bol_mono_x_s2.treefile > CUS.lineage_bm_s2.candidates.nwk
iqtree -s ../../CUS.mafft.trim.aln -m $MODEL -z CUS.lineage_bm_s2.candidates.nwk -n 0 -zb 10000 -au -T $T -seed 1 -pre CUS.lineage_bm_s2.AU -redo -quiet
echo "AU_DONE $(date)" >> timing.txt
