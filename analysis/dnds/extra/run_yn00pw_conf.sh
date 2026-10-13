#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../../.." && pwd)}"; export GMO_ROOT
source $(conda info --base)/etc/profile.d/conda.sh; conda activate PlantCad2
cd ${GMO_ROOT}/analysis/dnds/extra
for s in cus_Bol_conf cus_Bol_sp_conf acs_pKerr14a_conf nosODH_Bol_conf; do python3 yn00_pairwise_loose.py $s & done; wait
