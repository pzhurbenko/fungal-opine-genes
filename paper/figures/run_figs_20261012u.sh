#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../.." && pwd)}"; export GMO_ROOT
set -eo pipefail
source $(conda info --base)/etc/profile.d/conda.sh; conda activate PlantCad2
cd ${GMO_ROOT}/paper/figures
python make_fig1.py && echo FIG1_OK
python make_fig3.py && echo FIG3_OK
python make_figS_trees.py NOS && echo FIGS6_OK
