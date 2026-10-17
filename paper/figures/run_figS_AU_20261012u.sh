#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../.." && pwd)}"; export GMO_ROOT
set -eo pipefail
source $(conda info --base)/etc/profile.d/conda.sh; conda activate paper_build
cd ${GMO_ROOT}/paper/figures
python make_figS_AU.py && echo FIGS8_OK
