#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/.." && pwd)}"; export GMO_ROOT
# Rebuild figures (300 dpi PNG + 600 dpi TIFF) and supplementary tables (env PlantCad2).
set -e
P=${GMO_ROOT}/paper
cd $P/figures
for s in make_fig1.py make_fig2.py make_fig3.py make_fig4.py make_fig5_genera.py; do FIG_DPI=600 MPLBACKEND=Agg python $s > /dev/null; done
python export_tiff.py
for s in make_fig1.py make_fig2.py make_fig3.py make_fig4.py make_fig5_genera.py; do MPLBACKEND=Agg python $s > /dev/null; done
cd $P/supplementary && python build_supplementary.py > build_supplementary.log
echo REBUILD_DONE
