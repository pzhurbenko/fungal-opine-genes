#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/.." && pwd)}"; export GMO_ROOT
# Rebuild everything for the manuscript (env PlantCad2): figures (300 dpi PNG + 600 dpi TIFF), supplementary, manuscript md (marked + clean), HTML, DOCX.
set -e
P=${GMO_ROOT}/paper
cd $P/figures
for s in make_fig1.py make_fig2.py make_fig3.py make_fig4.py make_fig5_genera.py; do FIG_DPI=600 MPLBACKEND=Agg python $s > /dev/null; done
python export_tiff.py
for s in make_fig1.py make_fig2.py make_fig3.py make_fig4.py make_fig5_genera.py; do MPLBACKEND=Agg python $s > /dev/null; done
cd $P/supplementary && python build_supplementary.py > build_supplementary.log
# The manuscript text (draft_v2_*.md, references_v2.md) is not part of the public repository; the steps below run only if it is present.
if [ -f $P/draft_v2_methods.md ]; then
cd $P/manuscript_v2 && python build_manuscript_v2.py && python build_manuscript_v2.py --clean -o manuscript_v2_clean.md
cd $P && python build_paper_html.py || echo "build_paper_html exit $? (parser warnings are non-fatal)"
cd $P/manuscript_v2 && python build_docx.py -i manuscript_v2.md -o manuscript_v2.docx && python build_docx.py -i manuscript_v2_clean.md -o manuscript_v2_clean.docx
fi
echo REBUILD_DONE
