#!/bin/bash
# Pack the Supplementary Materials archive for submission (same file list as v15, 35 files).
# Run after build_supplementary.py (rebuild_all.sh). Usage: bash make_zip.sh [output.zip]
set -e
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../.." && pwd)}"; export GMO_ROOT
S=$GMO_ROOT/paper/supplementary
OUT=${1:-$GMO_ROOT/paper/submission/Supplementary_Materials_fungal_opine_genes.zip}
mkdir -p "$(dirname "$OUT")"
T=$(mktemp -d)
mkdir -p $T/figures $T/figures_pdf
cp $S/Supplementary_Materials.docx $S/Supplementary_Tables.xlsx $S/README_supplementary.md $S/Extended_Methods.md \
   $S/Data_S1_proteins_of_loci_of_bacterial_origin.faa $S/Data_S3_screened_assemblies.tsv $T/
cp $S/data/Supplementary_Data_trees.tar.gz $T/Data_S2_trees_alignments.tar.gz
cp $S/data/Data_S2b_AU_constraints.tar.gz $S/data/Data_S4_dNdS_codeml.tar.gz $T/
cp $S/figures/FigureS*.png $T/figures/
cp $S/figures/FigureS*.pdf $T/figures_pdf/
rm -f "$OUT"
(cd $T && zip -q -r -D "$OUT" Supplementary_Materials.docx Supplementary_Tables.xlsx README_supplementary.md Extended_Methods.md \
   Data_S1_proteins_of_loci_of_bacterial_origin.faa Data_S3_screened_assemblies.tsv Data_S2_trees_alignments.tar.gz \
   Data_S2b_AU_constraints.tar.gz Data_S4_dNdS_codeml.tar.gz figures figures_pdf)
rm -rf $T
unzip -l "$OUT" | tail -1
