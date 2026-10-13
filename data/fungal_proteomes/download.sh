#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../.." && pwd)}"; export GMO_ROOT
source $(conda info --base)/etc/profile.d/conda.sh; conda activate annot
cd ${GMO_ROOT}/data/fungal_proteomes
echo "START $(date)"
datasets download genome taxon fungi --annotated --include protein --dehydrated --filename fungi_annotated_protein.zip && \
unzip -o -q fungi_annotated_protein.zip -d fungi_annotated && \
datasets rehydrate --directory fungi_annotated --max-workers 10 && echo "REHYDRATED $(date)" && \
ls fungi_annotated/ncbi_dataset/data/*/protein.faa | wc -l && echo "DONE $(date)"
