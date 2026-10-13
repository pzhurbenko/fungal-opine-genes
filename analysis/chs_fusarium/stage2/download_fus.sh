#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../../.." && pwd)}"; export GMO_ROOT
source $(conda info --base)/etc/profile.d/conda.sh; conda activate annot
cd ${GMO_ROOT}/data/fungal_genomes/fusarium
echo "START $(date)"
datasets download genome accession --inputfile gca_list.txt --include genome,gff3 --dehydrated --filename fus.zip && \
unzip -o -q fus.zip -d fus && datasets rehydrate --directory fus --max-workers 10 && \
echo "genomes: $(ls fus/ncbi_dataset/data/*/*.fna | wc -l)" && echo "DONE $(date)" || echo "FAILED $(date)"
