#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../../.." && pwd)}"; export GMO_ROOT
set -eo pipefail
source $(conda info --base)/etc/profile.d/conda.sh; conda activate annot
cd ${GMO_ROOT}/analysis/fusarium_odh_island/rhypophila_flank
CAT=${GMO_ROOT}/analysis/refset/plasmids/plasmid_proteins_labelled.faa   # 7105 proteins, 42 plasmids
REP=${GMO_ROOT}/analysis/refset/step2/queries.faa                        # 389 representatives (T-DNA/opine/backbone/T1-T2)
diamond makedb --in $CAT -d plasmid_cat --quiet
diamond makedb --in $REP -d queries389 --quiet
F="6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore stitle"
diamond blastp -q flank50.faa -d plasmid_cat -o flank50_vs_plasmid.tsv --threads 8 --very-sensitive --evalue 1e-5 --max-target-seqs 5 --outfmt $F
diamond blastp -q flank50.faa -d queries389 -o flank50_vs_queries389.tsv --threads 8 --very-sensitive --evalue 1e-5 --max-target-seqs 5 --outfmt $F
diamond blastx -q flank50_windows.fna -d plasmid_cat -o windows_blastx_vs_plasmid.tsv --threads 8 --very-sensitive -F 15 --evalue 1e-5 --max-target-seqs 25 --outfmt $F
echo DMND_DONE
