#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../../../.." && pwd)}"; export GMO_ROOT
# Backbone-threshold sensitivity: collect contigs -> predict genes -> DIAMOND (blastp vs catalogue; blastx as in step-2 screen).
# Submit: bsub -K -q short -n 16 -R "span[hosts=1]" -R "rusage[mem=2000]" -J thr_bb -o work/%J.out -e work/%J.err "bash scripts/run.sh"
set -eo pipefail
source $(conda info --base)/etc/profile.d/conda.sh; conda activate annot
B=${GMO_ROOT}/analysis/thresholds/backbone; S=$B/scripts; W=$B/work
CAT=${GMO_ROOT}/analysis/refset/plasmids/plasmid_proteins_labelled.faa   # 7105 proteins, 42 plasmids
Q389=${GMO_ROOT}/analysis/refset/screen/step2/queries.dmnd             # 389 representatives used by the step-2 screen
cd $W
python -I $S/01_collect_contigs.py
PYTHONPATH=$W/pylib python $S/02_predict_genes.py
diamond makedb --in $CAT -d $W/plasmid_cat --quiet
# (1) requested: predicted proteins vs whole labelled catalogue
diamond blastp -q $W/proteins.faa -d $W/plasmid_cat -o $W/blastp_vs_catalog.tsv --threads 16 --very-sensitive \
  --evalue 1e-5 --max-target-seqs 25 --quiet --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore
# (2) exact replication of the step-2 screen (screen_one.sh): contig blastx vs 389 representatives
diamond blastx -q $W/contigs.fna -d $Q389 -o $W/blastx_vs_q389.tsv --threads 16 --very-sensitive -F 15 --range-culling --top 10 \
  --evalue 1e-5 --quiet --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qframe
python -I $S/03_summarize.py
# real-screen distributions (read-only on step2_loci.tsv); were submitted as separate bsub jobs
python -I $S/04_screen_distribution.py
python -I $S/05_screen_flag_sensitivity.py
bash $S/06_boletales_backbone.sh
echo RUN_DONE
