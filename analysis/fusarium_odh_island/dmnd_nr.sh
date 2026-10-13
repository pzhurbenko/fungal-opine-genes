: "${GMO_ROOT:=$(cd "$(dirname "$0")/../.." && pwd)}"; export GMO_ROOT
source $(conda info --base)/etc/profile.d/conda.sh; conda activate annot
cd ${GMO_ROOT}/analysis/fusarium_odh_island
date
diamond blastp -q island_queries.faa -d /shareDB/blastdb/nr/diamond/nr_10_2025.dmnd -o island_vs_nr.tsv --threads 64 -b 24 -c 1 --sensitive --max-target-seqs 1000 --evalue 1e-5 --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore stitle
date
