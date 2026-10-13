: "${GMO_ROOT:=$(cd "$(dirname "$0")/../.." && pwd)}"; export GMO_ROOT
source $(conda info --base)/etc/profile.d/conda.sh; conda activate annot
cd ${GMO_ROOT}/analysis/chs_fusarium
diamond blastp -q query5.faa -d /shareDB/blastdb/nr/diamond/nr_10_2025.dmnd -o diamond_nr_top100.tsv --threads 64 -b 24 -c 1 --max-target-seqs 100 --evalue 1e-5 --outfmt 6 qseqid sseqid pident length qlen slen qcovhsp evalue bitscore stitle
