# DIAMOND database of NCBI nr (the October 2025 release was used); build it with `diamond makedb --in nr.gz -d nr`
NR_DB=${NR_DB:-${GMO_ROOT:?set GMO_ROOT}/data/nr/nr.dmnd}
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../.." && pwd)}"; export GMO_ROOT
source $(conda info --base)/etc/profile.d/conda.sh; conda activate annot
cd ${GMO_ROOT}/analysis/pleurotus_gdpd
diamond blastp -q nr_queries.faa -d "$NR_DB" -o nr_hits.tsv --threads 64 -b 24 -c 1 --sensitive --max-target-seqs 500 --evalue 1e-5 --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore stitle
