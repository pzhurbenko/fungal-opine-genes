# DIAMOND database of NCBI nr (the October 2025 release was used); build it with `diamond makedb --in nr.gz -d nr`
NR_DB=${NR_DB:-${GMO_ROOT:?set GMO_ROOT}/data/nr/nr.dmnd}
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../.." && pwd)}"; export GMO_ROOT
source $(conda info --base)/etc/profile.d/conda.sh; conda activate annot
cd ${GMO_ROOT}/analysis/chs_fusarium
diamond blastp -q query5.faa -d "$NR_DB" -o diamond_nr_top100.tsv --threads 64 -b 24 -c 1 --max-target-seqs 100 --evalue 1e-5 --outfmt 6 qseqid sseqid pident length qlen slen qcovhsp evalue bitscore stitle
