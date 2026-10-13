#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../../../.." && pwd)}"; export GMO_ROOT
# Step 3 screening of ONE annotated fungal proteome: hmmsearch (family subgroup+OUTGROUP HMMs + 51 T-DNA lineage HMMs) and
# diamond blastp against all 7105 labelled plasmid proteins (T-DNA / opine_related / vir / tra / rep / mobile / other).
source $(conda info --base)/etc/profile.d/conda.sh; conda activate annot
D=${GMO_ROOT}/analysis/refset/screen
IDX=$(( ${1:-0} + LSB_JOBINDEX )); P=$(sed -n "${IDX}p" ${2:-$D/step3/proteomes.txt}); GCA=$(basename $(dirname $P))
O=$D/step3/out/$GCA; mkdir -p $O
[ -s $O/done ] && exit 0
hmmsearch --cpu 2 --noali -E 1e-3 --tblout $O/hmm.tbl $D/step3/all_screen.hmm $P > /dev/null
diamond blastp -q $P -d $D/plasmid_proteins.dmnd -o $O/dmnd.tsv --threads 2 --evalue 1e-5 --max-target-seqs 5 --quiet \
  --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore
echo "$GCA $(grep -vc '^#' $O/hmm.tbl) $(wc -l < $O/dmnd.tsv)" > $O/done
