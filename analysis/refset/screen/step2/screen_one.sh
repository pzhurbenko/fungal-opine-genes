#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../../../.." && pwd)}"; export GMO_ROOT
# Step 2 screening of ONE whole fungal assembly (unannotated genomes incl. Boletales):
#  (a) diamond blastx (frameshift-aware, -F 15) of the genome vs 389 representative plasmid/T-DNA/opine/family proteins (queries.faa);
#  (b) six-frame ORFs >= 60 aa (esl-translate) -> hmmsearch with 63 competitive HMMs (family subgroups+OUTGROUP, 51 T-DNA lineages).
# Usage: screen_one.sh <OFFSET> <list.tsv: GCA<TAB>fna<TAB>size>   (index = OFFSET + LSB_JOBINDEX)
source $(conda info --base)/etc/profile.d/conda.sh; conda activate annot
D=${GMO_ROOT}/analysis/refset/screen
IDX=$(( ${1:-0} + LSB_JOBINDEX )); LINE=$(sed -n "${IDX}p" ${2:-$D/step2/genomes.tsv})
GCA=$(echo "$LINE" | cut -f1); FNA=$(echo "$LINE" | cut -f2)
O=$D/step2/out/$GCA; mkdir -p $O
[ -s $O/done ] && exit 0
T0=$(date +%s)
diamond blastx -q $FNA -d $D/step2/queries.dmnd -o $O/dmnd.tsv --threads 2 --very-sensitive -F 15 --range-culling --top 10 \
  --evalue 1e-5 --quiet --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qframe
T1=$(date +%s)
esl-translate -l 60 $FNA | sed 's/ desc=.*//' > $O/orfs.faa
hmmsearch --cpu 2 --noali -E 1e-3 --tblout $O/hmm.tbl --domtblout $O/hmm.dom $D/step3/all_screen.hmm $O/orfs.faa > /dev/null
T2=$(date +%s)
# keep only ORFs with an HMM hit (for later extraction), drop the full six-frame file
grep -v '^#' $O/hmm.tbl | awk '{print $1}' | sort -u > $O/hit_orfs.txt
seqkit grep -f $O/hit_orfs.txt $O/orfs.faa > $O/hit_orfs.faa 2>/dev/null; rm -f $O/orfs.faa
echo "$GCA dmnd=$(wc -l < $O/dmnd.tsv) hmm=$(grep -vc '^#' $O/hmm.tbl) t_dmnd=$((T1-T0))s t_hmm=$((T2-T1))s" > $O/done
