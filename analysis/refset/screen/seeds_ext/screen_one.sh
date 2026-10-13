#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../../../.." && pwd)}"; export GMO_ROOT
# Seeds_ext screen: diamond blastp of annotated fungal proteomes vs the extended seed set (seeds_ext.faa: 11 seeds + T1 refs + T-DNA cluster members).
# Usage: screen_one.sh BLOCK_OFFSET [LIST]  -> processes proteomes lines (OFF+(LSB_JOBINDEX-1)*PER+1 .. +PER); PER=50 by default.
# Without LSB_JOBINDEX (local test): processes LIST lines given via env LINES="a b c" or whole LIST.
source $(conda info --base)/etc/profile.d/conda.sh; conda activate annot
D=${GMO_ROOT}/analysis/refset/screen; S=$D/seeds_ext
LIST=${2:-$D/step3/proteomes.txt}; PER=${PER:-50}; THR=${THR:-8}; DBSIZE=${DBSIZE:-2233656}
[ -s $S/seeds_ext.dmnd ] || diamond makedb --in $S/seeds_ext.faa -d $S/seeds_ext --quiet
mkdir -p $S/out
if [ "${LSB_JOBINDEX:-0}" -gt 0 ]; then
  A=$(( ${1:-0} + (LSB_JOBINDEX-1)*PER + 1 )); B=$(( A + PER - 1 )); LINES=$(sed -n "${A},${B}p" $LIST)
elif [ -z "$LINES" ]; then LINES=$(cat $LIST); fi
for P in $LINES; do
  GCA=$(basename $(dirname $P)); O=$S/out/$GCA
  [ -e $O.done ] && continue
  diamond blastp -q $P -d $S/seeds_ext.dmnd -o $O.tsv.tmp --threads $THR --very-sensitive -e 1e-3 --max-target-seqs $(grep -c ">" $S/seeds_ext.faa | awk "{print \$1+4}") \
    --dbsize $DBSIZE --quiet \
    --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore \
    && awk -F'\t' '$12>=50' $O.tsv.tmp > $O.tsv && rm -f $O.tmp $O.tsv.tmp && echo "$GCA $(wc -l < $O.tsv)" > $O.done
done
