#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../.." && pwd)}"; export GMO_ROOT
# susD check (2026-10-15): diamond blastp of all annotated fungal proteomes vs susD of pTiEU6 (ARU12249.1; the query file susD.faa also holds susL AAZ50401.1), same parameters as seeds_ext screen.
source $(conda info --base)/etc/profile.d/conda.sh; conda activate annot
D=${GMO_ROOT}/analysis/refset/screen; S=${GMO_ROOT}/analysis/susD_screen
LIST=${2:-$D/step3/proteomes.txt}; PER=${PER:-50}; THR=${THR:-8}; DBSIZE=${DBSIZE:-2233656}
mkdir -p $S/out
[ -s $S/susD.dmnd ] || diamond makedb --in $S/susD.faa -d $S/susD --quiet
if [ "${LSB_JOBINDEX:-0}" -gt 0 ]; then
  A=$(( ${1:-0} + (LSB_JOBINDEX-1)*PER + 1 )); B=$(( A + PER - 1 )); LINES=$(sed -n "${A},${B}p" $LIST)
elif [ -z "$LINES" ]; then LINES=$(cat $LIST); fi
for P in $LINES; do
  P=${P//\$\{GMO_ROOT\}/$GMO_ROOT}   # proteome list stores ${GMO_ROOT} literally
  GCA=$(basename $(dirname $P)); O=$S/out/$GCA
  [ -e $O.done ] && continue
  diamond blastp -q $P -d $S/susD.dmnd -o $O.tsv.tmp --threads $THR --very-sensitive -e 1e-3 --max-target-seqs 5 \
    --dbsize $DBSIZE --quiet \
    --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore \
    && awk -F'\t' '$12>=50' $O.tsv.tmp > $O.tsv && rm -f $O.tsv.tmp && echo "$GCA $(wc -l < $O.tsv)" > $O.done
done
