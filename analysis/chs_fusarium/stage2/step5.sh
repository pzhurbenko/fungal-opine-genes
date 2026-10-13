#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../../.." && pwd)}"; export GMO_ROOT
# +-50 kb around chs locus, diamond blastx vs plasmid/T-DNA queries (params as step2/screen_one.sh)
source $(conda info --base)/etc/profile.d/conda.sh; conda activate annot
S=${GMO_ROOT}/analysis/chs_fusarium/stage2; D=${GMO_ROOT}/data/fungal_genomes/fusarium/fus/ncbi_dataset/data
mkdir -p $S/step5; cd $S/step5
while read gca ctg s e name; do
  F=$(ls $D/$gca/*.fna); a=$((s-50000)); [ $a -lt 1 ] && a=1; b=$((e+50000))
  seqkit faidx $F $ctg:$a-$b 2>/dev/null | seqkit replace -p '.*' -r "${name}_${ctg}_${a}_$b" > $name.fa || continue
  diamond blastx -q $name.fa -d ${GMO_ROOT}/analysis/refset/screen/step2/queries.dmnd -o $name.tsv --threads 8 --very-sensitive -F 15 --range-culling --top 10 \
    --evalue 1e-5 --quiet --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qframe
done <<EOL
GCA_000151355.1 GG698909.1 900562 901779 vanettenii
GCA_000149955.2 CM000600.1 483558 484454 Fol4287
GCA_054098485.1 JBJLOX010000025.1 22466 23671 nematophilum1
GCA_054098565.1 JBRIOT010000006.1 3980459 3981664 nematophilum2
GCA_025433565.1 CM046522.1 934968 936182 sp_Ph1
EOL
