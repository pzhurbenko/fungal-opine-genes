: "${GMO_ROOT:=$(cd "$(dirname "$0")/../.." && pwd)}"; export GMO_ROOT
source $(conda info --base)/etc/profile.d/conda.sh; conda activate annot
W=${GMO_ROOT}/analysis/pleurotus_gdpd; cd $W
mkdir -p scan
for f in genomes/ncbi_dataset/data/GCA_*/*.fna genomes/ncbi_dataset/data/GCF_*/*.fna; do
  [ -s "$f" ] || continue
  a=$(basename $(dirname $f))
  miniprot -t 32 --gff -p 0.2 --outs=0.3 -N 5 "$f" q1.faa > scan/$a.gff 2> /dev/null
done
