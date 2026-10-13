#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../../../.." && pwd)}"; export GMO_ROOT
# Per-genome max backbone identity in Boletales assemblies of the step-2 screen (order from step2_loci.tsv col 5),
# and number of Boletales genomes with >=1 backbone-best locus at each threshold. Read-only on step2_loci.tsv.
L=${GMO_ROOT}/analysis/refset/screen/step2/step2_loci.tsv
O=${GMO_ROOT}/analysis/thresholds/backbone
awk -F'\t' 'NR>1 && $5=="Boletales"{g[$1]=$2; if($14 ~ /^(vir|tra|rep_partition|mobile)$/ && $16+0>m[$1]){m[$1]=$16+0; c[$1]=$6; l[$1]=$8; s[$1]=$14":"$15}}
  END{for(k in g) printf "%s\t%s\t%.1f\t%s\t%s\t%s\n",k,g[k],m[k]+0,c[k],l[k],s[k]}' $L | sort -t$'\t' -k3,3gr > $O/screen_boletales_genome_maxbb.tsv
for t in 30 35 40 45 50 55 60 65 70 75 80 85 90; do awk -F'\t' -v t=$t '$3>=t{n++} END{printf "%d\t%d\t%d\n",t,n+0,NR}' $O/screen_boletales_genome_maxbb.tsv; done > $O/screen_boletales_by_threshold.tsv
