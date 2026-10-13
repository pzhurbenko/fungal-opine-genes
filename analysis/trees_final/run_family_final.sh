#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../.." && pwd)}"; export GMO_ROOT
# Quick trees for one family: two alignments (MAFFT L-INS-i; hmmalign to a hypothesis-neutral Pfam profile, match columns only),
# ClipKIT smart-gap trimming, FastTree (LG+G) and IQ-TREE (LG+F+G4, UFBoot 1000, 32 threads) on each. Timing recorded in timing.txt.
# usage: run_family.sh <FAM>     (ACS|CUS|SUS|ORF358)
set -eo pipefail
FAM=$1; T=32
source $(conda info --base)/etc/profile.d/conda.sh
D=${GMO_ROOT}/analysis/trees_final/$FAM; cd $D
PFAM=${GMO_ROOT}/analysis/pfam/Pfam-A.hmm
declare -A PROF=( [ACS]=GDPD [CUS]=Cucumopine_C [SUS]=OCD_Mu_crystall [ORF358]=NONE )
tm() { echo "$FAM $1 $(( $(date +%s) - T0 ))s" >> timing.txt; }
# ---- 0. cut every sequence to the family-fold envelope (best domain hit over all round-2 subgroup + OUTGROUP HMMs, +-25 aa);
#         removes fused neighbours of chimeric models and unalignable tails; sequences without a hit or envelope < 100 aa are dropped
conda activate annot
HD=${GMO_ROOT}/analysis/loci_pipeline/round2/hmm
if [ "${NOCUT:-0}" = "1" ]; then cp ${FAM}_set.faa ${FAM}_cut.faa; echo "$FAM: NOCUT=1, full sequences used" >&2; else
cat $HD/${FAM}__*.hmm > ${FAM}.cut_profiles.hmm
hmmsearch --cpu $T --noali -E 1e-3 --domtblout ${FAM}.cut.domtbl ${FAM}.cut_profiles.hmm ${FAM}_set.faa > /dev/null
python3 - ${FAM} <<'PY'
import sys
F=sys.argv[1]; best={}
for l in open(f'{F}.cut.domtbl'):
    if l.startswith('#'): continue
    p=l.split(); q=p[0]; sc=float(p[13]); a,b=int(p[19]),int(p[20])
    if q not in best or sc>best[q][0]: best[q]=(sc,a,b)
recs=[]; h=None
for l in open(f'{F}_set.faa'):
    l=l.rstrip('\n')
    if l.startswith('>'): h=l[1:]; recs.append([h,''])
    elif l: recs[-1][1]+=l
out=[]; drop=[]
for h,s in recs:
    if h not in best: drop.append((h,'no_hit')); continue
    sc,a,b=best[h]; a=max(1,a-25); b=min(len(s),b+25)
    if b-a+1<100: drop.append((h,f'env{b-a+1}')); continue
    out.append((h,s[a-1:b]))
open(f'{F}_cut.faa','w').write(''.join(f'>{h}\n{s}\n' for h,s in out))
open(f'{F}_cut_dropped.tsv','w').write(''.join(f'{h}\t{r}\n' for h,r in drop))
L=sorted(len(s) for h,s in out); print(f'{F}: cut {len(out)} kept, {len(drop)} dropped; length median {L[len(L)//2]} min {L[0]} max {L[-1]}', file=sys.stderr)
PY
fi
# ---- 1. MAFFT L-INS-i
T0=$(date +%s); mafft --localpair --maxiterate 1000 --thread $T --quiet ${FAM}_cut.faa > ${FAM}.mafft.aln; tm mafft_linsi
# ---- 2. hmmalign to a profile
conda activate hmmer_pfam
T0=$(date +%s)
if [ "${PROF[$FAM]}" != "NONE" ]; then hmmfetch $PFAM ${PROF[$FAM]} > ${FAM}.profile.hmm
else hmmbuild --amino -n ${FAM}_set ${FAM}.profile.hmm ${FAM}.mafft.aln > /dev/null; fi   # ORF358: no Pfam family -> profile from the set itself
hmmalign --amino --trim -o ${FAM}.hmmalign.sto ${FAM}.profile.hmm ${FAM}_cut.faa
esl-alimask --rf-is-mask -o ${FAM}.hmmalign.match.sto ${FAM}.hmmalign.sto
esl-reformat -o ${FAM}.hmmalign.aln afa ${FAM}.hmmalign.match.sto
tm hmmalign
# ---- 3. trim + trees
conda activate phylo_env
for M in mafft hmmalign; do
  T0=$(date +%s); clipkit ${FAM}.$M.aln -m smart-gap -o ${FAM}.$M.trim0.aln -l > /dev/null
  # drop sequences that kept < 15% non-gap columns (fragments outside the profile -> IQ-TREE refuses all-gap sequences)
  python3 - ${FAM}.$M.trim0.aln ${FAM}.$M.trim.aln <<'PY'
import sys
recs=[]; h=None
for l in open(sys.argv[1]):
    l=l.rstrip('\n')
    if l.startswith('>'): h=l[1:]; recs.append([h,''])
    elif l: recs[-1][1]+=l
keep=[(h,s) for h,s in recs if s and (len(s)-s.count('-'))/len(s)>=0.15]
open(sys.argv[2],'w').write(''.join(f'>{h}\n{s}\n' for h,s in keep))
print(f'{sys.argv[2]}: kept {len(keep)}/{len(recs)} sequences', file=sys.stderr)
PY
  tm clipkit_$M
  T0=$(date +%s); FastTree -lg -gamma -quiet ${FAM}.$M.trim.aln > ${FAM}.$M.fasttree.nwk 2> ${FAM}.$M.fasttree.log; tm fasttree_$M
  T0=$(date +%s); iqtree -s ${FAM}.$M.trim.aln -m MFP -B 1000 -alrt 1000 -T $T -pre ${FAM}.$M.iqtree -redo -quiet > /dev/null 2>&1 || true; tm iqtree_$M
done
echo "$FAM DONE $(date)" >> timing.txt
# ---- 5. AU test of constrained topologies (MAFFT alignment): H_tdna = target fungi + T-DNA refs monophyletic; H_alt = target fungi + alternative donor/native monophyletic
conda activate annot
python3 ${GMO_ROOT}/analysis/trees_final/make_constraints.py $FAM ${FAM}.mafft.trim.aln
conda activate phylo_env
MODEL=$(grep -oP 'Best-fit model: \K\S+' ${FAM}.mafft.iqtree.iqtree | head -1); MODEL=${MODEL:-LG+F+G4}
for H in H_tdna H_alt; do
  [ -s ${FAM}.$H.constraint.nwk ] || continue
  T0=$(date +%s); iqtree -s ${FAM}.mafft.trim.aln -m $MODEL -g ${FAM}.$H.constraint.nwk -T $T -pre ${FAM}.$H -redo -quiet > /dev/null 2>&1 || true; tm iqtree_$H
done
cat ${FAM}.mafft.iqtree.treefile ${FAM}.H_tdna.treefile ${FAM}.H_alt.treefile 2>/dev/null > ${FAM}.candidates.nwk
T0=$(date +%s); iqtree -s ${FAM}.mafft.trim.aln -m $MODEL -z ${FAM}.candidates.nwk -n 0 -zb 10000 -au -T $T -pre ${FAM}.AU -redo -quiet > /dev/null 2>&1 || true; tm au_test
echo "$FAM FINAL_DONE $(date)" >> timing.txt
