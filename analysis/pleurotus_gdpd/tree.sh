: "${GMO_ROOT:=$(cd "$(dirname "$0")/../.." && pwd)}"; export GMO_ROOT
set -x
source $(conda info --base)/etc/profile.d/conda.sh
W=${GMO_ROOT}/analysis/pleurotus_gdpd; cd $W/tree
TF=${GMO_ROOT}/analysis/trees_final/ACS
cat $TF/ACS_cut.faa ../tree_add.faa > set.faa
conda activate hmmer_pfam
hmmalign --amino --trim -o al.sto $TF/ACS.profile.hmm set.faa
esl-alimask --rf-is-mask -o al.match.sto al.sto
esl-reformat -o al.aln afa al.match.sto
clipkit al.aln -m smart-gap -o al.trim0.aln -l > /dev/null
python3 - <<'PY'
recs=[]
for l in open('al.trim0.aln'):
    l=l.rstrip('\n')
    if l.startswith('>'): recs.append([l[1:],''])
    elif l: recs[-1][1]+=l
keep=[(h,s) for h,s in recs if s and (len(s)-s.count('-'))/len(s)>=0.15]
open('al.trim.aln','w').write(''.join(f'>{h}\n{s}\n' for h,s in keep))
print(len(keep),len(recs))
PY
conda activate phylo_env
iqtree -s al.trim.aln -m Q.PFAM+G4 -B 1000 -alrt 1000 -T 32 -pre pl_acs -redo -quiet > /dev/null 2>&1
iqtree -s al.trim.aln -m VT+I+R4 -B 1000 -alrt 1000 -T 32 -pre pl_acs_vt -redo -quiet > /dev/null 2>&1 || true
