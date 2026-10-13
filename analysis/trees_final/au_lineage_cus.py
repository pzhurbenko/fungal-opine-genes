#!/usr/bin/env python3
"""Per-lineage AU constraints for CUS (global make_constraints.py was uninformative: cus+mis merged).
Groups (labels in CUS.mafft.trim.aln, subgroup/organism from CUS_labels.tsv):
  CUSF = fungal cus (43); BOL = Boletales cus fungi (order Boletales via labels/loci table); REFCUS = REF_*_cus__ (21 T-DNA cus)
  STRAY4 = the 4 non-T-DNA bacterial refs that sit inside the cus clade in the ML tree
  ACT = fungal cus_actinobacteria (130) + REF_T2x_cus_actinobacteria (9); MIS = fungal mis (190) + REF_*_mis__ (8)
  BACT = all REF_T2x_cus_* (actino+beta+other, 26)
Hypotheses: cus_tdna (CUSF+REFCUS), cus_tdna4 (CUSF+REFCUS+STRAY4), cus_act (CUSF+ACT), cus_mis (CUSF+MIS), cus_bact (CUSF+BACT),
            cus_mono (CUSF alone monophyletic = one acquisition), bol_tdna (BOL+REFCUS), bol_act (BOL+ACT), bol_mis (BOL+MIS), bol_othercus (BOL + non-Boletales CUSF)
"""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import sys, os, csv
from Bio import Phylo
ALN=sys.argv[1]; OUT=os.path.dirname(os.path.abspath(ALN))
names=[l[1:].split()[0] for l in open(ALN) if l.startswith('>')]
lab={r['label']:r for r in csv.DictReader(open(f'{OUT}/CUS_labels.tsv'),delimiter='\t')}
loci={r['locus_id']:r for r in csv.DictReader(open(f'{GMO_ROOT}/analysis/loci_pipeline/round2/loci_table.tsv'),delimiter='\t')}
EXCL=['JBBPKX010002805']  # native fungal gene misassigned to cus (journal 2026-09-06); excluded from all groups, left unconstrained
cusf=[n for n in names if n.startswith('FUNGI_cus_') and 'actinobacteria' not in n and not any(e in n for e in EXCL)]
bol=[n for n in cusf if loci.get(lab[n]['locus'],{}).get('order')=='Boletales']
refcus=[n for n in names if n.startswith('REF_') and '_cus__' in n]
act=[n for n in names if 'cus_actinobacteria' in n]
mis=[n for n in names if n.startswith('FUNGI_mis_') or (n.startswith('REF_') and '_mis__' in n)]
bact=[n for n in names if n.startswith('REF_T2x_cus_')]
t=Phylo.read(f'{OUT}/CUS.mafft.iqtree.treefile','newick'); og=[l for l in t.get_terminals() if l.name.startswith('OG__')]; t.root_with_outgroup(og[0])
m=t.common_ancestor([l for l in t.get_terminals() if l.name in refcus])
stray4=[l.name for l in m.get_terminals() if l.name not in set(cusf)|set(refcus)]
print('STRAY inside cus clade:',stray4,file=sys.stderr)
print('cus fungi outside clade:',[n for n in cusf if n not in {l.name for l in m.get_terminals()}],file=sys.stderr)
H={'cus_tdna':cusf+refcus,'cus_tdna4':cusf+refcus+stray4,'cus_act':cusf+act,'cus_mis':cusf+mis,'cus_bact':cusf+bact,'cus_mono':cusf,
   'bol_tdna':bol+refcus,'bol_act':bol+act,'bol_mis':bol+mis,'bol_othercus':bol+[n for n in cusf if n not in bol]}
for h,grp in H.items():
    rest=[n for n in names if n not in set(grp)]
    open(f'{OUT}/CUS.{h}_x.constraint.nwk','w').write('((' + ','.join(grp) + '),' + ','.join(rest) + ');\n')
    print(h,'group',len(grp),'rest',len(rest),file=sys.stderr)
print('cusf',len(cusf),'bol',len(bol),'refcus',len(refcus),'act',len(act),'mis',len(mis),'bact',len(bact),file=sys.stderr)
