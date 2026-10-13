#!/usr/bin/env python3
"""Per-lineage AU constraints for ACS (replaces the uninformative global make_constraints.py for this family).
Groups (by label prefix in ACS.mafft.trim.aln):
  SIX   = the 6 intact Boletales ACS of the 'kerr14-like' lineage
  TDNA  = REF_*_acs__  (27 T-DNA acs refs, incl. Parasponia PON40589)
  KERR  = REF_*_kerr14__ + Parasponia PON40589
  ACC   = REF_*_gdpd_meso__ (3) + BACT_acc WP_183902790 (R. metallidurans) + WP_292444656 (Mesorhizobium) + FUNGI_gdpd_meso_* (60)
  NATV  = NATIVE_* + FUNGI_acs_* that are not SIX (native Agaricomycetes GDPD clade)
  MESO  = FUNGI_gdpd_meso_*  (59 Boletales, acc-locus lineage)
Hypotheses (each = one constrained ML tree):
  six_tdna: (SIX+TDNA) ; six_kerr: (SIX+KERR) ; six_acc: (SIX+ACC) ; six_natv: (SIX+NATV)
  meso_acc: (MESO+BACT_acc+REF gdpd_meso) ; meso_tdna: (MESO+TDNA)
"""
import sys, os
ALN=sys.argv[1]; OUT=os.path.dirname(os.path.abspath(ALN))
SIX=['JBBPKL010003317','JBBPKL010018021','JBBPNB010000203','JBBPNB010000771','JBBPOG010012756','JBBPRY010001778']
names=[l[1:].split()[0] for l in open(ALN) if l.startswith('>')]
six=[n for n in names if n.startswith('FUNGI_') and any(s in n for s in SIX)]
tdna=[n for n in names if n.startswith('REF_') and '_acs__' in n]
paras=[n for n in names if 'PON40589' in n]
kerr=[n for n in names if n.startswith('REF_') and '_kerr14__' in n]+paras
meso=[n for n in names if n.startswith('FUNGI_gdpd_meso')]
# acc-locus lineage = the actual sister group of the fungal gdpd_meso clade in the ML tree (BACT_acc as a whole is heterogeneous:
# torf6/acsx/'agrocinopine synthase' annotations sit with T-DNA acs, accF elsewhere)
bact=[n for n in names if (n.startswith('REF_') and '_gdpd_meso__' in n) or 'WP_183902790' in n or 'WP_292444656' in n]
natv=[n for n in names if n.startswith('NATIVE_') or (n.startswith('FUNGI_acs') and n not in six)]
H={'six_tdna':six+tdna,'six_kerr':six+kerr,'six_acc':six+bact+meso,'six_natv':six+natv,
   'meso_acc':meso+bact,'meso_tdna':meso+tdna}
for h,grp in H.items():
    rest=[n for n in names if n not in set(grp)]
    open(f'{OUT}/ACS.{h}.constraint.nwk','w').write('((' + ','.join(grp) + '),' + ','.join(rest) + ');\n')
    print(h,'group',len(grp),'rest',len(rest),file=sys.stderr)
print('six',len(six),'tdna',len(tdna),'kerr',len(kerr),'meso',len(meso),'bact',len(bact),'natv',len(natv),file=sys.stderr)
