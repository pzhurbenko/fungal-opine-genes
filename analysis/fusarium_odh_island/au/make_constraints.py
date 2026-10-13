#!/usr/bin/env python3
"""Build monophyly constraint trees (IQ-TREE -g) for the ODH island AU test.
Leaves are taken from ../tree/odh_tree.trim.aln; groups by name prefix:
  ISL_   - Fusarium poae / F. langsethiae island copies
  NAT_   - native fungal ODH (native_selected.tsv + Sporormia/Astrocystis from bact_selected.tsv)
  SEED_  - T-DNA seeds (nos x2, ocs, vis)
  FUNG_  - Rhypophila / Butyriboletus (nos-like fungi)
  BACT_Agrobacterium_larrymoorei_WP_309770759 - closest bacterial ODH (LARRY)
  BACT2_ - second A. larrymoorei ODH (WP_309770762), sits in the T-DNA part of ML tree
  BACT_  - other bacterial ODH
Each constraint: ((group),all other leaves); -> one bipartition, rest unconstrained.
"""
import os
D = os.path.dirname(os.path.abspath(__file__))
aln = os.path.join(D, '..', 'tree', 'odh_tree.trim.aln')
leaves = [l[1:].strip() for l in open(aln) if l.startswith('>')]
assert len(leaves) == len(set(leaves)) == 44, len(leaves)

LARRY = 'BACT_Agrobacterium_larrymoorei_WP_309770759'
P = lambda p: [x for x in leaves if x.startswith(p)]
ISL, NAT, SEED, FUNG, BACT2 = P('ISL_'), P('NAT_'), P('SEED_'), P('FUNG_'), P('BACT2_')
BACT = [x for x in P('BACT_') if x != LARRY]
NOS = P('SEED_nos')

H = {
    # name: (group, description)
    'H1_isl_native':  (ISL + NAT, 'ISL monophyletic with all native fungal ODH (fungal gene)'),
    'H1b_isl_sordario': (ISL + [x for x in NAT if 'Apiospora' in x or 'Astrocystis' in x], 'ISL + native Sordariomycetes ODH (Apiospora, Astrocystis): vertical-inheritance expectation for Fusarium'),
    'H1c_native_only': (NAT, 'control: native fungal ODH monophyletic (without ISL)'),
    'H2_isl_tdna':    (ISL + SEED + FUNG + BACT2, 'ISL inside ML T-DNA clade (seeds ocs/nos/vis + Rhypophila/Butyriboletus + larrymoorei WP_309770762)'),
    'H2s_isl_seeds':  (ISL + SEED, 'strict: ISL + T-DNA seeds ocs/nos/vis only'),
    'H3_isl_larry':   (ISL + [LARRY], 'control: ISL + A. larrymoorei WP_309770759 (ML clade)'),
    'H4_isl_bact':    (ISL + BACT, 'ISL + all other bacterial ODH, excluding both A. larrymoorei ODH'),
    'H5_isl_nos':     (ISL + NOS, 'ISL inside nos subfamily (nos seeds P00386 + Gillisella; no rhizobial NOSLIKE leaves in set)'),
}
with open(os.path.join(D, 'hypotheses.tsv'), 'w') as fo:
    fo.write('hypothesis\tn_in_group\tdescription\tgroup\n')
    for h, (g, desc) in H.items():
        assert len(set(g)) == len(g) and set(g) <= set(leaves)
        rest = [x for x in leaves if x not in g]
        open(os.path.join(D, h + '.nwk'), 'w').write('((' + ','.join(g) + '),' + ','.join(rest) + ');\n')
        fo.write(f'{h}\t{len(g)}\t{desc}\t{",".join(g)}\n')
        print(h, len(g))
