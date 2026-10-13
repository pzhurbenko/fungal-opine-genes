#!/usr/bin/env python3
# per-fold analysis of structure-based trees: prune to one fold, root on bacterial/native outgroup, test monophyly
import csv, sys
from ete3 import Tree
mp = {r[0]: r[1] for r in csv.reader(open('../map.tsv'), delimiter='\t')}
SIX = ['JBBPKL010003317', 'JBBPKL010018021', 'JBBPNB010000203', 'JBBPNB010000771', 'JBBPOG010012756', 'JBBPRY010001778']
def grp(n):
    if n.startswith('REF_'): return 'REF_' + n.split('__')[0].split('_', 2)[2]
    if n.startswith('FUNGI_'): return 'FUNGI_' + n.split('__')[0].split('_')[1]
    return n.split('__')[0]
FOLDS = {
 'SUS(OCD fold)': (lambda n: 'sus' in grp(n) or 'ocd' in grp(n), lambda n: grp(n) == 'REF_ocd_other_bacteria',
   [('fungal sus + T-DNA sus refs', lambda n: grp(n) in ('FUNGI_sus', 'REF_sus')), ('fungal sus only', lambda n: grp(n) == 'FUNGI_sus'), ('T-DNA sus refs only', lambda n: grp(n) == 'REF_sus')]),
 'ACS(GDPD fold)': (lambda n: not ('sus' in grp(n) or 'ocd' in grp(n)), lambda n: n.startswith('NATIVE_') or n.startswith('OG__'),
   [('gdpd_meso fungal', lambda n: grp(n) == 'FUNGI_gdpd'), ('gdpd_meso fungal + gdpd_meso refs + acc BACT', lambda n: grp(n) in ('FUNGI_gdpd', 'REF_gdpd_meso') or ('WP_183902790' in n or 'WP_292444656' in n)),
    ('six lineage', lambda n: any(s in n for s in SIX)), ('six + Parasponia PON40589', lambda n: any(s in n for s in SIX) or 'PON40589' in n), ('six + kerr14 + PON40589', lambda n: any(s in n for s in SIX) or 'PON40589' in n or grp(n) == 'REF_kerr14'),
    ('T-DNA acs refs (excl PON40589)', lambda n: grp(n) == 'REF_acs' and 'PON40589' not in n), ('T-DNA acs refs + six', lambda n: (grp(n) == 'REF_acs' and 'PON40589' not in n) or any(s in n for s in SIX))])}
for label, tf in [('struct_AA', 'struct_aa.treefile'), ('struct_3Di', 'struct_3di.treefile')]:
    try: T = Tree(tf, format=1)
    except Exception: print(label, 'not ready'); continue
    for l in T.iter_leaves(): l.add_feature('full', mp.get(l.name, l.name))
    for fold, (keep, rootsel, tests) in FOLDS.items():
        t = T.copy(); t.prune([l for l in t.iter_leaves() if keep(l.full)], preserve_branch_length=True)
        og = [l for l in t.iter_leaves() if rootsel(l.full)]
        if og:
            anc = t.get_common_ancestor(og) if len(og) > 1 else og[0]
            try: t.set_outgroup(anc)
            except Exception: t.set_outgroup(og[0])
        print(f'== {label} {fold}: leaves {len(t)}, outgroup {len(og)}')
        for name, sel in tests:
            leaves = [l for l in t.iter_leaves() if sel(l.full)]
            if len(leaves) < 2: continue
            m = t.get_common_ancestor(leaves); allm = m.get_leaves(); extra = [l.full[:55] for l in allm if not sel(l.full)]
            print(f'   {name}: n={len(leaves)} clade={len(allm)} mono={len(extra)==0} supp={m.name or "-"} intruders={len(extra)} {extra[:2]}')
