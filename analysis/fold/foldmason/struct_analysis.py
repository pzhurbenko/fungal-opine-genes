#!/usr/bin/env python3
# monophyly / clade composition of key groups in structure-based trees (AA of structural MSA; 3Di) vs sequence tree
import csv, sys
from ete3 import Tree
mp = {r[0]: r[1] for r in csv.reader(open('../map.tsv'), delimiter='\t')}
SIX = ['JBBPKL010003317', 'JBBPKL010018021', 'JBBPNB010000203', 'JBBPNB010000771', 'JBBPOG010012756', 'JBBPRY010001778']
def grp(n):
    if n.startswith('REF_'):
        return 'REF_' + n.split('__')[0].split('_', 2)[2]
    if n.startswith('FUNGI_'):
        return 'FUNGI_' + n.split('__')[0].split('_')[1]
    return n.split('__')[0]
tests = {
 'SUS: fungal sus + T-DNA sus refs': lambda n: grp(n) in ('FUNGI_sus', 'REF_sus'),
 'SUS: fungal sus only': lambda n: grp(n) == 'FUNGI_sus',
 'ACS: gdpd_meso fungal': lambda n: grp(n) == 'FUNGI_gdpd',
 'ACS: gdpd_meso + gdpd_meso refs': lambda n: grp(n) in ('FUNGI_gdpd', 'REF_gdpd_meso'),
 'ACS: six lineage': lambda n: any(s in n for s in SIX),
 'ACS: six + kerr14 ref': lambda n: any(s in n for s in SIX) or grp(n) == 'REF_kerr14',
 'ACS: T-DNA acs refs': lambda n: grp(n) == 'REF_acs',
}
for label, tf in [('struct_AA', 'struct_aa.treefile'), ('struct_3Di', 'struct_3di.treefile')]:
    try: t = Tree(tf, format=1)
    except Exception as e: print(label, 'missing', e); continue
    for l in t.iter_leaves(): l.add_feature('full', mp.get(l.name, l.name))
    og = [l for l in t.iter_leaves() if l.full.startswith('OG__')]
    if og: t.set_outgroup(og[0])
    print('==', label, 'leaves', len(t))
    for name, sel in tests.items():
        leaves = [l for l in t.iter_leaves() if sel(l.full)]
        if len(leaves) < 2: continue
        m = t.get_common_ancestor(leaves); allm = m.get_leaves()
        extra = [l.full[:60] for l in allm if not sel(l.full)]
        print(f'  {name}: n={len(leaves)} MRCA size={len(allm)} monophyletic={len(extra)==0} support={m.name} intruders={len(extra)} {extra[:3]}')
