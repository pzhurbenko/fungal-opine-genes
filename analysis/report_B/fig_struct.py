import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
#!/usr/bin/env python3
import os, csv
os.environ['QT_QPA_PLATFORM'] = 'offscreen'
from ete3 import Tree, TreeStyle, NodeStyle, TextFace
F = f'{GMO_ROOT}/analysis/fold'
mp = {r[0]: r[1] for r in csv.reader(open(f'{F}/map.tsv'), delimiter='\t')}
SIX = ['JBBPKL010003317', 'JBBPKL010018021', 'JBBPNB010000203', 'JBBPNB010000771', 'JBBPOG010012756', 'JBBPRY010001778']
def grp(n):
    if n.startswith('REF_'): return 'REF_' + n.split('__')[0].split('_', 2)[2]
    if n.startswith('FUNGI_'): return 'FUNGI_' + n.split('__')[0].split('_')[1]
    return n.split('__')[0]
COL = {'FUNGI_gdpd': '#1a9641', 'FUNGI_acs': '#7b3fa0', 'REF_acs': '#6a1b9a', 'REF_gdpd_meso': '#1b7837', 'REF_kerr14': '#c51b8a', 'BACT_acc': '#8c564b', 'OG': '#000000'}
t = Tree(f'{F}/foldmason/struct_aa.treefile', format=1)
for l in t.iter_leaves(): l.add_feature('full', mp.get(l.name, l.name))
keep = [l for l in t.iter_leaves() if not ('sus' in grp(l.full) or 'ocd' in grp(l.full))]
t.prune(keep, preserve_branch_length=True)
og = [l for l in t.iter_leaves() if l.full.startswith(('NATIVE_', 'OG__'))]
anc = t.get_common_ancestor(og)
t.set_outgroup(og[0] if anc is t else anc)
ts = TreeStyle(); ts.show_leaf_name = False; ts.show_scale = True; ts.branch_vertical_margin = 1; ts.scale = 400
for l in t.iter_leaves():
    g = grp(l.full); col = COL.get(g, '#33691e' if g.startswith('NATIVE') else '#444444')
    nm = l.full.replace('FUNGI_', '').replace('REF_', 'REF ')[:70]
    if any(s in l.full for s in SIX): nm = 'X ' + nm
    l.add_face(TextFace(nm, fsize=9, fgcolor=col), column=0)
    ns = NodeStyle(); ns['size'] = 0; l.set_style(ns)
for n in t.traverse():
    if not n.is_leaf():
        ns = NodeStyle(); ns['size'] = 0; n.set_style(ns)
for label, sel, color in [('acc-locus lineage (gdpd_meso): 48 fungal + refs', lambda n: grp(n) in ('FUNGI_gdpd', 'REF_gdpd_meso') or 'WP_183902790' in n or 'WP_292444656' in n, '#e5ffe5'),
                          ('six + Parasponia PON40589', lambda n: any(s in n for s in SIX) or 'PON40589' in n, '#f0e5ff'),
                          ('T-DNA acs refs', lambda n: grp(n) == 'REF_acs' and 'PON40589' not in n, '#f7f0ff')]:
    lv = [l for l in t.iter_leaves() if sel(l.full)]
    m = t.get_common_ancestor(lv); ns = NodeStyle(); ns['bgcolor'] = color; ns['size'] = 0; m.set_style(ns)
    m.add_face(TextFace(' ' + label, fsize=11, bold=True), column=0, position='branch-bottom')
for ext in ('png', 'svg', 'pdf'):
    t.render(f'fig/FB1_ACS_structural_tree.{ext}', tree_style=ts, w=2000, units='px', dpi=150)
print('done', len(t))
