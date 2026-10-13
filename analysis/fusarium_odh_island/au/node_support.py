#!/usr/bin/env python3
"""Pure-python: supports (SH-aLRT/UFBoot) of key bipartitions in an IQ-TREE newick tree.
usage: node_support.py tree.nwk  (default: ../tree/odh_tree.iqtree.treefile)"""
import os, re, sys
D = os.path.dirname(os.path.abspath(__file__))
f = sys.argv[1] if len(sys.argv) > 1 else os.path.join(D, '..', 'tree', 'odh_tree.iqtree.treefile')
s = open(f).read().strip().rstrip(';')

def parse(s):
    # returns list of (leafset, label) for internal nodes
    stack, nodes, i = [[]], [], 0
    while i < len(s):
        c = s[i]
        if c == '(':
            stack.append([]); i += 1
        elif c == ')':
            kids = stack.pop(); i += 1
            m = re.match(r"[^,():;]*", s[i:]); lab = m.group(0); i += len(lab)
            m2 = re.match(r":[^,();]*", s[i:]); i += len(m2.group(0)) if m2 else 0
            ls = frozenset().union(*kids); nodes.append((ls, lab)); stack[-1].append(ls)
        elif c == ',':
            i += 1
        else:
            m = re.match(r"([^,():;]+)(:[^,();]*)?", s[i:]); i += len(m.group(0))
            stack[-1].append(frozenset([m.group(1)]))
    return nodes

nodes = parse(s)
L = max(nodes, key=lambda x: len(x[0]))[0]
def sup(g):
    g = frozenset(g); c = L - g
    for ls, lab in nodes:
        if ls == g or ls == c: return lab or 'root'
    return 'NOT IN TREE'
P = lambda p: [x for x in L if x.startswith(p)]
ISL = P('ISL_'); LAR = ['BACT_Agrobacterium_larrymoorei_WP_309770759']
checks = [
    ('ISL monophyly (3 Fusarium)', ISL),
    ('F. poae pair', P('ISL_Fusarium_poae')),
    ('ISL + A. larrymoorei WP_309770759', ISL + LAR),
    ('nos seeds (P00386+Gillisella)', P('SEED_nos')),
    ('ocs + vis', P('SEED_ocs') + P('SEED_vis')),
    ('nos + Butyriboletus', P('SEED_nos') + P('FUNG_Buty')),
    ('Rhypophila', P('FUNG_Rhyp')),
    ('nos+Butyr+Rhypophila', P('SEED_nos') + P('FUNG_')),
    ('nos+FUNG+BACT2 larrymoorei', P('SEED_nos') + P('FUNG_') + P('BACT2')),
    ('T-DNA clade: seeds+FUNG+BACT2', P('SEED') + P('FUNG') + P('BACT2')),
    ('native fungal ODH (all NAT)', P('NAT_')),
    ('ISL+larry vs rest backbone = everything else', [x for x in L if x not in ISL + LAR]),
]
print('bipartition\tsupport_SHaLRT/UFBoot')
for n, g in checks: print(f'{n}\t{sup(g)}')
