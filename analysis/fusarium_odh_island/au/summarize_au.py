#!/usr/bin/env python3
"""Parse au_test.iqtree (USER TREES table) -> au_summary.tsv. Order of trees: au_trees.order.txt."""
import os, re
D = os.path.dirname(os.path.abspath(__file__))
order = [l.strip() for l in open(os.path.join(D, 'au_trees.order.txt')) if l.strip()]
desc = {'ML': 'unconstrained ML tree (Q.PFAM+G4)'}
for l in list(open(os.path.join(D, 'hypotheses.tsv')))[1:]:
    f = l.rstrip('\n').split('\t'); desc[f[0]] = f[2]
txt = open(os.path.join(D, 'au_test.iqtree')).read()
blk = txt.split('USER TREES')[1]
rows = []
for line in blk.splitlines():
    m = re.match(r'^\s*(\d+)\s+(-[\d.]+)\s+([\d.e-]+)\s+(.*)$', line)
    if not m: continue
    i, logl, dl, rest = m.groups()
    toks = rest.split()
    # bp-RELL [+-], p-KH [+-], p-SH [+-], c-ELW [+-], p-AU [+-]
    vals = [t for t in toks if t not in '+-']
    rows.append((int(i), logl, dl, vals[0], vals[1], vals[2], vals[3], vals[4]))
assert len(rows) == len(order), (len(rows), len(order))
with open(os.path.join(D, 'au_summary.tsv'), 'w') as fo:
    fo.write('hypothesis\tdescription\tlogL\tdeltaL\tbp-RELL\tp-KH\tp-SH\tc-ELW\tp-AU\tverdict_AU_0.05\n')
    for (i, logl, dl, bp, kh, sh, elw, au), h in zip(rows, order):
        v = 'rejected' if float(au) < 0.05 else 'not rejected'
        fo.write(f'{h}\t{desc[h]}\t{logl}\t{dl}\t{bp}\t{kh}\t{sh}\t{elw}\t{au}\t{v}\n')
print(open(os.path.join(D, 'au_summary.tsv')).read())
