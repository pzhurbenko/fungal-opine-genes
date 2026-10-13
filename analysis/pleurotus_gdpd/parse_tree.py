import sys
from Bio import Phylo
t=Phylo.read(sys.argv[1],"newick")
def cat(n):
    n=n.split("__")[0]
    if n.startswith("FUNGI_acs"): return "FUNGI_acs"
    if n.startswith("FUNGI_gdpd_meso"): return "FUNGI_gdpd_meso"
    if n.startswith("NATIVE"): return "NATIVE"
    if n.startswith("REF_T") : return n.split("_",2)[2] if n.count("_")>=2 else n
    return n
# root on Saccharomyces outgroup? keep as is; reroot on OG yeast
og=[c for c in t.get_terminals() if "Saccharomyces_cerevisiae" in c.name]
if og: t.root_with_outgroup(og[0])
q=[c for c in t.get_terminals() if c.name.startswith("QUERY")]
node=t.common_ancestor(q)
path=t.get_path(node)
anc=[t.root]+path
for depth,a in enumerate(reversed(anc[:-1])):
    if depth>=6: break
    tips=[c.name for c in a.get_terminals()]
    from collections import Counter
    cnt=Counter(cat(x) for x in tips)
    print(f"ancestor-{depth+1}: support={a.confidence or a.name} n={len(tips)} {dict(cnt)}")
    if len(tips)<=12: print("   ", tips)
print("query clade support:", node.confidence or node.name, [c.name for c in node.get_terminals()])
