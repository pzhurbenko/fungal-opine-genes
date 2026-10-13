#!/usr/bin/env python3
"""Check the fixed order of genes and lineages (paper v2).
(a) row order in Table1_funnel.tsv, TableS3b_fate.tsv (by gene), Table2a_lineages.tsv, TableS1_all_lineages.tsv (by lineage);
(b) enumerations of >=2 gene/lineage names (separated by ',', 'and', 'or', '/') and first cells of markdown tables in the v2 texts.
Enumerations on a line containing '[check order]' are reported separately and not counted as violations.
Usage: python check_order.py   (exit code 1 if violations)
"""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import csv, re, sys
P = f"{GMO_ROOT}/paper/"
GENE_ORDER = ["sus", "cus", "chs", "acs", "acs_pKerr14a", "nos", "ocs", "vis", "ags", "mas1", "mas2", "mis"]
LINEAGE_ORDER = ["sus-like (T-DNA type)", "cus-like (T-DNA type)", "chs-like (T-DNA type)", "acc-type GDPD", "pKerr14a-type acs-like",
                 "nos-subfamily opine dehydrogenase (Boletales)", "odh island: opine dehydrogenase",
                 "odh island: Xaa-Pro aminopeptidase", "odh island: third gene", "nos-subfamily ODH homolog (Rhypophila, not assigned)", "OCD-like, outside T-DNA sus clade", "actinobacterial-type", "mis lineage",
                 "native fungal GDPD (main clade)", "other fungal GDPD (position unresolved)", "other opine dehydrogenase (Boletales)", "not a homolog"]
viol = []
def rd(p): return list(csv.DictReader(open(p), delimiter="\t"))
def grank(g): gs = g.split(","); return (min(GENE_ORDER.index(x) for x in gs), len(gs))
def check_rows(name, keys):
    bad = [i for i in range(1, len(keys)) if keys[i] < keys[i - 1]]
    for i in bad: viol.append(f"{name}: row {i + 2} out of order")
    print(f"[table] {name}: {'OK' if not bad else str(len(bad)) + ' violations'}")
# ---------------- (a) tables
check_rows("Table1_funnel.tsv", [grank(r["gene"]) for r in rd(P + "tables/Table1_funnel.tsv")])
check_rows("TableS3b_fate.tsv", [(grank(r["gene"]), r["organism"]) for r in rd(P + "tables/TableS3b_fate.tsv")])
check_rows("Table2a_lineages.tsv", [LINEAGE_ORDER.index(r["group"]) for r in rd(P + "tables/Table2a_lineages.tsv")])  # column renamed 2026-10-12x
def s1rank(l):
    b = l.replace(" (not placed)", "")
    return LINEAGE_ORDER.index(b) if b in LINEAGE_ORDER else None
s1 = [s1rank(r["lineage"]) for r in rd(P + "tables/TableS1_all_lineages.tsv")]
# orphan '(not placed)' labels inherit the rank of the preceding row (they are anchored after a lineage of the same family)
k, last = [], 0
for x in s1: last = x if x is not None else last; k.append(last)
check_rows("TableS1_all_lineages.tsv", k)
# ---------------- (b) texts
I = r"[*_]{0,2}"  # markdown emphasis
LIN = [  # (regex, rank on the lineage scale)
    (rf"{I}sus{I}-like", 0), (rf"{I}cus{I}-like", 1), (rf"{I}chs{I}-like", 2), (r"acc-type GDPD", 3), (rf"pKerr14a-type {I}acs{I}-like", 4),
    (rf"{I}nos{I}-subfamily (?:opine dehydrogenase|ODH)?(?: \(?Boletales\)?| of Boletales)", 5), (rf"{I}nos{I}-subfamily (?:opine dehydrogenase|ODH)?(?: \(?Rhypophila\)?| of {I}Rhypophila{I})", 6),
    (rf"{I}nos{I}-subfamily (?:opine dehydrogenase|ODH|genes?|lineages?)", 5.5),
    (rf"{I}odh{I}[- ]island", 7), (r"OCD-like", 10), (r"actinobacterial-type", 11), (rf"{I}mis{I} lineage", 12), (r"native fungal GDPD", 13),
    (r"other fungal GDPD", 14), (r"other opine dehydrogenase", 15)]
GEN = [(rf"(?<![\w-]){I}acs{I} \(pKerr14a\)", 4), (rf"(?<![\w-]){I}acs{I}_pKerr14a", 4), (rf"pKerr14a {I}acs{I}(?![\w-])", 4)] + \
      [(rf"(?<![\w-]){I}{g}{I}(?![\w-]| \(pKerr14a\)| lineage|-like|-subfamily|-type| island|_)", i) for i, g in enumerate(GENE_ORDER) if g != "acs_pKerr14a"]
SEP = re.compile(r"^(?:\s*\([^()]{0,60}\))?[\s*_]*(?:,|;|/|\+|,? and|,? or|,? и)?[\s*_]*$")  # an item may carry a short parenthetical: "*sus* (succinamopine synthase), *cus* (...)"
def items(line):
    found, taken = [], []
    for kind, pats in (("L", LIN), ("G", GEN)):
        for pat, rk in pats:
            for m in re.finditer(pat, line):
                if any(m.start() < e and s < m.end() for s, e in taken): continue
                found.append((m.start(), m.end(), kind, rk, m.group())); taken.append((m.start(), m.end()))
    return sorted(found)
def enumerations(line):
    it, groups, cur = items(line), [], []
    for x in it:
        if cur and SEP.match(line[cur[-1][1]:x[0]]) and x[2] == cur[-1][2]: cur.append(x)
        else:
            if len(cur) >= 2: groups.append(cur)
            cur = [x]
    if len(cur) >= 2: groups.append(cur)
    return groups
flagged = []
for fn in ["draft_v2_results.md", "draft_v2_methods.md", "draft_v2_intro_additions.md"]:
    lines = open(P + fn).read().split("\n")
    for n, line in enumerate(lines, 1):
        for g in enumerations(line):
            rk = [x[3] for x in g]
            if any(rk[i] < rk[i - 1] for i in range(1, len(rk))):
                msg = f"{fn}:{n}: " + " | ".join(x[4] for x in g) + "   <<" + line[max(0, g[0][0] - 50):g[-1][1] + 30] + ">>"
                (flagged if "[check order]" in line else viol).append(msg)
    # markdown tables: first cells that name a lineage must follow LINEAGE_ORDER within one table
    tab = []
    for n, line in enumerate(lines + [""], 1):
        if line.startswith("|"):
            c = line.split("|")[1]; it = [x for x in items(c) if x[2] == "L"]
            if it: tab.append((n, it[0][3], it[0][4]))
        else:
            for i in range(1, len(tab)):
                if tab[i][1] < tab[i - 1][1]: viol.append(f"{fn}:{tab[i][0]}: table row '{tab[i][2]}' after '{tab[i - 1][2]}'")
            tab = []
print(f"[text] enumerations out of order: {sum(1 for v in viol if '.md:' in v)}; marked [check order]: {len(flagged)}")
for v in viol: print("VIOLATION", v)
for v in flagged: print("CHECK-ORDER", v)
print("violations:", len(viol))
sys.exit(1 if viol else 0)
