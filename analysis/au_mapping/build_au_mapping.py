"""Build the AU constraint mapping (Table S11 companion) and the au_constraints data folder.
Leaf sets are derived from the actual constraint Newick files; texts are curated below.
Run:  python3 -I build_au_mapping.py
"""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import sys, re, os, csv, shutil, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from parse import parse, leaves, groups
from comp import clades, mono, prune_clades

G = f"{GMO_ROOT}"
A = G + "/analysis"
SUP = G + "/paper/supplementary"
OUT = SUP + "/data/au_constraints"
TF = A + "/trees_final"
FI = A + "/fusarium_odh_island"
CH = A + "/chs_fusarium"

def rd_first_tree(p):
    return open(p).readline()

# ---------------------------------------------------------------- run definitions
# run code -> dict(fam folder, ML tree, run files, list of (n, id, constraint nwk, search prefix))
def cons(d, names):  # names: list of (n, cid, constraint_path, search_prefix)
    return names
RUNS = collections.OrderedDict()
RUNS["SUS"] = dict(fam="SUS", ml=TF + "/SUS/SUS.mafft.iqtree.treefile", d=TF + "/SUS",
    run=[TF + "/SUS/SUS.AU." + e for e in ("iqtree", "log", "trees")],
    items=[(1, "H_tdna", TF + "/SUS/SUS.H_tdna.constraint.nwk", TF + "/SUS/SUS.H_tdna"),
           (2, "H_alt", TF + "/SUS/SUS.H_alt.constraint.nwk", TF + "/SUS/SUS.H_alt")])
RUNS["SUS_nat"] = dict(fam="SUS", ml=TF + "/SUS/SUS.mafft.iqtree.treefile", d=TF + "/SUS/au_sus_native",
    run=[TF + "/SUS/au_sus_native/" + x for x in ("au_sus_native.iqtree", "au_sus_native.log", "au_sus_native.trees", "au_trees.nwk", "hyp_order.txt", "groups.txt", "run_au_sus_native.sh")],
    items=[(3, "sus_natv_core", TF + "/SUS/au_sus_native/SUS.sus_natv_core.constraint.nwk", TF + "/SUS/au_sus_native/SUS.sus_natv_core"),
           (4, "sus_natv_all", TF + "/SUS/au_sus_native/SUS.sus_natv_all.constraint.nwk", TF + "/SUS/au_sus_native/SUS.sus_natv_all")])
RUNS["ACS"] = dict(fam="GDPD", ml=TF + "/ACS/ACS.mafft.iqtree.treefile", d=TF + "/ACS",
    run=[TF + "/ACS/ACS.lineage.AU." + e for e in ("iqtree", "log", "trees", "order.txt")] + [TF + "/au_lineage_acs.sh"],
    items=[(i + 1, h, TF + "/ACS/ACS.%s.constraint.nwk" % h, TF + "/ACS/ACS." + h)
           for i, h in enumerate("six_tdna six_kerr six_acc six_natv meso_acc meso_tdna".split())])
RUNS["ACS_nat"] = dict(fam="GDPD", ml=TF + "/ACS/ACS.mafft.iqtree.treefile", d=TF + "/ACS/au_acc_native",
    run=[TF + "/ACS/au_acc_native/" + x for x in ("au_acc_native.iqtree", "au_acc_native.log", "au_acc_native.trees", "au_trees.nwk", "hyp_order.txt", "groups.txt", "run_au_acc_native.sh")],
    items=[(7, "meso_natv", TF + "/ACS/au_acc_native/ACS.meso_natv.constraint.nwk", TF + "/ACS/au_acc_native/ACS.meso_natv")])
RUNS["CUS"] = dict(fam="CUS", ml=TF + "/CUS/CUS.mafft.iqtree.treefile", d=TF + "/CUS",
    run=[TF + "/CUS/CUS.lineage_x.AU." + e for e in ("iqtree", "log", "trees", "order.txt")] + [TF + "/au_lineage_cus.sh"],
    items=[(i + 1, h, TF + "/CUS/CUS.%s.constraint.nwk" % h, TF + "/CUS/CUS." + h)
           for i, h in enumerate("cus_tdna_x cus_tdna4_x cus_act_x cus_mis_x cus_bact_x cus_mono_x bol_tdna_x bol_act_x bol_mis_x bol_othercus_x".split())])
RUNS["CUS_nb"] = dict(fam="CUS", ml=TF + "/CUS/CUS.mafft.iqtree.treefile", d=TF + "/CUS",
    run=[TF + "/CUS/CUS.lineage_nb.AU." + e for e in ("iqtree", "log", "trees", "order.txt")],
    items=[(i + 1, h, TF + "/CUS/CUS.%s.constraint.nwk" % h, TF + "/CUS/CUS." + h)
           for i, h in enumerate("nb_tdna_x nb_tdna4_x nb_act_x nb_mis_x nb_bact_x nb_mono_x".split())])
RUNS["CUS_bm"] = dict(fam="CUS", ml=TF + "/CUS/CUS.mafft.iqtree.treefile", d=TF + "/CUS/au_bolmono/run2",
    run=[TF + "/CUS/au_bolmono/run2/" + x for x in ("CUS.lineage_bm_s2.AU.iqtree", "CUS.lineage_bm_s2.AU.log", "CUS.lineage_bm_s2.AU.trees", "CUS.lineage_bm_s2.candidates.nwk", "run_au_bolmono_s2.sh")],
    items=[(5, "bol_mono_x", TF + "/CUS/au_bolmono/run2/CUS.bol_mono_x.constraint.nwk", TF + "/CUS/au_bolmono/run2/CUS.bol_mono_x_s2")])
RUNS["CHS"] = dict(fam="CHS", ml=None, mlfile=CH + "/au/au_trees.nwk", d=CH + "/au",
    run=[CH + "/au/" + x for x in ("au_test.iqtree", "au_test.log", "au_test.trees", "au_trees.nwk", "run_au.sh")],
    items=[(i + 1, h, CH + "/au/%s.nwk" % h, CH + "/au/" + h) for i, h in enumerate("chs_tdna chs_bactchs chs_ags chs_nat chs_bact chs_mono".split())])
_nest = "nest_seed fnest_seed nest_b139 fnest_b139 nest_b272871 fnest_b272871 nest_b272951 fnest_b272951 nest_b084 fnest_b084".split()
RUNS["CHS_nest"] = dict(fam="CHS", ml=None, mlfile=CH + "/au_nested/au_trees.nwk", d=CH + "/au_nested",
    run=[CH + "/au_nested/" + x for x in ("au_nested.iqtree", "au_nested.log", "au_nested.trees", "au_trees.nwk", "hyp_order.txt", "run_au_nested.sh")],
    items=[(i + 1, h, CH + "/au_nested/%s.nwk" % h, CH + "/au_nested/" + h) for i, h in enumerate(_nest)])
NC = FI + "/nos_check"
RUNS["NOS"] = dict(fam="NOS", ml=None, mlfile=NC + "/au_trees.nwk", d=NC,
    run=[NC + "/" + x for x in ("au_test.iqtree", "au_test.log", "au_test.trees", "au_trees.nwk", "run_au.sh")],
    items=[(1, "c1_bol_in_tdna", NC + "/c1_bol_in_tdna.nwk", NC + "/c1_bol_in_tdna"),
           (2, "c2_rh_in_tdna", NC + "/c2_rh_in_tdna.nwk", NC + "/c2_rh_in_tdna"),
           (3, "c3_fungal_mono", NC + "/c3_fungal_mono.nwk", NC + "/c3_fungal_mono")])
RUNS["NOS_nat"] = dict(fam="NOS", ml=None, mlfile=NC + "/au_native/au_trees.nwk", d=NC + "/au_native",
    run=[NC + "/au_native/" + x for x in ("au_native.iqtree", "au_native.log", "au_native.trees", "au_trees.nwk", "hyp_order.txt", "run_au_native.sh")],
    items=[(4, "c4_bol_nat", NC + "/au_native/c4_bol_nat.nwk", NC + "/au_native/c4_bol_nat"),
           (5, "c5_bol_rh_nat", NC + "/au_native/c5_bol_rh_nat.nwk", NC + "/au_native/c5_bol_rh_nat")])
RUNS["NOS_rh"] = dict(fam="NOS", ml=None, mlfile=NC + "/au_rhypo/au_trees.nwk", d=NC + "/au_rhypo",
    run=[NC + "/au_rhypo/" + x for x in ("au_rhypo.iqtree", "au_rhypo.log", "au_rhypo.trees", "au_trees.nwk", "hyp_order.txt", "groups.txt", "run_au_rhypo.sh")],
    items=[(6, "c6_rh_2sord", NC + "/au_rhypo/c6_rh_2sord.nwk", NC + "/au_rhypo/c6_rh_2sord"),
           (7, "c7_rh_sord", NC + "/au_rhypo/c7_rh_sord.nwk", NC + "/au_rhypo/c7_rh_sord"),
           (8, "c8_rh_nat", NC + "/au_rhypo/c8_rh_nat.nwk", NC + "/au_rhypo/c8_rh_nat")])
_odh = "H1_isl_native H1b_isl_sordario H1c_native_only H2_isl_tdna H2s_isl_seeds H3_isl_larry H4_isl_bact H5_isl_nos".split()
RUNS["ODH"] = dict(fam="ODH", ml=None, mlfile=FI + "/au/au_trees.nwk", d=FI + "/au",
    run=[FI + "/au/" + x for x in ("au_test.iqtree", "au_test.log", "au_test.trees", "au_trees.nwk", "au_trees.order.txt", "au_summary.tsv", "make_constraints.py", "run_au.sh")],
    items=[(i + 1, h, FI + "/au/%s.nwk" % h, FI + "/au/" + h) for i, h in enumerate(_odh)])

# ---------------------------------------------------------------- leaf-class helper
def key(l):
    if "__" in l:
        k = l.split("__")[0]
        return re.sub(r"_(intact|pseudogene_like|fragment|pseudogene)$", "", k)
    return "_".join(l.split("_")[:2]) if l.startswith("NOS_") else l.split("_")[0]

def grp(path):
    d = collections.OrderedDict(); cur = None
    for ln in open(path):
        ln = ln.rstrip("\n")
        if ln.endswith(":") or ln.startswith(("MESO", "NATV")):
            cur = ln; d[cur] = []
        elif ln.strip() and cur:
            d[cur].append(ln.strip())
    return d

# candidate sets (fungal genes under test) per run, for pruning; native sets for the explicit native-set check
def cand_of(run, all_leaves):
    if run.startswith("SUS"):
        g = grp(TF + "/SUS/au_sus_native/groups.txt")
        return frozenset(g[list(g)[0]])   # the 52 fungal sus-like genes of the sus-like clade
    if run.startswith("ACS"):
        six = [l for l in groups(TF + "/ACS/ACS.six_kerr.constraint.nwk")[1][0] if l.startswith("FUNGI_acs")]
        meso = [l for l in all_leaves if key(l) == "FUNGI_gdpd_meso"]
        return frozenset(six + meso)
    if run.startswith("CUS"):
        return frozenset(l for l in all_leaves if key(l).startswith("FUNGI_"))
    if run.startswith("CHS"):
        return frozenset(l for l in all_leaves if l.startswith("CHSlike"))
    if run.startswith("NOS"):
        return frozenset(l for l in all_leaves if l.startswith(("BOL", "FUNG")))
    if run == "ODH":
        return frozenset(l for l in all_leaves if l.startswith("ISL"))

def native_of(run, cons_leaves, all_leaves):
    """explicit 'native fungal' partner leaves inside the constraint, or None"""
    if run == "SUS_nat":
        g = grp(TF + "/SUS/au_sus_native/groups.txt"); ks = list(g)
        nat = set(g[ks[1]] + g[ks[2]] + (g[ks[3]] if "all" in cons_leaves[0] else []))
        return nat
    if run in ("ACS", "ACS_nat"):
        g = grp(TF + "/ACS/au_acc_native/groups.txt"); ks = list(g)
        return set(g[ks[1]])
    if run.startswith("NOS") or run == "ODH":
        return {l for l in cons_leaves[1] if l.startswith("NAT")}
    if run.startswith("CHS"):
        return {l for l in cons_leaves[1] if l.startswith("AGSnat")}
    return None

# ---------------------------------------------------------------- texts (curated)
# per cid: (biological_alternative, constraint_definition, interpretation, caveat)
NATCAV = ("The native fungal set forced into the clade is not itself monophyletic in the unconstrained ML tree (also not after the fungal genes under test are pruned from it); "
          "a rejection therefore shows that the fungal genes do not fall inside one clade together with all of these native genes, "
          "not by itself that they are of non-fungal origin.")
T = {}
def a(cid, *x): T[cid] = x
a("SUS-1", "The fungal sus-like genes are related to the T-DNA sus of agrobacterial plasmids (bacterial origin by transfer from a plasmid-bearing donor).",
  "Forced into one clade: 56 fungal OCD-fold leaves (52 sus-like genes + 4 fungal OCD-fold genes outside the sus-like clade) + 24 plasmid sus references (REF T1/T2 sus). All other leaves free.",
  "Not rejected: a single clade of the fungal sus-like genes and the plasmid T-DNA sus is as good as the ML tree. Compatible with a bacterial T-DNA origin of the fungal sus-like genes; it does not exclude other topologies and does not by itself show the direction of transfer.",
  "Positive (compatibility) test; non-rejection is not proof.")
a("SUS-2", "The fungal sus-like genes belong with the other (non-T-DNA) bacterial OCD-fold proteins rather than with the T-DNA sus.",
  "Forced into one clade: 56 fungal OCD-fold leaves + 14 other-bacterial OCD-fold proteins + 19 rhizobial chromosomal OCD-fold proteins + 10 outgroup proteins (5 fungal OCD-fold) + 2 native Sordariomycetes OCD-fold proteins. The 24 T-DNA sus references are left outside the clade.",
  "Rejected: the fungal sus-like genes do not group with the non-T-DNA OCD-fold proteins as one clade separate from the T-DNA sus; the T-DNA sus is the closer relative.",
  "Partner set contains fungal outgroup and native fungal proteins as well as bacterial ones; the test excludes this particular mixed grouping, not every possible non-T-DNA relationship.")
a("SUS_nat-3", "The fungal sus-like genes are descendants of native fungal OCD-fold genes (single native origin).",
  "Forced into one clade: 56 fungal OCD-fold leaves (52 sus-like + 4 outside the sus-like clade) + 2 native Sordariomycetes OCD-fold proteins (Pseudopestalotiopsis) = 58 leaves.",
  "Rejected: the fungal sus-like genes and the available native fungal OCD-fold homologs do not form one clade; evidence against a single native fungal origin of the sus-like genes.",
  NATCAV + " No control with the native set alone was run for this family.")
a("SUS_nat-4", "As SUS_nat-3, with all fungal OCD-fold proteins of the tree.",
  "Forced into one clade: the 58 leaves of SUS_nat-3 + the 5 fungal OCD-fold outgroup proteins = 63 leaves.",
  "Rejected: all fungal OCD-fold proteins of the tree do not form one clade; evidence against a single fungal ancestry of this set.",
  NATCAV + " No control with the native set alone was run for this family.")
a("ACS-1", "The pKerr14a-type acs-like fungal genes are related to the T-DNA acs.",
  "Forced into one clade: 6 fungal pKerr14a-type acs-like genes + 27 T-DNA acs references (8 REF T1 + 19 REF T2). Other leaves free.",
  "Rejected: the fungal pKerr14a-type genes do not group with the T-DNA acs clade; they are not nested in or sister to the T-DNA acs clade.",
  "Negative result for the T-DNA acs relationship; the genes are placed with the pKerr14a-type acs-like references instead (ACS-2).")
a("ACS-2", "The fungal pKerr14a-type genes are related to the acs-like genes of plasmid pKerr14a and the Parasponia cT-DNA fragment.",
  "Forced into one clade: 6 fungal pKerr14a-type genes + pKerr14a acs-like references (2: WP_280174108.1, CUX65485.1) + Parasponia PON40589 = 9 leaves.",
  "Not rejected: the fungal pKerr14a-type genes can form one clade with the pKerr14a acs-like references and the Parasponia protein. Compatible with a bacterial (pKerr14a-type) origin.",
  "Small reference set; non-rejection is not proof and does not give the direction.")
a("ACS-3", "The pKerr14a-type genes and the acc-type GDPD genes of fungi are one group (single GDPD lineage).",
  "Forced into one clade: 6 fungal pKerr14a-type genes + 60 fungal acc-type GDPD genes + 5 acc-locus bacterial references (2 REF T1, 1 REF T2, 2 acc-locus proteins) = 71 leaves.",
  "Not rejected: the two fungal GDPD groups together with the acc-locus bacterial references are compatible with one clade.",
  "Non-rejection (deltaL 62.7) shows only compatibility; it does not establish a common origin.")
a("ACS-4", "The pKerr14a-type genes are native fungal GDPD genes (vertical inheritance).",
  "Forced into one clade: 6 fungal pKerr14a-type genes + 138 native fungal GDPD (119 fungal GDPD from the acs-profile set + 19 native-fungal references of 8 fungal classes) = 144 leaves.",
  "Rejected: the fungal pKerr14a-type genes do not nest in a clade of the native fungal GDPD; evidence against a native fungal origin.",
  NATCAV + " No control with the native set alone was run.")
a("ACS-5", "The fungal acc-type GDPD genes group with the bacterial acc-locus references.",
  "Forced into one clade: 60 fungal acc-type GDPD + 5 acc-locus references = 65 leaves.",
  "Not rejected: the fungal acc-type genes can form one clade with the bacterial acc-locus references (the ML tree already has this grouping, deltaL = 0).",
  "Compatibility test only.")
a("ACS-6", "The acc-type GDPD genes are related to the T-DNA acs (same clade).",
  "Forced into one clade: 60 fungal acc-type GDPD + 27 T-DNA acs references = 87 leaves.",
  "Rejected: the acc-type genes do not group with the T-DNA acs; they are a separate bacterial-type lineage (acc locus).",
  "Negative test; separates the acc-type from the acs lineage.")
a("ACS_nat-7", "The fungal acc-type GDPD genes are native fungal GDPD genes (vertical inheritance).",
  "Forced into one clade: 60 fungal acc-type GDPD + 138 native fungal GDPD (119 fungal GDPD + 19 native references) = 198 leaves.",
  "Rejected: the acc-type GDPD genes do not fall in one clade with the native fungal GDPD; evidence against a native fungal origin.",
  NATCAV + " No control with the native set alone was run.")
CUSBM = "The Boletales cus-like genes"
a("CUS-1", "The fungal cus-like genes are related to the T-DNA cus of agrobacterial plasmids.",
  "Forced into one clade: 42 fungal cus-like (26 Boletales + 16 other fungi) + 21 plasmid cus references (8 REF T1 + 13 REF T2) = 63 leaves.",
  "Not rejected: the fungal cus-like genes can form one clade with the T-DNA cus; compatible with a single bacterial T-DNA origin of all fungal cus-like genes.",
  "Compatibility only; CUS-6 and CUS-10 reject a single fungal clade without the references, so the grouping depends on the references.")
a("CUS-2", "As CUS-1, with four bacterial references from the clade added.",
  "Forced into one clade: 42 fungal cus-like + 21 plasmid cus references + 4 other bacterial cus references (Trinickia, 2 Actinomycetota, Cohaesibacter) = 67 leaves.",
  "Not rejected (p = 1, deltaL = 0, the ML tree satisfies the constraint): the fungal cus-like genes lie within a bacterial cus clade that contains the T-DNA cus.",
  "Compatibility only.")
a("CUS-3", "The fungal cus-like genes belong with the actinobacterial-type fungal sequences (actinobacterial origin).",
  "Forced into one clade: 42 fungal cus-like + 130 fungal actinobacterial-type sequences + 9 actinobacterial cus references = 181 leaves.",
  "Rejected: the fungal cus-like genes do not form a clade with the actinobacterial-type sequences.",
  "Negative test.")
a("CUS-4", "The fungal cus-like genes belong with the mis lineage.",
  "Forced into one clade: 42 fungal cus-like + 190 fungal mis sequences + 8 mis references = 240 leaves.",
  "Rejected: the fungal cus-like genes do not form a clade with the mis lineage.",
  "Negative test.")
a("CUS-5", "The fungal cus-like genes belong with other (non-T-DNA) bacterial cus homologs.",
  "Forced into one clade: 42 fungal cus-like + 26 other bacterial cus references (9 actinobacterial, 7 betaproteobacterial, 10 other) = 68 leaves.",
  "Rejected (weakly, p = 0.0005): the fungal cus-like genes do not form one clade with these bacterial references.",
  "The bacterial references are not monophyletic in the ML tree; deltaL is large (1107) although p-AU is only 5e-4.")
a("CUS-6", "All fungal cus-like genes derive from a single acquisition (one fungal clade).",
  "Forced into one clade: the 42 fungal cus-like genes only (26 Boletales + 16 other fungi); no references.",
  "Rejected: the 42 fungal cus-like genes do not form one clade; at least two separate entries (acquisitions) into fungi, or a non-fungal position of some, are needed.",
  "Does not say where each fungal subset came from. Leaf set identical to CUS-10 (separate search; p differs slightly).")
a("CUS-7", "The Boletales cus-like genes are related to the T-DNA cus, separately from the other fungal cus-like genes.",
  "Forced into one clade: 26 Boletales cus-like + 21 plasmid cus references = 47 leaves; the 16 other fungal genes cannot be inside it.",
  "Rejected: Boletales cus-like genes do not form a clade with the T-DNA cus separate from the other fungal cus-like genes.",
  "Negative test.")
a("CUS-8", "The Boletales cus-like genes belong with the actinobacterial-type fungal sequences.",
  "Forced into one clade: 26 Boletales cus-like + 130 actinobacterial-type fungal sequences + 9 actinobacterial references = 165 leaves.",
  "Rejected.", "Negative test.")
a("CUS-9", "The Boletales cus-like genes belong with the mis lineage.",
  "Forced into one clade: 26 Boletales cus-like + 190 mis sequences + 8 mis references = 224 leaves.",
  "Rejected.", "Negative test.")
a("CUS-10", "The Boletales cus-like genes together with the remaining fungal cus-like genes form one fungal clade (references excluded).",
  "Forced into one clade: the same 42 fungal cus-like leaves as CUS-6 (26 Boletales + 16 other fungi); no references.",
  "Rejected: as CUS-6.",
  "Leaf set identical to CUS-6; separate search (p = 1.6e-05 vs 6.4e-05). Not an independent test.")
a("CUS_nb-1", "The fungal cus-like genes outside Boletales are related to the T-DNA cus, separately from the Boletales ones.",
  "Forced into one clade: 16 non-Boletales fungal cus-like + 21 plasmid cus references = 37 leaves.",
  "Rejected.", "Negative test.")
a("CUS_nb-2", "As CUS_nb-1, with four bacterial references from the clade added.",
  "Forced into one clade: 16 non-Boletales fungal cus-like + 21 plasmid cus references + 4 other bacterial cus references = 41 leaves.",
  "Rejected.", "Negative test.")
a("CUS_nb-3", "The non-Boletales fungal cus-like genes belong with the actinobacterial-type fungal sequences.",
  "Forced into one clade: 16 non-Boletales fungal cus-like + 130 actinobacterial-type fungal sequences + 9 actinobacterial references = 155 leaves.",
  "Rejected.", "Negative test.")
a("CUS_nb-4", "The non-Boletales fungal cus-like genes belong with the mis lineage.",
  "Forced into one clade: 16 non-Boletales fungal cus-like + 190 mis + 8 mis references = 214 leaves.",
  "Rejected.", "Negative test.")
a("CUS_nb-5", "The non-Boletales fungal cus-like genes belong with other bacterial cus homologs.",
  "Forced into one clade: 16 non-Boletales fungal cus-like + 26 other bacterial cus references = 42 leaves.",
  "Rejected.", "Negative test.")
a("CUS_nb-6", "All non-Boletales fungal cus-like genes derive from one acquisition.",
  "Forced into one clade: the 16 non-Boletales fungal cus-like genes only; no references.",
  "Rejected: the 16 genes do not form a clade among themselves (separate placements in the tree).",
  "Does not say where each came from.")
a("CUS_bm-5", "All Boletales cus-like genes derive from a single acquisition by Boletales.",
  "Forced into one clade: the 26 Boletales cus-like genes only; no references and no other fungi.",
  "Rejected: the 26 Boletales cus-like genes do not form a clade; they are intermingled with other fungal genes (more than one entry, or transfers within fungi).",
  "Separate user-tree set from CUS (the four Boletales constraints were re-tested in this run; p-values are comparable only within a run).")
a("CHS-1", "The fungal chs-like genes are descendants of the Ti-plasmid chs gene of the seed (pTi chrysopine synthase), as a clade apart from other Ti-plasmid chs genes.",
  "Forced into one clade: 11 fungal chs-like (10 Fusarium + 1 Dactylonectria) + the chs seed A0A2Z2PGA9 = 12 leaves. The four other Ti-plasmid chs genes cannot be inside it.",
  "Rejected: the fungal genes are not sister to the seed alone; in the ML tree the seed is inside a clade of Ti-plasmid chs genes.",
  "Tests a narrow placement; the fungal genes can still be nested within the Ti-plasmid chs clade (see CHS_nest).")
a("CHS-2", "The fungal chs-like genes are sister to the other four Ti-plasmid chs genes (without the seed).",
  "Forced into one clade: 11 fungal chs-like + 4 Ti-plasmid chs (WP_139103025.1, WP_272871282.1, WP_272951408.1, WP_084631803.1) = 15 leaves; the seed cannot be inside it.",
  "Rejected.", "Tests a narrow placement; see CHS_nest.")
a("CHS-3", "The fungal chs-like genes are related to the yeast ags (ags seed P27875).",
  "Forced into one clade: 11 fungal chs-like + ags seed P27875 = 12 leaves.",
  "Rejected.", "Single-leaf partner; tests the position of the ags seed only.")
a("CHS-4", "The fungal chs-like genes are descendants of native fungal ags homologs (vertical inheritance).",
  "Forced into one clade: 11 fungal chs-like + 19 native fungal ags homologs = 30 leaves.",
  "Rejected (p = 0.0017, deltaL 49.6): the fungal chs-like genes do not form one clade with the native fungal ags homologs; evidence against a native fungal origin.",
  "The native ags set is monophyletic in the ML tree (also when the fungal genes are pruned), so the native-set-non-monophyly caveat does not apply here; the rejection is modest in size (deltaL 49.6).")
a("CHS-5", "The fungal chs-like genes belong with other bacterial M24 aminopeptidases (not Ti-plasmid chs).",
  "Forced into one clade: 11 fungal chs-like + 22 other bacterial M24 aminopeptidases = 33 leaves.",
  "Rejected.", "Negative test.")
a("CHS-6", "All fungal chs-like genes (Fusarium + Dactylonectria) derive from one acquisition.",
  "Forced into one clade: the 11 fungal chs-like genes only; no references.",
  "Not rejected (p = 0.17, deltaL 4.6): a single fungal clade is as good as the ML tree; a single acquisition (or transfer within fungi) is compatible.",
  "Compatibility only; the direction and source are not tested by this constraint.")
for k, d, kk in (("seed", "the chs seed (pTi chrysopine synthase A0A2Z2PGA9)", "CHS_nest-1"), ):
    pass
_ti = [("seed", "the chs seed A0A2Z2PGA9"), ("b139", "Ti-plasmid chs WP_139103025.1"), ("b272871", "Ti-plasmid chs WP_272871282.1 (A. tumefaciens)"),
       ("b272951", "Ti-plasmid chs WP_272951408.1 (A. vitis)"), ("b084", "Ti-plasmid chs WP_084631803.1 (A. larrymoorei)")]
for i, (k, d) in enumerate(_ti):
    a("CHS_nest-%d" % (2 * i + 1), "The fungal chs-like genes are nested within the Ti-plasmid chs clade, as sister to %s." % d,
      "Forced into one clade: all 11 fungal chs-like + %s = 12 leaves." % d,
      "Rejected: the fungal genes (all 11) are not sister to this single Ti-plasmid chs gene.",
      "Narrow placement test of one position among several; rejection of one position does not reject nesting elsewhere in the clade.")
    a("CHS_nest-%d" % (2 * i + 2), "As the previous constraint, for the 10 Fusarium genes only (Dactylonectria free).",
      "Forced into one clade: 10 Fusarium chs-like + %s = 11 leaves; the Dactylonectria gene is left free." % d,
      "Rejected: the Fusarium genes are not sister to this single Ti-plasmid chs gene.",
      "Narrow placement test; see previous constraint.")
a("NOS-1", "The Boletales nos-subfamily opine dehydrogenase genes are related to the T-DNA nos.",
  "Forced into one clade: 6 Boletales nos-subfamily genes + the 9 members of the T-DNA nos clade of the ML tree (8 NOS T-DNA + 1 Agrobacterium WP_174048180). Other nos-like rhizobial paralogs free. 15 leaves.",
  "Rejected (p = 0.003, deltaL 35.6): the Boletales genes do not fit inside one clade with the T-DNA nos clade.",
  "Small-to-moderate deltaL (35.6); a likelihood loss, not a large one.")
a("NOS-2", "The Rhypophila nos-subfamily genes are related to the T-DNA nos.",
  "Forced into one clade: 3 Rhypophila nos-subfamily genes + the 9 members of the T-DNA nos clade = 12 leaves.",
  "Rejected (p = 0.009, deltaL 41.8).", "As NOS-1.")
a("NOS-3", "All fungal nos-subfamily genes (Boletales + Rhypophila) derive from a single acquisition.",
  "Forced into one clade: 6 Boletales + 3 Rhypophila nos-subfamily genes = 9 leaves; no references.",
  "Not rejected (p = 0.35, deltaL 5.0): a single fungal clade is compatible.",
  "Compatibility only; source not tested.")
a("NOS_nat-4", "The Boletales nos-subfamily genes are native fungal opine dehydrogenases (vertical inheritance).",
  "Forced into one clade: 6 Boletales nos-subfamily genes + 19 native fungal opine dehydrogenases = 25 leaves.",
  "Rejected: the Boletales genes and the native fungal ODH do not form one clade.",
  NATCAV + " No control with the native set alone was run for this family (compare ODH-3).")
a("NOS_nat-5", "All fungal nos-subfamily genes (Boletales + Rhypophila) derive from native fungal ODH (single native origin).",
  "Forced into one clade: 6 Boletales + 3 Rhypophila nos-subfamily genes + 19 native fungal opine dehydrogenases = 28 leaves.",
  "Rejected.", NATCAV)
a("NOS_rh-6", "The Rhypophila nos-subfamily genes are native fungal ODH, sister to their nearest Sordariomycetes relatives.",
  "Forced into one clade: 3 Rhypophila genes + the 2 native Sordariomycetes ODH adjacent to them in the ML tree (Cryphonectria parasitica, Gnomoniopsis smithogilvyi) = 5 leaves.",
  "Rejected (p = 0.002, deltaL 51.3): the Rhypophila genes do not group with their nearest native relatives, even with only these two relatives.",
  "The two native genes are monophyletic in the ML tree, so the native-set caveat does not apply here; but a two-leaf partner is a narrow test.")
a("NOS_rh-7", "The Rhypophila genes are native fungal ODH of Sordariomycetes (same class).",
  "Forced into one clade: 3 Rhypophila + 12 native Sordariomycetes ODH = 15 leaves.",
  "Rejected.", NATCAV)
a("NOS_rh-8", "The Rhypophila genes are native fungal ODH (any class).",
  "Forced into one clade: 3 Rhypophila + 19 native fungal ODH = 22 leaves.",
  "Rejected.", NATCAV)
a("ODH-1", "The odh-island genes of F. poae / F. langsethiae are native fungal opine dehydrogenases (vertical inheritance).",
  "Forced into one clade: 3 island genes (F. poae Fp133, F. poae DAOMC252244, F. langsethiae MFG217702) + 9 native fungal ODH = 12 leaves.",
  "Rejected (p = 3.5e-06, deltaL 305).",
  NATCAV + " The control ODH-3 (the 9 native ODH alone forced into one clade) is itself rejected (p = 0.0064, deltaL 36.5); only the difference in deltaL (about 269) is attributable to adding the island genes.")
a("ODH-2", "The island genes belong with their closest native relatives (Sordariomycetes ODH; vertical-inheritance expectation).",
  "Forced into one clade: 3 island genes + 2 native Sordariomycetes ODH (Apiospora, Astrocystis) = 5 leaves.",
  "Rejected (p = 0.0004, deltaL 467): island genes are not sister to the closest native Sordariomycetes ODH.",
  "The two native Sordariomycetes ODH are not monophyletic in the ML tree, so part of the penalty comes from forcing them together (see ODH-3). ")
a("ODH-3", "Control: the native fungal ODH form one clade (island genes free).",
  "Forced into one clade: the 9 native fungal ODH only = 9 leaves; island genes and all others free.",
  "Rejected (p = 0.0064, deltaL 36.5): the native fungal ODH of this tree are not one clade.",
  "Shows that the native set is itself not monophyletic in the ML tree; hence a rejected 'island + native' constraint cannot by itself prove HGT. Applies to ODH-1 and ODH-2, to NOS_nat-4/5, NOS_rh-7/8, ACS-4, ACS_nat-7, SUS_nat-3/4.")
a("ODH-4", "The island genes lie within the ML T-DNA ocs/nos/vis clade (donor-type position).",
  "Forced into one clade: 3 island genes + 4 seeds (nos P00386, nos Gillisella, ocs P0A395, vis A0A1S7UAT7) + 3 Rhypophila/Butyriboletus nos-subfamily genes + A. larrymoorei WP_309770762 = 11 leaves.",
  "Rejected (p = 2.5e-58, deltaL 125): the island genes are not inside the T-DNA ocs/nos/vis clade of the ML tree.",
  "Negative test; in the ML tree the island genes are sister to the A. larrymoorei WP_309770759 (ODH-6).")
a("ODH-5", "The island genes are related to the T-DNA ocs/nos/vis seeds only.",
  "Forced into one clade: 3 island genes + 4 seeds = 7 leaves.",
  "Rejected (p = 3.5e-57).", "Negative test.")
a("ODH-6", "Control: the island genes are sister to the A. larrymoorei ODH WP_309770759 (clade of the ML tree).",
  "Forced into one clade: 3 island genes + A. larrymoorei WP_309770759 = 4 leaves.",
  "Not rejected (p = 1, constraint equals the ML clade): the island genes group with the plasmid ODH of A. larrymoorei.",
  "Positive control and the main support for the bacterial affinity; does not give direction.")
a("ODH-7", "The island genes belong with other (non-larrymoorei) bacterial ODH.",
  "Forced into one clade: 3 island genes + 23 other bacterial ODH (both A. larrymoorei ODH excluded) = 26 leaves.",
  "Rejected (p = 2.1e-07, deltaL 283).", "Negative test.")
a("ODH-8", "The island genes lie inside the nos subfamily.",
  "Forced into one clade: 3 island genes + nos seed P00386 + Gillisella nos WP_070167540 = 5 leaves.",
  "Rejected (p = 1.4e-51).", "Negative test.")

# ---------------------------------------------------------------- build
s11 = list(csv.DictReader(open(SUP + "/TableS11_Approximately_unbiased_AU_topology_tests.tsv"), delimiter="\t"))
rows_by_run = []
run = -1
for r in s11:
    if r["tree_number"] == "0":
        run += 1
    else:
        rows_by_run.append((run, r))
ORDER = ["SUS", "SUS_nat", "CUS", "CUS_nb", "CUS_bm", "ACS", "ACS_nat", "CHS", "CHS_nest", "NOS", "NOS_nat", "NOS_rh", "ODH"]
sidx = {}
for run, r in rows_by_run:
    sidx[(ORDER[run], int(r["tree_number"]))] = r
assert len(sidx) == 60, len(sidx)

if os.path.isdir(OUT):
    shutil.rmtree(OUT)
os.makedirs(OUT)
mapping = []; manifest = []
def cp(src, dst):
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copy2(src, dst); return os.path.relpath(dst, OUT)

for rc, R in RUNS.items():
    mlt = parse(rd_first_tree(R["mlfile"]) if R.get("mlfile") else rd_first_tree(R["ml"]))
    al = frozenset(leaves(mlt)); cl = clades(mlt)
    cand = cand_of(rc, al)
    pal, pcl = prune_clades(mlt, cand)
    fam_dir = OUT + "/" + R["fam"] + "/" + rc
    runfiles = []
    for p in R["run"] + [R.get("mlfile") or R["ml"]]:
        if not os.path.exists(p): continue
        runfiles.append(cp(p, fam_dir + "/run/" + os.path.basename(p)))
    runfiles = sorted(set(runfiles))
    for n, name, cpath, pref in R["items"]:
        cid = "%s-%d" % (rc, n)
        allf, g = groups(cpath)
        gg = g[0]; assert len(g) == 1, (cid, len(g))
        files = [cp(cpath, fam_dir + "/constraint/" + os.path.basename(cpath))]
        for e in (".treefile", ".iqtree", ".log"):
            p = pref + e
            if os.path.exists(p):
                files.append(cp(p, fam_dir + "/search/" + os.path.basename(p)))
        cs = frozenset(gg)
        cd = frozenset(l for l in cs if l in cand); part = cs - cd
        nat = native_of(rc, (cpath, gg), al)
        if nat is not None:
            nat = frozenset(nat) & cs
        # monophyly in unconstrained ML tree
        part_mono = mono(part, pal, pcl) if len(part) > 1 else ("n.a." if not part else True)
        nat_mono = (mono(nat, pal, pcl) if nat and len(nat) > 1 else ("n.a." if not nat else True)) if nat is not None else "n.a."
        def yn(x): return x if isinstance(x, str) else ("yes" if x else "no")
        r11 = sidx[(rc, n)]
        bio, cdef, interp, cav = T[cid]
        if nat is None and len(part) > 1 and not part_mono:
            cav = cav + " The reference/partner set forced into the clade (%d leaves) is not itself a clade in the unconstrained ML tree (checked with the fungal genes under test pruned), so part of the likelihood loss may come from forcing it together." % len(part)
        mapping.append(dict(constraint_id=cid, family=r11["family"], tree_number=n,
            tested_topology_constraint=r11["tested_topology_constraint"], p_AU=r11["p_AU"], deltaL=r11["deltaL"], verdict=r11["verdict"],
            biological_alternative=bio, constraint_definition=cdef, constraint_leaves=len(gg), candidate_leaves=len(cd), partner_leaves=len(part),
            partner_set_monophyletic_in_ML_tree=yn(part_mono), native_fungal_set_in_constraint=("n.a." if nat is None or not nat else "%d leaves" % len(nat)),
            native_set_monophyletic_in_ML_tree=yn(nat_mono),
            interpretation=interp, caveat=cav, constraint_file=files[0]))
        manifest.append(dict(family=R["fam"], run=rc, constraint_id=cid, role="constraint", files=";".join(files)))
    manifest.append(dict(family=R["fam"], run=rc, constraint_id="(all constraints of the run)", role="AU test and ML/user trees", files=";".join(runfiles)))

cols = ["constraint_id", "family", "tree_number", "tested_topology_constraint", "p_AU", "deltaL", "verdict", "biological_alternative", "constraint_definition",
        "constraint_leaves", "candidate_leaves", "partner_leaves", "partner_set_monophyletic_in_ML_tree", "native_fungal_set_in_constraint",
        "native_set_monophyletic_in_ML_tree", "interpretation", "caveat", "constraint_file"]
with open(SUP + "/TableS11_constraint_mapping.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(mapping)
with open(OUT + "/MANIFEST.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, ["family", "run", "constraint_id", "role", "files"], delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(manifest)
print(len(mapping), "constraints;", len(manifest), "manifest rows")
for m in mapping:
    print(m["constraint_id"], m["verdict"], m["constraint_leaves"], m["candidate_leaves"], m["partner_leaves"], m["partner_set_monophyletic_in_ML_tree"], m["native_fungal_set_in_constraint"], m["native_set_monophyletic_in_ML_tree"])
