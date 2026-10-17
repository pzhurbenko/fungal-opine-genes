#!/usr/bin/env python3
"""Build Supplementary Tables (Biology/MDPI) from the finished analysis tables.

Re-runnable; reads only existing files, writes only into this directory.
Run with the paper_build python to also get Supplementary_Tables.xlsx (openpyxl).
"""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import csv, os, re, glob, html, collections

G = f"{GMO_ROOT}"
OUT = os.path.dirname(os.path.abspath(__file__))
A = G + "/analysis"

# ---------------------------------------------------------------- helpers
def rd(path, header=True, fields=None):
    with open(path, newline="") as fh:
        if header:
            return list(csv.DictReader(fh, delimiter="\t"))
        return [dict(zip(fields, r)) for r in csv.reader(fh, delimiter="\t")]

FAM = {"ACS": "GDPD", "CUS": "cucumopine synthase (Cus)", "SUS": "OCD fold (Sus)"}
def NOTE_FIX(v):
    v = v.replace("at the rule-B threshold", "at the identity threshold of T-DNA candidates")
    return re.sub(r"\s*\(20\d\d-\d\d-\d\d\)", "", v)
# 2026-10-13w: row-level notes added to the catalogue notes (S1b, hence S1a) and to S10b / S12
SPLIT_NOTE = "merged AUGUSTUS model split into two genes by HMM envelope (GDPD and mis); coordinates re-derived with miniprot"
MIS_NOTE = "AUGUSTUS had merged this gene with the nos-subfamily pseudogene 3.9 kb downstream (ODH_JBBPKX010000014.1); model re-derived as the single-exon ORF 13898-15079"
ODH_NOTE = "lies 3.9 kb downstream of the mis-group gene CUS_JBBPKX010000014_S2L2 (13898-15079), same strand; three nonsense mutations (Gln codons in the intact sisters)"
NOTE_ADD = {"ACS_JBBPKL010003317_L1": SPLIT_NOTE, "CUS_JBBPKL010003317_S2L1": SPLIT_NOTE, "ODH_JBBPKX010000014.1": ODH_NOTE, "CUS_JBBPKX010000014_S2L2": MIS_NOTE}
def add_note(old, new): return (old + "; " + new) if old else new
DETMAP = {"Table S1 (profile screen)": "Table S1b (profile screen)", "Table S12a": "Table S13a", "Table S14a": "Table S13b", "Table S11b": "Table S10b", "Table S3b": "Table S3"}
SUBG = {"kerr14": "pKerr14a-type acs-like", "gdpd_meso": "acc-type GDPD",
        "cus_actinobacteria": "actinobacterial-type", "sus": "sus-like", "cus": "cus-like",
        "mis": "mis group", "acs": "GDPD (acs profile), not placed on a tree", "OUTGROUP": "outgroup",
        "cus_other_bacteria": "other bacterial-type Cus",
        "cus_betaproteobacteria": "beta-proteobacterial-type Cus",
        "ocd_other_bacteria": "other bacterial OCD fold"}
STATUS = {"pseudogene_like": "pseudogene", "not_homolog": "not a homolog", "not_run": "not run"}

REPL = [
    (r'kerr14/acc-like \("six"\) lineage', "pKerr14a-type acs-like"),
    (r"acc-locus \(gdpd_meso\) clade", "acc-type GDPD"),
    (r"SUS T-DNA clade", "sus-like"), (r"CUS cus clade", "cus-like"), (r"\bCUS cus\b", "cus-like"),
    (r"ACS six lineage", "pKerr14a-type acs-like"), (r"ACS gdpd_meso vs acc refs", "acc-type GDPD vs acc-locus references"),
    (r"ACS gdpd_meso", "acc-type GDPD"),
    (r"gdpd_meso", "acc-type GDPD"), (r"cus_actinobacteria", "actinobacterial-type"),
    (r"pseudogene_like", "pseudogene"), (r"\bkerr14\b", "pKerr14a-type acs-like"),
    (r"FAM_ACS T1/T2", "GDPD family"), (r"FAM_ACS:", ""),
    (r"not in T-DNA/acc/six clades", "not in the T-DNA-type, acc-type or pKerr14a-type clades"),
    (r"the six tree-placed members", "the tree-placed pKerr14a-type acs-like members"),
    (r"mis lineage", "mis group"), (r"acs-HMM group", "GDPD (acs profile), not placed on a tree"),
    (r"cus-related small ORF", "other bacterial-type small ORF"), (r"cus-related ORF", "other bacterial-type ORF"),
]
REPL = [(re.compile(a), b) for a, b in REPL]

def clean(v):
    if not isinstance(v, str):
        return v
    for rx, b in REPL:
        v = rx.sub(b, v)
    return v

def disp(v):
    """display form written to files (internal logic keeps '(not placed)')"""
    if not isinstance(v, str):
        return v
    return clean(v).replace("(not placed)", "(profile only)")

def norm_lineage(v):
    # 2026-10-13h: displayed names of the odh-island genes (keys in the source tables are unchanged)
    v = v.replace("odh island: opine dehydrogenase", "odh island: odh").replace("odh island: Xaa-Pro aminopeptidase", "odh island: pepP").replace("odh island: third gene", "odh island: orf3")
    v = v.replace(" (T-DNA type)", "").replace("native fungal GDPD (main clade)", "native fungal GDPD")
    v = v.replace("other fungal GDPD (position unresolved)", "GDPD of unresolved position")
    v = v.replace("OCD-like, outside T-DNA sus clade", "OCD-like, outside the sus-like clade")
    v = v.replace("CUS", "Cus")
    return clean(v)

def num(x, nd=None):
    try:
        f = float(x)
    except (TypeError, ValueError):
        return None
    return f

TABLES = []   # (id, title, description, [(newcol, desc)], rows(list of dict newcol->val))

def proj(rows, spec):
    """spec: list of (src, new, desc[, fn]); returns rows with new keys, and column doc."""
    out = []
    for r in rows:
        o = {}
        for s in spec:
            src, new = s[0], s[1]
            v = r.get(src, "")
            if len(s) > 3:
                v = s[3](v)
            o[new] = v
        out.append(o)
    return out, [(s[1], s[2]) for s in spec]

def add(tid, title, desc, rows, spec, extra=None):
    rows, cols = proj(rows, spec)
    if extra:
        cols += extra
    TABLES.append((tid, title, desc, cols, rows))

fam_map = lambda v: FAM.get(v, v)
subg_map = lambda v: SUBG.get(v, v)
status_map = lambda v: STATUS.get(v, v)

# ---------------------------------------------------------------- S1
# S1 = unified catalogue of all lineages (paper/tables/TableS1_all_lineages.tsv, build_catalog.py: 2,791 profile-screen loci of the
# GDPD/Cus/Sus families + 73 loci of the targeted scans and blastx: chs-like, odh island, nos-subfamily ODH), joined with the
# per-locus details of the profile screen (TableS1_loci_catalog.tsv); journal 2026-10-12 suppl_sync_S1_S3b
# lineages of bacterial origin = main_lineage groups of Table 2a (same criterion as the column lineage_of_bacterial_origin of Table S18)
MAIN_RAW = {r["group"] for r in rd(G + "/paper/tables/Table2a_lineages.tsv") if r["main_lineage"] == "True"}
_cat = {r["locus_id"]: r for r in rd(G + "/paper/tables/TableS1_loci_catalog.tsv") if "ORF358" not in r["family"]}
s1 = []
for r in rd(G + "/paper/tables/TableS1_all_lineages.tsv"):
    c = _cat.get(r["locus_id"], {})
    row = dict(r)
    if r["locus_id"] in NOTE_ADD: row["note"] = add_note(r["note"], NOTE_ADD[r["locus_id"]])
    row["lineage_of_bacterial_origin"] = r["lineage"].replace(" (not placed)", "") in MAIN_RAW
    assert row["lineage_of_bacterial_origin"] == r["origin"].startswith(("T-DNA", "non-T-DNA")), r["locus_id"]
    row["basis"] = c.get("basis", "")
    row["hmm_subgroup"] = c.get("hmm_subgroup", "")
    row["division"] = c.get("division", ""); row["class"] = c.get("class", ""); row["family_tax"] = c.get("family_tax", "")
    row["species_cluster"] = c.get("species_cluster", ""); row["context_call"] = c.get("context_call", "")
    row["identity_T1T2"] = c.get("identity_T1T2", ""); row["species_source"] = c.get("species_source", "")
    s1.append(row)
assert len(s1) == 2864, len(s1)
s1_note = {}
PLEN_DESC = "Length of the fungal protein model, aa (for loci with a sequence in Data S1 it equals the length of that sequence; blank where no model could be built; for profile-screen loci with a gene model, the length of that model)"
# 2026-10-13af: Table S1d - sensitivity of the counts to the integration status (paper/tables/sensitivity_uncertain_loci.py)
_sens = rd(G + "/paper/tables/sensitivity_uncertain_loci.tsv")
add("S1d", "Sensitivity of the lineage counts to the integration status of the loci",
    "The counts of Table 2 (loci, intact loci, carrier species and assemblies per lineage) and the inputs of Figure 5 (genus rows and species-level clusters of Boletales with the lineage, and with an intact copy) "
    "computed from all 254 loci of the seven lineages (239 placed on the trees directly or through a representative, 13 fragments assigned by identity, 2 added from NCBI annotations after the trees were built; column all_loci) and from the loci whose integration in the fungal genome is confirmed or supported only (column confirmed_supported_only; Table S1a, column integration_status).",
    _sens, [("lineage", "lineage", "Lineage (prefixed by the table or figure the metric belongs to)"), ("metric", "metric", "Metric"),
            ("all_loci", "all_loci", "Value over all 254 loci of the seven lineages"), ("confirmed_supported_only", "confirmed_supported_only", "Value over the confirmed and supported loci only")])
add("S1b", "Complete catalogue of loci",
    "All loci of all lineages found in fungal genomes by the two routes of the screen: the profile screen (HMM/tblastn) of the "
    "GDPD, Cus and Sus families, and the targeted scans (chs-like loci, odh island, nos-subfamily opine dehydrogenases; "
    "details in Tables S10b, S13a and S13b). For each locus: family, assigned lineage and its origin, route of detection, evidence "
    "for the assignment, structural status and genomic coordinates. Loci marked 'not a homolog' are retained as a separate "
    "category (hits that failed the homology checks). The 254 loci of the seven lineages of bacterial origin (239 placed on the trees directly or through a representative, 13 fragments assigned by identity, 2 added from NCBI annotations after the trees were built) are listed separately, with additional columns, in Table S1a; "
    "the 152 fragments of the same lineages assigned by profile only are flagged here (column lineage_of_bacterial_origin = True, lineage with the suffix '(profile only)'). One row per locus.",
    s1, [
        ("locus_id", "locus_id", "Locus identifier (family prefix, contig, locus number)"),
        ("family", "family", "Protein family: GDPD, cucumopine synthase (Cus), OCD fold (Sus), M24 aminopeptidase (Chs/Ags) or opine dehydrogenase", fam_map),
        ("lineage", "lineage", "Final lineage assignment; '(profile only)' = assigned by HMM profile only, not placed on a tree", norm_lineage),
        ("origin", "origin", "Origin of the lineage as established in the paper (T-DNA gene; rhizobial locus outside the T-DNA; native fungal family; not analysed)"),
        ("lineage_of_bacterial_origin", "lineage_of_bacterial_origin", "True if the locus belongs to one of the seven lineages of bacterial origin (placed on the tree or assigned by profile only); the placed loci are also listed in Table S1a"),
        ("entry_route", "detection_route", "Route by which the locus was found: profile screen, targeted scan (miniprot), blastx of plasmid proteins, or protein-level screen"),
        ("source_table", "detailed_in", "Table with the per-locus details of the corresponding screen",
         lambda v: DETMAP.get(v, v)),
        ("basis", "assignment_basis", "Evidence for the lineage (profile-screen loci): phylogenetic tree, sequence similarity, or HMM profile only; blank for loci of the targeted scans (see detailed_in)",
         lambda v: {"hmm_only": "HMM only", "tree": "phylogenetic tree", "similarity": "sequence similarity"}.get(v, v)),
        ("hmm_subgroup", "hmm_profile_subgroup", "Best-matching HMM profile subgroup (the GDPD acs profile is not a tree-based lineage)", subg_map),
        ("status", "status", "Locus state: intact, pseudogene, fragment, remnant, ambiguous, not a homolog, not modelled (no gene model built)", status_map),
        ("organism", "organism", "Organism name (NCBI)"), ("species", "species", "Species name used in the paper"),
        ("division", "phylum", "Taxonomic phylum"), ("class", "class", "Taxonomic class"), ("order", "order", "Taxonomic order"),
        ("family_tax", "taxonomic_family", "Taxonomic family of the carrier fungus"), ("gca", "assembly", "Genome assembly accession"),
        ("species_cluster", "mash_species_cluster", "Genome cluster (mash) name; blank if not clustered"),
        ("contig", "contig", "Contig/scaffold accession"), ("contig_len", "contig_length_bp", "Contig length, bp"),
        ("start", "start", "Locus start on the contig (1-based)"), ("end", "end", "Locus end on the contig"),
        ("strand", "strand", "Strand"), ("n_exons", "n_exons", "Number of exons in the gene model"),
        ("protein_len", "protein_length_aa", PLEN_DESC),
        ("context_call", "genomic_context", "Flanking-gene context call (fungal, mixed, bacterial, short contig, no neighbours)",
         lambda v: v.replace("_", " ")),
        ("note", "note", "Free-text note", lambda v: NOTE_FIX(clean(v))),
        ("species_source", "species_name_source", "Where the species name comes from: NCBI binomial or mash cluster",
         lambda v: {"binomial": "NCBI binomial", "mash_cluster": "mash cluster"}.get(v, v)),
    ])
S1 = TABLES[-1][4]

# ---------------------------------------------------------------- S2
seeds = rd(A + "/refset/screen/seeds_ext/seeds_ext.tsv")
GRP = {"CUS/cus": "Cus (cus-like)", "CUS/mis": "Cus (mis lineage)", "SUS/sus": "Sus (sus-like)",
       "ACS/acs": "GDPD (T-DNA-type acs)", "ACS/kerr14": "GDPD (pKerr14a-type acs-like)"}
GENEFAM = {"cus": "Cus", "mis": "Cus", "sus": "Sus", "acs": "GDPD", "acs_pKerr14a": "GDPD"}
for r in seeds:
    o = r["origin"]
    m = re.match(r"T1 (\w+/\w+): T1\|(.*)$", o)
    if m:
        r["grp"] = GRP.get(m.group(1), m.group(1)); d = m.group(2)
        d = re.sub(r"^(Plant|Agro)_", "", d)
        r["desc"] = d.replace("_", " ")
        r["kind"] = "curated panel"
    else:
        r["grp"] = ""
        m2 = re.match(r"seeds11 seed: (.*)$", o)
        m3 = re.match(r"tdna cluster \S+: (.*)$", o)
        r["kind"] = "reference seed" if m2 else "T-DNA cluster member"
        r["desc"] = (m2 or m3).group(1) if (m2 or m3) else o
        r["desc"] = clean(r["desc"])
add("S2", "Seed proteins",
    "Query (seed) proteins used for the protein-level screen of fungal proteomes: bacterial Ti/Ri-plasmid proteins "
    "and plant cellular T-DNA (cT-DNA) proteins.",
    seeds, [
        ("gene", "gene", "Gene/protein class the seed represents"),
        ("seed_id", "seed_id", "Accession or identifier of the seed protein"),
        ("source", "source", "Seed source: bacterial or plant cT-DNA", lambda v: v.replace("plant_cTDNA", "plant cT-DNA")),
        ("grp", "seed_group", "Family/lineage group of the seed (blank for reference seeds of other genes)"),
        ("kind", "seed_set", "Seed set: curated panel, reference seed, or member of a T-DNA gene cluster"),
        ("desc", "description", "Organism / product / original record description"),
        ("length", "length_aa", "Seed length, aa"),
    ])

# ---------------------------------------------------------------- S3
RC = A + "/refset/screen/seeds_ext/reciprocal"
t1 = [r for r in rd(RC + "/table1_stats_dedup.tsv") if r["gene"].lower() != "orf358"]
# REMOVED 2026-10-13 (audit): Table S3a replaced by median/maximum identity columns in Table 1 (analysis/funnel/build_funnel.py)
# add("S3a", "Protein-level screen of fungal proteomes, by gene",
#     "Summary of the protein-level screen (seed proteins vs annotated fungal proteomes, GenBank/RefSeq duplicates of the same assembly removed): "
#     "number of proteomes with a hit and number of hits whose nearest relatives are non-fungal.",
#     t1, [
#         ("gene", "gene", "Seed gene/protein class"),
#         ("n_genomes_with_hit", "proteomes_with_hit", "Proteomes with at least one hit"),
#         ("median_pident", "median_identity_pct", "Median identity of the best hit to the seed, %"),
#         ("max_pident", "max_identity_pct", "Maximum identity of the best hit to the seed, %"),
#         ("n_tested", "proteins_tested", "Hits tested by reciprocal search"),
#         ("n_nonfungal", "with_nonfungal_nearest_relative", "Hits whose nearest relatives are non-fungal (index h > 0)"),
#         ("n_contam", "probable_contamination", "Of these, probable contamination (best bacterial hit identity >= 90%)"),
#         ("n_candidates", "candidates", "Hits remaining after excluding probable contamination"),
#         ("candidate_genera", "candidate_genera", "Genera of the candidates (number of hits)"),
#     ])
best = rd(A + "/refset/screen/seeds_ext/seeds_ext_best.tsv")
keep0 = {(r["gca"], r["protein"], r["gene"]) for r in best if r["dup_dropped"] == "0"}
drop1 = {(r["gca"], r["protein"], r["gene"]) for r in best if r["dup_dropped"] == "1"}
rec = rd(RC + "/reciprocal.tsv")
# identity to the nearest bacterial seed of the same gene (T-DNA candidates only; analysis/ctdna_plants_identity, journal 2026-10-12 tableS18)
pw_bact = {r["query"]: r["pident"] for r in rd(A + "/ctdna_plants_identity/pairwise_hits.tsv")
           if r["group"] == "fungal_ruleB" and r["status"] == "include" and r["pident"] != ""}
# fate of every query of the reverse search (analysis/funnel/build_funnel.py; 1,277 queries, 116 of them with non-fungal nearest relatives, as in Results 3.1 and Table 1)
def FATE5(f):
    if f.startswith("fungal nearest"): return "fungal nearest relatives (no further analysis)"
    if f.startswith("probable contamination"): return "probable contamination"
    if f.startswith("marginal"): return "member of a fungal family with a single better bacterial hit"
    if f.startswith("bacterial-type protein"): return "bacterial-type protein unrelated to rhizobia, not analysed further"
    na = "; not assigned to any lineage" if "not assigned" in f else ""
    if f.startswith("T-DNA-type candidate"): return "T-DNA candidate" + (" (nos homolog" + na + ")" if na else "")
    if f.startswith("other-rhizobial-locus candidate"): return "non-T-DNA (rhizobial) candidate" + (" (GDPD homolog" + na + ")" if na else "")
    return f
fate_of = {(r["gca"], r["protein"], r["gene"]): r["fate"] for r in rd(G + "/paper/tables/TableS3b_fate.tsv")}
sel, n_dup = [], 0
for r in rec:
    if r["gene"].lower() == "orf358":
        continue
    h = num(r["o_h"])
    k = (r["gca"], r["protein"], r["gene"])
    r["dup"] = "yes" if (k in drop1 and k not in keep0) else ""   # kept and flagged (journal 2026-10-12 suppl_sync_S1_S3b); previously excluded (112 rows)
    n_dup += r["dup"] == "yes"
    r["flag"] = "probable contamination" if (num(r["o_bact_pid"]) or 0) >= 90 else ""
    r["pid_bact_seed"] = pw_bact.get(r["protein"], "")
    r["fate"] = FATE5(fate_of.get(k, ""))
    r["nfn"] = "yes" if (h is not None and h > 0) else "no"
    sel.append(r)
assert len(sel) == 1277 and all(r["fate"] for r in sel), (len(sel), [r["protein"] for r in sel if not r["fate"]])
assert sum(r["nfn"] == "yes" for r in sel) == 116 and len({r["protein"] for r in sel}) == 1190
sel.sort(key=lambda r: (r["nfn"] != "yes", r["gene"], -(num(r["o_h"]) or 0), r["organism"]))
add("S3b", "Reverse search of the 1,277 query proteins against NCBI nr and the fate of each",
    "All queries of the reverse search of the protein-level screen: 1,277 queries (protein x seed-gene pairs) for 1,190 distinct proteins, because 83 proteins were detected by two or three seed genes. Each query was searched against NCBI nr and the best hits outside the protein's own order were compared (columns best_bacterial, best_plant and best_fungal_hit_outside_own_order). "
    "Index h = best bacterial or plant bit score minus best outside-order fungal bit score, as computed in the screen; the 116 queries with non-fungal closest relatives (h > 0; Results 3.1; column nonfungal_nearest = yes) are listed first, "
    "followed by the 1,161 queries whose closest relatives are fungal (h <= 0; nonfungal_nearest = no); each group is sorted by gene. The fate of each query in the funnel (column fate) is one of the five categories of Table 1 "
    "(probable contamination; member of a fungal family with a single better bacterial hit; bacterial-type protein unrelated to rhizobia, not analysed further; T-DNA candidate; non-T-DNA (rhizobial) candidate) "
    "or, for h <= 0, 'fungal nearest relatives (no further analysis)'. Queries from the GenBank/RefSeq duplicate of an assembly are kept and flagged (column duplicate_assembly).",
    sel, [
        ("gene", "gene", "Seed gene/protein class"), ("organism", "organism", "Organism"),
        ("gca", "assembly", "Genome assembly accession"), ("protein", "protein", "Protein accession"),
        ("pident_seed", "identity_to_seed_pct", "Identity to the nearest seed of the gene, %; the seeds include plant cT-DNA proteins, so this is not identity to the bacterial genes"),
        ("pid_bact_seed", "identity_to_bacterial_seed_pct", "Identity to the nearest bacterial T-DNA seed of the same gene (30 bacterial seeds; DIAMOND blastp), %; given only for the T-DNA candidates analysed in Tables S17a,b, blank otherwise"),
        ("cov_seed", "seed_coverage", "Fraction of the seed covered"),
        ("o_bact_org", "best_bacterial_hit_organism", "Best bacterial hit: organism"),
        ("o_bact_pid", "best_bacterial_hit_identity_pct", "Best bacterial hit: identity, %"),
        ("o_bact_bits", "best_bacterial_hit_bit_score", "Best bacterial hit: bit score"),
        ("o_plant_org", "best_plant_hit_organism", "Best plant hit: organism"),
        ("o_plant_pid", "best_plant_hit_identity_pct", "Best plant hit: identity, %"),
        ("o_plant_bits", "best_plant_hit_bit_score", "Best plant hit: bit score"),
        ("o_fung_org", "best_fungal_hit_outside_own_order_organism", "Best fungal hit outside the own order: organism"),
        ("o_fung_pid", "best_fungal_hit_outside_own_order_identity_pct", "Best fungal hit outside the own order: identity, %"),
        ("o_fung_bits", "best_fungal_hit_outside_own_order_bit_score", "Best fungal hit outside the own order: bit score"),
        ("o_h", "index_h", "Index h (see description)"),
        ("nfn", "nonfungal_nearest", "'yes' if index h > 0 (closest relatives outside the own order are non-fungal), otherwise 'no'"),
        ("flag", "flag", "'probable contamination' if best bacterial hit identity >= 90%"),
        ("dup", "duplicate_assembly", "'yes' if the assembly is the GenBank/RefSeq duplicate of another assembly of the screen (such proteins are counted in Table 1 but not in the proteome counts of Table 1 and Figure 1b)"),
        ("fate", "fate", "Fate of the query in the funnel (Table 1); 'fungal nearest relatives (no further analysis)' for h <= 0: probable contamination; member of a fungal family with a single better bacterial hit; bacterial-type protein unrelated to rhizobia, not analysed further; T-DNA candidate; non-T-DNA (rhizobial) candidate. Candidates not assigned to any lineage are marked and not counted among the lineages of Table 2"),
    ])

# ---------------------------------------------------------------- S4
gen = rd(A + "/boletales_mash_517/genomes.tsv", header=False, fields=["gca", "organism", "genus", "path"])
cl = rd(A + "/boletales_mash_517/clusters.tsv")
g2c = {}
for c in cl:
    for g in c["genomes"].split(";"):
        g2c[g] = c
acc1 = {r["accession"]: r for r in rd(G + "/data/fungal_genomes/boletales_all/accessions_20261006.tsv")}
accA = {r["GCA"]: r for r in rd(A + "/funnel/accessions.tsv")}
NOSL = "nos-subfamily opine dehydrogenase (Boletales)"
LINS = ["sus-like", "cus-like", "acc-type GDPD", "pKerr14a-type acs-like", NOSL]
PRI = {"intact": 3, "pseudogene": 2, "fragment": 1}
pres = collections.defaultdict(dict)
for r in S1:
    ln = r["lineage"]
    # 2026-10-13 (author decision): only loci placed on a tree are counted, for the nos-subfamily lineage as for the others
    if ln not in LINS:
        continue
    st = {"remnant": "fragment"}.get(r["status"], r["status"])
    if st not in PRI:
        st = "pseudogene" if st == "ambiguous" else st
    cur = pres[r["assembly"]].get(ln)
    if cur is None or PRI.get(st, 0) > PRI.get(cur, 0):
        pres[r["assembly"]][ln] = st
import json
NCBI = {}
with open(OUT + "/data/boletales517_ncbi_summary.jsonl") as fh:   # NCBI Datasets, main metadata source
    for ln_ in fh:
        j = json.loads(ln_)
        ai = j.get("assembly_info", {}); st_ = j.get("assembly_stats", {})
        bsm = ai.get("biosample", {}) or j.get("biosample", {}) or {}
        NCBI[j.get("accession", "")] = {
            "date": ai.get("release_date", ""), "level": ai.get("assembly_level", ""),
            "len": st_.get("total_sequence_length", ""), "nc": st_.get("number_of_contigs", ""),
            "n50": st_.get("contig_n50", ""), "bs": bsm.get("accession", ""),
            "bp": ai.get("bioproject_accession", ""), "sub": ai.get("submitter", ""),
            "tech": ai.get("sequencing_tech", ""),
            "wgs": (j.get("wgs_info", {}) or {}).get("wgs_project_accession", "")}
rows4 = []
for g in gen:
    a = g["gca"]; nb = NCBI.get(a, {}); c = g2c.get(a, {}); r1 = acc1.get(a, {}); rA = accA.get(a, {})
    sra = ";".join(re.findall(r"[SED]RR\d+", rA.get("used_in", "")))
    r = {"gca": a, "org": g["organism"], "genus": g["genus"], "cid": c.get("cluster_id", ""),
         "cname": c.get("cluster_name", ""), "conf": c.get("name_conflict", ""),
         "date": nb.get("date") or r1.get("release_date") or rA.get("release_date", ""),
         "level": nb.get("level") or r1.get("level") or rA.get("assembly_level", ""),
         "len": nb.get("len") or r1.get("total_len") or rA.get("length_bp", ""),
         "nc": nb.get("nc") or r1.get("n_contigs", ""),
         "n50": nb.get("n50") or r1.get("contig_n50") or rA.get("contig_N50_bp", ""),
         "bs": nb.get("bs") or rA.get("BioSample", ""), "bp": nb.get("bp") or rA.get("BioProject", ""),
         "sub": nb.get("sub", ""), "tech": nb.get("tech", ""), "wgs": nb.get("wgs", ""), "sra": sra}
    for ln in LINS:
        r[ln] = pres[a].get(ln, "-")
    rows4.append(r)
add("S4", "The 517 Boletales assemblies",
    "The 517 Boletales genome assemblies analysed, with genome cluster (mash), assembly metadata and presence of the five "
    "lineages found in Boletales (best state among the loci of the lineage in this assembly: intact > pseudogene > fragment; '-' = absent; "
    "all five lineages: loci placed on a tree). "
    "Assembly metadata (release date, level, length, contigs, N50, BioSample, BioProject, submitter, sequencing technology, WGS project) "
    "were retrieved from NCBI Datasets for all accessions (blank = not recorded by NCBI); SRA run only for assemblies used in read mapping.",
    rows4, [
        ("gca", "assembly", "Genome assembly accession"), ("org", "organism", "Organism name"),
        ("genus", "genus", "Genus"), ("cid", "genome_cluster_id", "Genome (mash) cluster number"),
        ("cname", "genome_cluster_name", "Accepted species name of the genome cluster"),
        ("conf", "cluster_name_conflict", "Conflicting NCBI names within the cluster, if any"),
        ("date", "release_date", "NCBI release date"), ("level", "assembly_level", "Assembly level"),
        ("len", "total_length_bp", "Total assembly length, bp"), ("nc", "n_contigs", "Number of contigs"),
        ("n50", "contig_N50_bp", "Contig N50, bp"), ("bs", "biosample", "BioSample accession (where available)"),
        ("bp", "bioproject", "BioProject accession (where available)"),
        ("sub", "submitter", "Submitter (NCBI)"), ("tech", "sequencing_technology", "Sequencing technology (NCBI, where recorded)"),
        ("wgs", "wgs_project", "WGS project accession (where recorded)"),
        ("sus-like", "sus-like", "Presence of sus-like (intact/pseudogene/fragment/-)"),
        ("cus-like", "cus-like", "Presence of cus-like"), ("acc-type GDPD", "acc-type GDPD", "Presence of acc-type GDPD"),
        ("pKerr14a-type acs-like", "pKerr14a-type acs-like", "Presence of pKerr14a-type acs-like"),
        (NOSL, "nos-subfamily ODH", "Presence of nos-subfamily opine dehydrogenase (ODH)"),
        ("sra", "sra_run", "SRA run used for read mapping (where applicable)"),
    ])

# ---------------------------------------------------------------- S5
HYP = {
 "SUS": ["fungal sus-like + T-DNA sus forced monophyletic",
         "fungal sus-like + alternative bacterial OCD-fold proteins forced monophyletic"],
 "ACS": ["pKerr14a-type acs-like + T-DNA acs forced monophyletic",
         "pKerr14a-type acs-like + pKerr14a acs references + Parasponia PON40589 forced monophyletic",
         "pKerr14a-type acs-like + acc-type GDPD forced monophyletic",
         "pKerr14a-type acs-like + native fungal GDPD forced monophyletic",
         "acc-type GDPD + acc-locus references forced monophyletic",
         "acc-type GDPD + T-DNA acs forced monophyletic"],
 "CUS": ["fungal cus-like + T-DNA cus forced monophyletic",
         "fungal cus-like + T-DNA cus + 4 bacterial references from the clade forced monophyletic",
         "fungal cus-like + actinobacterial-type forced monophyletic",
         "fungal cus-like + mis lineage forced monophyletic",
         "fungal cus-like + other bacterial homologs forced monophyletic",
         "all fungal cus-like in one clade (single acquisition)",
         "Boletales cus-like + T-DNA cus forced monophyletic, separate from other fungi",
         "Boletales cus-like + actinobacterial-type forced monophyletic",
         "Boletales cus-like + mis lineage forced monophyletic",
         "Boletales cus-like + remaining fungal cus-like (references excluded) forced monophyletic"]}
HYP["NOS"] = ["Boletales nos-subfamily genes + T-DNA nos forced monophyletic (excluding non-T-DNA rhizobial paralogs)",
              "Rhypophila nos-subfamily genes + T-DNA nos forced monophyletic (excluding non-T-DNA rhizobial paralogs)",
              "all fungal nos-subfamily genes (Boletales + Rhypophila) in one clade (single acquisition)"]
FAM["NOS"] = "opine dehydrogenase, nos subfamily (Boletales genes; Rhypophila homolog, not assigned to any lineage)"
# nos subfamily, second run with native fungal ODH constraints (2026-10-12, fusarium_odh_island/nos_check/au_native/RESULT.md; trees: ML, c1, c2, c3, c4, c5;
# c1-c3 re-tested in this set gave p = 0.0034, 0.0091, 0.32 — only ML, c4 and c5 are listed)
HYP["NOS_nat"] = ["(c1, listed above) Boletales nos-subfamily genes + T-DNA nos forced monophyletic", "(c2, listed above) Rhypophila nos-subfamily genes + T-DNA nos forced monophyletic",
                  "(c3, listed above) all fungal nos-subfamily genes in one clade",
                  "Boletales nos-subfamily genes + native fungal opine dehydrogenases (19 sequences) forced monophyletic (native origin of the Boletales genes)",
                  "Boletales + Rhypophila nos-subfamily genes + native fungal opine dehydrogenases forced monophyletic (single native origin of all fungal genes of the subfamily)"]
FAM["NOS_nat"] = FAM["NOS"]
SKIP = {"NOS_nat": {1, 2, 3}, "NOS_rh": {1, 2, 3, 4, 5}, "ACS_nat": {1, 2, 3, 4, 5, 6}}
# GDPD, second run (2026-10-12, trees_final/ACS/au_acc_native/RESULT.md): acc-type GDPD + native fungal GDPD; trees ML + the six constraints above (re-tested) + the new one
HYP["ACS_nat"] = HYP["ACS"] + ["acc-type GDPD + native fungal GDPD forced monophyletic (native origin of the acc-type genes)"]
FAM["ACS_nat"] = FAM["ACS"]
# Cus, Boletales monophyly run (2026-10-12, trees_final/CUS/au_bolmono/run2/RESULT.md): ML + bol_tdna_x, bol_act_x, bol_mis_x, bol_othercus_x (re-tested) + bol_mono_x
HYP["CUS_bm"] = HYP["CUS"][6:10] + ["all 26 Boletales cus-like sequences in one clade (single acquisition by Boletales)"]
FAM["CUS_bm"] = FAM["CUS"]
SKIP["CUS_bm"] = {1, 2, 3, 4}
# sus-like + native fungal OCD homologs (2026-10-12, trees_final/SUS/au_sus_native/RESULT.md; order hyp_order.txt: ML, H_tdna, H_alt, sus_natv_core, sus_natv_all; groups in groups.txt)
HYP["SUS_nat"] = ["(1) fungal sus-like + T-DNA sus", "(2) fungal sus-like + alternative bacterial OCD-fold proteins",
                  "fungal sus-like (52) + the native fungal OCD-fold homologs of the tree (2 native Sordariomycetes proteins and 4 fungal OCD-fold genes outside the sus-like clade) forced monophyletic",
                  "fungal sus-like + the native fungal OCD-fold homologs + the 5 fungal OCD-fold proteins of the outgroup (11 fungal sequences) forced monophyletic"]
FAM["SUS_nat"] = FAM["SUS"]
SKIP["SUS_nat"] = {1, 2}
# nos subfamily, third run: Rhypophila genes + native fungal ODH (2026-10-12, nos_check/au_rhypo/RESULT.md; trees: ML, c1-c5, c6, c7, c8; c1-c5 re-tested gave 0.0030, 0.0123, 0.34, 3.8e-44, 1.2e-07)
HYP["NOS_rh"] = ["(c1) Boletales + T-DNA nos", "(c2) Rhypophila + T-DNA nos", "(c3) Boletales + Rhypophila", "(c4) Boletales + native fungal ODH", "(c5) Boletales + Rhypophila + native fungal ODH",
                 "Rhypophila genes + the two native opine dehydrogenases of Sordariomycetes adjacent to the subfamily (Cryphonectria parasitica, Gnomoniopsis smithogilvyi) forced monophyletic",
                 "Rhypophila genes + all native opine dehydrogenases of Sordariomycetes (12 sequences) forced monophyletic",
                 "Rhypophila genes + all native fungal opine dehydrogenases (19 sequences) forced monophyletic"]
FAM["NOS_rh"] = FAM["NOS"]
# cus-like outside Boletales (16 sequences; same alignment/ML tree as CUS; constraints CUS.nb_*_x.constraint.nwk; order in CUS.lineage_nb.AU.order.txt)
HYP["CUS_nb"] = ["cus-like outside Boletales + T-DNA cus forced monophyletic, separate from Boletales cus-like",
                 "cus-like outside Boletales + T-DNA cus + 4 bacterial references from the clade forced monophyletic, separate from Boletales cus-like",
                 "cus-like outside Boletales + actinobacterial-type forced monophyletic",
                 "cus-like outside Boletales + mis lineage forced monophyletic",
                 "cus-like outside Boletales + other bacterial homologs forced monophyletic",
                 "all cus-like outside Boletales in one clade (single acquisition)"]
FAM["CUS_nb"] = FAM["CUS"]
# chs-like (Fusarium, Dactylonectria); order in analysis/chs_fusarium/au/run_au.sh: ML, chs_tdna, chs_bactchs, chs_ags, chs_nat, chs_bact, chs_mono (journal 2026-10-11 AU_chs)
HYP["CHS"] = ["fungal chs-like + T-DNA chs seed forced monophyletic, without the other Ti-plasmid chs genes",
              "fungal chs-like + the four other Ti-plasmid chs genes forced monophyletic, without the seed",
              "fungal chs-like + ags seed forced monophyletic",
              "fungal chs-like + native fungal ags homologs forced monophyletic",
              "fungal chs-like + other bacterial M24 aminopeptidases forced monophyletic",
              "all fungal chs-like (Fusarium + Dactylonectria) in one clade (single acquisition)"]
FAM["CHS"] = "chrysopine synthase (Chs; M24 aminopeptidase family)"
# chs-like, nested placements (2026-10-12, analysis/chs_fusarium/au_nested/RESULT.md; order hyp_order.txt): Ti-chs x = seed (pTi chrysopine synthase A0A2Z2PGA9),
# WP_139103025.1, WP_272871282.1 (A. tumefaciens), WP_272951408.1 (A. vitis), WP_084631803.1 (A. larrymoorei); nest_x = all 11 fungal + x; fnest_x = 10 Fusarium + x
_TI = [("the chs seed (pTi chrysopine synthase A0A2Z2PGA9)", "seed"), ("Ti-plasmid chs WP_139103025.1", "b139"),
       ("Ti-plasmid chs WP_272871282.1 (A. tumefaciens)", "b272871"), ("Ti-plasmid chs WP_272951408.1 (A. vitis)", "b272951"),
       ("Ti-plasmid chs WP_084631803.1 (A. larrymoorei)", "b084")]
HYP["CHS_nest"] = []
for _d, _k in _TI:
    HYP["CHS_nest"].append("all 11 fungal chs-like (Fusarium + Dactylonectria) + %s forced monophyletic (fungal genes nested within the Ti-plasmid chs clade next to %s)" % (_d, _k))
    HYP["CHS_nest"].append("the 10 Fusarium chs-like + %s forced monophyletic (Dactylonectria unconstrained)" % _d)
FAM["CHS_nest"] = FAM["CHS"]
# odh island of F. poae / F. langsethiae; order in analysis/fusarium_odh_island/au/au_trees.order.txt (descriptions from au_summary.tsv)
HYP["ODH"] = ["odh-island ODH + all native fungal ODH forced monophyletic",
              "odh-island ODH + native Sordariomycetes ODH (Apiospora, Astrocystis) forced monophyletic (vertical-inheritance expectation)",
              "control: native fungal ODH forced monophyletic (odh-island ODH unconstrained)",
              "odh-island ODH inside the ML T-DNA clade (ocs/nos/vis seeds + Rhypophila/Butyriboletus + A. larrymoorei WP_309770762) forced monophyletic",
              "odh-island ODH + T-DNA ocs/nos/vis seeds only forced monophyletic",
              "control: odh-island ODH + A. larrymoorei WP_309770759 forced monophyletic (clade of the ML tree)",
              "odh-island ODH + all other bacterial ODH (both A. larrymoorei ODH excluded) forced monophyletic",
              "odh-island ODH inside the nos subfamily (nos seed P00386 + Gillisella) forced monophyletic"]
FAM["ODH"] = "opine dehydrogenase, odh island (Fusarium)"
REFHYP = {"CUS_bm": "unconstrained maximum-likelihood tree (reference; separate test set for the monophyly of the Boletales cus-like genes; the four Boletales constraints above re-tested here gave p = 5.1e-07, 1.3e-56, 1.3e-46 and 3.2e-05)",
          "ACS_nat": "unconstrained maximum-likelihood tree (reference; separate test set with the native fungal GDPD constraint for the acc-type genes; constraints 1-6 above re-tested here gave p = 1.4e-31, 0.44, 0.16, 4.7e-61, 0.56 and 7.4e-11)",
          "CUS_nb": "unconstrained maximum-likelihood tree (reference; separate test set for cus-like outside Boletales)",
          "CHS_nest": "unconstrained maximum-likelihood tree (reference; separate test set: are the fungal chs-like genes nested within the Ti-plasmid chs clade?)",
          "SUS_nat": "unconstrained maximum-likelihood tree (reference; separate test set with the native fungal OCD-fold constraints for the sus-like genes; constraints 1-2 above re-tested here gave p = 0.26 and 5.9e-06)",
          "NOS_nat": "unconstrained maximum-likelihood tree (reference; separate test set with the native fungal opine dehydrogenase constraints; constraints 1-3 of the first set re-tested here gave p = 0.0034, 0.0091 and 0.32)",
          "NOS_rh": "unconstrained maximum-likelihood tree (reference; separate test set for the Rhypophila genes and the native fungal opine dehydrogenases; constraints 1-5 of the first two sets re-tested here gave p = 0.0030, 0.0123, 0.34, 3.8e-44 and 1.2e-07)"}
FILES = {"SUS": A + "/trees_final/SUS/SUS.AU.iqtree", "SUS_nat": A + "/trees_final/SUS/au_sus_native/au_sus_native.iqtree", "ACS": A + "/trees_final/ACS/ACS.lineage.AU.iqtree",
         "CUS": A + "/trees_final/CUS/CUS_AU_test_per_lineage.iqtree.txt",
         "CUS_nb": A + "/trees_final/CUS/CUS.lineage_nb.AU.iqtree",
         "CHS": A + "/chs_fusarium/au/au_test.iqtree",
         "CHS_nest": A + "/chs_fusarium/au_nested/au_nested.iqtree",
         "NOS": A + "/fusarium_odh_island/nos_check/au_test.iqtree",
         "NOS_nat": A + "/fusarium_odh_island/nos_check/au_native/au_native.iqtree",
         "NOS_rh": A + "/fusarium_odh_island/nos_check/au_rhypo/au_rhypo.iqtree",
         "ACS_nat": A + "/trees_final/ACS/au_acc_native/au_acc_native.iqtree",
         "CUS_bm": A + "/trees_final/CUS/au_bolmono/run2/CUS.lineage_bm_s2.AU.iqtree",
         "ODH": A + "/fusarium_odh_island/au/au_test.iqtree"}
rows5 = []
for f, fn in FILES.items():
    txt = open(fn).read()
    blk = txt.split("USER TREES")[1].split("deltaL  :")[0]
    trees = []
    for ln in blk.splitlines():
        m = re.match(r"\s*(\d+)\s+(-?[\d.]+)\s+([\d.e+-]+)\s+(.*)$", ln)
        if m:
            tok = m.group(4).split()
            trees.append((m.group(2), m.group(3), tok[-2]))
    assert len(trees) == len(HYP[f]) + 1, (f, len(trees))
    for i, (lg, dl, pau) in enumerate(trees):
        if i in SKIP.get(f, set()):
            continue
        if i == 0:
            hyp = REFHYP.get(f, "unconstrained maximum-likelihood tree (reference)")
        else:
            hyp = HYP[f][i - 1]
        p = float(pau)
        rows5.append({"fam": FAM[f], "n": i, "hyp": hyp, "lnl": lg, "dl": dl, "p": pau, "cid": f"{f}-{i}",
                      "v": "reference" if i == 0 else ("rejected" if p < 0.05 else "not rejected")})
add("S5", "Approximately unbiased (AU) topology tests",
    "AU tests of constrained topologies against the unconstrained ML tree, per family (IQ-TREE, 10,000 RELL resamplings). "
    "Verdict: rejected if p-AU < 0.05. The first row of each family is the unconstrained ML tree used as the reference "
    "(its deltaL is not necessarily 0 when a constrained tree has a higher likelihood). Within each family, Holm adjustment of the p-values "
    "over the listed constraints changes no verdict. Where a family has several reference rows, the rows form separate runs (separate user-tree sets), and p-values are comparable only within a run.",
    rows5, [
        ("fam", "family", "Protein family"), ("n", "tree_number", "Tree number in the IQ-TREE user-tree set (0 = unconstrained ML tree)"),
        ("hyp", "tested_topology_constraint", "Tested topology (constraint)"), ("lnl", "logL", "Log-likelihood of the tree"),
        ("dl", "deltaL", "Log-likelihood difference from the best tree in the set"), ("p", "p_AU", "p-value of the AU test"),
        ("v", "verdict", "rejected / not rejected (p-AU < 0.05)"),
        ("cid", "constraint_id", "Identifier of the constraint (run code and tree number); key to Table S11b and Data S2b")])
# 2026-10-13af: Table S11b - mapping of the AU constraints to the biological alternatives (reviewer point 3; TableS11_constraint_mapping.tsv by build_au_mapping.py)
_map = rd(OUT + "/TableS11_constraint_mapping.tsv")
assert len(_map) == 60, len(_map)
# 2026-10-13al: the wording of the interpretation of constraints with a non-monophyletic native set (and the ags seed description of CHS-3) now comes from build_au_mapping.py
assert all(r["interpretation"].startswith("not informative about a native origin") for r in _map if r["native_set_monophyletic_in_ML_tree"] == "no")
assert not any("yeast" in r["biological_alternative"] for r in _map)
_leg = open(OUT + "/TableS11_constraint_mapping_LEGEND.md", encoding="utf-8").read().strip()
_leg = re.sub(r"^\*\*Table S11b\.[^*]*\*\*\s*", "", _leg).replace("\n", " ")
_MAPDESC = {"constraint_id": "Identifier of the constraint (run code and tree number), as in Table S11 and Data S2b", "family": "Protein family", "tree_number": "Tree number in the IQ-TREE user-tree set",
            "tested_topology_constraint": "Tested topology (constraint), as in Table S11", "p_AU": "p-value of the AU test", "deltaL": "Log-likelihood difference from the best tree", "verdict": "rejected / not rejected (p-AU < 0.05)",
            "biological_alternative": "The biological alternative the constraint tests", "constraint_definition": "Leaf sets forced into one clade (all other leaves free)", "constraint_leaves": "Number of leaves in the forced clade",
            "candidate_leaves": "Number of fungal genes under test in the forced clade", "partner_leaves": "Number of partner leaves in the forced clade", "partner_set_monophyletic_in_ML_tree": "Whether the partner set is monophyletic in the unconstrained ML tree",
            "native_fungal_set_in_constraint": "Whether the constraint includes a native fungal set", "native_set_monophyletic_in_ML_tree": "Whether that native set is monophyletic in the unconstrained ML tree",
            "interpretation": "What rejection or non-rejection of the constraint allows one to conclude", "caveat": "Caveat on the interpretation (e.g. the native set is not itself monophyletic)", "constraint_file": "Constraint tree file in Data S2b"}
add("S5b", "Mapping of the AU constraints to the biological alternatives they test", _leg, _map, [(c, c, _MAPDESC.get(c, c.replace("_", " "))) for c in _map[0].keys()])

# ---------------------------------------------------------------- S6
cov = [r for r in rd(A + "/read_coverage/locus_coverage.tsv") if r["fam"] in FAM]
ln_of = {r["locus_id"]: r["lineage"] for r in S1}
for r in cov:
    r["lineage"] = ln_of.get(r["locus"], "")
add("S6", "Read depth at the loci",
    "Reads mapped back to the assemblies show whether each locus lies in the assembled genome: read depth at each GDPD/Cus/Sus locus relative to "
    "the genome and to its contig, and properties of the contig. Lineage from Table S1b.",
    cov, [
        ("gca", "assembly", "Genome assembly accession"), ("species", "species", "Species (genome label)"),
        ("locus", "locus_id", "Locus identifier (Table S1b)"), ("fam", "family", "Protein family", fam_map),
        ("lineage", "lineage", "Lineage (Table S1b)"),
        ("status", "status", "Locus state", status_map), ("contig", "contig", "Contig"), ("contig_len", "contig_length_bp", "Contig length, bp"),
        ("genome_median_depth", "genome_median_depth", "Median read depth of the genome"),
        ("contig_mean_depth", "contig_mean_depth", "Mean read depth of the contig"),
        ("contig_ratio", "contig_depth_ratio", "Contig depth / genome median depth"),
        ("locus_mean_depth", "locus_mean_depth", "Mean read depth over the locus"),
        ("locus_min_depth", "locus_min_depth", "Minimum read depth over the locus"),
        ("locus_ratio", "locus_depth_ratio", "Locus depth / genome median depth"),
        ("span_left", "reads_spanning_left_flank", "Read pairs spanning the left boundary of the locus"),
        ("span_right", "reads_spanning_right_flank", "Read pairs spanning the right boundary of the locus"),
        ("pp_frac_contig", "proper_pair_fraction_contig", "Fraction of properly paired reads on the contig"),
        ("pp_frac_genome", "proper_pair_fraction_genome", "Fraction of properly paired reads in the genome"),
        ("gc_contig", "gc_contig", "GC content of the contig"), ("gc_assembly", "gc_assembly", "GC content of the assembly"),
        ("mapq_contig", "mean_mapq_contig", "Mean mapping quality on the contig")])

# ---------------------------------------------------------------- S7
rep = [r for r in rd(A + "/indep_assemblies/contig_reproducibility.tsv") if "ORF358" not in r["loci"]]
add("S7", "Reproducibility of loci in independent assemblies",
    "The loci are reproduced in independent assemblies of the same species: for loci-bearing contigs of one assembly, the best matching sequence in an independent assembly of the same species and "
    "how much of the contig, the locus and its flanks is covered.",
    rep, [
        ("species", "species", "Species"), ("src_gca", "source_assembly", "Assembly containing the locus"),
        ("contig", "contig", "Contig in the source assembly"), ("contig_len", "contig_length_bp", "Contig length, bp"),
        ("loci", "loci", "Locus identifiers on the contig"), ("tgt_gca", "independent_assembly", "Independent assembly compared"),
        ("n_tgt_seqs", "n_matching_sequences", "Number of sequences in the independent assembly aligning to the contig"),
        ("best_tgt", "best_matching_sequence", "Best-matching sequence"), ("aligned_frac", "aligned_fraction", "Fraction of the contig aligned"),
        ("identity", "identity", "Identity of the alignment"), ("locus_covered_frac", "locus_covered_fraction", "Fraction of the locus covered"),
        ("left_flank_covered", "left_flank_covered_fraction", "Fraction of the left flank covered"),
        ("right_flank_covered", "right_flank_covered_fraction", "Fraction of the right flank covered")])

# ---------------------------------------------------------------- S8
isl = rd(A + "/island_order/island_gene_order.tsv")
for r in isl:
    k = r["kind"]
    if k == "neighbor":
        r["cls"], r["fam"], r["lin"], r["st"] = "neighboring gene", "", "", ""
    else:
        _, f, sg, st = k.split(":")
        r["cls"], r["fam"], r["lin"], r["st"] = "focal locus", FAM.get(f, f), SUBG.get(sg, sg), STATUS.get(st, st)
isl = [r for r in isl if "ORF358" not in r["kind"]]
for r in isl:   # gene identifiers g<n> are unique only within a contig: prefix with the contig
    if re.match(r"g\d+$", r["gene"]): r["gene"] = r["island_contig"] + ":" + r["gene"]
    r["note"] = MIS_NOTE if r["gene"] in ("T_CUS_JBBPKX010000014_S2L2", "JBBPKX010000014:g4") else ""
add("S8", "Gene content of the islands",
    "Genes on the contigs carrying the loci (focal loci and their annotated neighbours), with the best match to plasmid (T-DNA/Ti/Ri) proteins where found. "
    "Neighbouring genes are named <contig>:g<n> (the number is unique only within a contig).",
    isl, [
        ("island_contig", "contig", "Contig"), ("gca", "assembly", "Genome assembly accession"), ("gene", "gene", "Gene identifier (neighbouring genes: contig and gene number; focal loci: locus identifier)"),
        ("start", "start", "Start"), ("end", "end", "End"), ("strand", "strand", "Strand"), ("exons", "n_exons", "Number of exons"),
        ("cls", "gene_class", "Focal locus or neighbouring gene"), ("fam", "locus_family", "Family of the focal locus"),
        ("lin", "locus_lineage_or_profile", "Lineage/HMM subgroup of the focal locus"), ("st", "locus_status", "State of the focal locus"),
        ("best_plasmid_protein", "best_plasmid_protein", "Best-matching plasmid protein"), ("plasmid", "plasmid_record", "Plasmid/genome record of the match"),
        ("plasmid_gene", "plasmid_gene", "Gene name of the match"), ("product", "plasmid_product", "Product of the match"),
        ("in_TDNA", "match_in_T-DNA", "Whether the matched plasmid gene lies in a T-DNA (yes/no)"),
        ("plasmid_start", "plasmid_match_start", "Start of the match on the plasmid record"), ("plasmid_end", "plasmid_match_end", "End of the match on the plasmid record"),
        ("pident", "identity_pct", "Identity of the match, %"), ("qcov", "query_coverage", "Query coverage of the match"),
        ("evalue", "evalue", "E-value of the match"), ("note", "note", "Comment on the row")])

# ---------------------------------------------------------------- S9
dn = rd(A + "/dnds/dnds_summary.tsv")
add("S9a", "Selection analysis (dN/dS)",
    "Codon-model (PAML codeml/yn00) tests of selection on sequence sets of each family.",
    dn, [
        ("family", "family", "Protein family", fam_map), ("n_seqs", "n_sequences", "Number of sequences"),
        ("codon_sites", "codon_sites", "Codon sites in the alignment"), ("M0_omega", "omega_M0", "dN/dS (omega) under the one-ratio model M0"),
        ("M7_lnL", "lnL_M7", "Log-likelihood, model M7 (beta)"), ("M8_lnL", "lnL_M8", "Log-likelihood, model M8 (beta and omega)"),
        ("LRT_M7vsM8", "LRT_2dlnL_M7_vs_M8", "Likelihood-ratio statistic M7 vs M8"), ("p_M7vsM8", "p_M7_vs_M8", "p-value, M7 vs M8"),
        ("M8a_lnL", "lnL_M8a", "Log-likelihood, model M8a"), ("LRT_M8avsM8", "LRT_2dlnL_M8a_vs_M8", "Likelihood-ratio statistic M8a vs M8"),
        ("p_M8avsM8", "p_M8a_vs_M8", "p-value, M8a vs M8"),
        ("M8_positive_sites_BEB>0.95", "positively_selected_sites_BEB_gt_0.95", "Sites under positive selection (BEB posterior > 0.95) under M8"),
        ("yn00_median_dNdS", "pairwise_median_dNdS_yn00", "Median pairwise dN/dS (yn00)"), ("yn00_q90_dNdS", "pairwise_q90_dNdS_yn00", "90th percentile of pairwise dN/dS (yn00)"),
        ("yn00_pairs", "pairwise_comparisons_yn00", "Number of pairwise comparisons (yn00)")])
# Table S9b (conservation of reference-invariant positions) removed 2026-10-13q (restructuring plan).

# ---------------------------------------------------------------- S10
pl = rd(A + "/composition/per_locus.tsv")
for r in pl:
    r["lineage"] = norm_lineage(r["lineage"])
    r["cds_source"] = "extra genome-derived single-exon model" if r["cds_source"].startswith("extra2") else "gene model of the catalogue"
    r["bg_source"] = {"augustus_neighbors": "Augustus models of neighbouring genes", "ncbi_annotation_CDS": "NCBI annotation CDS",
                      "annotation_GFF_extra2": "annotation (GFF)"}.get(r["bg_source"], r["bg_source"])
# 2026-10-13o (audit A5): rows for the nos-subfamily ODH of Boletales (composition_nosODH.py: same columns as above), for the intact chs-like loci and for the intact
# odh-island genes of Fusarium (GC of the gene model and of the assembly from the scan tables; GC3 and background not computed for these -> blank)
NOSL = "nos-subfamily opine dehydrogenase (Boletales)"
for r in rd(A + "/composition/per_locus_nosODH.tsv"):
    r["lineage"] = NOSL
    r["cds_source"] = "miniprot gene model"
    r["bg_source"] = {"augustus_neighbors": "Augustus models of neighbouring genes"}.get(r["bg_source"], r["bg_source"])
    pl.append(r)
for src, lin_of in ((A + "/chs_fusarium/stage2/loci.tsv", lambda r: "chs-like"),
                    (A + "/fusarium_odh_island/scan/loci.tsv", lambda r: "odh island: " + {"odh": "odh", "pepP": "pepP", "g3": "orf3"}[r["gene"]])):
    for r in rd(src):
        if r["status"] != "intact": continue
        pl.append({"locus_id": r["locus_id"], "family": "M24 aminopeptidase (Chs/Ags)" if "chs" in src else "opine dehydrogenase / M24 aminopeptidase / hypothetical (odh island)",
                   "lineage": lin_of(r), "organism": r["organism"], "order": "Hypocreales", "gca": r["gca"],
                   "cds_len_nt": "", "cds_source": "miniprot gene model (targeted scan)", "gc": "%.4f" % (float(r["gene_gc"]) / 100),
                   "assembly_gc": "%.4f" % (float(r["assembly_gc"]) / 100), "bg_source": "not computed"})
for r in pl:
    r.setdefault("assembly_gc", "")
spec10 = [("locus_id", "locus_id", "Locus identifier"), ("family", "family", "Protein family", fam_map), ("lineage", "lineage", "Lineage"),
    ("organism", "organism", "Organism"), ("order", "order", "Taxonomic order"), ("gca", "assembly", "Genome assembly accession"),
    ("cds_len_nt", "cds_length_nt", "CDS length, nt"), ("cds_source", "cds_source", "Origin of the CDS sequence"),
    ("gc", "gc", "GC content of the CDS"), ("gc3", "gc3", "GC content at third codon positions"),
    ("assembly_gc", "assembly_gc", "GC content of the whole assembly (given for the chs-like and odh-island rows only; the background columns are empty for these rows)"),
    ("bg_source", "background_source", "Source of the background gene set (genes of the same genome)"),
    ("n_bg_genes", "n_background_genes", "Number of background genes"), ("bg_small_flag", "background_flag", "Background size flag"),
    ("bg_med_gc", "background_median_gc", "Background median GC"), ("bg_q1_gc", "background_q1_gc", "Background GC, first quartile"),
    ("bg_q3_gc", "background_q3_gc", "Background GC, third quartile"), ("bg_med_gc3", "background_median_gc3", "Background median GC3"),
    ("bg_q1_gc3", "background_q1_gc3", "Background GC3, first quartile"), ("bg_q3_gc3", "background_q3_gc3", "Background GC3, third quartile"),
    ("delta_gc3", "delta_gc3", "GC3 of the locus minus background median GC3"), ("gc3_pctile", "gc3_percentile", "Percentile of GC3 among background genes"),
    ("gc3_z", "gc3_z_score", "GC3 z-score against the background"), ("gc_pctile", "gc_percentile", "Percentile of GC among background genes"),
    ("codon_dist_1minusr_rscu", "codon_usage_distance", "Codon-usage distance (1 - r of RSCU) to the genome"),
    ("codon_dist_pctile", "codon_usage_distance_percentile", "Percentile of the codon-usage distance among length-matched background genes"),
    ("n_len_matched_bg", "n_length_matched_background", "Number of length-matched background genes"),
    ("bolpool_n_genes", "boletales_pool_n_genes", "Genes in the pooled Boletales reference set"),
    ("bolpool_med_gc3", "boletales_pool_median_gc3", "Median GC3 of the pooled Boletales set"),
    ("bolpool_delta_gc3", "boletales_pool_delta_gc3", "GC3 of the locus minus pooled Boletales median"),
    ("bolpool_gc3_pctile", "boletales_pool_gc3_percentile", "Percentile of GC3 in the pooled Boletales set"),
    ("bolpool_gc3_z", "boletales_pool_gc3_z_score", "GC3 z-score against the pooled Boletales set"),
    ("bolpool_codon_dist_pctile", "boletales_pool_codon_usage_distance_percentile", "Percentile of the codon-usage distance in the pooled set")]
add("S10a", "Nucleotide composition of the loci",
    "GC and GC3 of intact fungal locus CDSs compared with background genes of the same genome and with a pooled Boletales gene set (sus-like, cus-like, acc-type GDPD, pKerr14a-type and nos-subfamily loci of Boletales; two nos-subfamily pseudogenes included). For the intact chs-like and odh-island genes of Fusarium only the GC of the gene and of the assembly is given.",
    pl, spec10)
rb = rd(A + "/composition/reference_bacterial.tsv")
add("S10b", "Nucleotide composition of bacterial reference genes",
    "GC and GC3 of the bacterial plasmid reference CDSs.", rb, [
        ("ref_id", "reference_id", "Protein accession"), ("description", "description", "Gene and plasmid", clean),
        ("source_file", "source_record", "GenBank record of the CDS"), ("cds_len_nt", "cds_length_nt", "CDS length, nt"),
        ("gc", "gc", "GC content of the CDS"), ("gc3", "gc3", "GC content at third codon positions"), ("note", "note", "Product")])

# Table S11 (Dollo parsimony) removed 2026-10-13; later tables renumbered down by one (2026-10-13m).

# ---------------------------------------------------------------- S11
S2D = A + "/refset/screen/step2/boletales517_tdna_summary"
cg = rd(S2D + "/classic_genes.tsv")
cg = [r for r in cg if not r["gene"].startswith("rolC?") and r["gene"] != "orf13/orf14"]   # rolC? repeats the rolB row (same query); orf13/orf14 had no protein query
for r in cg:
    r["note"] = ""
    if r["gene"].startswith("rolB ("): r["gene"] = "rolB/rolC (AAK90962; one query serves both genes)"
    if r["gene"].startswith("acs family queries"):
        r["gene"] = "GDPD family (acs-family queries)"
        r["queries"] = "19 GDPD-family seed proteins of Table S2"
        r["note"] = "The hits are the acc-type GDPD and pKerr14a-type acs-like loci of Table 2, not canonical T-DNA acs"
    elif num(r["max_pident_best_locus"]) is not None and num(r["max_pident_best_locus"]) >= 45 and r["n_asm_strong_locus_best_subject"] == "0":
        r["note"] = "Best single-locus hit above 45% identity but covering less than half of the query: not counted as a strong hit"
add("S11a", "Oncogenes, hormone and opine synthesis genes of the T-DNA in Boletales assemblies",
    "Search of the 517 Boletales assemblies for the oncogenes, hormone synthesis genes and opine synthesis genes of the T-DNA (protein queries vs assemblies, blastx). "
    "A strong single-locus hit has an identity of at least 45% over at least half of the query. "
    "Plasmid maintenance and transfer genes (vir, tra, rep) are in Table S5b. The genes orf13/orf14 were not searched (no protein query).",
    cg, [
        ("gene", "gene", "Gene (query group) with query accession in brackets", clean), ("n_queries", "n_queries", "Number of query proteins"),
        ("queries", "queries", "Query protein accessions", clean),
        ("n_asm_strong_locus_best_subject", "assemblies_with_strong_locus_hit", "Assemblies with a strong single-locus hit"),
        ("max_pident_best_locus", "max_identity_pct_locus_hit", "Maximum identity of the best single-locus hit, %"),
        ("n_asm_chain_pid45_cov05", "assemblies_with_chained_hit", "Assemblies with a chained hit (identity >= 45%, query coverage >= 0.5)"),
        ("max_chain_pident_cov05", "max_identity_pct_chained_hit", "Maximum identity of chained hits with coverage >= 0.5, %"),
        ("note", "note", "Comment on the row")])
od = rd(S2D + "/opine_dehydrogenase_loci.tsv")
# lineage and gene-model status taken from the unified catalogue (paper/tables/build_catalog.py, incl. its NOS_OVERRIDE for ODH_JBBPWO010000281.1), so that S1 and S11b agree
cat_odh = {r["contig"]: r for r in rd(G + "/paper/tables/TableS1_all_lineages.tsv") if r["source_table"] == "Table S11b"}
for r in od:
    c = cat_odh[r["contig"]]
    r["lineage"], r["status"] = c["lineage"], c["status"]
    r["catalog_note"] = "" if c["note"].startswith("best plasmid protein") else c["note"]
    if "ODH_" + r["contig"] in NOTE_ADD: r["catalog_note"] = add_note(r["catalog_note"], NOTE_ADD["ODH_" + r["contig"]])
add("S11b", "Opine dehydrogenase-like loci in Boletales assemblies",
    "Loci in Boletales assemblies matching bacterial opine (nopaline) dehydrogenase queries.", od, [
        ("gca", "assembly", "Genome assembly accession"), ("organism", "organism", "Organism"), ("contig", "contig", "Contig"),
        ("contig_len", "contig_length_bp", "Contig length, bp"), ("start", "start", "Start"), ("end", "end", "End"), ("strand", "strand", "Strand"),
        ("subject", "query_protein", "Query protein (accession|source|record|organism|gene|product)", clean),
        ("pident", "identity_pct", "Identity, %"), ("subj_cov", "query_coverage", "Query coverage"),
        ("genome_is_carrier", "assembly_carries_lineage", "1 if the assembly carries a sus-like, cus-like, acc-type GDPD or pKerr14a-type acs-like locus"),
        ("contig_has_carrier_locus", "contig_carries_lineage_locus", "1 if the same contig carries such a locus"),
        ("carrier_lineages_on_contig", "lineages_on_contig", "Lineages of loci on the same contig", clean),
        ("lineage", "lineage", "Lineage in the catalogue (Table S1b); '(profile only)' = locus not placed on the nos tree"),
        ("status", "status", "Gene-model state (as in Table S1b)"), ("catalog_note", "note", "Reason for a manual reclassification, or other comment on the locus")])

# ---------------------------------------------------------------- S12
CH = A + "/chs_fusarium"
lc = rd(CH + "/stage2/loci.tsv")
add("S12a", "chs-like loci in Fusarium and relatives",
    "The 33 loci similar to the bacterial chs gene found in the 324 genomes of Fusarium and related genera (all are chs-like loci), with alignment to the reference proteins and gene-model state; one row per locus (per assembly in Table S7).",
    lc, [
        ("gca", "assembly", "Genome assembly accession"), ("organism", "organism", "Organism"), ("species", "species", "Species"),
        ("complex", "species_complex", "Species complex"), ("level", "assembly_level", "Assembly level"), ("locus_id", "locus_id", "Locus identifier"),
        ("contig", "contig", "Contig"), ("start", "start", "Start"), ("end", "end", "End"), ("strand", "strand", "Strand"),
        ("contig_len", "contig_length_bp", "Contig length, bp"), ("dist_left", "distance_to_contig_start_bp", "Distance to the contig start, bp"),
        ("dist_right", "distance_to_contig_end_bp", "Distance to the contig end, bp"),
        ("id_chs", "identity_to_chs", "Identity to the Ti-plasmid chs seed (A0A2Z2PGA9)"), ("cov_chs", "coverage_of_chs", "Coverage of the chs seed"),
        ("id_EEU", "identity_to_EEU41297.1", "Identity to the fungal reference EEU41297.1"), ("cov_EEU", "coverage_of_EEU41297.1", "Coverage of EEU41297.1"),
        ("id_KNB", "identity_to_KNB14377.1", "Identity to the reference protein KNB14377.1"),
        ("cov_KNB", "coverage_of_KNB14377.1", "Coverage of KNB14377.1"),
        ("status", "status", "intact / truncated / pseudogene"),
        ("primary_seed", "primary_seed", "Seed protein used for the gene model"), ("prim_qfrom", "seed_aln_start_aa", "Alignment start on the seed, aa"),
        ("prim_qto", "seed_aln_end_aa", "Alignment end on the seed, aa"), ("prim_len", "seed_length_aa", "Seed length, aa"),
        ("exons", "n_exons", "Exons in the model"), ("frameshifts", "frameshifts", "Frameshifts relative to the seed"),
        ("stops", "in_frame_stops", "In-frame stop codons"), ("missing", "missing_seed_segments", "Missing N-/C-terminal segments of the seed"),
        ("edge", "contig_edge", "Locus truncated by a contig edge"), ("gene_gc", "locus_gc_pct", "GC content of the locus, %"),
        ("assembly_gc", "assembly_gc_pct", "GC content of the assembly, %"),
        ("stop_pos_in_primary_seed", "stop_positions_in_seed_aa", "Positions of in-frame stops on the seed, aa"), ("note", "note", "Free-text note")])
# Table S12b (chs-like loci per species) replaced by the new Table S7 (324 assemblies), 2026-10-13q
runs = {r["run"]: r for r in rd(CH + "/expression/runs.tsv")}
smm = rd(CH + "/expression/summary.tsv")
for r in smm:
    r.update({k: v for k, v in runs.get(r["run"], {}).items() if k not in r})
add("S12c", "RNA-seq of Fusarium: expression of the chs-like locus",
    "RNA-seq runs analysed (library metadata) and expression of the chs-like locus relative to reference genes.", smm, [
        ("run", "sra_run", "SRA run accession"), ("species", "species_code", "Species code in the source analysis"),
        ("BioProject", "bioproject", "BioProject"), ("BioSample", "biosample", "BioSample"),
        ("organism_in_SRA", "organism", "Organism as recorded in SRA"), ("strain", "strain", "Strain"), ("condition", "condition", "Experimental condition"),
        ("platform", "platform", "Sequencing platform"), ("layout", "layout", "Library layout"), ("pairs_in_run", "read_pairs_in_run", "Read pairs in the run"),
        ("n_pairs_used", "read_pairs_used", "Read pairs used (subsample)"), ("mapped_pct", "mapped_pct", "Mapped reads, %"),
        ("reads_on_target", "reads_on_locus", "Reads on the locus"), ("target_TPM", "locus_TPM", "TPM of the locus"),
        ("target_reads_per_kb", "locus_reads_per_kb", "Reads per kb on the locus"),
        ("pct_all_genes_strictly_below", "pct_genes_with_lower_expression", "Percent of all genes with expression strictly below the locus"),
        ("target_TPM_over_median_expressed", "locus_TPM_over_median_expressed", "Locus TPM / median TPM of genes with reads"),
        ("actin_TPM", "actin_TPM", "TPM of actin"), ("EF1a_TPM", "EF1a_TPM", "TPM of EF1a"), ("GAPDH_TPM", "GAPDH_TPM", "TPM of GAPDH"),
        ("tubulinB_TPM", "tubulin_beta_TPM", "TPM of beta-tubulin"), ("lib_sense_frac_R1", "sense_fraction_R1", "Fraction of R1 reads in sense orientation"),
        ("target_breadth", "locus_breadth", "Fraction of the locus covered by reads")])

# ---------------------------------------------------------------- S13
ex = [r for r in rd(A + "/expression/expression_summary.tsv") if r["fam"] in FAM and "ORF358" not in r["locus"]]
for r in ex:
    r["lineage"] = ln_of.get(r["locus"], "")
add("S13", "RNA-seq of Boletaceae: expression of the loci",
    "Reads from public RNA-seq runs mapped to the loci of the corresponding assemblies; expression of each locus relative to the genes of the genome.",
    ex, [
        ("run", "sra_run", "SRA run accession"), ("gca", "assembly", "Genome assembly accession"), ("locus", "locus_id", "Locus identifier (Table S1b)"),
        ("fam", "family", "Protein family", fam_map), ("lineage", "lineage", "Lineage (Table S1b)"),
        ("status", "status", "Locus state", status_map), ("contig", "contig", "Contig"), ("start", "start", "Start"), ("end", "end", "End"),
        ("len_bp", "length_bp", "Locus length, bp"), ("reads", "reads", "Reads mapped to the locus"), ("reads_per_kb", "reads_per_kb", "Reads per kb"),
        ("genome_reads_per_kb_median_genes", "median_reads_per_kb_of_genes", "Median reads per kb of the genes of the genome"),
        ("percentile_among_genes", "percentile_among_genes", "Percentile of the locus among the genes of the genome")])

# ---------------------------------------------------------------- S14 (odh-like island of F. poae / F. langsethiae)
OI = A + "/fusarium_odh_island"
GN = {"odh": "odh (opine dehydrogenase family)", "pepP": "pepP (Xaa-Pro aminopeptidase, M24 family)", "g3": "orf3 (hypothetical protein)"}
gn_map = lambda v: GN.get(v, v)
g3_disp = lambda v: re.sub(r"^g3\b", "orf3", v)   # display only: g3 -> orf3 (source files unchanged)
lo = rd(OI + "/scan/loci.tsv")
add("S14a", "Loci of the odh-like island in Fusarium and relatives",
    "Loci similar to the three genes of the bacterial-type island (odh, pepP, orf3) found in the 324 genomes of Fusarium and related genera (per assembly in Table S7), "
    "with alignment to the reference proteins of F. poae Fp133 and gene-model state.", lo, [
        ("gca", "assembly", "Genome assembly accession"), ("organism", "organism", "Organism"), ("species", "species", "Species"),
        ("level", "assembly_level", "Assembly level"), ("gene", "island_gene", "Gene of the island", gn_map), ("seed", "reference_protein", "Reference protein of F. poae Fp133"),
        ("locus_id", "locus_id", "Locus identifier"), ("contig", "contig", "Contig"), ("start", "start", "Start"), ("end", "end", "End"), ("strand", "strand", "Strand"),
        ("contig_len", "contig_length_bp", "Contig length, bp"), ("dist_left", "distance_to_contig_start_bp", "Distance to the contig start, bp"),
        ("dist_right", "distance_to_contig_end_bp", "Distance to the contig end, bp"), ("identity", "identity_to_reference", "Identity to the reference protein"),
        ("coverage", "coverage_of_reference", "Coverage of the reference protein"), ("status", "status", "intact / truncated / pseudogene"),
        ("qfrom", "ref_aln_start_aa", "Alignment start on the reference, aa"), ("qto", "ref_aln_end_aa", "Alignment end on the reference, aa"), ("seed_len", "ref_length_aa", "Reference length, aa"),
        ("exons", "n_exons", "Exons in the model"), ("frameshifts", "frameshifts", "Frameshifts relative to the reference"), ("stops", "in_frame_stops", "In-frame stop codons"),
        ("missing", "missing_ref_segments", "Missing N-/C-terminal segments of the reference"), ("edge", "contig_edge", "Locus truncated by a contig edge"),
        ("gene_gc", "locus_gc_pct", "GC content of the locus, %"), ("assembly_gc", "assembly_gc_pct", "GC content of the assembly, %"),
        ("note", "note", "Free-text note")])
# Table S14b (odh island per species) replaced by the new Table S7 (324 assemblies), 2026-10-13q
rs = rd(OI + "/reverse_search_summary.tsv")
add("S14c", "Reverse search of the island genes against NCBI nr",
    "Best hits in NCBI nr (DIAMOND blastp, sensitive mode, up to 1,000 targets) for the three genes of the island in F. poae and F. langsethiae and for the two flanking fungal genes of F. poae Fp133. "
    "Hits from Fusarium and Neocosmospora are excluded from the 'best hit' columns. The three island genes have their best hits in three adjacent genes of one plasmid locus of Agrobacterium larrymoorei SORGH_AS_0974 (NZ_JAVIZC010000002.1: WP_309770756.1, WP_309770757.1/WP_309770758.1, WP_309770759.1).", rs, [
        ("query", "query_protein", "Query protein"), ("gene", "gene", "Gene and source", g3_disp), ("query_length_aa", "query_length_aa", "Query length, aa"),
        ("hits_total", "hits_total", "Hits in nr"), ("hits_Fusarium_Neocosmospora", "hits_Fusarium_Neocosmospora", "Hits from Fusarium/Neocosmospora"),
        ("best_non_Fusarium_hit", "best_hit", "Best hit outside Fusarium/Neocosmospora"), ("best_non_Fusarium_organism", "best_hit_organism", "Organism of the best hit"),
        ("best_non_Fusarium_identity", "best_hit_identity_pct", "Identity to the best hit, %"), ("best_non_Fusarium_aln_len", "best_hit_alignment_length_aa", "Alignment length, aa"),
        ("second_best_bacterial_identity", "second_bacterial_hit_identity_pct", "Identity to the second-best bacterial hit, %"),
        ("best_fungal_non_Fusarium_hit", "best_fungal_hit", "Best fungal hit outside Fusarium/Neocosmospora"), ("best_fungal_non_Fusarium_organism", "best_fungal_hit_organism", "Organism of the best fungal hit"),
        ("best_fungal_non_Fusarium_identity", "best_fungal_hit_identity_pct", "Identity to the best fungal hit, %"),
        ("top50_non_Fusarium_fungal", "fungal_among_top50", "Fungal proteins among the 50 best hits outside Fusarium/Neocosmospora"), ("top50_non_Fusarium_other", "non_fungal_among_top50", "Non-fungal proteins among the 50 best hits outside Fusarium/Neocosmospora")])
runs = {r["run"]: r for r in rd(OI + "/expression/runs_meta.tsv")}
smm = rd(OI + "/expression/summary.tsv")
for r in smm:
    r.update({k: v for k, v in runs.get(r["run"], {}).items() if k not in r})
add("S14d", "RNA-seq of Fusarium poae Fp133: expression of the island genes",
    "RNA-seq runs of the carrier strain F. poae Fp133 (library metadata) and expression of each gene of the island relative to the genes of the genome and to reference genes. "
    "The first 10 million reads of each run were mapped to the Fp133 assembly (GCA_046056315.1).", smm, [
        ("run", "sra_run", "SRA run accession"), ("BioProject", "bioproject", "BioProject"), ("BioSample", "biosample", "BioSample"), ("strain", "strain", "Strain"),
        ("condition", "condition", "Experimental condition"), ("platform", "platform", "Sequencing platform"), ("layout", "layout", "Library layout"),
        ("reads_in_run", "reads_in_run", "Reads in the run"), ("n_reads_used", "reads_used", "Reads used (subsample)"), ("mapped_pct", "mapped_pct", "Mapped reads, %"),
        ("gene", "island_gene", "Gene of the island and copy", g3_disp), ("locus_tag", "locus_tag", "Locus tag in the Fp133 annotation"),
        ("reads_on_target", "reads_on_gene", "Reads on the gene"), ("target_TPM", "gene_TPM", "TPM of the gene"), ("target_reads_per_kb", "gene_reads_per_kb", "Reads per kb on the gene"),
        ("pct_all_genes_strictly_below", "pct_genes_with_lower_expression", "Percent of all genes with expression strictly below the gene"),
        ("target_TPM_over_median_expressed", "gene_TPM_over_median_expressed", "Gene TPM / median TPM of genes with reads"),
        ("actin_TPM", "actin_TPM", "TPM of actin"), ("EF1a_TPM", "EF1a_TPM", "TPM of EF1a"), ("GAPDH_TPM", "GAPDH_TPM", "TPM of GAPDH"), ("tubulinB_TPM", "tubulin_beta_TPM", "TPM of beta-tubulin"),
        ("lib_sense_frac", "sense_fraction", "Fraction of reads in sense orientation over genes with >=200 reads (library strandedness)"),
        ("target_breadth", "gene_breadth", "Fraction of the gene covered by reads")])

# ---------------------------------------------------------------- fixed order of genes / lineages / families
# Same order as all tables/figures/texts of the paper (v2): genes with T-DNA-type lineages; genes with lineages from other rhizobial loci;
# genes without findings beyond rule A. Lineages: bacterial-origin first, then not analysed / native; '(not placed)' right after its lineage.
GENE_ORDER = ["sus", "cus", "chs", "acs", "acs_pKerr14a", "nos", "ocs", "vis", "ags", "mas1", "mas2", "mis"]
LINEAGE_ORDER = ["sus-like", "cus-like", "chs-like", "acc-type GDPD", "pKerr14a-type acs-like", "nos-subfamily opine dehydrogenase (Boletales)",
                 "odh island: odh", "odh island: pepP", "odh island: orf3",
                 "nos-subfamily ODH homolog (Rhypophila, not assigned)", "OCD-like, outside the sus-like clade", "OCD-like bacterial", "actinobacterial-type", "other bacterial-type Cus (beta-proteobacteria)",
                 "other bacterial-type Cus", "mis group", "other bacterial-type small ORF (cusrel181)", "rcorf12", "outgroup-HMM (Cus)", "native fungal GDPD",
                 "GDPD of unresolved position", "GDPD (acs profile), not placed on a tree", "no HMM subgroup", "outgroup-HMM", "other opine dehydrogenase (Boletales)", "not a homolog"]
FAMILY_ORDER = ["OCD fold (Sus)", "cucumopine synthase (Cus)", "GDPD"]  # sus, cus, acs (opine dehydrogenases after)
def gene_rank(g):
    gs = [x.strip() for x in str(g).split(",")]
    return (min(GENE_ORDER.index(x) if x in GENE_ORDER else len(GENE_ORDER) for x in gs), len(gs), str(g))
def lineage_rank(lin, fam=""):
    base = lin.replace(" (not placed)", "")
    if base == "outgroup-HMM" and fam.startswith("cucumopine"): base = "outgroup-HMM (Cus)"  # Cus outgroup-HMM fragments go with the Cus lineages
    return (LINEAGE_ORDER.index(base) if base in LINEAGE_ORDER else len(LINEAGE_ORDER), int(base != lin), lin)
def family_rank(f):
    return (FAMILY_ORDER.index(f) if f in FAMILY_ORDER else len(FAMILY_ORDER), f)
SORTS = {"S1b": lambda r: lineage_rank(r["lineage"], r["family"]),
         "S3a": lambda r: gene_rank(r["gene"]),
         "S3b": lambda r: (r["nonfungal_nearest"] != "yes", gene_rank(r["gene"]), r["organism"]),   # non-fungal-nearest rows first (2026-10-13al), then by gene, then by organism
         "S5": lambda r: family_rank(r["family"]),
         "S9a": lambda r: family_rank(r["family"]),
         "S10a": lambda r: lineage_rank(r["lineage"], r["family"])}
for i, (tid, title, desc, cols, rows) in enumerate(TABLES):
    if tid in SORTS:
        TABLES[i] = (tid, title, desc, cols, sorted(rows, key=SORTS[tid]))  # stable: previous order kept within a group
        unknown = {r["lineage"] for r in rows if "lineage" in r and lineage_rank(r["lineage"], r.get("family", ""))[0] == len(LINEAGE_ORDER)}
        if unknown: print("WARNING", tid, "lineages not in LINEAGE_ORDER:", unknown)
S1 = dict((t[0], t[4]) for t in TABLES)["S1b"]

# ---------------------------------------------------------------- S15 (sensitivity of thresholds; journal 2026-10-12 thresholds_sensitivity)
# Added after the sorting step: S15b is sorted here (contaminant contigs first, then fungal contigs by the fixed lineage order).
TH = A + "/thresholds"
no_six = lambda v: re.sub(r"\bsix\b", "6", v)
SECTION = {"probable contamination": "2.1", "rule A": "2.1", "rule B": "2.1", "rule C": "2.1", "plasmid backbone": "2.3",
           "chs-like locus": "2.2 (ii)", "pKerr14a-type fragment": "2.4"}
CLSN = {"B": "T-DNA candidate", "C": "non-T-DNA candidate", "single": "other bacterial-type protein", "marginal": "member of a family with a single better bacterial hit", "contamination": "probable contamination"}
def defang(v):
    """thresholds_summary.tsv uses internal rule names; display the terms of the paper (Table 1)"""
    v = re.sub(r"\b(B|C|single|marginal|contamination)\s?→\s?(B|C|single|marginal|contamination)\b", lambda m: CLSN[m.group(1)] + " → " + CLSN[m.group(2)], v)
    v = re.sub(r"\s→C\b", " → non-T-DNA candidate", v)
    v = re.sub(r"\s?→marginal\b", " → member of a family with a single better bacterial hit", v)
    v = v.replace("pass only via this rule", "pass only via this criterion")
    for a, b in (("rule A:", "bacterial-type protein:"), ("rule B:", "T-DNA candidate:"), ("rule C:", "non-T-DNA candidate:"),
                 ("plasmid backbone (vir/tra/rep) on the same contig", "vir/tra/rep genes on the same contig"),
                 ("(OR-rule)", "(alternative criterion)"), ("rule-A proteins", "bacterial-type proteins"), ("rule A", "bacterial-type"), ("rule B", "T-DNA candidate"), ("rule C", "non-T-DNA candidate"),
                 ("among A proteins", "among bacterial-type proteins"), ("marginal mostly", "members of families with a single better bacterial hit mostly"),
                 ("'backbone only in B. harrisonii'", "vir/tra/rep only in B. harrisonii"), ("backbone", "vir/tra/rep"), ("short_contig flag", "short-contig flag"),
                 ("acs-HMM remnants", "GDPD (acs profile) remnants"), ("single→", "other bacterial-type protein →")):
        v = v.replace(a, b)
    return clean(v)
s16a = rd(TH + "/thresholds_summary.tsv")
for r in s16a:
    r["methods_section"] = next((v for k, v in SECTION.items() if r["threshold"].startswith(k)), "")
    for k in ("threshold", "basis", "range_no_change", "what_changes_outside", "value"):
        r[k] = defang(r[k])
s16a_spec = [
    ("threshold", "threshold", "Threshold (classification step and statistic)"),
    ("methods_section", "methods_section", "Section of Materials and Methods where the threshold is used"),
    ("value", "value", "Value used in the paper"),
    ("basis", "basis", "Why the value is not arbitrary: a gap in the distribution of the statistic, or a range of values over which the result is unchanged; nearest proteins/loci on either side"),
    ("range_no_change", "range_no_change", "Range of threshold values over which the result (classification of proteins, loci or contigs) does not change"),
    ("what_changes_outside", "what_changes_outside", "What changes outside this range (proteins, loci or contigs that change class, with their values)", no_six)]
s16b = rd(TH + "/backbone/backbone_contigs.tsv")
for r in s16b:
    r["lineage"] = norm_lineage(r["lineage"])
    r["gc_minus_assembly"] = "%.1f" % (float(r["GC"]) - float(r["assembly_gc"])) if r.get("assembly_gc") else ""
s16b.sort(key=lambda r: (r["sample"] != "bacterial", lineage_rank(r["lineage"])[0] if r["sample"] != "bacterial" else 0))
s16b_spec = [
    ("sample", "sample", "Contig set: 'bacterial' = contigs carrying the 11 probable contaminants of the protein-level screen; 'fungal' = control contigs carrying loci of the catalogue"),
    ("contig", "contig", "Contig accession"),
    ("source", "source", "Protein (gene, best bacterial hit and its identity) or locus on the contig"),
    ("lineage", "lineage", "Class of the protein or lineage of the locus"),
    ("gca", "assembly", "Genome assembly accession"), ("organism", "organism", "Organism"),
    ("length", "contig_length_bp", "Contig length, bp"), ("GC", "contig_GC_pct", "GC content of the contig, %"),
    ("assembly_gc", "assembly_GC_pct", "GC content of the whole assembly, % (contaminant contigs only)"),
    ("gc_minus_assembly", "GC_difference_pct_points", "Contig GC minus assembly GC, percentage points"),
    ("genes_per_kb", "genes_per_kb", "Predicted genes per kb"),
    ("coding_frac", "coding_fraction", "Fraction of the contig covered by predicted coding sequences"),
    ("best_backbone_pid", "vir_tra_rep_best_identity_pct", "Best identity of a predicted protein to a plasmid maintenance and transfer protein (vir, tra, rep), %; 0 = no hit"),
    ("best_backbone_label", "vir_tra_rep_best_hit", "Best hit among vir, tra and rep proteins (category:gene:organism)"),
    ("n_backbone_hits_ge60", "vir_tra_rep_hits_ge60pct", "Predicted proteins with a vir, tra or rep hit at >=60% identity"),
    ("bx_best_backbone_pid", "vir_tra_rep_blastx_best_identity_pct", "Best identity of the contig (blastx, as in the genome-level screen) to a vir, tra or rep protein, %")]
for tid, title, desc, rows, spec in [
        ("S15a", "Classification thresholds and their sensitivity",
         "Every threshold used to classify proteins, loci and contigs, with the reason for its value and the range of values over which the result is unchanged "
         "(the classification was repeated over a grid of values; the distributions are shown in Figure S1). One row per threshold.", s16a, s16a_spec),
        ("S15b", "Contaminant and control contigs",
         "Contigs carrying the 11 probable contaminants of the protein-level screen (>=90% identity to a non-rhizobial bacterial protein) and 43 control contigs carrying loci of the catalogue: "
         "composition, gene density and plasmid maintenance and transfer genes (vir, tra, rep). None of the contaminant contigs carries Ti/Ri plasmid vir, tra or rep genes (best identity <=31.6%). "
         "The vir/tra/rep search of the Boletales assemblies is in Table S5b. One row per contig.", s16b, s16b_spec)]:
    pr, cols = proj(rows, spec)
    TABLES.append((tid, title, desc, cols, pr))

# S15c (2026-10-13o, audit A2): vir/tra/rep loci in the Boletales assemblies (provisional number; renumbering follows). Source: analysis/refset/screen/step2/boletales517_tdna_summary/vir_tra_rep_hits.tsv
vt = rd(S2D + "/vir_tra_rep_hits.tsv")
for r in vt:
    sp = r["subject"].split("|")   # accession|category|plasmid|organism|gene|product
    r["best_protein"] = "%s (%s, %s)" % (sp[5].replace("_", " "), sp[3].replace("_", " "), sp[0])
    r["cat_label"] = {"tra": "tra (conjugal transfer)", "vir": "vir (virulence)", "rep_partition": "rep (replication/partition)", "mobile": "mobile element"}.get(sp[1], sp[1])
    r["carrier"] = "yes" if r["genome_is_carrier"] == "1" else "no"
    r["carries"] = "yes" if r["contig_has_carrier_locus"] == "1" else "no"
add("S15c", "vir, tra and rep genes in the Boletales assemblies",
    "Contigs of the 517 Boletales assemblies carrying a locus that matches a Ti/Ri plasmid maintenance or transfer protein (vir, tra, rep) at >=60% identity "
    "(blastx of the assemblies against the plasmid protein catalogue; hits within 5 kb on one strand merged into a locus). Two assemblies have such contigs: "
    "Boletus harrisonii (a carrier assembly; 368-bp contig) and Pisolithus albus (not a carrier; 488-bp contig); neither contig carries a locus of the catalogue. One row per locus.",
    vt, [
        ("gca", "assembly", "Genome assembly accession"), ("organism", "organism", "Organism"),
        ("carrier", "assembly_carries_lineage", "'yes' if the assembly carries a sus-like, cus-like, acc-type GDPD or pKerr14a-type acs-like locus"),
        ("contig", "contig", "Contig"), ("contig_len", "contig_length_bp", "Contig length, bp"),
        ("carries", "contig_carries_lineage_locus", "'yes' if the same contig carries such a locus"),
        ("start", "start", "Start of the locus on the contig"), ("end", "end", "End of the locus on the contig"), ("strand", "strand", "Strand"),
        ("cat_label", "gene_class", "Class of the best-matching plasmid protein"),
        ("best_protein", "best_plasmid_protein", "Best-matching plasmid protein (product, organism, accession)"),
        ("pident", "identity_pct", "Identity, %"),
        ("subj_cov", "query_coverage", "Fraction of the plasmid protein covered by the alignment")])

# ---------------------------------------------------------------- S16 (dN/dS of the remaining lineages; journal 2026-10-12 dnds_extra)
SET17 = {"chs": "chs-like", "odh": "odh island: odh", "pepP": "odh island: pepP", "g3": "odh island: orf3",
         "nosODH_Bol": "nos-subfamily opine dehydrogenase (Boletales)", "cus_nonBol": "cus-like, outside Boletales (all intact copies)",
         "cus_nonBol_sp": "cus-like, outside Boletales (one copy per species)", "acs_pKerr14a": "pKerr14a-type acs-like",
         "cus_Bol": "cus-like, Boletales (all intact copies)", "cus_Bol_sp": "cus-like, Boletales (one copy per species)",
         "cus_Bol_conf": "cus-like, Boletales, loci with confirmed or supported integration only (7 genes)",
         "cus_Bol_sp_conf": "cus-like, Boletales, one gene per species, confirmed or supported only (4 genes)",
         "acs_pKerr14a_conf": "pKerr14a-type acs-like, confirmed or supported only (4 genes)",
         "nosODH_Bol_conf": "nos-subfamily ODH, Boletales, confirmed or supported only (3 genes; star tree, M0 only)"}
s17 = rd(A + "/dnds/extra/summary.tsv")
for r in s17: r["set_name"] = SET17.get(r["set"], r["set"])
_ord17 = {r["set"]: i for i, r in enumerate(s17)}
def _key17(r):   # lineage order; within a lineage the summary order, with each *_conf set directly after its parent set (2026-10-13al)
    par = r["set"][:-5] if r["set"].endswith("_conf") else r["set"]
    pn = SET17.get(par, par)
    return (lineage_rank(pn.split(",")[0].split(" (")[0] if par.startswith("cus") else pn), _ord17[par], r["set"] != par)
s17.sort(key=_key17)
s17_spec = [("set_name", "gene_set", "Gene set (intact copies; CDS of the Boletales nos-subfamily ODH from miniprot models)"), ("n_CDS", "n_CDS", "Number of CDS"),
            ("codons_strict", "codons_strict", "Codons in the alignment without gap columns"), ("codons_le20gap", "codons_le20pct_gaps", "Codons in the alignment with <=20% gaps per column"),
            ("M0_omega_strict", "M0_omega_strict", "codeml M0 omega, strict alignment"), ("M0_omega_le20gap", "M0_omega_le20pct_gaps", "codeml M0 omega, <=20% gaps alignment"),
            ("LRT_M7vsM8", "LRT_M7_vs_M8", "Likelihood-ratio statistic M7 vs M8"), ("p_M7vsM8", "p_M7_vs_M8", "p-value M7 vs M8"), ("p_M8avsM8", "p_M8a_vs_M8", "p-value M8a vs M8"),
            ("BEB_sites_P>0.95", "BEB_sites_P_gt_0.95", "Sites with Bayes empirical Bayes posterior > 0.95 under M8")]
for v, lab in (("strict", "strict alignment (yn00 removes all gap columns)"), ("le20gap", "pairwise gap deletion on the <=20% gaps alignment (yn00_pairwise_loose.py)")):
    s17_spec += [("yn00_%s_pairs_estimable" % v, "yn00_%s_pairs_estimable" % v, "yn00, %s: pairs with estimable dN/dS" % lab),
                 ("yn00_%s_pairs_used" % v, "yn00_%s_pairs_used" % v, "yn00, %s: pairs used" % lab),
                 ("yn00_%s_dNdS_median" % v, "yn00_%s_dNdS_median" % v, "yn00, %s: median pairwise dN/dS" % lab),
                 ("yn00_%s_dNdS_range" % v, "yn00_%s_dNdS_range" % v, "yn00, %s: range of pairwise dN/dS" % lab),
                 ("yn00_%s_dS_median" % v, "yn00_%s_dS_median" % v, "yn00, %s: median pairwise dS" % lab),
                 ("yn00_%s_dS_range" % v, "yn00_%s_dS_range" % v, "yn00, %s: range of pairwise dS" % lab),
                 ("yn00_%s_pairs_dS_gt2" % v, "yn00_%s_pairs_dS_gt2" % v, "yn00, %s: pairs with dS > 2 (saturation)" % lab),
                 ("yn00_%s_dNdS_median_dS_le2" % v, "yn00_%s_dNdS_median_dS_le2" % v, "yn00, %s: median dN/dS over pairs with dS <= 2" % lab)]
pr, cols = proj(s17, s17_spec)
TABLES.append(("S16", "Selection analysis of the remaining lineages (dN/dS)",
               "dN/dS of the intact copies of the lineages not covered by Table S9a: codeml M0, M7/M8, M8a and yn00 on codon alignments (MAFFT, IQ-TREE GTR+G tree), "
               "same procedure as Table S9a. One row per gene set.", cols, pr))

# ---------------------------------------------------------------- S17 (identity of plant cT-DNA and fungal rule-B proteins to bacterial T-DNA proteins; journal 2026-10-12 tableS18)
CT = A + "/ctdna_plants_identity"
GRP = {"plant_cTDNA": "plant cT-DNA", "fungal_ruleB": "fungal T-DNA candidate"}
s18a = rd(CT + "/summary.tsv"); s18b = rd(CT + "/pairwise_hits.tsv")
for r in s18a + s18b: r["group"] = GRP.get(r["group"], r["group"])
s18a.sort(key=lambda r: (gene_rank(r["gene"]), r["group"] != "plant cT-DNA"))
s18b.sort(key=lambda r: (gene_rank(r["gene"]), r["group"] != "plant cT-DNA"))
s18a_spec = [
    ("gene", "gene", "T-DNA gene"), ("group", "protein_group", "Protein group: plant cT-DNA (cellular T-DNA of naturally transgenic plants) or fungal T-DNA candidate (T-DNA-type fungal proteins)"),
    ("n", "n_proteins", "Number of proteins"),
    ("median_pident", "median_identity_pct", "Median identity to the nearest bacterial T-DNA seed of the same gene, %"),
    ("min", "min_identity_pct", "Minimum identity, %"), ("max", "max_identity_pct", "Maximum identity, %"),
    ("median_qcov", "median_query_coverage", "Median fraction of the protein covered by the alignment"),
    ("n_fulllen_scov70", "n_full_length", "Proteins covering at least 70% of the seed length"),
    ("median_pident_fulllen", "median_identity_full_length_pct", "Median identity of the full-length proteins, %"),
    ("min_fulllen", "min_identity_full_length_pct", "Minimum identity of the full-length proteins, %"),
    ("max_fulllen", "max_identity_full_length_pct", "Maximum identity of the full-length proteins, %"),
    ("n_unique_seq", "n_unique_sequences", "Number of distinct sequences"),
    ("median_pident_unique", "median_identity_unique_pct", "Median identity over distinct sequences, %"),
    ("median_pident_seed11", "median_identity_reference_seed_pct", "Median identity to the reference bacterial seed of the gene (canonical seed of the 11-seed set), %"),
    ("range_seed11", "range_identity_reference_seed_pct", "Range of identity to the reference seed, %"),
    ("median_pident_to_plant_cTDNA", "median_identity_to_plant_cTDNA_pct", "Fungal proteins only: median identity to the nearest plant cT-DNA protein of the same gene, %")]
s18b_spec = [
    ("gene", "gene", "T-DNA gene"), ("group", "protein_group", "Protein group: plant cT-DNA or fungal T-DNA candidate"),
    ("query", "protein", "Protein accession (or identifier of the miniprot model, Ntom_*)"),
    ("status", "status", "'include' = used in Table S17a; otherwise the reason for exclusion"),
    ("origin", "origin", "Source of the sequence"), ("qlen", "protein_length_aa", "Protein length, aa"),
    ("best_seed", "best_bacterial_seed", "Nearest bacterial seed of the same gene"),
    ("pident", "identity_to_bacterial_seed_pct", "Identity to that seed, %"),
    ("qcov", "query_coverage", "Fraction of the protein covered by the alignment"), ("scov", "seed_coverage", "Fraction of the seed covered by the alignment"),
    ("bits", "bit_score", "Bit score"), ("evalue", "e_value", "E-value"),
    ("seed11", "reference_seed", "Reference bacterial seed of the gene (11-seed set)"),
    ("pident_seed11", "identity_to_reference_seed_pct", "Identity to the reference seed, %"),
    ("qcov_seed11", "query_coverage_reference_seed", "Fraction of the protein covered by the alignment to the reference seed"),
    ("top_any_seed", "nearest_seed_any_gene", "Nearest bacterial seed of any gene (gene|accession)"),
    ("top_any_pident", "identity_to_nearest_seed_any_gene_pct", "Identity to the nearest seed of any gene, %"),
    ("best_plasmid", "nearest_plasmid_protein", "Nearest protein of the Ti/Ri plasmid catalogue (category|annotation)"),
    ("pident_plasmid", "identity_to_plasmid_protein_pct", "Identity to that protein, %"),
    ("qcov_plasmid", "query_coverage_plasmid_protein", "Fraction of the protein covered by the alignment to the plasmid protein"),
    ("best_plant_cTDNA", "nearest_plant_cTDNA_protein", "Fungal proteins only: nearest plant cT-DNA protein of the same gene"),
    ("pident_plant_cTDNA", "identity_to_plant_cTDNA_protein_pct", "Identity to the plant cT-DNA protein, %"),
    ("qcov_plant_cTDNA", "query_coverage_plant_cTDNA_protein", "Fraction of the protein covered by the alignment to the plant cT-DNA protein")]
for tid, title, desc, rows, spec in [
        ("S17a", "Identity of plant cT-DNA and fungal T-DNA-candidate proteins to bacterial T-DNA proteins: summary",
         "Summary by gene and protein group of the amino-acid identity to the nearest bacterial T-DNA seed of the same gene (30 bacterial seeds; DIAMOND blastp, very-sensitive). "
         "Plant cT-DNA proteins: seed-set members of plant origin, GenBank/Entrez entries and miniprot models in the Nicotiana tomentosiformis genome; fungal T-DNA candidates: proteins classified as T-DNA candidates (Table S3). "
         "Proteins excluded in Table S17b are not counted. One row per gene and group.", s18a, s18a_spec),
        ("S17b", "Identity of plant cT-DNA and fungal T-DNA-candidate proteins to bacterial T-DNA proteins: by protein",
         "Pairwise comparison of every plant cT-DNA and fungal T-DNA-candidate protein with the bacterial T-DNA seeds of the same gene, with the reference seed, with the Ti/Ri plasmid protein catalogue and, for fungal proteins, with the plant cT-DNA proteins. "
         "One row per protein.", s18b, s18b_spec)]:
    pr, cols = proj(rows, spec)
    TABLES.append((tid, title, desc, cols, pr))

norm_odh = lambda v: v.replace("odh island: opine dehydrogenase", "odh island: odh").replace("odh island: Xaa-Pro aminopeptidase", "odh island: pepP").replace("odh island: third gene", "odh island: orf3")  # 2026-10-13h display names
# ---------------------------------------------------------------- S18 (full per-group summary of the catalogue = former Table 2 of the main text; 2026-10-12x, review_critical 6.2)
t2a = rd(G + "/paper/tables/Table2a_lineages.tsv")
s19_spec = [("group", "group", "Tree-defined group of the catalogue (Table S1b). The first nine rows (main_lineage = True) are the seven lineages of bacterial origin (the odh island has three genes); the remaining rows are groups of native fungal families or groups whose origin was not analysed", norm_odh),
    ("family", "family", "Protein family"), ("origin", "origin", "Closest relatives / source of the group as established in Results 3.4 ('not analysed' for the groups outside the seven lineages)"),
    ("entry_route", "entry_route", "Route by which the loci entered the screen (Methods 2.2)"), ("n_loci", "n_loci", "Loci of the lineage in Table S1a (placed on a tree directly or through a representative, or assigned by identity, or added from NCBI annotations)"),
    ("intact", "intact", "Intact genes"), ("pseudogene", "pseudogene", "Pseudogenes (disrupted reading frame)"), ("truncated", "truncated", "Truncated copies (reading frame intact, coverage below 90%; chs-like loci of Fusarium)"),
    ("fragment", "fragment", "Fragments and remnants too short to be placed on a tree, assigned by identity (pKerr14a-type) or profile"), ("other", "other", "Loci with another status"),
    ("n_assemblies", "n_assemblies", "Assemblies with a placed locus"), ("n_species", "n_species", "Distinct NCBI species names"), ("n_orders", "n_orders", "Fungal orders"),
    ("Boletales", "assemblies_Boletales", "Carrier assemblies in Boletales"), ("Agaricales", "assemblies_Agaricales", "Carrier assemblies in Agaricales"), ("Hypocreales", "assemblies_Hypocreales", "Carrier assemblies in Hypocreales"),
    ("other_orders", "assemblies_other_orders", "Carrier assemblies in other orders (order and count)"), ("profile_only_fragments", "fragments_assigned_by_profile_only", "Loci detected only by the profile scan and not placed on a tree"),
    ("assemblies_incl_fragments", "assemblies_incl_fragments", "Assemblies with any locus of the group, including profile-only fragments"), ("main_lineage", "lineage_of_bacterial_origin", "True for the seven lineages of bacterial origin")]
add("S18", "Summary of the catalogue by group",
    "Per-group summary of the catalogue of loci (Table S1b): the seven lineages of bacterial origin (Table 2 of the main text) and the groups of the catalogue that belong to native fungal families or were not analysed further. "
    "One row per group.", t2a, s19_spec)

# ---------------------------------------------------------------- S1a (loci of bacterial origin) + Supplementary Data S1 (proteins); 2026-10-13d
# S1a = the rows of S1b flagged lineage_of_bacterial_origin; evidence_tables and protein_in_Data_S1 are collected from the tables built above
# (not from files on disk). Data S1 proteins: profile-screen models (loci_pipeline/round2/proteins), miniprot translations (##STA lines of
# mp_out/<GCA>/<seed>.gff) for chs-like and odh-island loci, nos_check/boletales_odh.faa for the nos-subfamily ODH of Boletales.
TT = dict((t[0], t[4]) for t in TABLES)
s1b_rows = TT["S1b"]
bact = [r for r in s1b_rows if r["lineage_of_bacterial_origin"]]
t2a_main = [r for r in t2a if r["main_lineage"] == "True"]
n_placed_s19 = sum(int(r["n_loci"]) for r in t2a_main); n_frag_s19 = sum(int(r["profile_only_fragments"]) for r in t2a_main)
assert len(bact) == sum(1 for r in s1b_rows if r["lineage_of_bacterial_origin"]) == n_placed_s19 + n_frag_s19, (len(bact), n_placed_s19, n_frag_s19)
assert all(r["origin"].startswith(("T-DNA", "non-T-DNA")) for r in bact)

def ckey(c): return re.sub(r"\.\d+$", "", c or "")
ev = collections.defaultdict(set); evc = collections.defaultdict(set)
for tid in ("S6", "S10a", "S12a", "S13", "S14a"):
    for r in TT[tid]: ev[r["locus_id"]].add(tid)
for r in TT["S7"]:
    evc[(r["source_assembly"], ckey(r["contig"]))].add("S7")
    for l in re.split(r"[;,\s|]+", r["loci"]):
        if l: ev[l].add("S7")
for tid, ga, cc in (("S8", "assembly", "contig"), ("S11b", "assembly", "contig"), ("S15b", "assembly", "contig")):
    for r in TT[tid]: evc[(r[ga], ckey(r[cc]))].add(tid)
ISL_GENE = {"odh island: odh": "odh", "odh island: pepP": "pepP", "odh island: orf3": "g3"}
s15c_genes = {("g3" if r["gene"].split()[0] == "orf3" else r["gene"].split()[0]) for r in TT["S14c"]}
FP133 = "GCA_046056315.1"   # assembly of the RNA-seq analysis (Table S14d)
s15d_keys = set()
for r in TT["S14d"]:
    m = re.match(r"(odh|pepP|orf3) \((CM\d+)", r["island_gene"])
    if m: s15d_keys.add((FP133, m.group(2), "g3" if m.group(1) == "orf3" else m.group(1)))
EVMAP = {"S6": "S8", "S10a": "S14", "S12a": "S13a", "S13": "S16", "S14a": "S13b", "S7": "S9", "S8": "S12", "S11b": "S10b", "S15b": "S5a", "S14c": "S13c", "S14d": "S16"}
def tid_key(t):
    m = re.match(r"S(\d+)(\w*)", t); return (int(m.group(1)), m.group(2))

# ---- proteins
def read_fa(path, key=lambda h: h.split("|")[0].split()[0]):
    d = {}; k = None
    for l in open(path):
        l = l.rstrip()
        if l.startswith(">"): k = key(l[1:]); d[k] = [l[1:], ""]
        elif k: d[k][1] += l.strip()
    return d
PRT = A + "/loci_pipeline/round2/proteins/"
prof_prot = {}
for f in ("ACS", "CUS", "SUS"):
    for suffix in ("_remnants_fragments", "_models"):   # models take precedence
        for k, (h, sq) in read_fa(PRT + f + suffix + ".faa").items():
            prof_prot[k] = (sq, "models" if suffix == "_models" else "remnants_fragments")
# 2026-10-13al: the two late cus-like loci (addendum; not on the trees) are NCBI-annotated proteins (KAO1061018.1, CAL9735043.1); sequences under the locus ids
for _k, (_h, _sq) in read_fa(A + "/dnds/extra/cus_nonBol/cus_nonBol.pep.faa").items():
    if _k in ("CUS_JBZWEJ010000015_A1", "CUS_OZ188520_A1"): prof_prot[_k] = (_sq, "ncbi_annotation")
assert sum(1 for v in prof_prot.values() if v[1] == "ncbi_annotation") == 2
# 2026-10-13w: JBBPKL010003317 - the merged AUGUSTUS model (657 aa) is replaced by the two proteins split by the HMM envelope (251 aa GDPD, 203 aa mis);
# coordinates re-derived with miniprot (analysis/loci_pipeline/checks/split.gff)
SPLIT_SEQ = {}
for _l in open(A + "/curated/split_models/JBBPKL010003317_ACS_CUS_split_by_HMM_envelope.faa"):
    _l = _l.rstrip()
    if _l.startswith(">"): _k = _l[1:].split("__")[1]; SPLIT_SEQ[_k] = ""
    else: SPLIT_SEQ[_k] += _l
_l2 = [x.rstrip() for x in open(A + "/curated/split_models/JBBPKX010000014_CUS_S2L2_single_exon_ORF_13898-15079.faa")]
SPLIT_SEQ["CUS_JBBPKX010000014_S2L2"] = "".join(_l2[1:])   # 2026-10-13w: AUGUSTUS model (616 aa) merged this gene with the nos-subfamily pseudogene; single-exon ORF 13898-15079 (-)
assert sorted(SPLIT_SEQ) == ["ACS_JBBPKL010003317_L1", "CUS_JBBPKL010003317_S2L1", "CUS_JBBPKX010000014_S2L2"] and [len(v) for _, v in sorted(SPLIT_SEQ.items())] == [251, 203, 393]
for _k, _v in SPLIT_SEQ.items(): prof_prot[_k] = (_v, "models")
def mp_translation(gff, contig, start, end, strand):
    sta = ""
    with open(gff) as fh:
        for line in fh:
            if line.startswith("##STA"): sta = line.rstrip("\n").split("\t")[1].strip(); continue
            if line[0] == "#": continue
            f = line.split("\t")
            if len(f) > 8 and f[2] == "mRNA" and f[0] == contig and f[3] == str(start) and f[4] == str(end) and f[6] == strand:
                return sta
    return None
mp_prot = {}
for tid, base in (("S12a", A + "/chs_fusarium/stage2/mp_out/"), ("S14a", A + "/fusarium_odh_island/scan/mp_out/")):
    seedcol = "primary_seed" if tid == "S12a" else "reference_protein"
    for r in TT[tid]:
        sq = mp_translation("%s%s/%s.gff" % (base, r["assembly"], r[seedcol]), r["contig"], r["start"], r["end"], r["strand"])
        if sq: mp_prot[r["locus_id"]] = (sq.rstrip("*"), "miniprot")
nos_prot = {}
for k, (h, sq) in read_fa(A + "/fusarium_odh_island/nos_check/boletales_odh.faa", key=lambda h: h).items():
    m = re.search(r"_(\d{9}\.\d)_(\S+)$", k)
    nos_prot["ODH_" + m.group(2)] = ("GCA_" + m.group(1), sq)

s1a_all = []   # all 406 loci of the seven lineages (placed + not placed); S1a = the placed ones only (2026-10-13e)
fasta_of = {}; src_of = {}; mism = []; len_mism = []; len_changed = collections.Counter()
for r in bact:
    o = dict(r)
    o["placed_on_tree"] = "no" if r["lineage"].endswith("(not placed)") else "yes"
    o["lineage_clean"] = r["lineage"].replace(" (not placed)", "")
    e = set(ev.get(r["locus_id"], ())) | set(evc.get((r["assembly"], ckey(r["contig"])), ()))
    g = ISL_GENE.get(o["lineage_clean"])
    if g:
        if g in s15c_genes: e.add("S14c")
        if (r["assembly"], ckey(r["contig"]), g) in s15d_keys: e.add("S14d")
    o["evidence_tables"] = "; ".join(sorted({EVMAP.get(t, t) for t in e}, key=tid_key))
    sq = None
    if r["detection_route"].startswith("profile screen"):
        if r["locus_id"] in prof_prot: sq, src = prof_prot[r["locus_id"]]
    elif r["locus_id"] in mp_prot: sq, src = mp_prot[r["locus_id"]]
    elif r["locus_id"] in nos_prot:
        if nos_prot[r["locus_id"]][0] == r["assembly"]: sq, src = nos_prot[r["locus_id"]][1], "nos_check"
        else: mism.append(r["locus_id"])
    if sq:
        src_of[r["locus_id"]] = src
        fasta_of[r["locus_id"]] = ">%s lineage=%s status=%s species=%s assembly=%s coords=%s:%s-%s(%s)\n%s\n" % (
            r["locus_id"], o["lineage_clean"], r["status"], r["species"], r["assembly"], r["contig"], r["start"], r["end"], r["strand"],
            "\n".join(sq.rstrip("*")[i:i + 80] for i in range(0, len(sq.rstrip("*")), 80)))  # terminal stop codon dropped, internal stops of pseudogenes kept
        # protein_length_aa = length of the fungal protein model (sequence without the terminal stop); replaces the aligned length of the bacterial seed
        # (miniprot loci) and the blank (nos-ODH); in S1b too (r is the S1b row itself)
        newlen = str(len(sq.rstrip("*")))
        if r["protein_length_aa"] != newlen:
            if r["detection_route"].startswith("profile screen"): len_mism.append((r["locus_id"], r["protein_length_aa"], newlen))
            len_changed[(o["lineage_clean"], o["placed_on_tree"])] += 1
        r["protein_length_aa"] = o["protein_length_aa"] = newlen
    o["protein_in_Data_S1"] = "yes" if sq else "no"
    s1a_all.append(o)
assert not mism, mism
print("protein_length_aa: changed for", sum(len_changed.values()), "loci (lineage, placed):", sorted(len_changed.items()))
# all other S1b rows with a gene model in ACS/CUS/SUS_models.faa: same rule (length without the terminal stop; the 'len=' of the faa header counts it) (2026-10-13f)
models_seq = {}
for f in ("ACS", "CUS", "SUS"):
    for k, (h, sq) in read_fa(PRT + f + "_models.faa").items(): models_seq[k] = sq
models_seq.update(SPLIT_SEQ)   # 2026-10-13w
_cat_ids = {r["locus_id"] for r in s1b_rows}
_nb_new = _nb_m1 = 0
for r in s1b_rows:
    if r["lineage_of_bacterial_origin"] or r["locus_id"] not in models_seq: continue
    nl = str(len(models_seq[r["locus_id"]].rstrip("*")))
    if r["protein_length_aa"] != nl:
        _nb_new += 1
        if r["protein_length_aa"] and int(r["protein_length_aa"]) - int(nl) == 1: _nb_m1 += 1
        r["protein_length_aa"] = nl
_missing = [k for k in models_seq if k not in _cat_ids]
print("protein_length_aa (S1b, other groups with a model in ACS/CUS/SUS_models.faa): new length for", _nb_new, "rows, of them -1:", _nb_m1)
print("models in ACS/CUS/SUS_models.faa without a row in the catalogue:", len(_missing), _missing[:10])
print("protein_length_aa: profile-screen loci whose model length differs from the previous value:", len(len_mism), len_mism[:15])
s1a = [o for o in s1a_all if o["placed_on_tree"] == "yes"]
s1a.sort(key=lambda r: (lineage_rank(r["lineage"], r["family"]), r["species"], r["contig"], num(r["start"]) or 0))
assert len(s1a) == n_placed_s19 == 254, (len(s1a), n_placed_s19)
assert len(s1a_all) - len(s1a) == n_frag_s19 == 152, (len(s1a_all), n_frag_s19)
fasta = []; src_cnt = collections.Counter(); noprot = collections.Counter()
for o in s1a:
    if o["locus_id"] in fasta_of:
        fasta.append(fasta_of[o["locus_id"]]); src_cnt[src_of[o["locus_id"]]] += 1
    else:
        noprot[(o["lineage_clean"] if o["detection_route"].startswith("blastx") else o["detection_route"][:14], o["status"])] += 1
DATA_S1 = "Data_S1_proteins_of_loci_of_bacterial_origin.faa"
with open(os.path.join(OUT, DATA_S1), "w") as fh:
    fh.write("".join(fasta))
n_prot = len(fasta); n_noprot = len(s1a) - n_prot
assert n_prot == sum(1 for r in s1a if r["protein_in_Data_S1"] == "yes") == sum(src_cnt.values())
print("Data S1:", n_prot, "sequences", dict(src_cnt), "; without protein:", n_noprot, dict(noprot))
ln_np = collections.Counter(r["lineage_clean"] for r in s1a)
print("S1a loci by lineage:", sorted(ln_np.items()))
n_frag = len(s1a_all) - len(s1a)
# 2026-10-13af: integration status per locus (paper/tables/integration_status.py; reviewer point 4)
INTST = {r["locus_id"]: r["status"] for r in rd(G + "/paper/tables/integration_status.tsv")}
assert set(INTST) == {o["locus_id"] for o in s1a}, "integration_status.tsv out of date: run paper/tables/integration_status.py"
for o in s1a: o["integration_status"] = INTST[o["locus_id"]]
print("S1a integration status:", collections.Counter(INTST.values()))
add("S1a", "Loci of bacterial origin",
    "The %d loci of the seven lineages of bacterial origin (239 placed on the phylogenetic trees directly or through a representative, 13 fragments assigned by identity, 2 added from NCBI annotations after the trees were built; named after their closest bacterial relatives) - one row per locus, "
    "with coordinates, structural status, the route of detection and the supplementary tables in which the locus appears. "
    "The %d additional fragments and remnants of the same lineages assigned by HMM profile only (not placed on a tree) are listed in Table S1b "
    "(lineage_of_bacterial_origin = True, lineage with the suffix '(profile only)'). Protein sequences in Data S1 (%d of the %d loci)." % (len(s1a), n_frag, n_prot, len(s1a)) +
    " Integration status (column integration_status), assigned by fixed rules (Methods 2.5) "
    "from the per-neighbour classification of the flanking genes (nr hits; hits to the carrier genus alone do not count), Table S8 (read depth), Table S9 (independent assemblies), Tables S13a, S13b and the NCBI assembly reports: "
    "confirmed = at least one of (i) a gene classified as fungal on both sides of the locus, "
    "(ii) the locus and both flanks covered (>= 0.5) by an independent assembly of the species (Table S9), "
    "(iii) the locus on a sequence listed as an assembled-molecule chromosome in the NCBI assembly report, >= 1 kb from both ends, or "
    "(iv) read-depth ratio 0.5-2 with junction reads on both flanks (Table S8); "
    "uncertain (not assessable; applied next) = a contig shorter than the locus plus 2 kb, or a locus touching a contig end; "
    "supported = none of the above and at least one of a fungal gene on one side only, a contig >= 100 kb with the locus away from its ends and a depth ratio not outside 0.5-2, "
    "a depth ratio 0.5-2 without junction reads on both flanks, or a gene of the same lineage present in two or more assemblies of the species without verification of the flanks (Tables S13a, S13b); "
    "the remaining loci (a depth ratio outside 0.5-2 without other support, or no classified neighbours and no other support) are uncertain.",
    s1a, [
        ("locus_id", "locus_id", "Locus identifier (as in Table S1b and in the headers of Data S1)"),
        ("lineage_clean", "lineage", "Lineage of bacterial origin (seven lineages; the odh island has three genes)"),
        ("origin", "origin", "Origin of the lineage as established in the paper (T-DNA gene or rhizobial locus outside the T-DNA)"),
        ("status", "status", "Locus state: intact, pseudogene, truncated, fragment, remnant or ambiguous"),
        ("species", "species", "Species name used in the paper"), ("order", "order", "Taxonomic order"),
        ("taxonomic_family", "taxonomic_family", "Taxonomic family of the carrier fungus"), ("assembly", "assembly", "Genome assembly accession"),
        ("contig", "contig", "Contig/scaffold accession"), ("contig_length_bp", "contig_length_bp", "Contig length, bp"),
        ("start", "start", "Locus start on the contig (1-based)"), ("end", "end", "Locus end on the contig"), ("strand", "strand", "Strand"),
        ("n_exons", "n_exons", "Number of exons in the gene model"), ("protein_length_aa", "protein_length_aa", PLEN_DESC),
        ("genomic_context", "genomic_context", "Flanking-gene context call (fungal, mixed, bacterial, short contig, no neighbours)"),
        ("integration_status", "integration_status", "Integration status of the locus in the fungal genome: confirmed, supported or uncertain (rules in the table legend)"),
        ("detection_route", "detection_route", "Route by which the locus was found: profile screen, targeted scan (miniprot), blastx of plasmid proteins, or protein-level screen"),
        ("detailed_in", "detailed_in", "Table with the per-locus details of the corresponding screen"),
        ("evidence_tables", "evidence_tables", "Supplementary tables in which the locus appears: by locus_id (S8, S14, S13a, S13b, S16), by contig and assembly (S9, S12, S10b, S5a; contig-level, "
         "so also listed for other loci on the same contig), for the odh island by gene (S13c, reverse search of the island genes; S16, RNA-seq of the Fp133 copies); blank = none"),
        ("protein_in_Data_S1", "protein_in_Data_S1", "yes if the translation of the locus is in Data S1 (no gene model, hence no protein, for the others)"),
        ("note", "note", "Free-text note"),
    ])
_s1a = TABLES.pop()
TABLES.insert([t[0] for t in TABLES].index("S1b"), _s1a)

# ---------------------------------------------------------------- renumbering and merges (2026-10-13q, suppl_restructure_plan.md section B)
# internal ids above are the old ones; the final numbering is assigned here. New <- old:
# S1a<-S1a S1b<-S1b S1c<-S18 S2 S3<-S3b S4<-S15a S5a<-S15b S5b<-S15c S6<-S4 S7<-new(324 assemblies; replaces S12b+S14b) S8<-S6 S9<-S7
# S10a<-S11a S10b<-S11b S11<-S5 S12<-S8 S13a<-S12a S13b<-S14a S13c<-S14c S14<-S10a+S10b S15<-S9a+S16 S16<-S13+S12c+S14d S17a,b<-S17a,b; S9b deleted
OLD = {t[0]: t for t in TABLES}
FINAL = []
def moved(nid, oid, title=None, desc=None):
    t = OLD[oid]
    FINAL.append((nid, title or t[1], desc or t[2], t[3], t[4]))

# ---- S7 (new): the 324 assemblies of Fusarium, Neocosmospora and Dactylonectria
f324 = rd(G + "/paper/tables/fusarium_324_assemblies.tsv")
assert len(f324) == 324, len(f324)
S7_SPEC = [
    ("assembly", "assembly", "Genome assembly accession"), ("organism", "organism", "Organism name (NCBI)"), ("species", "species", "Species name used in the paper"),
    ("group", "group", "Group of the assembly: F. oxysporum, other Fusarium, FSSC/Neocosmospora (Fusarium solani species complex) or Dactylonectria"),
    ("assembly_level", "assembly_level", "Assembly level"),
    ("chs_like_status", "chs_like_status", "State of the chs-like gene in the assembly: intact, truncated, pseudogene or absent (the best state if there are several loci)"),
    ("chs_like_loci", "chs_like_loci", "Number of chs-like loci (Table S13a)"),
    ("odh_status", "odh_status", "State of the odh gene of the odh island (intact, truncated, pseudogene or absent)"), ("odh_loci", "odh_loci", "Number of odh loci (Table S13b)"),
    ("pepP_status", "pepP_status", "State of the pepP gene of the odh island"), ("pepP_loci", "pepP_loci", "Number of pepP loci"),
    ("orf3_status", "orf3_status", "State of the orf3 gene of the odh island"), ("orf3_loci", "orf3_loci", "Number of orf3 loci"),
    ("total_loci", "total_loci", "Total number of loci of the four genes in the assembly")]
pr, cols = proj(f324, S7_SPEC)

# ---- S14: composition, loci + bacterial references
S14_SPEC = [
    ("entry_type", "entry_type", "'fungal locus' (intact CDS of a locus) or 'bacterial reference' (CDS of a Ti/Ri plasmid reference gene)"),
    ("id", "id", "Locus identifier (fungal loci) or protein accession (bacterial references)"),
    ("family", "family", "Protein family (fungal loci)"), ("lineage", "lineage", "Lineage (fungal loci)"),
    ("description", "description", "Bacterial references: gene, plasmid, GenBank record of the CDS and product"),
    ("organism", "organism", "Organism (fungal loci)"), ("order", "order", "Taxonomic order (fungal loci)"), ("assembly", "assembly", "Genome assembly accession (fungal loci)"),
    ("cds_length_nt", "cds_length_nt", "CDS length, nt"), ("cds_source", "cds_source", "Origin of the CDS sequence (fungal loci)"),
    ("gc", "gc", "GC content of the CDS"), ("gc3", "gc3", "GC content at third codon positions (not computed for the intact chs-like and odh-island rows)"),
    ("assembly_gc", "assembly_gc", "GC content of the whole assembly (given for the chs-like and odh-island rows only; the background columns are empty for these rows)"),
    ("background_source", "background_source", "Source of the background gene set (genes of the same genome)"),
    ("n_background_genes", "n_background_genes", "Number of background genes"),
    ("background_median_gc", "background_median_gc", "Background median GC"), ("background_median_gc3", "background_median_gc3", "Background median GC3"),
    ("delta_gc3", "delta_gc3", "GC3 of the locus minus background median GC3"), ("gc3_percentile", "gc3_percentile", "Percentile of GC3 among background genes"),
    ("gc_percentile", "gc_percentile", "Percentile of GC among background genes"),
    ("codon_usage_distance_percentile", "codon_usage_distance_percentile", "Percentile of the codon-usage distance (1 - r of RSCU to the genome) among length-matched background genes"),
    ("boletales_pool_delta_gc3", "boletales_pool_delta_gc3", "GC3 of the locus minus the median GC3 of the pooled Boletales gene set"),
    ("boletales_pool_gc3_percentile", "boletales_pool_gc3_percentile", "Percentile of GC3 in the pooled Boletales set"),
    ("boletales_pool_codon_usage_distance_percentile", "boletales_pool_codon_usage_distance_percentile", "Percentile of the codon-usage distance in the pooled Boletales set")]
comp = []
for r in OLD["S10a"][4]:
    o = dict(r); o["entry_type"] = "fungal locus"; o["id"] = r["locus_id"]; o["description"] = ""; comp.append(o)
for r in OLD["S10b"][4]:
    o = {k: "" for k in [c[0] for c in S14_SPEC]}
    o.update(entry_type="bacterial reference", id=r["reference_id"], cds_length_nt=r["cds_length_nt"], gc=r["gc"], gc3=r["gc3"],
             description=clean("; ".join(x for x in (r["description"], "GenBank record " + r["source_record"].replace(".gb", ""), r["note"]) if x)))
    comp.append(o)
n_fung = sum(1 for r in comp if r["entry_type"] == "fungal locus")
assert n_fung == len(OLD["S10a"][4])
pr14, cols14 = proj(comp, S14_SPEC)

# ---- S15: dN/dS of all lineages
d9 = OLD["S9a"][4]; d16 = OLD["S16"][4]
def loose_sf(fam):
    # <=20%-gap codeml M0 run of the sus-like / cus-like (all intact copies) sets (analysis/dnds/loose; same CDS, alignment and tree as the strict run, checked 2026-10-13)
    import re as _re
    w = _re.search(r"omega \(dN/dS\) =\s*([\d.]+)", open(A + "/dnds/loose/%s.M0.out" % fam).read()).group(1)
    return str(int(open(A + "/dnds/loose/%s.codon.phy" % fam).readline().split()[1]) // 3), w
LOOSE15 = {"OCD fold (Sus)": loose_sf("SUS"), "cucumopine synthase (Cus)": loose_sf("CUS")}
def dn_rows():
    out = []
    for r in d9:
        out.append({"gene_set": {"OCD fold (Sus)": "sus-like", "cucumopine synthase (Cus)": "cus-like (all intact copies)"}.get(r["family"], r["family"]),
                    "n": r["n_sequences"], "cs": r["codon_sites"], "cl": LOOSE15[r["family"]][0], "m0s": r["omega_M0"], "m0l": LOOSE15[r["family"]][1], "lrt": r["LRT_2dlnL_M7_vs_M8"], "p78": r["p_M7_vs_M8"],
                    "p8a": r["p_M8a_vs_M8"], "beb": r["positively_selected_sites_BEB_gt_0.95"], "yn": r["pairwise_median_dNdS_yn00"], "yn2": "", "gt2": "", "np": r["pairwise_comparisons_yn00"]})
    for r in d16:
        out.append({"gene_set": r["gene_set"], "n": r["n_CDS"], "cs": r["codons_strict"], "cl": r["codons_le20pct_gaps"], "m0s": r["M0_omega_strict"], "m0l": r["M0_omega_le20pct_gaps"],
                    "lrt": r["LRT_M7_vs_M8"], "p78": r["p_M7_vs_M8"], "p8a": r["p_M8a_vs_M8"], "beb": r["BEB_sites_P_gt_0.95"], "yn": r["yn00_strict_dNdS_median"],
                    "yn2": r["yn00_strict_dNdS_median_dS_le2"], "gt2": r["yn00_strict_pairs_dS_gt2"], "np": r["yn00_strict_pairs_used"]})
    return out
dn_all = dn_rows()
assert len(dn_all) == 16, len(dn_all)
S15_SPEC = [
    ("gene_set", "gene_set", "Gene set (intact copies; CDS of the Boletales nos-subfamily ODH from miniprot models)"), ("n", "n_sequences", "Number of CDS"),
    ("cs", "codons_strict", "Codons in the alignment without gap columns"), ("cl", "codons_le20pct_gaps", "Codons in the alignment with <=20% gaps per column"),
    ("m0s", "omega_M0_strict", "codeml M0 omega, strict alignment"), ("m0l", "omega_M0_le20pct_gaps", "codeml M0 omega, <=20% gaps alignment"),
    ("lrt", "LRT_M7_vs_M8", "Likelihood-ratio statistic, model M7 vs M8"), ("p78", "p_M7_vs_M8", "p-value, M7 vs M8"), ("p8a", "p_M8a_vs_M8", "p-value, M8a vs M8"),
    ("beb", "BEB_sites_P_gt_0.95", "Sites under positive selection (Bayes empirical Bayes posterior > 0.95) under M8"),
    ("np", "yn00_pairs_used", "Pairwise comparisons used in yn00 (strict alignment)"),
    ("yn", "yn00_median_dNdS", "Median pairwise dN/dS (yn00, strict alignment)"),
    ("yn2", "yn00_median_dNdS_dS_le2", "Median pairwise dN/dS over pairs with dS <= 2 (yn00, strict alignment; blank for sus-like and cus-like, not computed)"),
    ("gt2", "yn00_pairs_dS_gt2", "Pairs with dS > 2 (saturation; blank for sus-like and cus-like, not computed)")]
pr15, cols15 = proj(dn_all, S15_SPEC)

# ---- S16: RNA-seq at the loci (three datasets)
rows16 = []
for r in OLD["S13"][4]:
    rows16.append({"dataset": "Boletaceae loci", "sra_run": r["sra_run"], "assembly": r["assembly"], "locus_id": r["locus_id"], "family": r["family"], "lineage": r["lineage"],
                   "status": r["status"], "contig": r["contig"], "start": r["start"], "end": r["end"], "length_bp": r["length_bp"], "reads_on_locus": r["reads"],
                   "reads_per_kb": r["reads_per_kb"], "median_reads_per_kb_of_genes": r["median_reads_per_kb_of_genes"], "percentile_among_genes": r["percentile_among_genes"]})
for r in OLD["S12c"][4]:
    rows16.append({"dataset": "Fusarium chs-like locus", "sra_run": r["sra_run"], "organism": r["organism"], "strain": r["strain"], "bioproject": r["bioproject"], "biosample": r["biosample"],
                   "condition": r["condition"], "platform": r["platform"], "layout": r["layout"], "read_pairs_in_run": r["read_pairs_in_run"], "read_pairs_used": r["read_pairs_used"],
                   "mapped_pct": r["mapped_pct"], "lineage": "chs-like", "reads_on_locus": r["reads_on_locus"], "locus_TPM": r["locus_TPM"], "reads_per_kb": r["locus_reads_per_kb"],
                   "percentile_among_genes": r["pct_genes_with_lower_expression"], "locus_TPM_over_median_expressed": r["locus_TPM_over_median_expressed"],
                   "actin_TPM": r["actin_TPM"], "EF1a_TPM": r["EF1a_TPM"], "GAPDH_TPM": r["GAPDH_TPM"], "tubulin_beta_TPM": r["tubulin_beta_TPM"],
                   "sense_fraction": r["sense_fraction_R1"], "locus_breadth": r["locus_breadth"]})
for r in OLD["S14d"][4]:
    rows16.append({"dataset": "F. poae odh island", "sra_run": r["sra_run"], "assembly": FP133, "organism": "Fusarium poae", "strain": r["strain"], "bioproject": r["bioproject"], "biosample": r["biosample"],
                   "condition": r["condition"], "platform": r["platform"], "layout": r["layout"], "reads_in_run": r["reads_in_run"], "reads_used": r["reads_used"], "mapped_pct": r["mapped_pct"],
                   "locus_id": r["island_gene"], "locus_tag": r["locus_tag"], "lineage": "odh island: " + g3_disp(r["island_gene"]).split()[0],
                   "reads_on_locus": r["reads_on_gene"], "locus_TPM": r["gene_TPM"], "reads_per_kb": r["gene_reads_per_kb"], "percentile_among_genes": r["pct_genes_with_lower_expression"],
                   "locus_TPM_over_median_expressed": r["gene_TPM_over_median_expressed"], "actin_TPM": r["actin_TPM"], "EF1a_TPM": r["EF1a_TPM"], "GAPDH_TPM": r["GAPDH_TPM"],
                   "tubulin_beta_TPM": r["tubulin_beta_TPM"], "sense_fraction": r["sense_fraction"], "locus_breadth": r["gene_breadth"]})
S16_SPEC = [(k, k, d) for k, d in [
    ("dataset", "Dataset: 'Boletaceae loci' (public RNA-seq runs mapped to the loci of the corresponding assemblies), 'Fusarium chs-like locus' (reference strains F. vanettenii 77-13-4 and F. oxysporum f. sp. lycopersici 4287) or 'F. poae odh island' (strain Fp133)"),
    ("sra_run", "SRA run accession"), ("organism", "Organism as recorded in SRA (Fusarium datasets)"), ("strain", "Strain (Fusarium datasets)"),
    ("bioproject", "BioProject (Fusarium datasets)"), ("biosample", "BioSample (Fusarium datasets)"), ("condition", "Experimental condition (Fusarium datasets)"),
    ("platform", "Sequencing platform (Fusarium datasets)"), ("layout", "Library layout (Fusarium datasets)"),
    ("read_pairs_in_run", "Read pairs in the run (chs-like dataset, paired-end)"), ("read_pairs_used", "Read pairs used, subsample (chs-like dataset)"),
    ("reads_in_run", "Reads in the run (odh-island dataset)"), ("reads_used", "Reads used, subsample (odh-island dataset)"),
    ("mapped_pct", "Mapped reads, % (Fusarium datasets)"),
    ("assembly", "Genome assembly accession (Boletaceae and odh-island datasets)"), ("locus_id", "Locus identifier (Table S1b) or, for the odh island, gene and copy"),
    ("locus_tag", "Locus tag in the Fp133 annotation (odh-island dataset)"), ("family", "Protein family (Boletaceae dataset)"), ("lineage", "Lineage"),
    ("status", "Locus state (Boletaceae dataset)"), ("contig", "Contig (Boletaceae dataset)"), ("start", "Locus start (Boletaceae dataset)"), ("end", "Locus end (Boletaceae dataset)"),
    ("length_bp", "Locus length, bp (Boletaceae dataset)"), ("reads_on_locus", "Reads mapped to the locus or gene"), ("reads_per_kb", "Reads per kb on the locus or gene"),
    ("median_reads_per_kb_of_genes", "Median reads per kb of the genes of the genome (Boletaceae dataset)"),
    ("locus_TPM", "TPM of the locus or gene (Fusarium datasets)"),
    ("percentile_among_genes", "Boletaceae dataset: percentile of the locus among the genes of the genome; Fusarium datasets: percent of all genes with expression strictly below the locus or gene"),
    ("locus_TPM_over_median_expressed", "Locus TPM / median TPM of genes with reads (Fusarium datasets)"),
    ("actin_TPM", "TPM of actin (Fusarium datasets)"), ("EF1a_TPM", "TPM of EF1a (Fusarium datasets)"), ("GAPDH_TPM", "TPM of GAPDH (Fusarium datasets)"),
    ("tubulin_beta_TPM", "TPM of beta-tubulin (Fusarium datasets)"),
    ("sense_fraction", "Fraction of reads in sense orientation: R1 reads for the chs-like dataset; reads over genes with >=200 reads (library strandedness) for the odh-island dataset; the two are not comparable"),
    ("locus_breadth", "Fraction of the locus or gene covered by reads (Fusarium datasets)")]]
pr16 = [{k: r.get(k, "") for k, _, _ in S16_SPEC} for r in rows16]
assert len(pr16) == len(OLD["S13"][4]) + len(OLD["S12c"][4]) + len(OLD["S14d"][4])

# ---- assemble in the new order
moved("S1a", "S1a"); moved("S1b", "S1b"); moved("S1c", "S18", title="Summary of the catalogue by group"); moved("S1d", "S1d")
moved("S2", "S2"); moved("S3", "S3b"); moved("S4", "S15a"); moved("S5a", "S15b")
moved("S5b", "S15c"); moved("S6", "S4")
FINAL.append(("S7", "The 324 assemblies of Fusarium, Neocosmospora and Dactylonectria",
    "The 324 genome assemblies of Fusarium, Neocosmospora and Dactylonectria screened for the chs-like gene and for the three genes of the odh island (odh, pepP, orf3): "
    "one row per assembly with the state of each gene and the number of loci (%d assemblies carry a chs-like locus, %d carry genes of the odh island; loci in Tables S13a and S13b)."
    % (sum(r["chs_like_status"] != "absent" for r in f324), sum(any(r[k] != "absent" for k in ("odh_status", "pepP_status", "orf3_status")) for r in f324)), cols, pr))
moved("S8", "S6"); moved("S9", "S7"); moved("S10a", "S11a"); moved("S10b", "S11b"); moved("S11", "S5"); moved("S11b", "S5b"); moved("S12", "S8")
moved("S13a", "S12a"); moved("S13b", "S14a"); moved("S13c", "S14c")
FINAL.append(("S14", "Nucleotide composition of the loci and of bacterial reference genes",
    "GC and GC3 of the intact fungal locus CDSs (sus-like, cus-like, acc-type GDPD, pKerr14a-type and nos-subfamily loci of Boletales, two nos-subfamily pseudogenes included) compared with background genes of the same genome and with a pooled Boletales gene set; "
    "for the intact chs-like and odh-island genes of Fusarium only the GC of the gene and of the assembly is given. The last rows (entry_type = bacterial reference) give GC and GC3 of the CDS of the Ti/Ri plasmid reference genes. "
    "One row per locus or reference gene.", cols14, pr14))
FINAL.append(("S15", "dN/dS of the intact genes of all lineages",
    "Codon-model tests of selection (PAML codeml M0, M7/M8, M8a; pairwise yn00) on the intact copies of each lineage with at least two copies of a gene, on codon alignments (MAFFT, IQ-TREE tree). "
    "For the sus-like and cus-like rows (all intact copies) the <=20%-gap alignment was analysed with codeml M0 only (same sequences and tree as the strict run); the remaining rows were analysed on both alignments with all tests. "
    "cus-like is given for all intact copies, and separately for Boletales and for outside Boletales, each for all intact copies and for one copy per species (longest CDS of each species for the Boletales set). The four sets marked 'confirmed or supported only' exclude the loci whose integration status is uncertain (Table S1a); for the nos-subfamily set of three genes M0 was run on a star tree and no site models were fitted. One row per gene set.", cols15, pr15))
FINAL.append(("S16", "RNA-seq at the loci: Boletaceae, Fusarium chs-like locus and Fusarium poae odh island",
    "Expression of the loci in public RNA-seq runs, in three datasets (column dataset): Boletaceae (reads mapped to the loci of the corresponding assemblies), the chs-like locus in the Fusarium reference strains and the genes of the odh island in F. poae Fp133 "
    "(first 10 million reads of each run mapped to the Fp133 assembly). Columns that apply to only some datasets are marked; blank = not applicable. One row per run and locus.", [(c, d) for _, c, d in S16_SPEC], pr16))
moved("S17a", "S17a"); moved("S17b", "S17b")
# ---------------------------------------------------------------- S18 (2026-10-14, author): the per-lineage evidence matrix, formerly Table 3 of the main text
s18m = rd(G + "/paper/tables/TableS18_evidence_matrix.tsv")
cols18 = [("Test", "Row: criterion, test or summary statistic (see legend)")] + [(c, "Evidence for the %s lineage" % c) for c in list(s18m[0].keys())[1:]]
FINAL.append(("S18", "Evidence for each lineage by the three criteria of horizontal transfer and the tests against contamination",
    "The seven lineages of bacterial origin (cus-like split into Boletales and other orders): evidence by the three criteria of the Introduction (rows 1-5), the prerequisite that contamination is excluded (rows 6-7), the fate of the genes (row 8) and the closest sampled relatives (row 9); Table 3 of the main text gives the same lineages in words. "
    "Clusters with: the bacterial genes with which the fungal genes group in the tree of the family, with the support of the node that joins them (SH-aLRT/UFBoot, %; the support of the fungal clade itself is in Figure 3 of the main text); the alternative placements rejected by the AU test, and which constraints are informative, are in Tables S11 and S11b. For the pKerr14a-type lineage the grouping node is weakly supported and placement with the acc-type GDPD is not rejected (p = 0.415 and 0.132), so criterion 1 is met with low confidence. "
    "Patchy distribution: species with the gene among the species screened (358 species-level clusters of Boletales; 90 species of Fusarium and relatives; Tables S6, S7). Genomic context: a single gene among fungal genes or an island with genes of other acquired lineages (Table S12). Intact genes with introns: genes with at least one intron over all intact genes of the lineage (Table S1a). "
    "GC: median difference of the genes' GC from the carrier genome in percentage points (Boletales, from the median of the neighbouring genes; Fusarium, from the assembly) and the number of genes above the 95th percentile of the genome's genes (Table S14); met, at least half of the genes above that percentile or a difference above 5 points; partly, fewer than half; not met otherwise. "
    "Contamination excluded by: the tests passed by at least one locus of the lineage: F, fungal genes on both flanks; A, the locus at the same place with the same flanks in an independent assembly of the species (Table S9); S, a gene of the lineage present in several assemblies of the species, flanks not compared; C, the locus on an assembled chromosome; R, read depth equal to the genome and reads spanning the junctions (Table S8); P, no plasmid vir, tra or rep genes on the contig (Table S5b); a letter is omitted when the test was not done. "
    "Loci with integration confirmed / supported / not assessable: per-locus status (Table S1a, column integration_status; rules in Methods 2.5 of the main text). omega, dN/dS of the M0 model (Table S15; in parentheses, pairs with dS <= 2); transcription from RNA-seq (Table S16); n.d., not determined. "
    "Closest sampled relatives and locus context: affinities, not established sources of transfer. Numbers of loci, assemblies and species are in Table 2 of the main text.", cols18, s18m))

TABLES = FINAL
TNAME = {"S1a": "Table S1a"}
# ---------------------------------------------------------------- write
NAMES = {}
for tid, title, desc, cols, rows in TABLES:
    fn = "Table%s_%s.tsv" % (tid, re.sub(r"[^A-Za-z0-9]+", "_", title).strip("_")[:50])
    NAMES[tid] = fn
    names = [c for c, _ in cols]
    with open(os.path.join(OUT, fn), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(names)
        for r in rows:
            w.writerow([disp(r[c]) if isinstance(r[c], str) else r[c] for c in names])

# ---------------------------------------------------------------- Data S3: all screened assemblies
boletales_asm = rd(A + "/boletales_mash_517/genomes.tsv", header=False, fields=["gca", "organism", "genus", "path"])
other_asm = rd(G + "/paper/tables/other_orders_1516_assemblies.tsv")
ds3 = [{"accession": r["gca"], "organism": r["organism"], "order": "Boletales", "class": "Agaricomycetes", "set": "Boletales (Table S6)", "note": ""} for r in boletales_asm]
ds3 += [{"accession": r["accession"], "organism": r["organism"], "order": r["order"], "class": r["class"], "set": "other orders", "note": r["note"]} for r in other_asm]
ds3 += [{"accession": r["assembly"], "organism": r["organism"], "order": "Hypocreales", "class": "Sordariomycetes", "set": "Fusarium and related genera (Table S7)", "note": ""} for r in f324]
assert len(boletales_asm) == 517 and len(other_asm) == 1386 and len(ds3) == 2227, (len(boletales_asm), len(other_asm), len(ds3))
_acc = collections.Counter(r["accession"] for r in ds3)
assert len(_acc) == 2215 and sum(1 for v in _acc.values() if v > 1) == 12, (len(_acc), sum(1 for v in _acc.values() if v > 1))   # 12 assemblies of Fusarium and relatives are in both the 'other orders' and the 324 sets
for r in ds3:
    if _acc[r["accession"]] > 1: r["note"] = (r["note"] + "; " if r["note"] else "") + "also listed in the other set (Fusarium and related genera / other orders)"
DATA_S3 = "Data_S3_screened_assemblies.tsv"
with open(os.path.join(OUT, DATA_S3), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["accession", "organism", "order", "class", "set", "note"])
    for r in ds3: w.writerow([r[k] for k in ("accession", "organism", "order", "class", "set", "note")])

FIGS = [("S1", "Distributions of the statistics underlying the thresholds", "figures/FigureS1_thresholds_distributions.png",
         "The 11 thresholds of the classification of proteins, loci and contigs are shown in 12 panels (the vir/tra/rep threshold has two panels, k and l). Each panel gives the distribution of the underlying statistic "
         "(grey, the comparison group; colour, the group passing the threshold), the threshold (dashed line) and the range of threshold values over which the result is unchanged (grey band). "
         "(a) probable contamination; (b, c) bacterial-type protein; (d-f) T-DNA candidate; (g) non-T-DNA (rhizobial) candidate; (h, i) chs-like loci; (j) fragments assigned to the pKerr14a-type lineage; (k, l) vir/tra/rep genes on contigs. "
         "The numerical results of the sensitivity analysis are in Table S4."),
        ("S2", "Full tree of the OCD-fold family (Sus)", "figures/FigureS2_tree_sus.png",
         "The fungal sus-like genes lie within the clade of T-DNA sus genes. Full maximum-likelihood tree of the OCD-fold family (125 sequences, 434 sites, Q.PFAM+R5), from which panel (a) of Figure 3 was collapsed. "
         "Leaf colours: T-DNA sus clade (fungal sus-like genes with the Ti/Ri plasmid and plant cT-DNA references), Rhizobiales chromosomal OCD, other bacterial OCD, native fungal homologs, outgroup; key lineages in bold. "
         "Node labels, SH-aLRT/UFBoot (%), shown when UFBoot is at least 50; rooted as in Figure 3; scale bar, substitutions per site. Tree, alignment and IQ-TREE report in Data S2 (sus.*)."),
        ("S3", "Full tree of the cucumopine synthase family (Cus)", "figures/FigureS3_tree_cus.png",
         "The fungal cus-like genes lie within the T-DNA cus clade, interleaved with the references. Full maximum-likelihood tree of the cucumopine synthase family (436 sequences, 810 sites, JTT+R7), from which panel (b) of Figure 3 was collapsed: "
         "T-DNA cus clade (fungal cus-like genes of Boletales and of other orders with the Ti/Ri plasmid, plant cT-DNA and four other bacterial references), actinobacterial-type clade, mis group, other bacterial homologs, outgroup. "
         "In panel (b) of Figure 3 one native fungal gene detected by the cus profile (Boletus austroedulis; excluded from the analysis) is omitted from the schematic; it is labelled here as \"native fungal gene, excluded\". Supports, rooting and scale as in Figure S2; because of the size of the tree the labels are small and the vector files (PDF, in figures_pdf/) should be used for reading. Data S2: cus.*."),
        ("S4", "Full tree of the M24 aminopeptidase family (Chs, Ags)", "figures/FigureS4_tree_chs.png",
         "The Fusarium chs-like genes are sister to the Ti-plasmid chs clade and are not nested within it. Full maximum-likelihood tree of the M24 aminopeptidase family (58 sequences, 687 sites, Q.PFAM+I+G4), from which panel (c) of Figure 3 was collapsed: "
         "Fusarium and Dactylonectria chs-like genes, Ti-plasmid chs clade (the chs of pTi14-818c and other agrobacterial chs), native fungal Ags homologs, other bacterial M24 aminopeptidases; rooted on the ags seed. Supports and scale as in Figure S2. Data S2: chs.*."),
        ("S5", "Full tree of the GDPD family", "figures/FigureS5_tree_gdpd.png",
         "The acc-type and pKerr14a-type GDPD genes group with the GDPD of rhizobial acc-like loci and with the acs of pKerr14a, respectively, and not with the T-DNA acs clade. Full maximum-likelihood tree of the GDPD family (259 sequences, 989 sites, JTT+R5), from which panel (d) of Figure 3 was collapsed: "
         "T-DNA acs clade, acc-type GDPD with the bacterial acc-locus references, pKerr14a-type acs-like with the acs of pKerr14a and the Parasponia cT-DNA fragment, native fungal GDPD, fungal GDPD of unresolved position, outgroup. Supports, rooting and scale as in Figure S2. Data S2: gdpd.*."),
        ("S6", "Full tree of the opine dehydrogenase family, nos subfamily", "figures/FigureS6_tree_nos.png",
         "The nos-subfamily genes of Boletales form a clade sister to the bacterial genes of the subfamily; the Rhypophila homolog is shown but not assigned to any lineage. Full maximum-likelihood tree of the opine dehydrogenase family around the nos subfamily (64 sequences, 390 sites, Q.PFAM+R4), from which panel (e) of Figure 3 was collapsed: "
         "nos-subfamily genes of Boletales, the Rhypophila nos homolog, T-DNA nos clade, rhizobial nos-like genes, odh-island genes, native fungal and other bacterial opine dehydrogenases; rooted on the ocs/vis seeds. Supports and scale as in Figure S2. Data S2: nos.*."),
        ("S7", "Full tree of the opine dehydrogenase family, odh island", "figures/FigureS7_tree_odh.png",
         "The odh-island genes of Fusarium are sister to the plasmid opine dehydrogenase of Agrobacterium larrymoorei. Full maximum-likelihood tree of the opine dehydrogenase family around the odh island (44 sequences, 423 sites, Q.PFAM+G4), from which panel (f) of Figure 3 was collapsed: "
         "odh-island genes of F. poae and F. langsethiae, the plasmid opine dehydrogenase of A. larrymoorei (donor locus), nos subfamily, native fungal and other bacterial opine dehydrogenases; rooted on the ocs/vis seeds. Supports and scale as in Figure S2. Data S2: odh.*."),
        ("S8", "Approximately unbiased (AU) tests of the lineage hypotheses", "figures/FigureS8_AU_tests.png",
         "AU tests of the topological constraints for all lineages, drawn from Table S11. Each row is one constrained tree; the dot shows p-AU on a logarithmic scale (dashed line, p = 0.05); squares mark the unconstrained maximum-likelihood tree of each run; "
         "dlnL, log-likelihood difference from the best tree. Yellow, constraint not rejected; black, rejected. Hypotheses are worded as in Table S11 (Methods 2.6)."),
        ("S9", "Structural phylogeny of the GDPD family", "figures/FigureS9_structural_gdpd.png",
         "A tree built from predicted protein structures reproduces the sequence-based GDPD lineages. Maximum-likelihood tree built from the structural alignment (FoldMason, 3Di and amino-acid alphabets) of ColabFold models of the GDPD family (Methods 2.6). Shaded: the T-DNA acs references, the acc-type lineage (48 fungal genes with the acc-locus references) and the pKerr14a-type lineage (fungal genes with the acs of pKerr14a and the Parasponia cT-DNA fragment). "
         "Leaf colours: fungal loci (intact, pseudogene-like), T-DNA and plant cT-DNA references, bacterial references, native fungal homologs, outgroup."),
        ("S10", "All opine islands of Boletaceae", "figures/FigureS10_islands.png",
         "All opine islands of Boletaceae: the 13 contigs (7 assemblies of 6 species or species-level taxa) carrying two or more loci of bacterial type at distinct positions (Table S12), drawn to a common scale as in Figure 4a; "
         "each window spans the loci with 12 kb on each side (within the gene set of Table S12). Genes of the lineages of bacterial origin (sus-like, cus-like, acc-type GDPD, nos-subfamily opine dehydrogenase) are coloured as in Figure 4; loci of the nos-subfamily opine dehydrogenase lineage, found by blastx and therefore absent from Table S12, are included (Table S1b; contigs with such a locus and no other locus of bacterial type are not shown); "
         "other genes of bacterial type (mis group; other bacterial-type ORF) and the GDPD gene of unresolved position are shown in separate colours (legend); fungal neighbour genes grey; transposases and integrases black; dashed outline, pseudogene, fragment or remnant. "
         "Lineages and gene states are those of Table S1b. Arrows show the strand. Repeats are not drawn, because repeat annotation exists only for the contigs of Figure 4a. "
         "Contigs with a single locus are not shown. On JBBPKX010000014 the nos-subfamily gene (pseudogene) lies on the same strand as the intact mis-group gene, 3.9 kb downstream of it."),
        ("S11", "Pairwise dN/dS of the intact genes of each lineage", "figures/FigureS11_dnds.png",
         "Pairwise dN/dS (yn00) between intact genes of each lineage: (a) sus-like genes of Boletales (31 genes, 464 pairs); (b) all cus-like genes (27 genes, 347 pairs); "
         "(c) cus-like genes outside Boletales (20 genes; 59 pairs with dS <= 2 shown, 120 of 179 pairs excluded as saturated); (d) chs-like genes of Fusarium (6 genes, 13 pairs); "
         "(e) pKerr14a-type acs-like genes (6 genes, 15 pairs; M0 omega for the strict 99-codon and the 219-codon alignment); (f) nos-subfamily opine dehydrogenase genes of Boletales (4 genes, 5 pairs); "
         "(g) odh-island genes of Fusarium (5 loci, 7 pairs; unreliable, dS about 0.02). Histograms (a-d) or individual pairs (e-g); solid coloured line, median pairwise dN/dS; "
         "dotted line, omega of the M0 model (codeml); dashed black line, dN/dS = 1. The acc-type GDPD lineage (3 intact of 60 loci) was not tested. Values in Table S15."),
        ("S12", "The strain-specific insertion on contig JBBPRY010001778 of Boletus edulis", "figures/FigureS12_bedulis_insertion.png",
         "The strain-specific insertion on contig JBBPRY010001778 of Boletus edulis (7.4 kb; Figure 4b). Top, read depth along the contig (smoothed over 200 bp; dashed line, genome median 91x); the repeat-insertion junction at 1,313 bp (dotted line) is spanned by 19 reads (9 read pairs), which supports the continuity of the contig across this junction. "
         "Bottom, the gene model of the pKerr14a-type acs-like gene (five exons, purple), flanking repeats (hatched), and the 5.9-kb segment (1,313-7,202 bp; red) absent from five other assemblies of the species, with the flanks present in other assemblies (green). Read depth of all loci in Table S8; reproducibility in independent assemblies in Table S9."),
        ("S13", "Both chromosomal copies of the odh island of Fusarium poae Fp133", "figures/FigureS13_fp133_island_two_copies.png",
         "Both chromosomal copies of the odh island of F. poae Fp133 compared with the homologous locus of the A. larrymoorei plasmid NZ_JAVIZC010000002.1 (Figure 4c shows the chromosome 2 copy only). Top, plasmid locus (110.9-137.3 kb), where the cassette is bounded by IS3-, IS110- and IS66-family transposases; "
         "middle, chromosome 2 (CM100255.1, 9,450.2-9,476.6 kb; 346.9 kb to the chromosome end); bottom, chromosome 3 (CM100256.1, 8,410.4-8,436.8 kb, shown as the reverse complement; 217.6 kb to the chromosome end). Ribbons connect the homologous genes orf3, pepP and odh; dashed outline, pseudogene or fragment; arrows indicate strand. Scale bar, 5 kb.")]

# ---------------------------------------------------------------- Data S4: codeml / yn00 inputs and outputs of the dN/dS sets of Table S15 (2026-10-13al)
import tarfile, shutil
DATA_S4 = "Data_S4_dNdS_codeml.tar.gz"
DN = A + "/dnds"
# (directory, Table S15 gene set, file prefix, has_tree_and_codeml)
S4_SETS = [(DN, "sus-like", "SUS", True), (DN, "cus-like (all intact copies)", "CUS", True),
           (DN + "/extra/cus_nonBol", "cus-like, outside Boletales (all intact copies)", "cus_nonBol", True), (DN + "/extra/cus_nonBol_sp", "cus-like, outside Boletales (one copy per species)", "cus_nonBol_sp", True),
           (DN + "/extra/cus_Bol", "cus-like, Boletales (all intact copies)", "cus_Bol", True), (DN + "/extra/cus_Bol_sp", "cus-like, Boletales (one copy per species)", "cus_Bol_sp", True),
           (DN + "/extra/chs", "chs-like", "chs", True), (DN + "/extra/acs_pKerr14a", "pKerr14a-type acs-like", "acs_pKerr14a", True),
           (DN + "/extra/nosODH_Bol", "nos-subfamily opine dehydrogenase (Boletales)", "nosODH_Bol", True), (DN + "/extra/odh", "odh island: odh", "odh", True),
           (DN + "/extra/pepP", "odh island: pepP", "pepP", False), (DN + "/extra/g3", "odh island: orf3", "g3", False)] + [
          (DN + "/extra/" + k, SET17[k], k, True) for k in ("cus_Bol_conf", "cus_Bol_sp_conf", "acs_pKerr14a_conf", "nosODH_Bol_conf")]
S4_FILES = [("codon.phy", "codon alignment (PHYLIP, input of codeml and yn00)"), ("codon.fna", "codon alignment (FASTA)"), ("tree.nwk", "tree given to codeml (Newick)"),
            ("M0.ctl", "codeml control file, model M0 (one omega)"), ("M0.out", "codeml output, model M0"),
            ("M7M8.ctl", "codeml control file, site models M7 and M8"), ("M7M8.out", "codeml output, site models M7 and M8 (LRT, BEB sites)"),
            ("M8a.ctl", "codeml control file, site model M8a"), ("M8a.out", "codeml output, site model M8a"),
            ("yn00.ctl", "yn00 control file (pairwise dN/dS)"), ("yn00.out", "yn00 output (pairwise dN, dS, omega)"), ("tree.iqtree", "IQ-TREE report of the tree")]
s4_man = []
S4_TMP = os.path.join(OUT, "data", "_s4_tmp")
shutil.rmtree(S4_TMP, ignore_errors=True)
for d0, setname, pre, _t in S4_SETS:
    for sub, tag in (("", "strict"), ("/loose", "le20pct_gaps")):
        for suf, dsc in S4_FILES:
            f = "%s%s/%s.%s" % (d0, sub, pre, suf)
            if not os.path.exists(f): continue
            rel = "%s/%s/%s.%s" % (pre, tag, pre, suf)
            os.makedirs(os.path.dirname(os.path.join(S4_TMP, rel)), exist_ok=True)
            shutil.copy(f, os.path.join(S4_TMP, rel))
            s4_man.append((setname, rel, ("<=20%%-gap alignment: " if tag != "strict" else "strict alignment: ") + dsc))
assert len(s4_man) > 100 and {m[0] for m in s4_man} == {x[1] for x in S4_SETS}
with open(os.path.join(S4_TMP, "MANIFEST.tsv"), "w") as fh:
    fh.write("set\tfile\tdescription\n" + "".join("%s\t%s\t%s\n" % m for m in s4_man))
os.makedirs(os.path.join(OUT, "data"), exist_ok=True)
with tarfile.open(os.path.join(OUT, "data", DATA_S4), "w:gz") as tf:
    for fn in ["MANIFEST.tsv"] + sorted(m[1] for m in s4_man):
        ti = tf.gettarinfo(os.path.join(S4_TMP, fn), arcname=fn); ti.uid = ti.gid = 0; ti.uname = ti.gname = ""; ti.mtime = 0
        with open(os.path.join(S4_TMP, fn), "rb") as fh: tf.addfile(ti, fh)
shutil.rmtree(S4_TMP)
_s4mb = os.path.getsize(os.path.join(OUT, "data", DATA_S4)) / 1e6
print("Data S4: %d files, %.2f MB" % (len(s4_man) + 1, _s4mb)); assert _s4mb < 20
# ---------------------------------------------------------------- Extended Methods: concise submission version (markdown) -> Extended_Methods.docx
EXT_METHODS_SRC = OUT + "/src_extended_methods/Extended_Methods_submission.md"

def new_docx_document():
    """A4 document with 2.5 cm margins and Palatino 10 pt (shared by Supplementary_Materials.docx and Extended_Methods.docx)"""
    from docx import Document
    from docx.shared import Cm, Pt
    doc = Document()
    sec = doc.sections[0]; sec.page_width = Cm(21); sec.page_height = Cm(29.7)
    for m in ("left_margin", "right_margin", "top_margin", "bottom_margin"): setattr(sec, m, Cm(2.5))
    st = doc.styles["Normal"]; st.font.name = "Palatino Linotype"; st.font.size = Pt(10)
    return doc

_INLINE_RX = re.compile(r"(\*\*(?:[^*]|\*[^*]+\*)+\*\*|\*[^*\s][^*]*?\*|`[^`]+`)")
def add_inline_runs(p, text):
    """add text to paragraph p, rendering **bold**, *italic* and `code` (backticks removed)"""
    from docx.shared import Pt
    for part in _INLINE_RX.split(text):
        if not part: continue
        if part.startswith("**") and part.endswith("**") and len(part) > 4:
            for sub in _INLINE_RX.split(part[2:-2]):
                if not sub: continue
                if sub.startswith("*") and sub.endswith("*") and len(sub) > 2 and not sub.startswith("**"):
                    r = p.add_run(sub[1:-1]); r.italic = True
                else:
                    r = p.add_run(sub)
                r.bold = True
        elif part.startswith("*") and part.endswith("*") and len(part) > 2:
            p.add_run(part[1:-1]).italic = True
        elif part.startswith("`") and part.endswith("`"):
            r = p.add_run(part[1:-1]); r.font.name = "Courier New"; r.font.size = Pt(9)
        else:
            p.add_run(part)

def md_to_docx(md_path, docx_path):
    """Convert a simple markdown file to DOCX: # headings, paragraphs, '-'/'*' bullet lists, numbered lists ('1. ', kept as written),
    pipe tables (as Word tables, bold header row), **bold**, *italic*, `code`. Blank lines separate paragraphs."""
    from docx.shared import Cm
    with open(md_path, encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    doc = new_docx_document()
    para = []
    def flush():
        if para:
            add_inline_runs(doc.add_paragraph(), " ".join(x.strip() for x in para))
            para.clear()
    i = 0
    while i < len(lines):
        ln = lines[i]; st = ln.strip()
        m_head = re.match(r"^(#{1,6})\s+(.*)$", st)
        m_bul = re.match(r"^[-*+]\s+(.*)$", st)
        m_num = re.match(r"^(\d+)[.)]\s+(.*)$", st)
        if not st:
            flush()
        elif m_head:
            flush()
            h = doc.add_heading(level=min(len(m_head.group(1)), 4))
            add_inline_runs(h, m_head.group(2).strip())
        elif st.startswith("|"):
            flush()
            rows = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                cells = [c.strip() for c in lines[i].strip().strip("|").split("|")]
                if not all(re.fullmatch(r":?-{2,}:?", c) for c in cells if c) or not any(cells):
                    rows.append(cells)
                i += 1
            i -= 1
            ncol = max(len(r) for r in rows)
            tb = doc.add_table(rows=len(rows), cols=ncol); tb.style = "Table Grid"
            for ri, r in enumerate(rows):
                for ci in range(ncol):
                    cell = tb.cell(ri, ci); cp = cell.paragraphs[0]
                    txt = r[ci] if ci < len(r) else ""
                    if ri == 0: txt = "**" + txt.replace("**", "") + "**" if txt else txt
                    add_inline_runs(cp, txt)
            doc.add_paragraph()
        elif m_bul and not st.startswith("**"):
            flush()
            add_inline_runs(doc.add_paragraph(style="List Bullet"), m_bul.group(1))
        elif m_num:
            flush()
            p = doc.add_paragraph()
            p.paragraph_format.left_indent = Cm(0.9); p.paragraph_format.first_line_indent = Cm(-0.9)
            add_inline_runs(p, m_num.group(1) + ".\t" + m_num.group(2))
        elif re.fullmatch(r"[-*_]{3,}", st):
            flush()
        else:
            para.append(ln)
        i += 1
    flush()
    doc.save(docx_path)

if not os.path.exists(EXT_METHODS_SRC):
    raise FileNotFoundError("Extended Methods source not found: " + EXT_METHODS_SRC + " (the concise submission version of the Extended Methods must be placed there)")
md_to_docx(EXT_METHODS_SRC, os.path.join(OUT, "Extended_Methods.docx"))
print("Extended_Methods.docx written")
DATA = [("S1", "Proteins of the loci of bacterial origin", DATA_S1),
        ("S2", "Protein-family trees, alignments and IQ-TREE reports", "Data_S2_trees_alignments.tar.gz"),
        ("S2b", "Constraint trees and IQ-TREE outputs of the AU tests", "Data_S2b_AU_constraints.tar.gz"),
        ("S3", "Genome assemblies screened", DATA_S3),
        ("S4", "Codon alignments, trees and codeml/yn00 files of the dN/dS analyses", DATA_S4)]
DATA_DESC = {
 "S1": "Predicted protein sequences of the loci of Table S1a (the %d loci of the seven lineages: 239 placed on the trees directly or through a representative, 13 fragments assigned by identity, 2 added from NCBI annotations after the trees were built; the fragments assigned by profile only are not included). "
       "Header: `>locus_id lineage= status= species= assembly= coords=contig:start-end(strand)`; stop codons of pseudogenes are kept as `*`. "
       "Sources: gene models of the profile screen (%d; %s), the NCBI protein annotations of the two cus-like loci added after the trees were built (%d), miniprot translations of the chs-like and odh-island loci of the targeted scans (%d), and miniprot models of the nos-subfamily opine dehydrogenase of Boletales (%d). "
       "%s FASTA file (%d sequences)." % (
           len(s1a), src_cnt["models"] + src_cnt["remnants_fragments"], "%d models of intact genes and pseudogenes, %d translations of fragments and remnants" % (src_cnt["models"], src_cnt["remnants_fragments"]),
           src_cnt["ncbi_annotation"], src_cnt["miniprot"], src_cnt["nos_check"],
           ("%d loci have no sequence (column protein_in_Data_S1 = no of Table S1a; no gene model could be built): %s." % (n_noprot, "; ".join("%s, status %s: %d" % (k[0], k[1], v) for k, v in sorted(noprot.items()))) if n_noprot else "Every locus of Table S1a has a sequence."), n_prot),
 "S2": "For each of the protein-family trees of Figure 3 and Figures S2-S7 (sus, cus, chs, gdpd, nos, odh): the maximum-likelihood tree with SH-aLRT/UFBoot supports (`<family>.treefile`), the UFBoot consensus tree (`<family>.contree`), "
       "the trimmed protein alignment given to IQ-TREE (`<family>.aln.faa`) and the IQ-TREE report with model and run settings (`<family>.iqtree`). MANIFEST.tsv lists for each file the family, the number of sequences and sites and the substitution model. Compressed archive (tar.gz).",
 "S2b": "Constraint trees (Newick) and IQ-TREE outputs (.treefile, .iqtree, .log) of the 60 constrained maximum-likelihood searches behind Table S11, with the IQ-TREE AU reports (.iqtree, .log, .trees) and the user-tree sets of each run; "
        "MANIFEST.tsv gives for each constraint_id (as in Tables S11 and S11b) the files. Layout <family>/<run>/{constraint,search,run}. Compressed archive (tar.gz, 2.1 MB).",
 "S4": "For each of the 16 gene sets of Table S15: the codon alignment (`.codon.phy`, `.codon.fna`), the tree given to codeml (`.tree.nwk`), the codeml control files and outputs (M0; site models M7/M8 and M8a where run) and the yn00 control file and output (pairwise dN/dS), in two sub-folders per set: `strict` (alignment used for the main results) and `le20pct_gaps` (alignment keeping columns with up to 20% gaps; M0 and yn00 only). "
       "The sets odh island: pepP and orf3 have only yn00 files (2 and 3 sequences, no tree). The set nos-subfamily ODH, Boletales, confirmed or supported only has M0 on a star tree and no site models. MANIFEST.tsv lists set, file and description. Compressed archive (tar.gz).",
 "S3": "The 2,215 distinct genome assemblies screened, 2,227 rows because 12 assemblies of Fusarium relatives belong to both the 1,386 and the 324 set: the 517 assemblies of Boletales (Table S6), 1,386 assemblies of other fungal orders (accession, organism, order, class; a note where the order was inferred) and the 324 assemblies of Fusarium, Neocosmospora and Dactylonectria (Table S7), "
       "in one table with the column set identifying the group (12 assemblies of Fusarium-related genera were included in both the screen of other orders and the Fusarium set and appear twice; 2,215 distinct assemblies). Tab-separated file."}

def table_head(tid):
    return "Table %s" % tid

with open(os.path.join(OUT, "README_supplementary.md"), "w") as fh:
    fh.write("# Supplementary Materials: Tables, Figures and Data\n\nTables S1a-S18 are in `Supplementary_Tables.xlsx` (one sheet per table; each sheet starts with the table title and legend, then the column header and the data; blank = not available). "
             "This archive contains: `Supplementary_Tables.xlsx` (Tables S1a-S18); `Supplementary_Materials.docx` (the legends of all figures, tables and data files); `Extended_Methods.docx` (Extended Methods, the full text of the Methods); `figures/` (Figures S1-S13 as PNG) and `figures_pdf/` (the same figures as PDF); the data files Data S1, S2, S2b, S3 and S4, described at the end of this file; and this file. "
             "Nomenclature: families OCD fold (Sus), cucumopine synthase (Cus), GDPD, M24 aminopeptidase (Chs/Ags), opine dehydrogenase; "
             "lineages of bacterial origin (Table 2): sus-like, cus-like, chs-like, acc-type GDPD, pKerr14a-type acs-like, nos-subfamily opine dehydrogenase (Boletales), odh island; "
             "other groups of the catalogue (Table S1c): mis group, actinobacterial-type, native fungal GDPD, GDPD of unresolved position, OCD-like outside the sus-like clade.\n\n## Supplementary Tables\n\n")
    for tid, title, desc, cols, rows in TABLES:
        fh.write("### Table %s. %s\n\nSheet `%s` of Supplementary_Tables.xlsx (%d rows, %d columns)\n\n%s\n\n| Column | Description |\n|---|---|\n" %
                 (tid, title, tid, len(rows), len(cols), disp(desc)))
        for c, d in cols:
            fh.write("| %s | %s |\n" % (c, disp(d)))
        fh.write("\n")
    fh.write("## Supplementary Figures\n\n")
    for fid, ftitle, ffile, fcap in FIGS:
        fh.write("### Figure %s. %s\n\nFiles: `%s` (PDF: `%s`)\n\n%s\n\n" % (fid, ftitle, ffile, ffile.replace("figures/", "figures_pdf/").replace(".png", ".pdf"), fcap))
    fh.write("## Supplementary Data\n\n")
    for did, dtitle, dfile in DATA:
        fh.write("### Data %s. %s\n\nFile: `%s`\n\n%s\n\n" % (did, dtitle, dfile, DATA_DESC[did]))

# ---------------------------------------------------------------- xlsx (README sheet + one sheet per table; row 1 title, row 2 legend, row 3 header)
try:
    import openpyxl
    from openpyxl.styles import Font
    wb = openpyxl.Workbook(); ws = wb.active; ws.title = "README"
    ws.append(["Supplementary Tables S1a-S18. Each sheet starts with the table title (row 1) and legend (row 2), then the column header and the data. Figures S1-S13 and Data S1-S4: see Supplementary_Materials.docx."])
    ws.append([]); ws.append(["Table", "Title", "Legend", "Rows"])
    for tid, title, desc, cols, rows in TABLES:
        ws.append([tid, title, disp(desc), len(rows)])
    ws.append([]); ws.append(["Table", "Column", "Description"])
    for tid, title, desc, cols, rows in TABLES:
        for c, d in cols:
            ws.append([tid, c, disp(d)])
    for tid, title, desc, cols, rows in TABLES:
        ws = wb.create_sheet(tid)
        names = [c for c, _ in cols]
        ws.append(["Table %s. %s" % (tid, title)]); ws["A1"].font = Font(bold=True)
        ws.append([disp(desc)])
        ws.append(names)
        for c in ws[3]: c.font = Font(bold=True)
        for r in rows:
            ws.append([disp(r[c]) if isinstance(r[c], str) else r[c] for c in names])
        ws.freeze_panes = "A4"
    wb.save(os.path.join(OUT, "Supplementary_Tables.xlsx"))
    print("xlsx written")
except ImportError:
    print("openpyxl not available: xlsx skipped")

# ---------------------------------------------------------------- Supplementary_Materials.md / .docx
TITLE = "Opine metabolism genes horizontally transferred into fungal genomes"
md = ["# Supplementary Materials for: %s" % TITLE, "",
      "This document contains the legends of Figures S1-S13, Tables S1a-S18 and Data S1-S4 and S2b. The tables are in Supplementary_Tables.xlsx (one sheet per table), the figures in the folders figures/ (PNG) and figures_pdf/ (PDF), the data files as named below. The full text of the Methods is given in Extended Methods (Extended_Methods.docx).", "", "## Figures", ""]
for fid, ftitle, ffile, fcap in FIGS:
    md += ["![Figure S%s](%s)" % (fid[1:], ffile), "", "**Figure %s. %s.** %s" % (fid, ftitle, fcap), ""]
md += ["## Tables", ""]
for tid, title, desc, cols, rows in TABLES:
    md += ["**Table %s. %s.** %s Columns: %s. See Supplementary_Tables.xlsx, sheet %s." % (tid, title, disp(desc), "; ".join("%s, %s" % (c, disp(d)) for c, d in cols), tid), ""]
md += ["## Data", ""]
for did, dtitle, dfile in DATA:
    md += ["**Data %s. %s.** File: %s. %s" % (did, dtitle, dfile, DATA_DESC[did]), ""]
open(os.path.join(OUT, "Supplementary_Materials.md"), "w").write("\n".join(md) + "\n")
try:
    from docx.shared import Cm
    doc = new_docx_document()
    def par(text, bold_lead=None):
        p = doc.add_paragraph()
        if bold_lead:
            r = p.add_run(bold_lead + " "); r.bold = True
        p.add_run(text.replace("`", ""))
        return p
    doc.add_heading("Supplementary Materials for: %s" % TITLE, level=1)
    par(md[2])
    doc.add_heading("Figures", level=2)
    for fid, ftitle, ffile, fcap in FIGS:
        doc.add_picture(os.path.join(OUT, ffile), width=Cm(16))
        par(fcap, "Figure %s. %s." % (fid, ftitle))
    doc.add_heading("Tables", level=2)
    for tid, title, desc, cols, rows in TABLES:
        par(disp(desc), "Table %s. %s." % (tid, title))
        par("; ".join("%s, %s" % (c, disp(d)) for c, d in cols) + ".", "Columns:")
        par("See Supplementary_Tables.xlsx, sheet %s." % tid)
    doc.add_heading("Data", level=2)
    for did, dtitle, dfile in DATA:
        par("File: %s. %s" % (dfile, DATA_DESC[did]), "Data %s. %s." % (did, dtitle))
    doc.save(os.path.join(OUT, "Supplementary_Materials.docx"))
    print("docx written")
except ImportError:
    print("python-docx not available: docx skipped")

# ---------------------------------------------------------------- paragraph for the article
SHORT_DATA = {"S1": "protein sequences of the loci of bacterial origin", "S2": "protein-family trees, alignments and IQ-TREE reports", "S3": "genome assemblies screened", "S2b": "constraint trees and IQ-TREE outputs of the AU tests", "S4": "codon alignments, trees and codeml/yn00 files of the dN/dS analyses"}
def low(t): return t if t[:2].isupper() or t.startswith(("chs", "vir")) else t[0].lower() + t[1:]
PARA_T = {"S16": "RNA-seq at the loci of Boletaceae and Fusarium"}
items = ["Figure S%s: %s" % (f[0][1:], low(f[1])) for f in FIGS] + ["Table S%s: %s" % (t[0][1:], PARA_T.get(t[0]) or low(t[1]).replace(": ", ", ")) for t in TABLES] + ["Data S%s: %s" % (d[0][1:], SHORT_DATA[d[0]]) for d in DATA]
para = "**Supplementary Materials:** The following supporting information can be downloaded at: https://www.mdpi.com/xxx/s1, " + "; ".join(items) + "."
os.makedirs(G + "/paper/submission", exist_ok=True)
open(G + "/paper/submission/supplementary_materials_paragraph.md", "w").write(para + "\n")

# ---------------------------------------------------------------- checks
BAD = [("cyrillic", re.compile("[\u0400-\u04ff]")), ("gdpd_meso", re.compile("gdpd_meso")), ("kerr14", re.compile("kerr14")),
       ("six", re.compile(r"\bsix\b")), ("Tree_clade", re.compile("Tree_clade")),
       ("Acs1", re.compile(r"\bAcs\d")), ("Cus1", re.compile(r"\bCus\d")), ("ORF358", re.compile("ORF358|orf358")),
       ("T1|/T2|", re.compile(r"T[123]\|")), ("tier T1/T2/T3", re.compile(r"\bT[123]\b")),
       ("pseudogene_like", re.compile("pseudogene_like")), ("cus_actinobacteria", re.compile("cus_actinobacteria")),
       ("backbone", re.compile("backbone")), ("rule A/B/C", re.compile(r"\brule[ -]?[ABC]\b")), ("acs-HMM", re.compile("acs-HMM")),
       ("not placed", re.compile(r"\(not placed\)")), ("mis lineage", re.compile("mis lineage")), ("species_code", re.compile("species_code")),
       ("path", re.compile(r"analysis/|paper/|\.py\b|\.tsv\b")), ("date", re.compile(r"2026-\d\d-\d\d")), ("Supplementary Table/Figure", re.compile("Supplementary (Table|Figure) ?S")),
       ("1,516/2,357", re.compile("1,516|2,357"))]
viol = []
files = sorted(glob.glob(OUT + "/Table*.tsv")) + [OUT + "/README_supplementary.md", OUT + "/Supplementary_Materials.md"]
for f in files:
    for i, line in enumerate(open(f), 1):
        for nm, rx in BAD:
            if nm == "path" and not f.endswith(".tsv"): continue
            if nm == "date" and f.endswith(".tsv") and "Boletales_assemblies" in f: continue
            m = rx.search(line)
            if m:
                viol.append((os.path.basename(f), i, nm, line.strip()[:120]))
print("violations:", len(viol))
cnt = collections.Counter((v[0], v[2]) for v in viol)
for k, n in sorted(cnt.items()):
    print("  ", k, n)
for v in viol[:40]:
    print("   e.g.", v)
for tid, title, desc, cols, rows in TABLES:
    print(NAMES[tid], len(rows), "rows", len(cols), "cols")
