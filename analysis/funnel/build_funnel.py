#!/usr/bin/env python3
"""Funnel of the protein-level screen (Table 1 rebuilt) and fate of every protein with non-fungal nearest relatives (Table S3b rebuilt).
Input: refset/screen/seeds_ext/reciprocal/reciprocal.tsv (reverse search vs nr), seeds_ext_best.tsv (best hit per genome and gene).
Rules (all on hits OUTSIDE the query's own order):
  contamination : best bacterial hit >= 90% identity
  A  bacterial-type protein : >=35 of the 50 best hits are non-fungal AND bit-score margin (best bacterial/plant - best fungal) >= 50
  B  T-DNA-type candidate   : A and identity to the T-DNA seed >= 35% (coverage >= 0.4) and >= 0.9 x identity to the best bacterial hit; B and C are propagated within gene x genus
  C  other-rhizobial-locus candidate : A, not B, best bacterial hit from Rhizobiales (Hyphomicrobiales) with >= 50% identity
Outputs: paper/tables/Table1_funnel.tsv, paper/tables/TableS3b_fate.tsv, analysis/funnel/funnel_numbers.txt
"""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import pandas as pd, re, collections
R = f"{GMO_ROOT}/analysis/refset/screen/seeds_ext/"
OUT = f"{GMO_ROOT}/paper/tables/"
RHIZ = re.compile(r"^(Agrobacterium|Rhizobium|Allorhizobium|Neorhizobium|Pararhizobium|Sinorhizobium|Ensifer|Mesorhizobium|Bradyrhizobium|Devosia|Labrys|Bosea|Roseibium|"
                  r"Phyllobacterium|Ochrobactrum|Brucella|Methylobacterium|Aminobacter|Shinella|Gillisella|Rhizobiaceae|Hyphomicrobiales|Rhizobiales|Xanthobacter|Azorhizobium|"
                  r"Nitrobacter|Rhodopseudomonas|Pseudorhizobium|Chelativorans|Pelagibacterium|Rhizobium/Agrobacterium|Afipia|Microvirga|Hoeflea|Aquamicrobium|Nitratireductor)", re.I)
# fixed gene order shared by all tables/figures of the paper (v2): genes with T-DNA-type lineages; genes with lineages from other rhizobial loci; genes without findings beyond rule A
GENES = ["sus", "cus", "chs", "acs", "acs_pKerr14a", "nos", "ocs", "vis", "ags", "mas1", "mas2", "mis"]
def gene_rank(g):  # multi-gene labels ("acs,acs_pKerr14a") go by their first gene in GENES order, after the single-gene label
    gs = [x.strip() for x in str(g).split(",")]
    return (min(GENES.index(x) for x in gs), len(gs))
r = pd.read_csv(R + "reciprocal/reciprocal.tsv", sep="\t")
best = pd.read_csv(R + "seeds_ext_best.tsv", sep="\t")
best = best[best.dup_dropped == 0]  # one assembly per GenBank/RefSeq pair, as in Fig. 1b
r["nonfungal_nearest"] = r.o_h > 0
r["contamination"] = r.nonfungal_nearest & (r.o_bact_pid >= 90)
r["ruleA"] = r.nonfungal_nearest & ~r.contamination & (r.o_top50_bact >= 35) & (r.o_h >= 50)
r["genus"] = r.organism.str.split().str[0]
r["ruleB"] = r.ruleA & (r.pident_seed >= 35) & (r.cov_seed >= 0.4) & (r.pident_seed >= 0.9 * r.o_bact_pid)
r["ruleC"] = r.ruleA & ~r.ruleB & (r.o_bact_pid >= 50) & r.o_bact_org.fillna("").apply(lambda o: bool(RHIZ.match(o)))
# rules B and C are applied per gene x genus: if one protein of a group passes, all bacterial-type proteins of the group follow it (partial/annotation variants of one lineage)
gB = r[r.ruleB].groupby(["gene", "genus"]).size().index; gC = r[r.ruleC].groupby(["gene", "genus"]).size().index
key = list(zip(r.gene, r.genus))
r["ruleB"] = r.ruleA & pd.Series([k in gB for k in key], index=r.index)
r["ruleC"] = r.ruleA & ~r.ruleB & pd.Series([k in gC for k in key], index=r.index)
def fate(x):
    if not x.nonfungal_nearest: return "fungal nearest relatives"
    if x.contamination: return "probable contamination"
    if not x.ruleA: return "marginal: member of a fungal family with a single better bacterial hit"
    if x.ruleB: return "T-DNA-type candidate (rule B): genome-level screen"
    if x.ruleC: return "other-rhizobial-locus candidate (rule C): targeted analysis"
    return "bacterial-type protein of an unrelated family, single species: not analysed further"
r["fate"] = r.apply(fate, axis=1)
# Post-hoc fate of the Pleurotus GDPD protein (rule C): in the GDPD tree it falls into none of the lineages (journal 2026-10-12 pleurotus_verdict)
FATE_OVERRIDE = {("Pleurotus", "acs,acs_pKerr14a"): "other-rhizobial-locus candidate (rule C): GDPD homolog; not assigned to any lineage (not counted)",
                 ("Rhypophila", "nos"): "T-DNA-type candidate (rule B): nos homolog; not assigned to any lineage (not counted)"}  # Rhypophila: base of the nos subfamily, rule-B threshold (2026-10-12)
r["fate"] = [FATE_OVERRIDE.get((g, ge), f) if (rc or rb) else f for g, ge, f, rc, rb in zip(r.genus, r.gene, r.fate, r.ruleC, r.ruleB)]
# ---- per-gene funnel (a protein detected with several seeds is counted under each gene it was detected with, as in the screen)
# Table 1 columns (2026-10-12aa, user decision: no rule labels, no 'tested'/'bacterial type' columns): candidates taken forward with hosts grouped
# as "Phylum: Order (genera)"; everything else with non-fungal nearest relatives = "not pursued" (contamination + marginal + bacterial-type non-rhizobial).
PHYLUM = {"Agaricomycetes": "Basidiomycota", "Tremellomycetes": "Basidiomycota", "Ustilaginomycetes": "Basidiomycota", "Pucciniomycetes": "Basidiomycota",
          "Sordariomycetes": "Ascomycota", "Saccharomycetes": "Ascomycota", "Dothideomycetes": "Ascomycota", "Eurotiomycetes": "Ascomycota",
          "Leotiomycetes": "Ascomycota", "Pezizomycetes": "Ascomycota", "Lecanoromycetes": "Ascomycota", "Orbiliomycetes": "Ascomycota",
          "Schizosaccharomycetes": "Ascomycota", "Taphrinomycetes": "Ascomycota", "Mortierellomycetes": "Mortierellomycota", "Mucoromycetes": "Mucoromycota"}
unknown_classes = set()
def hosts(sub):
    """'n: Phylum: Order (*Genus*, ...); Order (...); Phylum: ...' ordered by phylum (Basidiomycota first), order, genus"""
    n = len(sub)
    if n == 0: return "0"
    grp = collections.OrderedDict()
    for cl, od, ge in zip(sub.own_class, sub.own_order, sub.genus):
        ph = PHYLUM.get(cl)
        if ph is None: unknown_classes.add(cl); ph = cl
        grp.setdefault((ph, od), []).append(ge)
    bylab = collections.OrderedDict()  # phylum once per group: "Phylum: Order (genera); Order (genera)"
    for (ph, od), gs in sorted(grp.items(), key=lambda kv: (kv[0][0] != "Basidiomycota", kv[0][0], kv[0][1])):
        seen = []
        for g in gs:
            if g not in seen: seen.append(g)
        bylab.setdefault(ph, []).append(f"{od} ({', '.join('*' + g + '*' for g in seen)})")
    return f"{n}: " + "; ".join(f"{ph}: " + "; ".join(ods) for ph, ods in bylab.items())
rows = []
for g in GENES:
    sub = r[r.gene.str.split(",").apply(lambda s: g in s)]
    b = best[best.gene == g]
    nf = int(sub.nonfungal_nearest.sum()); nB = int(sub.ruleB.sum()); nC = int(sub.ruleC.sum())
    # 2026-10-12ab (user): forward/reverse search made explicit in the headers, no colons; 'not pursued' dropped (composition in the caption)
    # 2026-10-13o (audit A4): median and maximum identity of the best hit to the seed (one best hit per proteome, GenBank/RefSeq duplicates removed; replaces Table S3a)
    rows.append({"gene": g, "proteomes with a forward-search hit": len(b),
                 "median identity of best hits to the seed, %": round(float(b.pident.median()), 1), "maximum identity of best hits to the seed, %": round(float(b.pident.max()), 1),
                 "proteins checked by reverse search": len(sub),
                 "proteins with non-fungal closest relatives": nf, "T-DNA candidates": hosts(sub[sub.ruleB]), "non-T-DNA candidates": hosts(sub[sub.ruleC])})
pd.DataFrame(rows).to_csv(OUT + "Table1_funnel.tsv", sep="\t", index=False)
if unknown_classes: print("WARNING classes without phylum mapping:", unknown_classes)
# ---- S3b: every query of the reverse search (2026-10-13al: all 1,277 rows; column nonfungal_nearest marks the 116 with h > 0)
r["nonfungal_nearest_yn"] = r.nonfungal_nearest.map({True: "yes", False: "no"})
cols = ["protein", "gca", "organism", "own_class", "own_order", "gene", "pident_seed", "cov_seed", "o_bact_org", "o_bact_pid", "o_bact_bits", "o_plant_org", "o_plant_pid",
        "o_fung_org", "o_fung_pid", "o_fung_bits", "o_h", "nonfungal_nearest_yn", "o_top50_bact", "o_top50_plant", "o_top50_fung", "n_fung_genera_all", "contamination", "ruleA", "ruleB", "ruleC", "fate"]
s3 = r.copy(); s3["_rank"] = s3.gene.map(gene_rank)
s3 = s3.sort_values(["_rank", "organism", "fate"], kind="mergesort")[cols].rename(columns={"nonfungal_nearest_yn": "nonfungal_nearest"})  # by gene (GENES order), then by genus/organism
s3.to_csv(OUT + "TableS3b_fate.tsv", sep="\t", index=False)
with open(f"{GMO_ROOT}/analysis/funnel/funnel_numbers.txt", "w") as f:
    f.write(f"proteins tested (unique): {len(r)}\nnon-fungal nearest: {int(r.nonfungal_nearest.sum())}\nprobable contamination: {int(r.contamination.sum())}\n"
            f"rule A: {int(r.ruleA.sum())}\nrule B: {int(r.ruleB.sum())}\nrule C: {int(r.ruleC.sum())}\n"
            f"A not B/C: {int((r.ruleA & ~r.ruleB & ~r.ruleC).sum())}\nmarginal: {int((r.nonfungal_nearest & ~r.contamination & ~r.ruleA).sum())}\n\n")
    for k, v in r.groupby("fate").size().items(): f.write(f"{v:4d}  {k}\n")
    f.write("\nrule B by gene/genus:\n" + r[r.ruleB].groupby(["gene", "genus"]).size().to_string() + "\n")
    f.write("\nrule C by gene/genus:\n" + r[r.ruleC].groupby(["gene", "genus"]).size().to_string() + "\n")
    f.write("\nA, not pursued, by gene/genus:\n" + r[r.ruleA & ~r.ruleB & ~r.ruleC].groupby(["gene", "genus"]).size().to_string() + "\n")
print(open(f"{GMO_ROOT}/analysis/funnel/funnel_numbers.txt").read())
