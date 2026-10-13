#!/usr/bin/env python3
"""Unified catalogue of loci of ALL lineages found in fungal genomes (Table S1 v2 = TableS1_all_lineages.tsv) and its per-lineage summary (Table2a_lineages.tsv).
Sources (all existing tables; nothing recomputed):
  profile screen (HMM/blastx, 517 Boletales + 1,386 other assemblies): paper/tables/TableS1_loci_catalog.tsv (families GDPD, Cus, Sus)
  targeted scan of 324 Fusarium/Neocosmospora/Dactylonectria genomes (miniprot): analysis/chs_fusarium/stage2/loci.tsv (chs-like), analysis/fusarium_odh_island/scan/loci.tsv (odh island)
  blastx of 388 plasmid proteins on the 517 Boletales assemblies: analysis/refset/screen/step2/boletales517_tdna_summary/opine_dehydrogenase_loci.tsv (nos-subfamily ODH)
  annotated proteomes (protein-level screen): Rhypophila nos-subfamily proteins (reverse search table)
"""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import csv, collections, re
G = f"{GMO_ROOT}"; T = G + "/paper/tables/"; A = G + "/analysis/"
def rd(p): return list(csv.DictReader(open(p), delimiter="\t"))
COLS = ["lineage", "family", "origin", "entry_route", "locus_id", "gca", "organism", "species", "order", "contig", "contig_len", "start", "end", "strand", "status", "n_exons", "protein_len", "note", "source_table"]
ORIGIN = {"sus-like (T-DNA type)": "T-DNA (sus)", "cus-like (T-DNA type)": "T-DNA (cus)", "acc-type GDPD": "non-T-DNA rhizobial locus (acc-like catabolic locus)",
          "pKerr14a-type acs-like": "non-T-DNA rhizobial locus (Ti-type plasmid pKerr14a)", "mis lineage": "mis family (not analysed)",
          "actinobacterial-type": "non-rhizobial bacterial homologs (not analysed)", "other bacterial-type CUS": "bacterial homologs (not analysed)",
          "native fungal GDPD (main clade)": "native fungal family", "other fungal GDPD (position unresolved)": "native fungal GDPD family (not analysed)",
          "OCD-like, outside T-DNA sus clade": "OCD-fold homologs outside the sus-like clade (not analysed)", "not a homolog": "not a homolog of the lineages", "T-DNA acs clade": "T-DNA (acs)", "nos-subfamily ODH homolog (Rhypophila, not assigned)": "not assigned to any lineage (not counted)"}
FAM = {"ACS": "GDPD", "CUS": "cucumopine synthase (Cus)", "SUS": "OCD fold (Sus)"}
STATUS = {"pseudogene_like": "pseudogene", "fragment": "fragment", "remnant": "remnant", "intact": "intact", "ambiguous": "ambiguous", "not_homolog": "not a homolog", "not_run": "not modelled"}
out = []
# ---- 1. profile screen
for r in rd(T + "TableS1_loci_catalog.tsv"):
    lin = r["lineage_final"]; base = lin.replace(" (not placed)", "")
    out.append(dict(lineage=lin, family=FAM[r["family"]], origin=ORIGIN.get(base, "not placed on a tree (not analysed)"), entry_route="profile screen (HMM/blastx) of 517 Boletales + 1,386 other assemblies",
                    locus_id=r["locus_id"], gca=r["gca"], organism=r["organism"].replace("_", " "), species=r["species"], order=r["order"], contig=r["contig"], contig_len=r["contig_len"],
                    start=r["start"], end=r["end"], strand=r["strand"], status=STATUS.get(r["status"], r["status"]), n_exons=r["n_exons"], protein_len=r["protein_len"], note=r["note"], source_table="Table S1 (profile screen)"))
# ---- 2. chs-like loci (Fusarium scan)
for r in rd(A + "chs_fusarium/stage2/loci.tsv"):
    if r["cls"] != "chs_like": continue
    out.append(dict(lineage="chs-like (T-DNA type)", family="M24 aminopeptidase (Chs/Ags)", origin="T-DNA (chs)", entry_route="targeted scan (miniprot) of 324 Fusarium-related genomes",
                    locus_id=r["locus_id"], gca=r["gca"], organism=r["organism"], species=r["species"], order="Hypocreales", contig=r["contig"], contig_len=r["contig_len"],
                    start=r["start"], end=r["end"], strand=r["strand"], status=r["status"], n_exons=r["exons"], protein_len=int(r["prim_qto"]) - int(r["prim_qfrom"]) + 1, note=r["note"], source_table="Table S12a"))
# ---- 3. odh island (Fusarium scan)
GN = {"odh": "odh island: opine dehydrogenase", "pepP": "odh island: Xaa-Pro aminopeptidase", "g3": "odh island: third gene"}
for r in rd(A + "fusarium_odh_island/scan/loci.tsv"):
    out.append(dict(lineage=GN[r["gene"]], family="opine dehydrogenase (odh) / M24 aminopeptidase / hypothetical", origin="non-T-DNA rhizobial locus (opine-catabolic plasmid of Agrobacterium larrymoorei)",
                    entry_route="targeted scan (miniprot) of 324 Fusarium-related genomes", locus_id=r["locus_id"], gca=r["gca"], organism=r["organism"], species=r["species"], order="Hypocreales",
                    contig=r["contig"], contig_len=r["contig_len"], start=r["start"], end=r["end"], strand=r["strand"], status=r["status"], n_exons=r["exons"],
                    protein_len=int(r["qto"]) - int(r["qfrom"]) + 1, note=r["note"], source_table="Table S14a"))
# ---- 4. nos-subfamily ODH of Boletales (blastx of plasmid proteins) ; gene-model state from nos_check (miniprot with KAG8213183.1)
state = {("GCA_038089435.1", "JBBPKX010000014.1"): "pseudogene", ("GCA_038089435.1", "JBBPKX010000046.1"): "pseudogene"}
for r in rd(A + "refset/screen/step2/boletales517_tdna_summary/opine_dehydrogenase_loci.tsv"):
    nosrel = "nos" in r["subject"]
    lin = "nos-subfamily opine dehydrogenase (Boletales)" if nosrel else "other opine dehydrogenase (Boletales)"
    st = state.get((r["gca"], r["contig"]), "intact" if nosrel and float(r["pident"]) > 36 else "not modelled")
    if not nosrel: st = "not modelled"
    out.append(dict(lineage=lin, family="opine dehydrogenase (nos subfamily)" if nosrel else "opine dehydrogenase", origin="non-T-DNA rhizobial locus (nos subfamily, outside T-DNA nos)" if nosrel else "bacterial homologs (not analysed)",
                    entry_route="blastx of 388 plasmid proteins on 517 Boletales assemblies", locus_id=f"ODH_{r['contig']}", gca=r["gca"], organism=r["organism"], species=" ".join(r["organism"].split()[:2]),
                    order="Boletales", contig=r["contig"], contig_len=r["contig_len"], start=r["start"], end=r["end"], strand=r["strand"], status=st, n_exons="1" if nosrel else "", protein_len="",
                    note=f"best plasmid protein {r['subject'].split('|')[0]} ({r['subject'].split('|')[4]}), {r['pident']}% identity", source_table="Table S11b"))
# manual overrides of nos-ODH loci (journal 2026-10-12 nosODH_count): ODH_JBBPWO010000281.1 (Boletaceae sp. TRTC157229) - coordinates are blastx HSP bounds, not a CDS;
# translation has 6 in-frame stops (3 broken ORFs in one frame), 35-37% identity to the other 6 nos-ODH loci of Boletales (76-100% among themselves), no miniprot model,
# not in the nos tree / AU tests -> pseudogene, not placed. The same override is applied to Table S11b in paper/supplementary/build_supplementary.py.
NOS_OVERRIDE = {"ODH_JBBPWO010000281.1": dict(lineage="nos-subfamily opine dehydrogenase (Boletales) (not placed)", status="pseudogene",
                                              note="6 in-frame stops in the blastx HSP frame; 35% identity to the Boletales nos-ODH genes; not in the nos tree")}
for r in out:
    r.update(NOS_OVERRIDE.get(r["locus_id"], {}))
# ---- 5. Rhypophila nos-subfamily genes (annotated proteomes)
RH = [("KAM7187585.1", "GCA_054392455.1", "Rhypophila sp. PSN 637", "JBTAFM010000166.1", 96792, 35649, 36989, "-", 2), ("KAK4208163.1", "GCA_033297075.1", "Rhypophila decipiens", "MU858255.1", 61824, 24006, "", "+", 2),
      ("KAM7213959.1", "GCA_054392445.1", "Rhypophila decipiens", "JBTAFL010000151.1", 93318, 50674, "", "+", 2)]
for p, g, o, c, cl, s, e, st, ex in RH:
    out.append(dict(lineage="nos-subfamily ODH homolog (Rhypophila, not assigned)", family="opine dehydrogenase (nos subfamily)", origin="not assigned to any lineage (not counted)",
                    entry_route="protein-level screen of annotated proteomes (reverse search)", locus_id=p, gca=g, organism=o, species=" ".join(o.split()[:2]), order="Sordariales", contig=c, contig_len=cl,
                    start=s, end=e, strand=st, status="intact", n_exons=ex, protein_len="", note="NCBI protein " + p + "; nos homolog at the rule-B threshold (33-36% to nos), base of the nos subfamily; not assigned to any lineage (2026-10-12)", source_table="Table S3b"))
# ---- fixed lineage order shared by all tables/figures/texts of the paper (v2): bacterial-origin lineages (T-DNA type, then other rhizobial loci), then lineages not analysed / native
ORDER = ["sus-like (T-DNA type)", "cus-like (T-DNA type)", "chs-like (T-DNA type)", "acc-type GDPD", "pKerr14a-type acs-like", "nos-subfamily opine dehydrogenase (Boletales)",
         "odh island: opine dehydrogenase", "odh island: Xaa-Pro aminopeptidase", "odh island: third gene"]
TAIL = ["nos-subfamily ODH homolog (Rhypophila, not assigned)", "OCD-like, outside T-DNA sus clade", "actinobacterial-type", "mis lineage", "native fungal GDPD (main clade)", "other fungal GDPD (position unresolved)",
        "other opine dehydrogenase (Boletales)", "not a homolog"]
LINEAGE_ORDER = ORDER + TAIL
# profile-only fragments without a lineage of their own are placed right after the closest lineage of the same family
ANCHOR = {"OCD-like bacterial": "OCD-like, outside T-DNA sus clade", "other bacterial-type CUS": "actinobacterial-type", "other bacterial-type CUS (beta-proteobacteria)": "actinobacterial-type",
          "cus-related small ORF (cusrel181)": "mis lineage", "rcorf12": "mis lineage", "acs-HMM group": "other fungal GDPD (position unresolved)",
          "no HMM subgroup": "other fungal GDPD (position unresolved)"}
def lineage_rank(r):
    lin = r["lineage"]; base = lin.replace(" (not placed)", ""); np_ = lin != base
    if base in LINEAGE_ORDER: return (LINEAGE_ORDER.index(base), int(np_), "")
    if base == "outgroup-HMM": base = "mis lineage" if r["family"].startswith("cucumopine") else "other fungal GDPD (position unresolved)"
    else: base = ANCHOR.get(base)
    if base is None: return (len(LINEAGE_ORDER), 2, lin)  # unknown label: at the very end
    return (LINEAGE_ORDER.index(base), 2, lin)
out = [r for _, r in sorted(enumerate(out), key=lambda ir: (lineage_rank(ir[1]), ir[0]))]  # stable: original order within a lineage
with open(T + "TableS1_all_lineages.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, COLS, delimiter="\t"); w.writeheader(); w.writerows(out)
# ---- summary per lineage (bacterial-origin lineages first)
rows = []
present = {r["lineage"] for r in out if not r["lineage"].endswith("(not placed)")}
for lin in [l for l in LINEAGE_ORDER if l in present] + sorted(present - set(LINEAGE_ORDER)):
    rs = [r for r in out if r["lineage"] == lin]
    if not rs: continue
    if lin == "nos-subfamily ODH homolog (Rhypophila, not assigned)": continue  # 2026-10-12v: not a lineage; loci stay in Table S1, not in Table 2
    np_ = [r for r in out if r["lineage"] == lin + " (not placed)"]  # profile-only fragments not placed on a tree: counted separately
    st = collections.Counter(r["status"] for r in rs)
    oc = collections.defaultdict(set)
    for r in rs: oc[r["order"] or "NA"].add(r["species"])
    og = collections.defaultdict(set)  # 2026-10-12v: carrier assemblies (gca) by order, placed loci only
    for r in rs: og[r["order"] or "NA"].add(r["gca"])
    other = "; ".join(f"{k} {len(v)}" for k, v in sorted(og.items(), key=lambda kv: (-len(kv[1]), kv[0])) if k not in ("Boletales", "Agaricales", "Hypocreales"))
    rows.append(dict(lineage=lin, family=rs[0]["family"], origin=rs[0]["origin"], entry_route=rs[0]["entry_route"], n_loci=len(rs), intact=st["intact"], pseudogene=st["pseudogene"],
                     truncated_fragment=st["truncated"] + st["fragment"] + st["remnant"], other=len(rs) - st["intact"] - st["pseudogene"] - st["truncated"] - st["fragment"] - st["remnant"],
                     n_genomes=len({r["gca"] for r in rs}), n_species=len({r["species"] for r in rs}), n_orders=len(oc),
                     Boletales=len(og["Boletales"]), Agaricales=len(og["Agaricales"]), Hypocreales=len(og["Hypocreales"]), other_orders=other or "0", profile_only_fragments=len(np_), genomes_incl_fragments=len({r["gca"] for r in rs + np_}), main_lineage=lin in ORDER))
# 2026-10-12x (review_terms): column names — "group" (tree-defined group; the first seven rows are the lineages of bacterial origin), "n_assemblies"
for r in rows:
    r["group"] = r.pop("lineage"); r["n_assemblies"] = r.pop("n_genomes"); r["assemblies_incl_fragments"] = r.pop("genomes_incl_fragments")
    r["truncated"] = sum(1 for x in out if x["lineage"] == r["group"] and x["status"] == "truncated")
    r["fragment"] = r.pop("truncated_fragment") - r["truncated"]  # fragments / remnants: too short to be placed in the tree (assigned by identity or profile)
rows = [{k: r[k] for k in ["group", "family", "origin", "entry_route", "n_loci", "intact", "pseudogene", "truncated", "fragment", "other", "n_assemblies", "n_species", "n_orders",
                            "Boletales", "Agaricales", "Hypocreales", "other_orders", "profile_only_fragments", "assemblies_incl_fragments", "main_lineage"]} for r in rows]
with open(T + "Table2a_lineages.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
# ---- Table 2 of the main text (2026-10-12x, review_critical 6.2): the seven lineages of bacterial origin, 8 columns; the full summary (Table2a) goes to the supplement (Table S18)
SRC = {"sus-like (T-DNA type)": "T-DNA *sus*", "cus-like (T-DNA type)": "T-DNA *cus*", "chs-like (T-DNA type)": "Ti-plasmid *chs*",
       "acc-type GDPD": "acc-like opine catabolic locus of rhizobia", "pKerr14a-type acs-like": "*acs* of the Ti plasmid pKerr14a; cT-DNA fragment of *Parasponia*",
       "nos-subfamily opine dehydrogenase (Boletales)": "rhizobial *nos* subfamily, outside T-DNA *nos*",
       "odh island: opine dehydrogenase": "opine-catabolic plasmid of *A. larrymoorei*", "odh island: Xaa-Pro aminopeptidase": "opine-catabolic plasmid of *A. larrymoorei*", "odh island: third gene": "opine-catabolic plasmid of *A. larrymoorei*"}
# 2026-10-13h (comments block 1): entry route in words
ROUTE = {"profile screen (HMM/blastx) of 517 Boletales + 1,386 other assemblies": "HMM profiles of the family vs. assemblies", "targeted scan (miniprot) of 324 Fusarium-related genomes": "miniprot scan of *Fusarium* assemblies",
         "blastx of 388 plasmid proteins on 517 Boletales assemblies": "blastx of plasmid proteins vs. Boletales assemblies"}
# identity of intact genes to the nearest reference, % (range; median) - values from Results 3.4 and Table S14c; pKerr14a-type not reported in the text (n.d.)
IDENT = {"sus-like (T-DNA type)": "38–68; median 46", "cus-like (T-DNA type)": "46–65; median 57", "chs-like (T-DNA type)": "59–60 (chrysopine synthase)",
         "acc-type GDPD": "60–66 (GDPD of *Mesorhizobium*, *Neorhizobium*)", "pKerr14a-type acs-like": "n.d.",
         "nos-subfamily opine dehydrogenase (Boletales)": "36–40 (nopaline synthase)", "odh island: opine dehydrogenase": "78–79", "odh island: Xaa-Pro aminopeptidase": "87", "odh island: third gene": "70"}
NAME = {"nos-subfamily opine dehydrogenase (Boletales)": "nos-subfamily opine dehydrogenase (Boletales)", "odh island: opine dehydrogenase": "odh island: *odh*",
        "odh island: Xaa-Pro aminopeptidase": "odh island: *pepP*", "odh island: third gene": "odh island: *orf3*"}
main = []
for r in rows:
    if not r["main_lineage"]: continue
    byo = f'{r["Boletales"]} / {r["Agaricales"]} / {r["Hypocreales"]} / {r["other_orders"]}'
    main.append({"lineage": NAME.get(r["group"], r["group"]), "family": r["family"], "source (closest bacterial relatives)": SRC[r["group"]], "entry route": ROUTE[r["entry_route"]], "identity to nearest reference, % (intact genes; range, median)": IDENT[r["group"]],
                 "loci": r["n_loci"], "intact / pseudogene / truncated / fragment": f'{r["intact"]} / {r["pseudogene"]} / {r["truncated"]} / {r["fragment"]}',
                 "assemblies: Boletales / Agaricales / Hypocreales / other orders": byo, "species": r["n_species"]})
with open(T + "Table2_main.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, list(main[0]), delimiter="\t"); w.writeheader(); w.writerows(main)
print(len(out), "loci in unified catalogue")
for r in rows: print(f"{'*' if r['main_lineage'] else ' '} {r['group'][:46]:46s} loci {r['n_loci']:4d} intact {r['intact']:4d} pseudo {r['pseudogene']:4d} trunc {r['truncated']:3d} frag {r['fragment']:3d} +frag(profile) {r['profile_only_fragments']:3d} assemblies {r['n_assemblies']:4d}(+{r['assemblies_incl_fragments']-r['n_assemblies']}) species {r['n_species']:4d} orders {r['n_orders']}")
