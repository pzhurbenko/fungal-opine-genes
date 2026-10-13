#!/usr/bin/env python3
"""thresholds_summary.tsv: one row per threshold of the paper (supplementary table). Numbers come from
funnel_stability.tsv / funnel_oat.tsv (sensitivity_funnel.py), chs/, pkerr/, backbone/ (see those folders' scripts)."""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import pandas as pd
T = f"{GMO_ROOT}/analysis/thresholds/"
st = pd.read_csv(T + "funnel_stability.tsv", sep="\t").set_index("param")
rng = lambda k, u="": f"{st.loc[k, 'stable_from']:g}–{st.loc[k, 'stable_to']:g}{u}"
F = "analysis/thresholds/sensitivity_funnel.py; funnel_stability.tsv; funnel_oat.tsv; funnel_grid.tsv; funnel_flips.tsv; margins_116.tsv; analysis/funnel/build_funnel.py"
rows = [
 ("probable contamination: best bacterial hit identity", "2.2 protein screen (build_funnel.py)", "≥90%", "gap in distribution: 79.1% (F. langsethiae ocs, rule C) → 90.5% (Mortierella mas2); no value in between among 116 proteins",
  rng("contam", "%"), "≤79%: F. langsethiae ocs (79.1) and nos (78.0) move C→contamination (+ more below 75%: 14 at 75%); >90.4%: Mortierella mas2 (90.5) contamination→single (10 at 95%)",
  F + "; contamination_gap.tsv"),
 ("rule A: non-fungal among 50 best hits", "2.2 protein screen, rule A", "≥35/50", "stability range (bimodal: rule-A proteins 37–50, marginal mostly <35)",
  rng("A_hits"), ">37: Melanogaster sus (37/50) B→marginal; at 45: 5 proteins change (Melanogaster sus, Pleurotus acs ×2→marginal, Butyriboletus nos, Kwoniella vis)", F),
 ("rule A: bit-score margin bacterial/plant − fungal", "2.2 protein screen, rule A", "≥50 bits", "gap in distribution: no protein with margin 39–63 bits (38 Ambispora → 64 Leucoagaricus); Beauveria nos/ocs/vis margins 5–17 bits, far below",
  rng("A_margin", " bits"), "<39: Ambispora chs (38) marginal→single; ≤32: + Mycena nos/vis (32); ≤28: Aspergillus chs RAL15815 →C; ≤23 Trichosporon mas1; ≤16 Beauveria nos/vis (single species); >64: Leucoagaricus acs (64) single→marginal; ≥70–80: Rhypophila nos, Butyriboletus nos B→marginal", F),
 ("rule B: identity to T-DNA seed", "2.2 protein screen, rule B", "≥35%", "stability range (narrow): lower edge Leucoagaricus acs 30.1%, upper edge Rhypophila nos 36.0%",
  rng("B_ident", "%"), "≤30%: Leucoagaricus acs single→B; >35.9%: Rhypophila nos (3 proteins, one at 36.0%, two propagated) B→single; 45%: also Dactylonectria chs (44.0), Butyriboletus nos (40.0), Melanogaster sus (41.5)", F),
 ("rule B: coverage of T-DNA seed", "2.2 protein screen, rule B", "≥0.4", "stability range; Fusarium chs proteins with cov 0.29 pass via gene×genus propagation",
  rng("B_cov"), ">0.49: Melanogaster sus (cov 0.497) B→single; nothing changes down to 0", F),
 ("rule B: seed identity ≥ ratio × best bacterial identity", "2.2 protein screen, rule B", "≥0.9", "stability range (gap among A proteins 0.56 → 0.91)",
  rng("B_ratio"), ">0.91: Rhypophila nos (ratio 0.911) B→single; at 1.0 also Dactylonectria chs, Amanita & Hebeloma cus (B→C), 6 changes; <0.56: Cutaneotrichosporon mas1 single→B", F),
 ("rule C: best Rhizobiales hit identity", "2.2 protein screen, rule C", "≥50%", "stability range: Kwoniella vis (Bosea, 42.5%) below, Pleurotus acs (Ensifer, 51.0%) above",
  rng("C_ident", "%"), "<42.5%: Kwoniella vis single→C; >51%: Pleurotus pulmonarius acs (2 proteins) C→single; ≥57.5: Butyriboletus acs; ≥64.8 Boletus reticuloceps acs", F),
 ("plasmid backbone (vir/tra/rep) on the same contig", "2.3 contig filter (aggregate_step2.py CONTAM_PID=60)", "≥60%",
  "stability range for the stated result; no gap in genome-wide distribution (427 contigs ≥60%, 1401 ≥50%, 156 ≥70%)",
  "contig flag of family loci: 0 at any 40–80%; genome flag (≥3 contigs): 21 loci at 60–75%; Boletales 'backbone only in B. harrisonii': 55.2–89.6%",
  "genome flag: 23 at 55%, 27 at 50%, 15 at 80%; Boletales: at ≤55% aff. Porphyrellus (55.1) + six assemblies 50–52.5%; short_contig flag strongly threshold-dependent (237 at 40, 77 at 60, 9 at 70). None of the 11 contaminant contigs carries backbone at any threshold (max 31.6%)",
  "analysis/thresholds/backbone/ (backbone_contigs.tsv, backbone_thresholds.tsv, screen_flag_by_threshold.tsv, screen_boletales_*.tsv, scripts/); analysis/refset/screen/step2/aggregate_step2.py"),
 ("chs-like locus: identity to chrysopine synthase seed", "2.x chs_fusarium (parse_loci.py:109)", "≥45%",
  "gap in distribution: native homologues ≤37.4% (≤33.3% at cov≥0.25) → chs lineage ≥44.0% (proteins); miniprot loci 53.9–65.4% (no native homologues reported by miniprot)",
  "loci: 0–58.5% (33 loci); proteins: 37.5–44.0%", "loci: 60% −D. estremocensis L1 (58.5); >65.4% −D. macrodidyma; proteins: 35% +Venturia (36.3), Funneliformis (37.4) (both cov<0.25); 45% vs 40%: −D. estremocensis KAH7120428.1 (44.0, mis-modelled, 58.5% by miniprot)",
  "analysis/thresholds/chs/ (chs_identity.tsv, chs_thresholds.tsv, chs_sensitivity.py); analysis/chs_fusarium/stage2/parse_loci.py"),
 ("chs-like locus: identity to F. vanettenii EEU41297 (OR-rule)", "2.x chs_fusarium (parse_loci.py:109)", "≥70%",
  "gap: no locus between 65.7% and 75.2%; 18 F. oxysporum loci (16 truncated 75.2–76.6%, 2 pseudogenes 87.0%) pass only via this rule",
  "0–75.2% (33 loci)", "80–85%: −16 truncated F. oxysporum (17 loci left); 90%: 15. Lowering below 70% adds nothing",
  "analysis/thresholds/chs/ (vanettenii_identity.tsv, vanettenii_thresholds.tsv)"),
 ("pKerr14a-type fragment → intact member", "2.x / Table 1 note (build_table1.py:145; journal 2026-10-06)", "≥60%",
  "gap in distribution: unassigned ≤49.4% → assigned ≥62.7% (top hit; ≥67.2% by max identity)", "49.5–62.7% (top hit); 49.5–67.2% (max identity)",
  "45%: +2 acs-HMM remnants (JBBPNQ010000397, JBBPWX010062892); 40%: +CEMI01314509; 65%: −JBBPNQ010000794 (51 aa); 70%: −B. eberwhitei JBBPNO, −B. harrisonii JBBPRV (10 left); nearest intact member not unique (top-2 differ by 0–5 points)",
  "analysis/thresholds/pkerr/ (fragment_identity.tsv, fragment_thresholds.tsv, pkerr_fragments.py); paper/tables/build_table1.py"),
]
cols = ["threshold", "used_in", "value", "basis", "range_no_change", "what_changes_outside", "source_files"]
pd.DataFrame(rows, columns=cols).to_csv(T + "thresholds_summary.tsv", sep="\t", index=False)
print(pd.DataFrame(rows, columns=cols)[["threshold", "value", "range_no_change"]].to_string(index=False))
