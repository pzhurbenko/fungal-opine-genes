# AU test: Rhypophila nos-subfamily genes + native fungal ODH (2026-10-12)

Alignment `../nos_tree.trim.aln` (64 leaves), model Q.PFAM+R4 (as in `../run_au.sh`, `../au_native/run_au_native.sh`), constrained search with `-g`,
then `iqtree -z au_trees.nwk -n 0 -zb 10000 -au` (ML + earlier c1–c5 + new c6–c8). LSF 712191, `run_au_rhypo.sh`, order in `hyp_order.txt`, composition in `groups.txt`.

Groups: RH = 3 Rhypophila genes (FUNG_*); 2SORD = NAT_Cryphonectria_parasitica_XP_040775628 + NAT_Gnomoniopsis_smithogilvyi_KAJ4397655 (closest native ODH in the tree);
SORD = 12 native Sordariomycetes ODH (Cryphonectria, Hirsutella, Purpureocillium ×2, Claviceps ×3, Diaporthaceae sp., Gnomoniopsis, Thozetella, Apiospora, Astrocystis); NAT = all 19 native fungal ODH.

Position in the ML tree (root ocs/vis): (RH, (BOL, (T-DNA nos + rhizobial nos-like)))) 100/100 — Rhypophila at the base of the subfamily; in the UFBoot consensus 2SORD is the next branch outward (100), in the ML tree their position is unresolved (36.2/46).

| # | hypothesis | constraint | logL | ΔL | p-AU |
|---|---|---|---|---|---|
| 1 | ML | — | −21721.11 | 0 | 1.00 |
| 2 | c1_bol_in_tdna | BOL + T-DNA nos monophyletic | −21756.74 | 35.6 | 0.0030 |
| 3 | c2_rh_in_tdna | RH + T-DNA nos monophyletic | −21762.95 | 41.8 | 0.0123 |
| 4 | c3_fungal_mono | BOL + RH monophyletic | −21726.11 | 5.0 | 0.34 (not rejected) |
| 5 | c4_bol_nat | BOL + NAT | −21911.33 | 190.2 | 3.8e-44 |
| 6 | c5_bol_rh_nat | BOL + RH + NAT | −21939.24 | 218.1 | 1.2e-07 |
| 7 | **c6_rh_2sord** | **RH + 2SORD monophyletic** | −21772.44 | 51.3 | **0.0020 (rejected)** |
| 8 | **c7_rh_sord** | **RH + 12 native Sordariomycetes ODH** | −22102.27 | 381.2 | **1.1e-37 (rejected)** |
| 9 | **c8_rh_nat** | **RH + all 19 native fungal ODH** | −21875.36 | 154.3 | **2.1e-46 (rejected)** |

Note: p-AU for c1–c5 in this set of 9 trees differ slightly from the earlier runs (c1 0.003/0.0034, c2 0.009/0.0091, c3 0.35/0.32, c4 3.4e-05, c5 1.0e-26) — conclusions are the same; Table S5 keeps the earlier values for c1–c5.

Conclusion (fact): the Rhypophila genes group neither with the closest native Sordariomycete ODH nor with native ODH in general; they sit at the base of the nos subfamily, outside T-DNA nos (0.009) and outside the Boletales genes. By the user's decision (2026-10-12) they are NOT counted as a lineage: three single genes across all assemblies of the genus, 33–36 % identity at the rule-B threshold, position at the subfamily boundary; described by facts as a borderline rule-B hit, "not assigned to any lineage (not counted)".
