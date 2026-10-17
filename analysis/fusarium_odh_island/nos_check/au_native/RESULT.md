# AU test: nos subfamily — Boletales nos-ODH + native fungal ODH (2026-10-12)

Alignment `../nos_tree.trim.aln` (64 leaves), model Q.PFAM+R4 (as in `../run_au.sh`), constrained search with `-g`, then
`iqtree -z au_trees.nwk -n 0 -zb 10000 -au` (ML + c1–c3 from `../` + new c4, c5). LSF 711483, `run_au_native.sh`.

Groups (prefixes in `nos_tree_set.faa`): BOL = 6 Boletales nos-ODH (Boletus ×4, Butyriboletus ×2); FUNG = 3 Rhypophila;
NAT = 19 native fungal ODH (NAT_*); T-DNA nos = NOS_TDNA_* (8); rhizobial nos-like = NOSLIKE_* (10).

| # | hypothesis | constraint | logL | ΔL | p-AU |
|---|---|---|---|---|---|
| 1 | ML | — | −21721.11 | 0 | 1.00 |
| 2 | c1_bol_in_tdna | BOL + NOS_TDNA + NOSLIKE_A.tumefaciens_WP_174048180 monophyletic | −21756.74 | 35.6 | 0.0034 |
| 3 | c2_rh_in_tdna | FUNG (Rhypophila) + T-DNA nos monophyletic | −21762.95 | 41.8 | 0.0091 |
| 4 | c3_fungal_mono | BOL + FUNG monophyletic | −21726.11 | 5.0 | 0.32 (not rejected) |
| 5 | **c4_bol_nat** | **BOL + NAT (19 native fungal ODH) monophyletic** | −21911.33 | 190.2 | **3.4e-05 (rejected)** |
| 6 | c5_bol_rh_nat | BOL + FUNG + NAT monophyletic ("all fungal ODH share one native source") | −21939.24 | 218.1 | 1.0e-26 (rejected) |

Note: p-AU for c1–c3 differ slightly from `../au_test.iqtree` (0.0030/0.0089/0.345) because the test uses a different tree set (6 instead of 4);
conclusions are the same. Conclusion: the Boletales nos-ODH cannot be native fungal ODH — vertical origin from fungal ODH is rejected;
origin from the rhizobial/Agrobacterium nos subfamily remains the only non-rejected alternative
(apart from BOL+Rhypophila monophyly, which does not contradict it).
