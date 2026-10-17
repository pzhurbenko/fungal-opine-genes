# AU test: are the fungal chs-like genes nested within the Ti-plasmid chs clade? (2026-10-12)

Alignment `../chs_tree.trim.aln` (58 leaves), model Q.PFAM+I+G4 (BIC, as in `../au/run_au.sh`), constrained search with `-g`,
then `iqtree -z au_trees.nwk -n 0 -zb 10000 -au` (ML + 10 constraints). LSF 711482, `run_au_nested.sh`, order in `hyp_order.txt`.

ML topology (root SEED_ags_P27875): (Fusarium ×10 [99.9/100], (Dactylonectria, (B084, (B272951, (B272871, (SEED_chs, B139)))))) —
Dactylonectria is sister to all five Ti-chs, the Fusarium clade is sister to the node (Dactylonectria + Ti-chs). No fungal gene lies inside the Ti-chs clade.

Ti-chs (5): SEED_chs_A0A2Z2PGA9 (seed, chrysopine synthase pTi), BACT_WP_139103025.1 (Rhizobium/Agrobacterium group), BACT_WP_272871282.1
(A. tumefaciens), BACT_WP_272951408.1 (A. vitis), BACT_WP_084631803.1 (A. larrymoorei). Fungal chs-like: 11 CHSlike_* (10 Fusarium + Dactylonectria estremocensis).

Constraints: nest_x = all 11 fungal + Ti-chs x monophyletic (fungi nested in Ti-chs next to x); fnest_x = 10 Fusarium + x (Dactylonectria free).

| # | hypothesis | logL | ΔL | p-AU |
|---|---|---|---|---|
| 1 | ML (fungi sister to the Ti-chs clade) | −24638.05 | 0 | 1.00 |
| 2 | nest_seed (11 fungal + SEED_chs) | −24771.31 | 133.3 | 1.9e-04 |
| 3 | fnest_seed (10 Fusarium + SEED_chs) | −24785.33 | 147.3 | 8.4e-62 |
| 4 | nest_b139 (+ WP_139103025.1) | −24771.31 | 133.3 | 1.9e-04 |
| 5 | fnest_b139 | −24785.33 | 147.3 | 4.3e-60 |
| 6 | nest_b272871 (+ A. tumefaciens WP_272871282.1) | −24742.25 | 104.2 | 1.2e-04 |
| 7 | fnest_b272871 | −24758.09 | 120.0 | 5.2e-58 |
| 8 | nest_b272951 (+ A. vitis WP_272951408.1) | −24753.76 | 115.7 | 2.3e-05 |
| 9 | fnest_b272951 | −24771.84 | 133.8 | 2.6e-40 |
| 10 | nest_b084 (+ A. larrymoorei WP_084631803.1) | −24713.84 | 75.8 | 3.8e-06 |
| 11 | fnest_b084 | −24732.97 | 94.9 | 3.7e-05 |

All 10 "nested" placements are rejected (p-AU ≤ 1.9e-04; Holm correction for 10 tests does not change the conclusion, max p·10 = 1.9e-03).
The smallest ΔL is for the placement next to A. larrymoorei (75.8). nest_seed has the same logL as the earlier chs_tdna (`../au/au_test.iqtree`, ΔL 133.26) — the same constraint.

**Interpretation: "sister, not nested"** — the fungal chs-like genes are the closest relatives of the Ti-plasmid chs but are not nested in their clade;
the donor (T-DNA chs or a related plasmid gene outside the sample) cannot be determined from the topology. Combined with the earlier tests
(`../au/RESULT` in au_test.iqtree): fungi + seed chs without Ti-chs 2.0e-48; fungi + 4 Ti-chs without seed 1.5e-72; with native AGS 1.7e-03;
with seed ags 8.8e-04; with other bacterial M24 4.1e-15; fungi monophyletic 0.167 (not rejected).
