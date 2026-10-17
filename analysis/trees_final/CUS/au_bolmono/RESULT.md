# AU test: all cus-like Boletales are monophyletic (2026-10-12)

Alignment `../CUS.mafft.trim.aln` (436 leaves, 810 columns), model JTT+R7 (as in `../../au_lineage_cus.sh`), constrained search with
`-g`, then `iqtree -z ... -n 0 -zb 10000 -au`. Comparison set: ML + four earlier Boletales constraints (`../CUS.bol_*_x.treefile`)
+ the new bol_mono_x. Constraint: 26 cus-like Boletales (`bol26.txt`; group BOL from au_lineage_cus.py, by `order == Boletales` in loci_table,
without JBBPKX010002805) form a single clade; the remaining 410 leaves are free.

Runs: run1 (LSF 711484, -T 32, seed 1) hung after iteration 10 for >70 min with no progress in the log at full CPU load — killed at 04:28;
result is **run2/** (LSF 712190, -T 16, seed 2, -ninit 20): search 3 320 s, AU `run2/CUS.lineage_bm_s2.AU.iqtree`.

| # | hypothesis | logL | ΔL | p-AU (this set) | p-AU in lineage_x (earlier) |
|---|---|---|---|---|---|
| 1 | ML | −79268.36 | 0 | 1.00 | 0.53 |
| 2 | bol_tdna_x (BOL + 21 REF T-DNA cus) | −79784.94 | 516.6 | 5.1e-07 | 1.3e-03 |
| 3 | bol_act_x (BOL + actinobacterial cus) | −79909.52 | 641.2 | 1.3e-56 | 5e-59 |
| 4 | bol_mis_x (BOL + mis) | −80428.11 | 1159.7 | 1.3e-46 | 3.3e-77 |
| 5 | bol_othercus_x (BOL + 16 cus-like outside Boletales) | −79400.45 | 132.1 | 3.2e-05 | 1.6e-05 |
| 6 | **bol_mono_x (26 cus-like Boletales monophyletic)** | −79609.51 | 341.2 | **4.9e-06 (rejected)** | — |

Conclusion: cus-like Boletales are not monophyletic — in the ML tree they are scattered across the T-DNA cus clade, interleaved with 21 REF, 16 fungal genes outside Boletales and
plant cT-DNA (journal AU_cus_nb_interpretation); the hypothesis of a single cus acquisition by the Boletales ancestor (with subsequent losses) is rejected,
which is consistent with ≥2 independent entries into Boletales by ML/UFBoot.
