# AU test: acc-type GDPD + native fungal GDPD (2026-10-12)

Alignment `../ACS.mafft.trim.aln` (259 leaves), model JTT+R5 (BIC, as in `../ACS.lineage.AU`), constrained search with `-g` (2 574 s, 32 threads),
then `iqtree -z au_trees.nwk -n 0 -zb 10000 -au` (ML + 6 earlier GDPD constraints + the new one). LSF 712220, `run_au_acc_native.sh`, order in `hyp_order.txt`, group composition in `groups.txt`.

Groups (by leaf labels, as in `../../au_lineage_acs.py`): MESO = 60 FUNGI_gdpd_meso (acc-type GDPD); NATV = 138 (NATIVE_* + FUNGI_acs_* except the six pKerr14a-type).

| # | hypothesis | logL | ΔL | p-AU |
|---|---|---|---|---|
| 1 | ML | −77045.09 | 0 | 1.00 |
| 2 | six_tdna (pKerr14a-type + T-DNA acs) | −78767.86 | 1722.8 | 1.4e-31 |
| 3 | six_kerr (pKerr14a-type + pKerr14a acs + Parasponia) | −77060.36 | 15.3 | 0.44 |
| 4 | six_acc (pKerr14a-type + acc-type + acc-locus refs) | −77100.34 | 55.2 | 0.16 |
| 5 | six_natv (pKerr14a-type + native GDPD) | −80303.87 | 3258.8 | 4.7e-61 |
| 6 | meso_acc (acc-type + acc-locus refs) | −77045.91 | 0.8 | 0.56 |
| 7 | meso_tdna (acc-type + T-DNA acs) | −79161.32 | 2116.2 | 7.4e-11 |
| 8 | **meso_natv (acc-type + native fungal GDPD, 138)** | −80542.68 | **3497.6** | **1.1e-95 (rejected)** |

Conclusion: a native origin of the acc-type GDPD is rejected; the earlier verdicts (items 2–7) are reproduced (p in the earlier run: 1.4e-60, 0.42, 0.13, 2.4e-06, 1.0, 5.9e-47).
