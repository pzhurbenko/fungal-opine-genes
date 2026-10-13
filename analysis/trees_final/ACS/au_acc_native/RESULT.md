# AU-тест: acc-type GDPD + нативные грибные GDPD (2026-10-12)

Выравнивание `../ACS.mafft.trim.aln` (259 листьев), модель JTT+R5 (BIC, как `../ACS.lineage.AU`), поиск с ограничением `-g` (2 574 с, 32 потока),
затем `iqtree -z au_trees.nwk -n 0 -zb 10000 -au` (ML + 6 прежних ограничений GDPD + новое). LSF 712220, `run_au_acc_native.sh`, порядок — `hyp_order.txt`, состав групп — `groups.txt`.

Группы (по меткам листьев, как в `../../au_lineage_acs.py`): MESO = 60 FUNGI_gdpd_meso (acc-type GDPD); NATV = 138 (NATIVE_* + FUNGI_acs_* кроме шести pKerr14a-type).

| # | гипотеза | logL | ΔL | p-AU |
|---|---|---|---|---|
| 1 | ML | −77045.09 | 0 | 1.00 |
| 2 | six_tdna (pKerr14a-type + T-DNA acs) | −78767.86 | 1722.8 | 1.4e-31 |
| 3 | six_kerr (pKerr14a-type + pKerr14a acs + Parasponia) | −77060.36 | 15.3 | 0.44 |
| 4 | six_acc (pKerr14a-type + acc-type + acc-locus refs) | −77100.34 | 55.2 | 0.16 |
| 5 | six_natv (pKerr14a-type + native GDPD) | −80303.87 | 3258.8 | 4.7e-61 |
| 6 | meso_acc (acc-type + acc-locus refs) | −77045.91 | 0.8 | 0.56 |
| 7 | meso_tdna (acc-type + T-DNA acs) | −79161.32 | 2116.2 | 7.4e-11 |
| 8 | **meso_natv (acc-type + native fungal GDPD, 138)** | −80542.68 | **3497.6** | **1.1e-95 (отвергнуто)** |

Вывод: нативное происхождение acc-type GDPD отвергнуто; прежние вердикты (п. 2–7) воспроизведены (p в прежнем прогоне: 1.4e-60, 0.42, 0.13, 2.4e-06, 1.0, 5.9e-47).
