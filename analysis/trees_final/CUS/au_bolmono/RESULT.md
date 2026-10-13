# AU-тест: все cus-like Boletales монофилетичны (2026-10-12)

Выравнивание `../CUS.mafft.trim.aln` (436 листьев, 810 колонок), модель JTT+R7 (как `../../au_lineage_cus.sh`), поиск с ограничением
`-g`, затем `iqtree -z ... -n 0 -zb 10000 -au`. Набор сравнения: ML + четыре прежних Boletales-ограничения (`../CUS.bol_*_x.treefile`)
+ новое bol_mono_x. Ограничение: 26 cus-like Boletales (`bol26.txt`; группа BOL из au_lineage_cus.py, по `order == Boletales` в loci_table,
без JBBPKX010002805) образуют одну кладу; остальные 410 листьев свободны.

Запуски: run1 (LSF 711484, -T 32, seed 1) завис после итерации 10 на >70 мин без прогресса в логе при полной загрузке CPU — убит 04:28;
результат — **run2/** (LSF 712190, -T 16, seed 2, -ninit 20): поиск 3 320 с, AU `run2/CUS.lineage_bm_s2.AU.iqtree`.

| # | гипотеза | logL | ΔL | p-AU (этот набор) | p-AU в lineage_x (ранее) |
|---|---|---|---|---|---|
| 1 | ML | −79268.36 | 0 | 1.00 | 0.53 |
| 2 | bol_tdna_x (BOL + 21 REF T-DNA cus) | −79784.94 | 516.6 | 5.1e-07 | 1.3e-03 |
| 3 | bol_act_x (BOL + актинобактериальные cus) | −79909.52 | 641.2 | 1.3e-56 | 5e-59 |
| 4 | bol_mis_x (BOL + mis) | −80428.11 | 1159.7 | 1.3e-46 | 3.3e-77 |
| 5 | bol_othercus_x (BOL + 16 cus-like вне Boletales) | −79400.45 | 132.1 | 3.2e-05 | 1.6e-05 |
| 6 | **bol_mono_x (26 cus-like Boletales монофилетичны)** | −79609.51 | 341.2 | **4.9e-06 (отвергнуто)** | — |

Вывод: cus-like Boletales не монофилетичны — в ML они разбросаны по T-DNA-кладе cus вперемешку с 21 REF, 16 грибными вне Boletales и
растительными cT-DNA (журнал AU_cus_nb_interpretation); гипотеза одного приобретения cus предком Boletales (с последующими потерями) отвергнута,
что согласуется с ≥2 независимыми входами в Boletales по ML/UFBoot.
