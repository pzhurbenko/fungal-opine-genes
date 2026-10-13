# AU-тест: грибные sus-like + нативные грибные OCD-гомологи монофилетичны (2026-10-12)

Выравнивание `../SUS.mafft.trim.aln` (125 листьев), модель Q.PFAM+R5 (BIC, как `../SUS.AU`), поиск с ограничением `-g` (-T 16, -seed 2),
затем `iqtree -z au_trees.nwk -n 0 -zb 10000 -au` (ML + прежние H_tdna/H_alt + два новых ограничения). LSF au_sus_native (cn349),
скрипт `run_au_sus_native.sh`, порядок `hyp_order.txt`, состав групп по меткам листьев — `groups.txt`. Время: поиски 162 с + 156 с, AU 4 с.

Группы: SUS = 52 FUNGI_sus внутри T-DNA sus-клады (36 intact + 16 pseudogene-like; клада = наименьшая, содержащая 24 REF_*_sus);
NATIVE = 2 NATIVE_Sordariomycetes (Pseudopestalotiopsis chinensis, P. theae); fungOCD = 4 FUNGI_sus вне sus-клады
(Echinosphaeria canescens JAOXNW010100945, Austroboletus olivaceoglutinosus JBBPSN010000183, Boletus speciosus QFEG01000233,
Retiboletus ornatipes QFEH01000351); ogF = 5 грибных OCD аутгруппы (Imleria badia, Aureobasidium melanogenum, Podila verticillata,
Pestalotiopsis fici, Schizosaccharomyces osmophilus).

| # | гипотеза | ограничение | logL | ΔL | p-AU |
|---|---|---|---|---|---|
| 1 | ML | — | −42857.01 | 0 | 1.00 |
| 2 | H_tdna | SUS + 24 T-DNA/cT-DNA sus референсов монофилетичны | −42871.63 | 14.6 | 0.26 (не отвергнуто) |
| 3 | H_alt | SUS + все альтернативы (OG, NATIVE, REF chromosomal/bacteria) | −43033.26 | 176.3 | 5.9e-06 |
| 4 | **sus_natv_core** | **SUS + NATIVE + fungOCD (6 нативных грибных OCD) монофилетичны** | −43609.88 | 752.9 | **9.5e-08 (отвергнуто)** |
| 5 | **sus_natv_all** | **SUS + NATIVE + fungOCD + ogF (11 грибных OCD) монофилетичны** | −43566.86 | 709.8 | **9.9e-67 (отвергнуто)** |

Вывод: нативное (грибное) происхождение sus-like отвергнуто для обоих составов «нативных» (p ≤ 9.5e-08); положение внутри T-DNA sus-клады
не отвергнуто (0.26), как прежде. Симметрично тестам GDPD (acc, pKerr14a + native) и nos (Boletales + native).
