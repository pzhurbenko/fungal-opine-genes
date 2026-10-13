# AU-тест: вложены ли грибные chs-like в кладу Ti-плазмидных chs? (2026-10-12)

Выравнивание `../chs_tree.trim.aln` (58 листьев), модель Q.PFAM+I+G4 (BIC, как `../au/run_au.sh`), поиск с ограничением `-g`,
затем `iqtree -z au_trees.nwk -n 0 -zb 10000 -au` (ML + 10 ограничений). LSF 711482, `run_au_nested.sh`, порядок — `hyp_order.txt`.

ML-топология (корень SEED_ags_P27875): (Fusarium ×10 [99.9/100], (Dactylonectria, (B084, (B272951, (B272871, (SEED_chs, B139)))))) —
Dactylonectria сестринская всем пяти Ti-chs, Fusarium-клада сестринская узлу (Dactylonectria + Ti-chs). Ни один грибной ген внутри клады Ti-chs.

Ti-chs (5): SEED_chs_A0A2Z2PGA9 (seed, chrysopine synthase pTi), BACT_WP_139103025.1 (Rhizobium/Agrobacterium group), BACT_WP_272871282.1
(A. tumefaciens), BACT_WP_272951408.1 (A. vitis), BACT_WP_084631803.1 (A. larrymoorei). Грибные chs-like: 11 CHSlike_* (10 Fusarium + Dactylonectria estremocensis).

Ограничения: nest_x = все 11 грибных + Ti-chs x монофилетичны (грибы вложены в Ti-chs рядом с x); fnest_x = 10 Fusarium + x (Dactylonectria свободна).

| # | гипотеза | logL | ΔL | p-AU |
|---|---|---|---|---|
| 1 | ML (грибы сестринские кладе Ti-chs) | −24638.05 | 0 | 1.00 |
| 2 | nest_seed (11 грибных + SEED_chs) | −24771.31 | 133.3 | 1.9e-04 |
| 3 | fnest_seed (10 Fusarium + SEED_chs) | −24785.33 | 147.3 | 8.4e-62 |
| 4 | nest_b139 (+ WP_139103025.1) | −24771.31 | 133.3 | 1.9e-04 |
| 5 | fnest_b139 | −24785.33 | 147.3 | 4.3e-60 |
| 6 | nest_b272871 (+ A. tumefaciens WP_272871282.1) | −24742.25 | 104.2 | 1.2e-04 |
| 7 | fnest_b272871 | −24758.09 | 120.0 | 5.2e-58 |
| 8 | nest_b272951 (+ A. vitis WP_272951408.1) | −24753.76 | 115.7 | 2.3e-05 |
| 9 | fnest_b272951 | −24771.84 | 133.8 | 2.6e-40 |
| 10 | nest_b084 (+ A. larrymoorei WP_084631803.1) | −24713.84 | 75.8 | 3.8e-06 |
| 11 | fnest_b084 | −24732.97 | 94.9 | 3.7e-05 |

Все 10 «вложенных» положений отвергнуты (p-AU ≤ 1.9e-04; Holm-поправка на 10 тестов не меняет вывода, max p·10 = 1.9e-03).
Наименьшая ΔL — у положения рядом с A. larrymoorei (75.8). nest_seed совпадает по logL с прежним chs_tdna (`../au/au_test.iqtree`, ΔL 133.26) — то же ограничение.

**Интерпретация: «sister, not nested»** — грибные chs-like являются ближайшими родственниками Ti-плазмидных chs, но не вложены в их кладу;
донор (T-DNA chs либо родственный плазмидный ген вне выборки) не определяется по топологии. Сочетание с прежними тестами
(`../au/RESULT` в au_test.iqtree): грибы + seed chs без Ti-chs 2.0e-48; грибы + 4 Ti-chs без seed 1.5e-72; с нативными AGS 1.7e-03;
с seed ags 8.8e-04; с прочими бактериальными M24 4.1e-15; грибы монофилетичны 0.167 (не отвергнуто).
