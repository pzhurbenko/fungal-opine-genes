# AU test: fungal sus-like + native fungal OCD homologs are monophyletic (2026-10-12)

Alignment `../SUS.mafft.trim.aln` (125 leaves), model Q.PFAM+R5 (BIC, as in `../SUS.AU`), constrained search with `-g` (-T 16, -seed 2),
then `iqtree -z au_trees.nwk -n 0 -zb 10000 -au` (ML + earlier H_tdna/H_alt + two new constraints). LSF au_sus_native (cn349),
script `run_au_sus_native.sh`, order in `hyp_order.txt`, group composition by leaf labels in `groups.txt`. Time: searches 162 s + 156 s, AU 4 s.

Groups: SUS = 52 FUNGI_sus inside the T-DNA sus clade (36 intact + 16 pseudogene-like; clade = smallest one containing the 24 REF_*_sus);
NATIVE = 2 NATIVE_Sordariomycetes (Pseudopestalotiopsis chinensis, P. theae); fungOCD = 4 FUNGI_sus outside the sus clade
(Echinosphaeria canescens JAOXNW010100945, Austroboletus olivaceoglutinosus JBBPSN010000183, Boletus speciosus QFEG01000233,
Retiboletus ornatipes QFEH01000351); ogF = 5 fungal OCD outgroup (Imleria badia, Aureobasidium melanogenum, Podila verticillata,
Pestalotiopsis fici, Schizosaccharomyces osmophilus).

| # | hypothesis | constraint | logL | ΔL | p-AU |
|---|---|---|---|---|---|
| 1 | ML | — | −42857.01 | 0 | 1.00 |
| 2 | H_tdna | SUS + 24 T-DNA/cT-DNA sus references monophyletic | −42871.63 | 14.6 | 0.26 (not rejected) |
| 3 | H_alt | SUS + all alternatives (OG, NATIVE, REF chromosomal/bacteria) | −43033.26 | 176.3 | 5.9e-06 |
| 4 | **sus_natv_core** | **SUS + NATIVE + fungOCD (6 native fungal OCD) monophyletic** | −43609.88 | 752.9 | **9.5e-08 (rejected)** |
| 5 | **sus_natv_all** | **SUS + NATIVE + fungOCD + ogF (11 fungal OCD) monophyletic** | −43566.86 | 709.8 | **9.9e-67 (rejected)** |

Conclusion: a native (fungal) origin of the sus-like genes is rejected for both "native" compositions (p ≤ 9.5e-08); placement inside the T-DNA sus clade
is not rejected (0.26), as before. Symmetric to the GDPD (acc, pKerr14a + native) and nos (Boletales + native) tests.
