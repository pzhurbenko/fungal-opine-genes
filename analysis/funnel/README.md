# carriers_by_order.tsv (2026-10-13, analysis "a", review point M11)
Inputs: paper/supplementary/TableS3_Proteins_with_non_fungal_closest_relatives_and_the.tsv (116 proteins);
  analysis/refset/screen/step3/genome_meta.tsv (6,847 assemblies, gca + taxid); NCBI taxdump in analysis/refset/screen/seeds_ext/reciprocal/taxdump/.
Universe: GenBank/RefSeq pairs merged by numeric accession (GCA preferred, then highest version) = 6,123 proteomes, as in aggregate_seeds_ext.py.
Order and phylum: taxid of each assembly walked up the taxdump lineage; no order rank = "unassigned (<class/phylum>)".
Candidate proteome: >=1 protein with fate "T-DNA candidate" or "non-T-DNA (rhizobial) candidate" and duplicate_assembly != yes (43 proteins, 42 proteomes).
share_pct = 100 * n_with_candidates / n_proteomes of that order; orders with zero candidates are not listed.
n_TDNA / n_nonTDNA count proteomes by fate (a proteome with both would count in both columns).
genera_with_candidates = genus of the organism name, semicolon-separated. Script: analysis/funnel/carriers_by_order.py (run with python3 -I).
