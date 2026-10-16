# Numerical thresholds of the classification

This table lists every numerical threshold used to classify hits and loci, the range of values over which the result of the classification does not change, and the section of the Methods where the threshold is applied.
It is a copy of Table S4 (`paper/supplementary/TableS4_Classification_thresholds_and_their_sensitivity.tsv`), which also gives the basis of each threshold and what changes outside the range.

| Step / quantity | Threshold | Range over which the result is unchanged | Where used (Methods section) |
|---|---|---|---|
| probable contamination: best bacterial hit identity | ≥90% | 79.1–90.4% | 2.1 |
| bacterial-type protein: non-fungal among 50 best hits | ≥35/50 | 0–37 | 2.1 |
| bacterial-type protein: bit-score margin bacterial/plant − fungal | ≥50 bits | 39–64 bits | 2.1 |
| T-DNA candidate: identity to T-DNA seed | ≥35% | 30.1–35.9% | 2.1 |
| T-DNA candidate: coverage of T-DNA seed | ≥0.4 | 0–0.49 | 2.1 |
| T-DNA candidate: seed identity ≥ ratio × best bacterial identity | ≥0.9 | 0.56–0.91 | 2.1 |
| non-T-DNA candidate: best Rhizobiales hit identity | ≥50% | 42.5–50.9% | 2.1 |
| vir/tra/rep genes on the same contig | ≥60% | contig flag of family loci: 0 at any 40–80%; genome flag (≥3 contigs): 21 loci at 60–75%; Boletales vir/tra/rep only in B. harrisonii: 55.2–89.6% | 2.3 |
| chs-like locus: identity to chrysopine synthase seed | ≥45% | loci: 0–58.5% (33 loci); proteins: 37.5–44.0% | 2.2 (ii) |
| chs-like locus: identity to F. vanettenii EEU41297 (alternative criterion) | ≥70% | 0–75.2% (33 loci) | 2.2 (ii) |
| pKerr14a-type fragment → intact member | ≥60% | 49.5–62.7% (top hit); 49.5–67.2% (max identity) | 2.4 |
