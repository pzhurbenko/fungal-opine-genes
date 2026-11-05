# Analyses added for the revision of 4 November 2026

The two M0 sensitivity calculations and four WGS mappings below were completed. The original phylogenetic and expression analyses were retained.

| Analysis | Main result | Files |
|---|---|---|
| SUS without unassessable loci | 15 sequences; strict / ≤20%-gap ω = 0.21949 / 0.22671 | `dnds/SUS_conf/` |
| Non-Boletales cus without unassessable loci | 19 sequences; ω = 0.09842 / 0.09965 | `dnds/cus_nonBol_conf/` |
| CLIB189 cus-like WGS | Full coverage; 128.6× locus / 127× genome; 72/73 edge-spanning reads at MAPQ ≥20 | `wgs/ERR12742700/` |
| Fol4287 chs-like WGS | Full coverage; 13.9× / 13×; 1/4 edge-spanning reads | `wgs/SRR305351/` |
| Fp133 WGS | All 7 loci covered, 0.77–1.05 relative Illumina depth; 65/47 Nanopore molecules span chr2/chr3 blocks and both fungal anchors | `wgs/SRR23927955/`, `wgs/SRR23927956/` |
| Plasmid provenance | 42 records; candidate intervals in 23, 302 proteins from 22; all sequences match source CDSs | `plasmids/` |

Both starting ω values (0.4, 1.5) converged for each new M0 fit. Only M0 was fitted for these two new subsets. WGS summaries separate primary MAPQ ≥0 and ≥20 alignments; gene-edge continuity and whole-block linkage are separate tests. Run provenance files give exact public read accessions, assembly accessions, FASTA SHA256 values, FASTQ URLs and MD5 values. Fp133 and CLIB189 share their assembly BioSamples; Fol4287 is a separate deposit of the same named strain. The status of CUS_OZ188520_A1 changes from supported to confirmed by the existing D2 rule. Counts: 95 confirmed, 87 supported, 72 unassessable, 254 total; 182 assessable loci include 72 intact copies.

## Rerunning

Set `GMO_CONDA_ENVS` to the directory containing the `annot`, `phylo_env`, `minimap2_env` and `PlantCad2` environments, or adjust the executable paths in the scripts. Python needs Biopython, numpy and pysam as applicable. The two dN/dS jobs can be run as `python dnds/run_filtered.py SUS_conf` and `python dnds/run_filtered.py cus_nonBol_conf`; original CDS/protein inputs are under the existing `analysis/dnds/` tree. Each job uses two threads. The recorded input audit is from the integration categories at the time of fitting; CLIB189's later change from supported to confirmed does not change subset membership.

For each run, obtain the reference assembly at the path in `wgs/map_reads.py`, verify its SHA256 against `provenance.json`, run `python wgs/download_ranges.py RUN`, then `python wgs/map_reads.py RUN` on a node with 16 CPUs and about 16 GB RAM. Mapping calls `summarize_reads.py`. The raw reads and BAMs are large and are excluded from the repository. The spanning-read CIGARs are in `wgs/SRR23927956/spanning_alignments.tsv`. `anchor_search.json` records the completed DIAMOND search and database snapshot; `anchors_vs_nr.tsv` retains all 800 reported hits. The anchor search establishes positive fungal homology and is not an absence test.

`python plasmids/audit_plasmids.py` verifies all 302 proteins against the compressed source GenBank records. Legacy sequence IDs remain available in the original catalogue and Table S19b; corrected provenance does not alter the amino-acid inputs. Do not substitute the new headers into old analysis files without applying the mapping.

## Supplementary files

`supplement_tables/` contains the audited exports for the revised/new sheets. `paper/supplementary/revision_exports.py` integrates these with the existing table builder and appends the completed M0 outputs to Data S4. The original supplementary rebuild command therefore retains S19–S20 and Data S5. `Data_S5_revision_evidence.tar.gz` holds the compact read/provenance evidence and a checksum manifest. Raw output files for the new codeml fits are deliberately included despite the repository-wide `*.out` exclusion.
