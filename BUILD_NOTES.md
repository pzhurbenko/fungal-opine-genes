# BUILD_NOTES

Assembled 2026-10-13 from the working project (`GMO/`) following `paper/submission/repo_inventory.md`. Originals were
copied, never moved or edited; this directory (`GMO/repo/`) holds copies only. No `git init`, nothing pushed.

## What was copied from where (relative to the working project root)

| In repo | Source | Rule |
|---|---|---|
| `pipeline/` | `pipeline/` | all files except `__pycache__`; `scripts/common.py` and `run_pipeline.sh` got a small edit (below); `README.md` translated to English |
| `analysis/<dir>` | `analysis/<dir>` | in-paper folders of inventory section 2: refset/{plasmids,hmm_v0,neighbourhood,screen/{seeds_ext,step2,step3,tdna_hmm}} + refset/step2/queries*.faa + outgroup json, loci_pipeline/{addenda,guides,registry,round2 (tables only, plus `loci/ACS_JBBPRY010001778_L1/augustus.gff3`)}, trees_final, chs_fusarium, fusarium_odh_island, pleurotus_gdpd, fold, dnds, composition, conservation, islands, island_order, synteny/islands_hq, indep_assemblies, read_coverage, expression, boletales_mash_517, thresholds, funnel, ctdna_plants_identity, curated/*.fasta (report_A, report_B, clustering and conservation were removed on 2026-10-16, see the last section) |
| file filter | | code (`.py .sh .R .yaml .md`) and small data (`.tsv .txt .faa .fasta .aln .json .jsonl .nwk .treefile .iqtree .fna .gff .gff3 .csv .ctl .phy .nex`) up to 400 KB; explicit larger inclusions: `all_screen.hmm` (8.9 MB), `funnel_grid.tsv`, `all_blast_stats.tsv`, `seeds_ext_best.tsv`, `reciprocal.tsv`, `round2/loci_table.tsv`, `step2_family_calls.tsv`, `registry/loci.tsv`, `step3/proteomes.txt`, `step3/genome_meta.tsv`, `plasmids/cds_table.tsv`, `islands_genes.tsv`, `donor/NZ_JAVIZC010000002.1.gb`, `boletales517_ncbi_summary.jsonl`, `TableS1b`, `TableS1_*` catalogue inputs, `Data_S1` |
| `data/` | `data/tsv/merged/all_blast_stats.tsv`, `data/tsv/rep_seqs/*.tsv`, `data/fungal_genomes/{gca_list.txt,fusarium/gca_list.txt,boletales_all/accessions_20261006.tsv}`, `data/fungal_proteomes/download.sh` | accession lists and pipeline inputs only |
| `paper/figures`, `paper/tables`, `paper/supplementary` | `paper/...` | build scripts, `common.py`, `export_tiff.py`, small inputs, supplementary tables S1-S17, Data S1/S3, `supplementary/data/trees`; no PNG/PDF/SVG/TIFF, no docx/xlsx |
| `paper/manuscript_v2/` | `paper/manuscript_v2/` | (removed 2026-10-16: manuscript-text tooling is not part of the repository) |
| `paper/rebuild_all.sh` | `paper/rebuild_all.sh` | figures and supplementary only (manuscript steps removed 2026-10-16) |
| `docs/methods_extended.md` | `paper/manuscript_v2/methods_extended.md` (working project) | |
| `env/` | `conda env export` of annot, phylo_env, mashtree_env, minimap2_env, colabfold, paper_build, orthofinder, hmmer_pfam | `<env>.yml` = `--from-history`, `<env>.pinned.yml` = `--no-builds` |
| `LICENSE`, `LICENSE-DATA`, `.gitignore`, `README.md`, `gmo_paths.py` | written new | |

Excluded: inventory section 5, donor_meta, donor_search, generax, alien_scan, distribution, orf358_origin, accF, acc_test,
boletales_mash (140 genomes), dollo_* / make_figures_517 / make_figS_mash517_v2 / run_mashtree.sh / cluster_changes in
boletales_mash_517, seeds11, `analysis/trees`, `analysis/loci`, `my_code`, all `*.bak_*`, `bak_*/`, LSF `*.out`/`*.err`,
`__pycache__`, `tmp*/`, `quarantine/`, `mp_out/`, `round2/hmm*` (HMM tables), `.dmnd`, `.gb` (except one donor record), BAM/reads,
`.hmm` files in `tdna_hmm` and `hmm_v0` (rebuild with `hmmbuild` from the `.aln` files).

## Path replacements (copies only)

* the absolute path of the working project on the cluster -> `GMO_ROOT`.
  * Python (`from gmo_paths import GMO_ROOT`, inserted after the docstring together with a `sys.path` line that points to the
    repo root): the string literals were turned into f-strings (`{GMO_ROOT}`); literals in comments/docstring-like statements became `$GMO_ROOT`.
    30 scripts. `analysis/report_A/build_report.py` (uses py3.12+ f-string syntax) was edited by hand in the same way.
  * Shell (`: "${GMO_ROOT:=$(cd "$(dirname "$0")/<..>" && pwd)}"; export GMO_ROOT` inserted after the shebang): 20 scripts, literal -> `${GMO_ROOT}`.
    The default uses `$0`; when a script is submitted through `bsub` from another location, export `GMO_ROOT` first.
  * YAML / md / txt / tsv: literal -> `${GMO_ROOT}` (`$GMO_ROOT` in comments). Path columns in a few tables (assembly locations) therefore
    contain the text `${GMO_ROOT}` and must be expanded before use.
* the user-specific conda installation path -> `$(conda info --base)/etc/profile.d/conda.sh` (shell), `$CONDA_BASE/...` (docs).
  the user-specific `<conda base>/envs/<env>/bin/<tool>` -> `$(conda info --base)/envs/<env>/bin/<tool>` (shell) or the bare command name (Python:
  `diamond` in `ctdna_plants_identity/run_identity.py`; `python` in docstring "Run:" lines; the matching env must be active).
* Two non-mechanical edits for the pipeline: `pipeline/scripts/common.py` expands `${GMO_ROOT}` in `project_dir` (and sets a default
  `GMO_ROOT`), `pipeline/run_pipeline.sh` expands variables in values read from the config.
* Check: a grep of the repository for the cluster home prefix and the user name returns nothing. All `.py` files parse (py_compile with the `annot` interpreter), all `.sh` pass `bash -n`.

## Left as is

* Russian comments inside code and in some tables/docs (only `pipeline/README.md` was translated).
* Scripts were not run against data; nothing was executed beyond syntax checks.

## Orphans and open points

* Outputs with no generating script (inventory section 6): `reciprocal/table1_stats_dedup.tsv`, chs/odh `reverse_search_summary.tsv`,
  odh `expression/runs_meta.tsv`, chs `expression/summary.tsv`, CUS `lineage_nb` AU constraints, MMseqs2 clustering and `tdna_hmm` building,
  `boletales517_ncbi_summary.jsonl` (an `datasets summary genome accession` dump), the data-trees tarball. They are shipped as data; scripts still to be written or documented.
* (moved to `analysis/curated/quarantine.tsv` on 2026-10-16) `analysis/quarantine/quarantine.tsv` (512 bytes, read by the pipeline config) was copied as the single exception to the "no quarantine/" rule; remove it and the
  `quarantine_tsv` key if undesired.
* `analysis/pfam/Pfam-A.hmm` and an nr DIAMOND database are referenced by scripts but not shipped; `README.md` says where they must be.
* `refset/step2/queries.faa` count (389 vs 388 in Methods) still to be checked, as in the inventory.
* The manuscript text (`draft_v2_*.md`, `references_v2.md`, reference map) is deliberately not included; `paper/rebuild_all.sh` skips the manuscript step without it.
* Total is ~69 MB / ~1,660 files (above the ~50 MB target): `trees_final` (19 MB, AU-test constrained-run reports) and `refset/screen/step3` (10 MB) dominate. Candidates to prune: constrained-run `*.iqtree` reports in `trees_final/*/`, `hmmalign` variants (1.8 MB).
* GitHub owner `pzhurbenko`, repository `fungal-opine-genes`; citation line in README is a placeholder.
* ClipKIT 2.7.0 is a pip install (not in the conda export); BLAST+ version is not named in the Methods.

## Refresh 2026-10-13 (evening), after the review response

All sources listed above were re-copied (overwrite) with the same filters and the same path replacement (hand-edited files
`pipeline/scripts/common.py`, `pipeline/run_pipeline.sh`, `pipeline/README.md`, `analysis/report_A/build_report.py`, `paper/rebuild_all.sh` were kept as they were,
their sources being unchanged). Added:
* `paper/supplementary/`: rebuilt `build_supplementary.py` and the whole current `*.tsv` set, new Tables S1d, S11b, `TableS11_constraint_mapping.tsv` + legend,
  `data/Data_S2b_AU_constraints.tar.gz` and `data/au_constraints/MANIFEST.tsv` (the unpacked 12 MB constraint tree is not duplicated; it is inside the tarball).
* `paper/tables/`: `integration_status.py/.tsv`, `integration_status_summary.tsv`, `sensitivity_uncertain_loci.py/.tsv`, `intact_but_uncertain_loci.tsv`, updated catalogue and Tables 1-2.
* `paper/manuscript_v2/`: `build_docx_mdpi.py`, `reference_map.tsv`, updated builders; `paper/figures/`: `make_fig3_phylogram.py`, `make_figS10/12/13`, `Fig5_subfamilies_Tremble2024.tsv`.
* `analysis/funnel/carriers_by_order.py/.tsv`; `analysis/au_mapping/` (from `_au_mapping_tmp`: `build_au_mapping.py`, `parse.py`, `comp.py`);
  `analysis/pkerr14a_borders/` (the `.gb` record is not shipped; the script log is `borders_pkerr.log.txt`); `analysis/ensifer_sus_check/` (no `ncbi/` downloads);
  `analysis/data_checks_20261013/` (notes, scripts, small outputs); `analysis/dnds/extra/` with `cus_Bol/`, `cus_Bol_sp/`, `prep_cusBol.py`.
* Not shipped (as before): `PATCH_*.md` working notes, per-locus `round2/loci` and `round2/proteins` trees, mash/Dollo files in `boletales_mash_517`.
* Two scripts (`analysis/expression/summarize.py`, `analysis/dnds/extra/yn00_pairwise_loose.py`) have f-strings the tokenizer-based replacement could not rewrite; edited by hand to `{GMO_ROOT}`.

## Refresh 2026-10-13 night

Only the `paper/` build scripts and tables were re-copied (same filters, same path replacement); nothing else in the repo was touched. Changed:
* `paper/tables/`: `Table1_funnel.tsv` (the two identity columns merged into one), `build_catalog.py` (simplified IDENT dict) and its output `Table2_main.tsv`.
* `paper/manuscript_v2/build_manuscript_v2.py`: author block, affiliations and author contributions filled in.
* `paper/figures/`: `make_fig5_genera.py` (Figure 5 rework) and `Fig2_genus_counts.tsv`.
* Manuscript text (`draft_v2_*.md`, `manuscript_v2*.md`) and PNG figures are still not shipped, as before.

## Refresh 2026-10-13 (late), response to the review of v15

Same filters and path replacement (script kept outside the repo). Changed or added:
* `analysis/refset/screen/seeds_ext/reciprocal/`: new `scripts/02b_taxonomy_fix.py` (mandatory after `02_taxonomy.py`): organism names of nr titles are taken
  from the trailing bracket(s) only (`orgs_from_stitle`, also used by `03_reciprocal.py`), and names that match the taxdump only after truncation,
  single-word names and multi-kingdom homonyms are re-resolved by the taxid of the hit protein (NCBI esummary; cache `eutils_cache.json`,
  `name_taxid.json`, log `taxonomy_fix_log.tsv`, dropped description fragments `nonorganism_brackets.json`). Fixes the bacterial genus *Bosea*
  (now *Allobosea*) being read as the plant genus. Rebuilt `reciprocal.tsv`, `summary.tsv`, `h_gt0_order.tsv`, `taxonomy_stats.json`;
  the funnel (`analysis/funnel/funnel_numbers.txt`, `paper/tables/Table1_funnel.tsv`) is unchanged, `TableS3b_fate.tsv` differs only in best-hit organism names.
* `paper/tables/integration_status.py/.tsv/_summary.tsv`: the 2026-10-13 classifier (F2 from per-neighbour classifications, R2 -> supported, CHR by sequence role);
  the previous repo copy was the v14 version.
* `paper/supplementary/`: `build_supplementary.py` (Table S1a legend now states these rules; Table S3 is the full reverse search, 1,277 rows, replacing the old 116-row file),
  rebuilt tables, `Data_S1`, `data/Data_S4_dNdS_codeml.tar.gz`, new `make_zip.sh` (submission archive, 35 files).
* `analysis/dnds/extra/`: the four `*_conf` sets (loci with confirmed or supported integration only), `prep_conf.py`, `run_conf.sh`, `run_yn00pw_conf.sh`, updated `summary.tsv`.
* `docs/methods_extended.md`: synchronised with manuscript v16 (section numbers 2.1-2.11, integration-status rules, filtered dN/dS sets, references [M1]-[M26]).
* `paper/manuscript_v2/build_manuscript_v2.py`: current builder. Manuscript text and PNG figures are still not shipped.

## Cleanup 2026-10-16 (repository limited to what is needed to reproduce the published analyses)

Removed: `analysis/report_A`, `analysis/report_B`, `analysis/quarantine`, `analysis/data_checks_20261013`, `analysis/synteny`, `analysis/clustering`,
`analysis/conservation`, `paper/manuscript_v2`, `paper/build_paper_html.py`. Files that remaining scripts read were relocated first:
* tables of the carriers (`carriers_tree_based.tsv`, `carriers_summary.tsv`, `accessions.tsv`) -> `analysis/funnel/`;
  `island_JAHRGV010000066_genes.tsv` -> `analysis/islands/`; `bedulis_contig_depth.tsv` -> `analysis/read_coverage/`;
  `CUS_AU_test_per_lineage.iqtree.txt` -> `analysis/trees_final/CUS/`; `quarantine.tsv` -> `analysis/curated/` (the `quarantine_tsv` key of the pipeline configs follows);
  `split.gff` and `assembly_reports/` (copied from the working project) -> `analysis/loci_pipeline/checks/`.
* byte-identical copies already present elsewhere are used instead of duplicates: trees and leaf labels -> `analysis/trees_final/<FAM>/<FAM>.mafft.iqtree.treefile`, `<FAM>_labels.tsv`;
  SUS/ACS AU reports -> `trees_final/SUS/SUS.AU.iqtree`, `trees_final/ACS/ACS.lineage.AU.iqtree`; `island_gene_order.tsv` -> `analysis/island_order/`;
  `dnds_summary.tsv` -> `analysis/dnds/`; `expression_summary.tsv` -> `analysis/expression/`; `locus_coverage.tsv` -> `analysis/read_coverage/`;
  `contig_reproducibility.tsv` -> `analysis/indep_assemblies/`.
* `build_supplementary.py` takes `Extended_Methods.md` from `docs/methods_extended.md`.
Added: `analysis/susD_screen/`, `paper/figures/make_graphical_abstract_v2.py`, `docs/thresholds.md`.
