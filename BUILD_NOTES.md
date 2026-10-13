# BUILD_NOTES

Assembled 2026-10-13 from the working project (`GMO/`) following `paper/submission/repo_inventory.md`. Originals were
copied, never moved or edited; this directory (`GMO/repo/`) holds copies only. No `git init`, nothing pushed.

## What was copied from where (relative to the working project root)

| In repo | Source | Rule |
|---|---|---|
| `pipeline/` | `pipeline/` | all files except `__pycache__`; `scripts/common.py` and `run_pipeline.sh` got a small edit (below); `README.md` translated to English |
| `analysis/<dir>` | `analysis/<dir>` | in-paper folders of inventory section 2: refset/{plasmids,hmm_v0,neighbourhood,screen/{seeds_ext,step2,step3,tdna_hmm}} + refset/step2/queries*.faa + outgroup json, loci_pipeline/{addenda,guides,registry,round2 (tables only, plus `loci/ACS_JBBPRY010001778_L1/augustus.gff3`)}, trees_final, chs_fusarium, fusarium_odh_island, pleurotus_gdpd, fold, dnds, composition, conservation, islands, island_order, synteny/islands_hq, indep_assemblies, read_coverage, expression, boletales_mash_517, thresholds, funnel, ctdna_plants_identity, report_A, report_B, clustering/scripts, curated/*.fasta |
| file filter | | code (`.py .sh .R .yaml .md`) and small data (`.tsv .txt .faa .fasta .aln .json .jsonl .nwk .treefile .iqtree .fna .gff .gff3 .csv .ctl .phy .nex`) up to 400 KB; explicit larger inclusions: `all_screen.hmm` (8.9 MB), `funnel_grid.tsv`, `all_blast_stats.tsv`, `seeds_ext_best.tsv`, `reciprocal.tsv`, `round2/loci_table.tsv`, `step2_family_calls.tsv`, `registry/loci.tsv`, `step3/proteomes.txt`, `step3/genome_meta.tsv`, `plasmids/cds_table.tsv`, `islands_genes.tsv`, `donor/NZ_JAVIZC010000002.1.gb`, `boletales517_ncbi_summary.jsonl`, `TableS1b`, `TableS1_*` catalogue inputs, `Data_S1` |
| `data/` | `data/tsv/merged/all_blast_stats.tsv`, `data/tsv/rep_seqs/*.tsv`, `data/fungal_genomes/{gca_list.txt,fusarium/gca_list.txt,boletales_all/accessions_20261006.tsv}`, `data/fungal_proteomes/download.sh` | accession lists and pipeline inputs only |
| `paper/figures`, `paper/tables`, `paper/supplementary` | `paper/...` | build scripts, `common.py`, `export_tiff.py`, small inputs, supplementary tables S1-S17, Data S1/S3, `supplementary/data/trees`; no PNG/PDF/SVG/TIFF, no docx/xlsx |
| `paper/manuscript_v2/` | `paper/manuscript_v2/` | `build_manuscript_v2.py`, `build_docx.py`, `build_references.py` only; `paper/build_paper_html.py` |
| `paper/rebuild_all.sh` | `paper/rebuild_all.sh` | manuscript steps wrapped in `if [ -f $P/draft_v2_methods.md ]` because the manuscript text is not shipped |
| `docs/methods_extended.md` | `paper/manuscript_v2/methods_extended.md` | |
| `env/` | `conda env export` of annot, phylo_env, mashtree_env, minimap2_env, colabfold, PlantCad2, orthofinder, hmmer_pfam | `<env>.yml` = `--from-history`, `<env>.pinned.yml` = `--no-builds` |
| `LICENSE`, `LICENSE-DATA`, `.gitignore`, `README.md`, `gmo_paths.py` | written new | |

Excluded: inventory section 5, donor_meta, donor_search, generax, alien_scan, distribution, orf358_origin, accF, acc_test,
boletales_mash (140 genomes), dollo_* / make_figures_517 / make_figS_mash517_v2 / run_mashtree.sh / cluster_changes in
boletales_mash_517, seeds11, `analysis/trees`, `analysis/loci`, `my_code`, all `*.bak_*`, `bak_*/`, LSF `*.out`/`*.err`,
`__pycache__`, `tmp*/`, `quarantine/`, `mp_out/`, `round2/hmm*` (HMM tables), `.dmnd`, `.gb` (except one donor record), BAM/reads,
`.hmm` files in `tdna_hmm` and `hmm_v0` (rebuild with `hmmbuild` from the `.aln` files).

## Path replacements (copies only)

* the absolute path of the working project (GMO) and of the Claude project directory (GMO_matweeva) on the cluster -> `GMO_ROOT`.
  * Python (`from gmo_paths import GMO_ROOT`, inserted after the docstring together with a `sys.path` line that points to the
    repo root): the string literals were turned into f-strings (`{GMO_ROOT}`); literals in comments/docstring-like statements became `$GMO_ROOT`.
    30 scripts. `analysis/report_A/build_report.py` (uses py3.12+ f-string syntax) was edited by hand in the same way.
  * Shell (`: "${GMO_ROOT:=$(cd "$(dirname "$0")/<..>" && pwd)}"; export GMO_ROOT` inserted after the shebang): 20 scripts, literal -> `${GMO_ROOT}`.
    The default uses `$0`; when a script is submitted through `bsub` from another location, export `GMO_ROOT` first.
  * YAML / md / txt / tsv: literal -> `${GMO_ROOT}` (`$GMO_ROOT` in comments). Path columns in a few tables (assembly locations) therefore
    contain the text `${GMO_ROOT}` and must be expanded before use.
* the user-specific miniforge `etc/profile.d/conda.sh` -> `$(conda info --base)/etc/profile.d/conda.sh` (shell), `$CONDA_BASE/...` (docs).
  the user-specific `miniforge3/envs/<env>/bin/<tool>` -> `$(conda info --base)/envs/<env>/bin/<tool>` (shell) or the bare command name (Python:
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
* `analysis/quarantine/quarantine.tsv` (512 bytes, read by the pipeline config) was copied as the single exception to the "no quarantine/" rule; remove it and the
  `quarantine_tsv` key if undesired.
* `analysis/pfam/Pfam-A.hmm` and an nr DIAMOND database are referenced by scripts but not shipped; `README.md` says where they must be.
* `refset/step2/queries.faa` count (389 vs 388 in Methods) still to be checked, as in the inventory.
* The manuscript text (`draft_v2_*.md`, `references_v2.md`, reference map) is deliberately not included; `paper/rebuild_all.sh` skips the manuscript step without it.
* Total is ~69 MB / ~1,660 files (above the ~50 MB target): `trees_final` (19 MB, AU-test constrained-run reports) and `refset/screen/step3` (10 MB) dominate. Candidates to prune: constrained-run `*.iqtree` reports in `trees_final/*/`, `hmmalign` variants (1.8 MB).
* GitHub owner `pzhurbenko`, repository `fungal-opine-genes`; citation line in README is a placeholder.
* ClipKIT 2.7.0 is a pip install (not in the conda export); BLAST+ version is not named in the Methods.
