# Opine metabolism genes horizontally transferred into fungal genomes

Code, small input data and derived tables behind the article "Opine metabolism genes horizontally transferred into fungal genomes" (Biology, MDPI, submitted 2026).
This repository holds the code, small input data and derived tables behind a screen of 6,123 fungal proteomes and 2,215 fungal assemblies, which yielded
candidate loci that fall into seven independent lineages of bacterial origin (gene families named ACS, CUS and SUS after the
T-DNA opine synthase genes they resemble, plus the chs-like and odh-like loci of Fusarium and relatives). The pipeline builds gene models on the real contigs, classifies the genomic context, and the paper scripts rebuild
all figures, tables and supplementary files.

**Citation:** to be added on publication.

## Repository layout

| Path | Content |
|---|---|
| `pipeline/` | Locus pipeline (registry of tblastn hits -> contigs -> miniprot/AUGUSTUS gene models -> HMM validation -> diamond-vs-nr context -> `loci_table.tsv`): `run_pipeline.sh`, `config*.yaml`, `scripts/` (`01`-`08`, `common.py`), `README.md` |
| `analysis/` | One folder per analysis step used in the paper (same relative structure as in the working project): `refset/` (seeds, plasmid catalogue, protein screen, step 2/3 screens of other orders), `funnel/`, `loci_pipeline/` (registry, guides, round-2 table, `checks/` = miniprot re-derivation `split.gff` and NCBI assembly reports used by the integration-status rules), `trees_final/` (alignments, ML trees, AU tests), `chs_fusarium/`, `fusarium_odh_island/`, `pleurotus_gdpd/`, `fold/` (ColabFold, FoldMason), `dnds/`, `composition/`, `islands/`, `island_order/`, `indep_assemblies/`, `read_coverage/`, `expression/`, `boletales_mash_517/`, `thresholds/`, `ctdna_plants_identity/`, `susD_screen/` (susD/susL check), `au_mapping/`, `pkerr14a_borders/`, `ensifer_sus_check/`, `curated/` (curated query sets and `quarantine.tsv`) |
| `paper/figures/` | Figure scripts (`make_fig*.py`, `make_figS*.py`, `common.py`), `export_tiff.py`, small plot inputs; `make_graphical_abstract_v2.py` (graphical abstract; uses the Nimbus Sans fonts of `/usr/share/fonts/urw-base35` if present, otherwise DejaVu Sans) |
| `paper/tables/` | `build_table1.py`, `build_catalog.py`, `check_order.py` and the main tables / catalogue inputs |
| `paper/supplementary/` | `build_supplementary.py`, supplementary tables S1-S17, Data S1 (proteins of the loci), Data S2b (AU constraints, tar.gz), Data S3 (screened assemblies), tree files |
| `paper/rebuild_all.sh` | Rebuilds figures (300 dpi PNG and 600 dpi TIFF) and the supplementary tables |
| `docs/methods_extended.md` | Extended Methods: full, parameter-level description of the methods (section numbers follow Section 2 of the article) |
| `docs/thresholds.md` | Table of every numerical threshold of the classification, with the range over which the result is unchanged (copy of Table S4) |
| `analysis/susD_screen/` | DIAMOND screen of the annotated proteomes for *Agrobacterium* susD (ARU12249.1) and susL (AAZ50401.1): `screen_one.sh`, `best_hits.tsv`, `SUMMARY.txt`, `README.md` |
| `data/` | Small inputs only (BLAST hit table, representative sequences, assembly accession lists, proteome download script) |
| `env/` | Conda environments: `<env>.yml` (`--from-history`) and `<env>.pinned.yml` (full versions, no builds) |
| `gmo_paths.py` | Defines `GMO_ROOT` |
| `BUILD_NOTES.md` | How this repository was assembled, path replacements, open points |

## Project root (`GMO_ROOT`)

All scripts take the project root from the environment variable `GMO_ROOT` (shell: `${GMO_ROOT}`; Python:
`from gmo_paths import GMO_ROOT`). If it is not set, the root of this repository is used.

```bash
export GMO_ROOT=/path/to/fungal-opine-genes     # optional; default = repository root
```

Large inputs that are not in the repository are expected under `$GMO_ROOT/data/` (genome assemblies and proteomes downloaded from
NCBI, reads from SRA, a DIAMOND database of NCBI nr) and `$GMO_ROOT/analysis/pfam/Pfam-A.hmm` (Pfam-A, used by the domain scans
and tree scripts). Some tab-separated tables contain `$GMO_ROOT/...` in path columns (local genome locations); resolve them with
`envsubst` or by `sed "s|\$GMO_ROOT|$GMO_ROOT|g"` before using them.

## Environments and tools

Create the environments with `conda env create -f env/<env>.yml` (unpinned) or `env/<env>.pinned.yml` (exact versions).

| Tool | Version | Environment |
|---|---|---|
| DIAMOND | 2.2.6 | `annot` |
| HMMER | 3.4 | `annot`, `hmmer_pfam` |
| miniprot | 0.18 | `annot` |
| AUGUSTUS | 3.1 | `annot` |
| MAFFT | 7.526 | `annot`, `phylo_env` |
| ClipKIT | 2.7.0 | `pip install clipkit==2.7.0` (installed with pip, not captured in the conda files) |
| IQ-TREE | 3.1.1 | `phylo_env` |
| minimap2 | 2.30 | `minimap2_env` |
| SAMtools | 1.24 | `annot` |
| ColabFold | 1.6.2 | `colabfold` |
| FoldMason | 4.dd3c235 | `phylo_env` |
| MMseqs2 | 18.8cc5c | `orthofinder` |
| Mash / mashtree | 2.3 / 1.2.0 | `mashtree_env` |
| PAML | 4.10.7 | `mashtree_env` |
| Python stack (Python 3.10, matplotlib, pandas, Biopython, python-docx, openpyxl; used only for figures, tables and the supplementary files) | see file | `paper_build` |

## Running the locus pipeline

The pipeline and several analysis scripts submit jobs to an LSF cluster (`bsub`); on another system, run the commands that each script wraps directly, with the same environments.

```bash
cd $GMO_ROOT/pipeline
conda activate annot
./run_pipeline.sh registry && ./run_pipeline.sh fetch && ./run_pipeline.sh guides
./run_pipeline.sh annotate 1 && ./run_pipeline.sh validate 1 && ./run_pipeline.sh diamond 1 \
  && ./run_pipeline.sh classify 1 && ./run_pipeline.sh summarize 1
```
The cluster steps use LSF (`bsub`); see `pipeline/README.md` and `pipeline/config.yaml` (queue, resources, nr database path).
The individual analyses in `analysis/` are shell/Python scripts that were run on the same cluster; each folder reads its inputs
from the paths given at the top of its scripts.

## Regenerating figures, tables and supplementary files

```bash
conda activate paper_build
bash $GMO_ROOT/paper/rebuild_all.sh
```
This runs the figure scripts (Figures 1-5, 600 dpi TIFF export) and `paper/supplementary/build_supplementary.py` (Tables S1-S18,
Supplementary Materials). Table 1 and the loci catalogue
are produced by `paper/tables/build_table1.py` and `build_catalog.py`; the funnel (Table 1 numbers) by `analysis/funnel/build_funnel.py`.

## Which folder produces what

Section numbers are those of the Methods (`docs/methods_extended.md`). Scripts are listed relative to the folder; "–" means that no table or figure of the paper is derived directly from that step, or that this could not be established from the scripts.

| Methods section | Folder(s) | Main script(s) | Table / figure of the paper |
|---|---|---|---|
| 2.1 Seed proteins and protein-level screen | `analysis/refset/screen/seeds_ext/`, `analysis/refset/screen/step3/`, `analysis/funnel/`, `analysis/thresholds/`, `analysis/ctdna_plants_identity/`, `analysis/susD_screen/` | `seeds_ext/build_seeds_ext.py`, `screen_one.sh`, `aggregate_seeds_ext.py`, `reciprocal/scripts/`, `funnel/build_funnel.py`, `funnel/carriers_by_order.py`, `thresholds/sensitivity_funnel.py`, `thresholds/build_summary.py`, `ctdna_plants_identity/run_identity.py`, `susD_screen/screen_one.sh` | Table 1 (`funnel/build_funnel.py`), Figure 2 (`paper/figures/make_fig2.py`), Tables S2, S3, S4, S17a/b, Figure S1 |
| 2.2 Genome-level screen | `analysis/refset/screen/step2/`, `analysis/refset/screen/tdna_hmm/`, `analysis/chs_fusarium/stage2/`, `analysis/fusarium_odh_island/scan/`, `analysis/pleurotus_gdpd/`, `data/` | `step2/aggregate_step2.py`, `step2/summarize_step2.py`, `chs_fusarium/stage2/parse_loci.py`, `fusarium_odh_island/scan/parse_loci.py` | Figure 1 (`paper/figures/make_fig1.py`), Tables S6, S7, S10a, S10b, S13a, S13b |
| 2.3 Reference plasmid catalogue and contamination control | `analysis/refset/plasmids/`, `analysis/pkerr14a_borders/`, `analysis/refset/screen/step2/` | `plasmids/build_plasmid_catalog.py`, `pkerr14a_borders/borders_pkerr.py` | Tables S5a, S5b |
| 2.4 Gene models and the catalogue of loci | `pipeline/`, `analysis/loci_pipeline/`, `analysis/curated/`, `paper/tables/` | `pipeline/run_pipeline.sh`, `pipeline/scripts/01`-`08`, `paper/tables/build_table1.py`, `paper/tables/build_catalog.py` | Table 2 (`build_catalog.py`), Tables S1a-S1c (via `paper/supplementary/build_supplementary.py`), Data S1 |
| 2.5 Tests against bacterial contamination | `analysis/loci_pipeline/` (neighbour classification), `analysis/ensifer_sus_check/`, `analysis/fusarium_odh_island/`, `analysis/indep_assemblies/`, `analysis/read_coverage/` | `pipeline/scripts/06_classify_neighbors.py`, `indep_assemblies/run_indep.sh`, `indep_assemblies/summarize.py`, `read_coverage/map_one.sh`, `read_coverage/summarize.py`, `paper/tables/integration_status.py` | Tables S1d, S8, S9, S13c, S18; Figure S12 (`make_figS12_bedulis.py`) |
| 2.6 Phylogenetic trees, lineage assignment and AU tests | `analysis/trees_final/`, `analysis/au_mapping/`, `analysis/chs_fusarium/`, `analysis/fusarium_odh_island/`, `analysis/pleurotus_gdpd/`, `analysis/fold/` | `trees_final/run_family_final.sh`, `trees_final/make_constraints.py`, `trees_final/au_test.sh`, `trees_final/au_lineage_*.py`, `au_mapping/build_au_mapping.py`, `fold/run_fold.sh` | Figure 3 (`make_fig3.py`, `make_fig3_phylogram.py`), Tables S11, S11b, Data S2/S2b, Figures S2-S8 (`make_figS_trees.py`, `make_figS_AU.py`), Figure S9 |
| 2.7 Genomic context and islands | `analysis/islands/`, `analysis/island_order/`, `analysis/fusarium_odh_island/` | `islands/islands.py`, `island_order/island_order.py`, `island_order/analyze.py` | Figure 4 (`make_fig4.py`), Table S12, Figures S10, S13 (`make_figS10_islands.py`, `make_figS13_fp133.py`) |
| 2.8 Compositional analysis | `analysis/composition/` | `scripts/composition.py`, `scripts/composition_nosODH.py`, `scripts/fig_top.py` | Table S14 |
| 2.9 Test for selection (dN/dS) | `analysis/dnds/` | `prep.py`, `run_dnds.sh`, `parse_paml.py`, `extra/` | Table S15, Figure S11 (`make_fig5.py`), Data S4 |
| 2.10 Expression analysis (RNA-seq) | `analysis/expression/`, `analysis/chs_fusarium/expression/`, `analysis/fusarium_odh_island/expression/` | `expression/summarize.py`, `expression/map_rna.sh` | Table S16 |
| 2.11 Species clusters of Boletales | `analysis/boletales_mash_517/` | `run_mash.sh`, `clusters.py`, `make_fig2_genus.py` | Figure 5 (`make_fig5_genera.py`), Table S6 |
| (summary of all lineages) | `paper/tables/`, `paper/figures/` | `paper/tables/build_catalog.py`, `paper/figures/make_graphical_abstract_v2.py` | Table S18 (`paper/tables/TableS18_evidence_matrix.tsv`), graphical abstract; Table 3 (-) |

## Data availability

The assemblies and proteomes analysed are public NCBI records; their accessions are listed in
`paper/supplementary/Data_S3_screened_assemblies.tsv` (Data S3) and, for subsets, in `data/fungal_genomes/` and
`analysis/refset/screen/`. Proteins of all loci of bacterial origin are in `paper/supplementary/Data_S1_proteins_of_loci_of_bacterial_origin.faa`.
Further tables: Table S1d (sensitivity of the lineage counts to the integration status of the loci; computed by
`paper/tables/integration_status.py` and `sensitivity_uncertain_loci.py`, outputs `integration_status*.tsv`), Table S11b (mapping of each
AU constraint to the biological alternative it tests; `analysis/au_mapping/`, `paper/supplementary/TableS11_constraint_mapping.tsv` with legend) and
Data S2b (constraint files, alignments and reports of the AU tests: `paper/supplementary/data/Data_S2b_AU_constraints.tar.gz`, file list in
`paper/supplementary/data/au_constraints/MANIFEST.tsv`). Additional checks from the review response are in `analysis/pkerr14a_borders/`,
`analysis/ensifer_sus_check/`, `analysis/loci_pipeline/checks/` and `analysis/dnds/extra/` (dN/dS of cus-like Boletales genes, `cus_Bol*`).
Raw genomes, reads and the NCBI nr database are not redistributed.

## Licences

Code: MIT (`LICENSE`). Data and derived tables, alignments, trees and sequence sets: CC BY 4.0 (`LICENSE-DATA`).
