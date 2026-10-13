# fungal-opine-genes

Opine metabolism genes of rhizobial (Agrobacterium/Rhizobium T-DNA) origin in fungal genomes. This repository holds the code,
small input data and derived tables behind a screen of 6,123 fungal proteomes and 2,215 fungal assemblies, which yielded
candidate loci that fall into seven independent lineages of bacterial origin (gene families named ACS, CUS and SUS after the
T-DNA opine synthase genes they resemble, plus the chs-like and odh-like loci of Fusarium and relatives). The pipeline builds gene models on the real contigs, classifies the genomic context, and the paper scripts rebuild
all figures, tables and supplementary files.

**Citation:** [TO FILL: citation]

## Repository layout

| Path | Content |
|---|---|
| `pipeline/` | Locus pipeline (registry of tblastn hits -> contigs -> miniprot/AUGUSTUS gene models -> HMM validation -> diamond-vs-nr context -> `loci_table.tsv`): `run_pipeline.sh`, `config*.yaml`, `scripts/` (`01`-`08`, `common.py`), `README.md` |
| `analysis/` | One folder per analysis step used in the paper (same relative structure as in the working project): `refset/` (seeds, plasmid catalogue, protein screen, step 2/3 screens of other orders), `funnel/`, `loci_pipeline/` (registry, guides, round-2 table), `trees_final/` (alignments, ML trees, AU tests), `chs_fusarium/`, `fusarium_odh_island/`, `pleurotus_gdpd/`, `fold/` (ColabFold, FoldMason), `dnds/`, `conservation/`, `composition/`, `islands/`, `island_order/`, `synteny/islands_hq/`, `indep_assemblies/`, `read_coverage/`, `expression/`, `boletales_mash_517/`, `thresholds/`, `ctdna_plants_identity/`, `report_A/`, `report_B/`, `clustering/scripts/`, `curated/` |
| `paper/figures/` | Figure scripts (`make_fig*.py`, `make_figS*.py`, `common.py`), `export_tiff.py`, small plot inputs |
| `paper/tables/` | `build_table1.py`, `build_catalog.py`, `check_order.py` and the main tables / catalogue inputs |
| `paper/supplementary/` | `build_supplementary.py`, supplementary tables S1-S17, Data S1 (proteins of the loci), Data S2b (AU constraints, tar.gz), Data S3 (screened assemblies), tree files |
| `paper/rebuild_all.sh` | Rebuilds figures (300 dpi PNG and 600 dpi TIFF), supplementary tables and, if the manuscript text is present, the manuscript files |
| `docs/methods_extended.md` | Extended Methods |
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
| Python plotting/build stack (pandas 2.2.3, matplotlib 3.9.2, Biopython 1.84, ete3 3.1.3, openpyxl 3.1.5, python-docx 1.2.0) | see file | `PlantCad2` |

## Running the locus pipeline

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
conda activate PlantCad2
bash $GMO_ROOT/paper/rebuild_all.sh
```
This runs the figure scripts (Figures 1-5, 600 dpi TIFF export), `paper/supplementary/build_supplementary.py` (Tables S1-S17,
Supplementary Materials) and, only if the manuscript sources are present, the manuscript build. Table 1 and the loci catalogue
are produced by `paper/tables/build_table1.py` and `build_catalog.py`; the funnel (Table 1 numbers) by `analysis/funnel/build_funnel.py`.

## Data availability

The assemblies and proteomes analysed are public NCBI records; their accessions are listed in
`paper/supplementary/Data_S3_screened_assemblies.tsv` (Data S3) and, for subsets, in `data/fungal_genomes/` and
`analysis/refset/screen/`. Proteins of all loci of bacterial origin are in `paper/supplementary/Data_S1_proteins_of_loci_of_bacterial_origin.faa`.
Further tables: Table S1d (sensitivity of the lineage counts to the integration status of the loci; computed by
`paper/tables/integration_status.py` and `sensitivity_uncertain_loci.py`, outputs `integration_status*.tsv`), Table S11b (mapping of each
AU constraint to the biological alternative it tests; `analysis/au_mapping/`, `paper/supplementary/TableS11_constraint_mapping.tsv` with legend) and
Data S2b (constraint files, alignments and reports of the AU tests: `paper/supplementary/data/Data_S2b_AU_constraints.tar.gz`, file list in
`paper/supplementary/data/au_constraints/MANIFEST.tsv`). Additional checks from the review response are in `analysis/pkerr14a_borders/`,
`analysis/ensifer_sus_check/`, `analysis/data_checks_20261013/` and `analysis/dnds/extra/` (dN/dS of cus-like Boletales genes, `cus_Bol*`).
Raw genomes, reads and the NCBI nr database are not redistributed.

## Licences

Code: MIT (`LICENSE`). Data and derived tables, alignments, trees and sequence sets: CC BY 4.0 (`LICENSE-DATA`).
