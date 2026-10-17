# GMO loci pipeline — gene models and genomic context for candidate T-DNA genes in fungi

Version 1, 2026-09-05.

## Why

The original fungal candidate "proteins" (ACS, CUS, SUS, ORF358) are **concatenations of tblastn HSPs**: translations of
contig pieces that aligned to a reference protein. This is fine for screening and clustering, but:

* the sequence is truncated to the shape of the query, hence biased towards similarity to the bacterial/plant protein
  (towards the HGT hypothesis);
* there is no DNA, so GC/codon usage, introns, dN/dS and neighbouring genes cannot be computed;
* a truncated HSP cannot distinguish a pseudogene from a divergent gene;
* about 2% of the concatenations carry duplications at HSP junctions, and some contain stop codons of unknown origin.

The pipeline replaces these proxies with real gene models on the real contigs and at the same time collects what is needed for
horizontal-transfer criteria: introns, neighbouring genes and their origin, GC profile, island borders.

## Key ideas (read before running)

1. **Annotate all loci, not only representatives.** 2007 WGS contigs instead of 204. Computation is cheap, downloading takes
   an hour or two. Rows flagged `Was_deleted=yes` ("Low homology", "Too short") are processed too, with a flag: filtering is
   done after the models and HMMs, not on HSPs.
2. **Locus = family x contig x HSP cluster** closer than 5 kb. One contig can carry several loci and several families.
3. **A contig <= 1 Mb is downloaded entirely**, a longer one as a +-100 kb window. The gene model is built in a +-15 kb window;
   neighbours are searched adaptively (15 -> 50 -> 100 kb -> contig end) until three genes with >= 2 exons are collected on each
   side. The distance at which this succeeded is itself a measure of "insertion into an island vs. into an ordinary region".
4. **Three tiers of guide proteins for miniprot; the tier is recorded in every model:**
   * T1 — the original Agro_/Plant_ references (after quarantine);
   * T2 — agro hits from `rep_seqs` (column Subj_Accession_agro), keeping only those >= 40% identical to T1; this removes
     chromosomal paralogues of the same fold (ornithine cyclodeaminases, reductases, TIR proteins); table `guides/<FAM>_T2_review.tsv`;
   * T3 — full-length fungal models that passed the check in the previous round (intact only, HMM pass, >= 30% to T1/T2, not in a
     bacterial context) plus NCBI-annotated proteins that passed the HMM.
   Fungal guides are allowed because miniprot uses them only to build the gene structure; each model proves its homology to
   T-DNA separately (item 6). Two rounds are enough.
   **Reference cleaning** (`guides.exclude_T1`, 2026-09-05): chromosomal paralogues and fragments were removed from T1 (SUS: 4 OCD;
   ACS: Sphagnum, a Parasponia fragment, 3 unassignable GDPD/PI-PLC singletons; CUS: a 138-aa DUF3830 protein). They are kept in
   `<FAM>_T1_excluded.faa` as the outgroup pool for trees. Check: `scripts/check_T1_homogeneity.py` -> `<FAM>_T1_check.tsv`;
   diamond of the references against nr is in `refcheck/`. Effect: the ACS HMM stopped admitting 145 yeast Gde1p (native GDPD of S. cerevisiae).
   **Subgroups** (`guides.subgroups`): ACS = acs / gdpd_meso / kerr14; CUS = cus / mis / rcorf12 / cusrel181; SUS, ORF358 — one each.
   One HMM per subgroup; a model gets its `subgroup` from the best HMM. gdpd_meso is a small bacterial GDPD subfamily
   (Mesorhizobium/Neorhizobium), closest to Boletales-ACS (58-65%); the classical pTi acs is only 25-29% identical to them.
5. **Two independent predictors.** miniprot (splice-aware protein alignment) gives the gene structure, frameshifts and stops
   explicitly. AUGUSTUS ab initio with fungal parameters per class (+ CDSpart hints from miniprot) gives an independent model and
   all neighbours. T-DNA genes of bacterial origin may be **intronless**; AUGUSTUS with fungal parameters tends to cut them with
   false introns, so when the two disagree we trust miniprot for T1/T2.
6. **Homology is checked independently of the guide:** one HMM per subgroup (MAFFT -> hmmbuild) **plus a competing outgroup HMM**
   (`guides.outgroups`: native fungal paralogues — yeast Gde1p, Aspergillus GDPD — and chromosomal bacterial GDPD/OCD of
   Agrobacterium/Rhizobium). Every model is run against all of them; `subgroup` = best HMM; if the best is OUTGROUP, the status is
   `not_homolog`. This removes the need to narrow the HMMs: 147 yeast Gde1p go to the outgroup, Boletales-ACS stay in gdpd_meso.
   Statuses: `intact` / `pseudogene_like` (frameshift or stop) / `fragment` (< 60% of the T1 median) / `not_homolog` (OUTGROUP is
   better, or HMM E > 1e-5 at identity < 25%, or an AUGUSTUS model without HMM) / `no_model`.
   **Fallback model**: miniprot finds no seeds at ~40% identity (CUS of Epicoccum with HMM E 1e-70 would have remained without a model);
   if miniprot is empty, the model is the AUGUSTUS gene overlapping the HSP (`model_source=augustus`), and it must pass the HMM.
   `identity` is to the guide that built the model (for T3 often 1.00 to itself — uninformative); look at `identity_T1T2`.
7. **Context = the main contamination control.** Neighbours are classified by the best diamond hits against nr (a shared database
   for the cluster, without taxonomy -> organism from the header -> kingdom via NCBI Taxonomy, cached). `context_call`:
   fungal_context / mixed / bacterial_context / no_neighbors / short_contig. Plus GC of the context vs GC of the assembly.
8. **Diamond against nr runs as one job for all proteins of the round:** the time is determined by the pass over the 371 GB
   database, not by the number of queries. `-b 24 -n 64` -> 20-30 min instead of 2 h with `-b 6 -n 16`.

## Layout

```
pipeline/
  config.yaml            all parameters and paths (the only thing to edit for another machine; project_dir = ${GMO_ROOT})
  run_pipeline.sh        stage driver (see the file header)
  scripts/
    common.py            fasta/tsv, rate-limited E-utilities with retries, config
    01_build_registry.py all_blast_stats.tsv -> registry/loci.tsv, contigs.tsv, protein_subjects.tsv
    02_fetch_contigs.py  contigs (efetch; bypass of the search index for fresh records) + assemblies.tsv (GCA, annotated?, GC)
    03_build_guides.py   guides/<FAM>_T1.faa, _T2.faa (+review), _T3.faa (from the previous round), _guides.faa
    04_annotate_locus.py one locus: miniprot + AUGUSTUS + neighbours + GC -> round<N>/loci/<locus_id>/locus.json
    05_validate_models.py family HMMs, hmmsearch, statuses -> round<N>/validation.tsv
    06_classify_neighbors.py collect -> diamond (LSF) -> parse: neighbors_class.tsv, locus_context.tsv
    07_summarize.py      round<N>/loci_table.tsv, proteins/<FAM>_models.{faa,fna}, guides/<FAM>_T3.faa
    08_add_step2_loci.py adds loci found in the step-2 screen of other orders to the registry
analysis/loci_pipeline/  (out_dir) registry/ contigs/ guides/ round1/ round2/ logs/
```

## How to run

```bash
export GMO_ROOT=/path/to/project            # default: repository root (see ../gmo_paths.py)
source $(conda info --base)/etc/profile.d/conda.sh; conda activate annot   # env with all tools (see config; ../env/annot.yml)
export NCBI_API_KEY=...                     # optional, speeds up NCBI 3x
cd $GMO_ROOT/pipeline
./run_pipeline.sh registry                  # seconds
./run_pipeline.sh fetch                     # ~1-2 h, resumable (skips what is downloaded)
./run_pipeline.sh guides                    # minutes; REVIEW guides/<FAM>_T2_review.tsv
./run_pipeline.sh annotate 1                # LSF job array, waits for completion
./run_pipeline.sh validate 1
./run_pipeline.sh diamond 1                 # LSF, one job
./run_pipeline.sh classify 1
./run_pipeline.sh summarize 1               # -> round1/loci_table.tsv ; guides/<FAM>_T3.faa
./run_pipeline.sh guides && ./run_pipeline.sh annotate 2 && ... summarize 2     # round 2 with T3
```
Test on a subset: `fetch --only contigs.txt`, `annotate 1 --only loci.txt --local` (no LSF, on the login node).

## Output to look at

`round<N>/loci_table.tsv` — one row per locus. Key columns: `status`, `context_call`, `n_exons`, `frameshift`,
`stop_codons`, `protein_len`, `len_frac_of_T1_median`, `guide_tier`, `identity`, `hmm_evalue`, `hmm_other_family_evalue`
(whether the family was mixed up), `augustus_target_exons` (agreement with miniprot), `nb_fungi/nb_bacteria/nb_TE`,
`window_left_bp/right_bp` (how far the window had to be extended), `gc_context`, `gc_delta_vs_assembly`, `other_fams_on_contig`
(clusters of T-DNA-like genes), `was_deleted` (the original user mark).

`round<N>/proteins/<FAM>_models.faa/.fna` — new curated sets of full-length proteins and CDS instead of HSP concatenations:
the input for alignments, trees and codon analysis.

## Known pitfalls

* Fresh WGS records (JB...) are not in the esearch/esummary index -> `02` downloads them directly by efetch on `ACC.1`.
* `bsub -R rusage[mem=X]` — X is per slot; for `-n 16` it is x16. Do not set 60000.
* An LSF job array does not accept `bsub -K`; the driver waits for the array with `bwait -w "ended(<id>)"` (fallback — polling bjobs).
  Right after bsub `bjobs <id>` may return empty — that does not mean "finished".
* miniprot: the `##STA` line (translation) comes BEFORE its mRNA line. The script accounts for it.
* miniprot finds no seeds at identity < 25-30%: "no model from T1" is not "not a homologue". The HMM and the tree decide.
* The diamond nr database was built without `--taxonmap`: the taxon comes only from the hit header.
* AUGUSTUS `--species` as a list: `ls $AUGUSTUS_CONFIG_PATH/species`; there is no own model for Dothideomycetes, botrytis_cinerea is used.
* Scripts for bsub must be on a shared disk (not in /tmp of the login node).

## What the pipeline does NOT do (the next steps of the paper)

Trees with Agro/Plant anchors and native fungal paralogues; synteny with related genomes; read coverage from SRA; expression;
dN/dS; structural alignment (ESMFold/FoldMason) for ORF358 and divergent ACS. The input for all of these is
`proteins/<FAM>_models.*` and `loci_table.tsv`.
