# Ensifer sp. NM-2 sus-like OCD gene (WP_107028680.1): genuine rhizobial gene or contamination?

Date: 2026-10-13. Review point M12 / manuscript 4.5. Analysis "b" of the post-review plan.
All records in `ncbi/`, script `scripts/analyze_context.py`, tables in `results/`.

## 1. Identity of the sequence

| item | value |
|---|---|
| Tree leaf (paper/supplementary/data/trees/sus.treefile) | `REF_T1_sus__Agro_WP_107028680_1_hypothetical_protein_Ensifer_sp_NM-2` |
| Protein | WP_107028680.1 (= PSS59732.1), 361 aa, PGAP name "Rossmann-fold NAD(P)-binding domain-containing protein" (older: hypothetical protein) |
| Organism | *Ensifer* sp. NM-2 (taxid 2109730), Rhizobiaceae; submitters' intended name "Ensifer ninsenella sp. nov." (Lund et al., Aarhus University, unpublished) |
| Assembly | GCF_003024455.1 / GCA_003024455.1 (ASM302445v1), released 2018-03-28, contig level, 223 contigs, 8.07 Mb, N50 275 kb, SPAdes 3.11, Illumina MiSeq, 49.7x; NCBI taxonomy check: OK |
| BioSample | SAMN08661303: strain NM-2, isolation source = soil, sample type = soil, Denmark, collected 2017-02-02 (no host given) |
| Source contig | NZ_PVZR01000036.1 (= PVZR01000036.1), "NODE_36_length_60479_cov_22.693233", 60,479 bp, 51 CDS + 1 tRNA |
| Gene coordinates | 45,794..46,879, minus strand, 1,086 bp = 361 codons + stop, single CDS, no introns, locus_tag C6558_RS36880 |

IPG (identical protein group) lists only this one genome: no identical protein in any other assembly, bacterial or fungal.

## 2. Genomic context on the contig (±8 genes; full list in results/all_cds_on_contig.tsv)

| rank | protein | strand | aa | GC% | product (PGAP) |
|---|---|---|---|---|---|
| -8 | WP_107028672.1 | - | 389 | 60.6 | FMNH2-dependent alkanesulfonate monooxygenase |
| -7 | (pseudo) | + | – | 52.1 | dienelactone hydrolase-related enzyme |
| -6 | WP_107028675.1 | - | 439 | 56.8 | sensor histidine kinase |
| -5 | WP_107028676.1 | - | 224 | 57.2 | response regulator transcription factor |
| (tRNA) | | - | | | tRNA-Ile (40,531..40,606) |
| -4 | WP_245449331.1 | - | 298 | 53.0 | alpha/beta fold hydrolase |
| -3 | WP_107028677.1 | + | 299 | 56.4 | LysR family transcriptional regulator |
| -2 | WP_107028678.1 | - | 213 | 54.4 | hydrolase |
| -1 | WP_107028679.1 | + | 186 | 56.3 | FMN reductase |
| **0** | **WP_107028680.1** | **-** | **361** | **58.5** | **sus-like OCD (this gene)** |
| +1 | WP_107028681.1 | + | 454 | 55.9 | hypothetical protein |
| +2 | WP_107028682.1 | + | 244 | 54.0 | hypothetical protein |
| +3 | WP_107028696.1 | - | 297 | 62.2 | LysR family transcriptional regulator |
| +4 | WP_107028683.1 | + | 283 | 58.8 | transporter substrate-binding domain-containing protein |
| +5 | WP_107028684.1 | + | 101 | 60.8 | (2Fe-2S)-binding protein |
| +6 | WP_107028685.1 | + | 495 | 66.5 | FAD/NAD(P)-dependent oxidoreductase |
| +7 | WP_107028686.1 | + | 388 | 63.2 | NAD(P)/FAD-dependent oxidoreductase |
| +8 | WP_107028687.1 | + | 237 | 62.5 | ABC transporter permease (start of an ABC permease/ATPase operon, +8..+10) |

Whole contig: 51 CDS in 60.5 kb (one gene per 1.19 kb, coding-dense, no introns anywhere, all single-part CDS); contents are typical rhizobial accessory-genome material: two complete ABC transporter operons (branched-chain amino acid type at 25–31 kb; +8..+10), two-component system, LysR/TetR/HTH regulators, SDR and flavin oxidoreductases, Dyp peroxidase, a tRNA-Ile, and three IS elements (IS110 and IS5 pseudo-transposases at 15 kb and 22 kb, intact IS256 transposase at the contig end 59.6–60.5 kb). Genes +5..+7 and +8..+10 overlap by 1–4 bp (translational coupling) — a bacterial operon signature. Nothing on the contig looks eukaryotic (no multi-exon models, no eukaryote-specific products), and there are no plasmid replication/partition (repABC) or conjugation genes, so the contig is most likely chromosomal or megaplasmid accessory DNA, not a Ti/Ri-type plasmid.

## 3. GC content and read coverage (results/gc_summary.tsv, results/contig_coverage.tsv)

| item | GC% |
|---|---|
| gene (1,086 bp) | 58.5 (GC3 = 77.6 %, strongly bacterial-high-GC codon usage) |
| ±20 kb window (25,794..60,479; truncated by contig end) | 59.0 |
| whole contig NZ_PVZR01000036.1 | 58.5 |
| whole assembly GCF_003024455.1 | 61.0 (median of 49 contigs ≥10 kb: 60.8) |

SPAdes k-mer coverage from contig headers: this contig 22.7; median of contigs ≥10 kb 21.7, length-weighted mean 20.1, range 15.5–30.5. The contig is at exactly the genome's modal coverage — a contaminating organism would almost always show a different depth. GC of the contig (58.5 %) is 2.5 points below the genome mean, typical for a rhizobial accessory/island region (it carries IS elements), and far from the Boletales sus-like carriers (contigs ~44–50 %) or Mortierella (51.2 % over the 200-kb window of contig JANATF010000012 held locally).

## 4. Homologues (project tables + fresh remote BLASTP vs nr)

From the project's own searches (analysis/refset/neighbourhood/seeds_nr500.tsv, analysis/donor_meta/results/queries_v2_self.tsv, globdb_v2.tsv):

- Closest sequence of any kind: the *Mortierella alpina* sus-like pseudogene SUS_JANATF010000012_L1 (GCA_029449245.1, contig JANATF010000012 of 1.69 Mb, fungal context, 1 exon, 345 aa) — **60.9 % identity over 345 aa**. In the published sus tree Ensifer NM-2 and Mortierella are sisters (both leaves 0.237 subst/site, joint stem 0.858), inside the group "Boletales + Mortierella + Ensifer" (92.7/94).
- Boletales sus-like genes: 38–46 % (e.g. Porphyrellus sp. 42.9 %, Boletellus deceptivus 43.8 %).
- Closest bacterial proteins in nr: *Rhizobium rhizogenes* Ri-plasmid Sus-like WP_174158420.1 41.7 %, Rhizobium/Agrobacterium WP_034521016.1 43.0 %, *Sinorhizobium meliloti* hypothetical proteins 37 %; plant cT-DNA Sus (Ipomoea, Diospyros, Euphorbia) 40–46 %; T-DNA succinamopine synthase ARU12438 36.7 %. GlobDB isolates: best 43.9 % (GCF_027460065, a rhizobium). No bacterial sequence above 44 %.
- No fungal or bacterial sequence is near-identical (>90 %) to the Ensifer gene; it is a singleton divergent sus-like gene.

Fresh remote BLASTP vs nr (NCBI URL API, RID in `ncbi/blast_rid.txt`, submitted 2026-10-13 06:29) was still queued after 25 min and was not waited for; the project's own nr search (seeds_nr500.tsv, 73 hits, same query) already gives the nr picture above. To collect it later:
`curl -s "https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=$(cut -d= -f2 ncbi/blast_rid.txt)&FORMAT_TYPE=Tabular" -o ncbi/blastp_nr_WP_107028680.txt` (RIDs expire after ~36 h).

## 5. Verdict

The gene WP_107028680.1 is a genuine bacterial gene of the *Ensifer* sp. NM-2 genome, not fungal contamination: it sits in the middle of a 60-kb coding-dense contig (51 intron-less CDS, overlapping operons, tRNA, IS elements, two-component system, ABC transporters) whose every other gene is an ordinary rhizobial gene, the contig has the same SPAdes coverage as the rest of the assembly (22.7 vs median 21.7), its GC (58.5 %, GC3 77.6 %) matches an Ensifer accessory region rather than any fungal carrier (44–50 %), and no fungal sequence in the project or in nr is more than 61 % identical to it — a contaminant would be near-identical to its source. The reverse scenario (the fungal genes being Ensifer contamination) is equally excluded for the same reason (≤61 % identity, and the fungal loci are on megabase fungal contigs). What the data do show is that the Ensifer gene is a divergent sus-like OCD whose nearest relative is the *Mortierella alpina* soil-fungus pseudogene (61 %), not the Ri/Ti-plasmid Sus proteins (≤44 %); both organisms are soil-derived. Uncertainty: the assembly is contig-level with no reads examined here, the region carries IS elements so it may be a mobile/accessory island rather than core chromosome, and with a single copy in one strain (IPG: one genome) we cannot say whether this sus-like gene is widespread in Ensifer; the Ensifer–Mortierella pairing in the tree rests on a long shared stem (0.86) and should be described as "the closest known bacterial relative of the fungal sus-like group" rather than as its source.
