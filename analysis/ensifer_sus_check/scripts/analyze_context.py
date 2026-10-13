#!/usr/bin/env python3
"""Context check for Ensifer sp. NM-2 WP_107028680.1 (sus-like OCD) on contig NZ_PVZR01000036.1.
Usage: python3 -I analyze_context.py <contig.gb> <genome.fna.gz> <out_dir>
Outputs: neighbours.tsv, gc_summary.tsv, contig_coverage.tsv (from SPAdes headers).
"""
import sys, gzip, re, os, statistics
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

gb, fna, out = sys.argv[1:4]
os.makedirs(out, exist_ok=True)
TARGET = "WP_107028680.1"

rec = SeqIO.read(gb, "genbank")
L = len(rec.seq)
cds = [f for f in rec.features if f.type == "CDS"]
cds.sort(key=lambda f: int(f.location.start))
idx = next(i for i, f in enumerate(cds) if f.qualifiers.get("protein_id", [""])[0] == TARGET)
tgt = cds[idx]
ts, te = int(tgt.location.start), int(tgt.location.end)

def has_intron(f):
    return len(f.location.parts) > 1

with open(f"{out}/neighbours.tsv", "w") as fh:
    fh.write("rank\tlocus_tag\tprotein_id\tstart\tend\tstrand\tlen_aa\tparts\tgc_pct\tproduct\n")
    lo, hi = max(0, idx - 8), min(len(cds), idx + 9)
    for i in range(lo, hi):
        f = cds[i]
        s, e = int(f.location.start), int(f.location.end)
        prod = f.qualifiers.get("product", ["?"])[0]
        pid = f.qualifiers.get("protein_id", ["-"])[0]
        lt = f.qualifiers.get("locus_tag", ["-"])[0]
        aa = len(f.qualifiers.get("translation", [""])[0])
        gc = 100 * gc_fraction(rec.seq[s:e])
        tag = "TARGET" if i == idx else str(i - idx)
        fh.write(f"{tag}\t{lt}\t{pid}\t{s+1}\t{e}\t{f.location.strand}\t{aa}\t{len(f.location.parts)}\t{gc:.1f}\t{prod}\n")

w0, w1 = max(0, ts - 20000), min(L, te + 20000)
gc_gene = 100 * gc_fraction(rec.seq[ts:te])
gc_win = 100 * gc_fraction(rec.seq[w0:w1])
gc_contig = 100 * gc_fraction(rec.seq)
# GC3 of the gene (codon third positions)
g = tgt.extract(rec.seq)
third = str(g)[2::3]
gc3 = 100 * sum(c in "GC" for c in third) / len(third)

# whole genome and per-contig coverage from SPAdes headers
tot_gc = tot_len = 0
cov = []
with gzip.open(fna, "rt") as fh:
    for r in SeqIO.parse(fh, "fasta"):
        tot_gc += sum(r.seq.count(b) for b in "GCgc"); tot_len += len(r.seq)
        m = re.search(r"NODE_(\d+)_length_(\d+)_cov_([\d.]+)", r.description)
        if m:
            cov.append((r.id, int(m.group(2)), float(m.group(3)), 100 * gc_fraction(r.seq)))
gc_genome = 100 * tot_gc / tot_len
cov.sort(key=lambda x: -x[1])
with open(f"{out}/contig_coverage.tsv", "w") as fh:
    fh.write("contig\tlength\tspades_cov\tgc_pct\n")
    for c in cov:
        fh.write(f"{c[0]}\t{c[1]}\t{c[2]:.2f}\t{c[3]:.1f}\n")
big = [c for c in cov if c[1] >= 10000]
covs = [c[2] for c in big]
lenw = sum(c[1] for c in big)
cov_wmean = sum(c[1] * c[2] for c in big) / lenw
this = next(c for c in cov if c[0] == rec.id)
gcs_big = [c[3] for c in big]

with open(f"{out}/gc_summary.tsv", "w") as fh:
    fh.write("item\tvalue\n")
    fh.write(f"contig\t{rec.id}\ncontig_length_bp\t{L}\ncds_on_contig\t{len(cds)}\n")
    fh.write(f"target_coords\t{ts+1}..{te} strand {tgt.location.strand}\ntarget_len_bp\t{te-ts}\ntarget_parts\t{len(tgt.location.parts)}\n")
    fh.write(f"gc_gene_pct\t{gc_gene:.1f}\ngc3_gene_pct\t{gc3:.1f}\ngc_window_20kb_pct\t{gc_win:.1f}\ngc_contig_pct\t{gc_contig:.1f}\ngc_genome_pct\t{gc_genome:.1f}\n")
    fh.write(f"window_coords\t{w0+1}..{w1}\n")
    fh.write(f"contig_spades_cov\t{this[2]:.2f}\ncontigs_ge10kb\t{len(big)}\ncov_median_ge10kb\t{statistics.median(covs):.2f}\ncov_lenweighted_mean_ge10kb\t{cov_wmean:.2f}\n")
    fh.write(f"cov_min_max_ge10kb\t{min(covs):.2f}..{max(covs):.2f}\ngc_median_contigs_ge10kb\t{statistics.median(gcs_big):.1f}\n")
    fh.write(f"genome_total_bp\t{tot_len}\n")
print(open(f"{out}/gc_summary.tsv").read())
print(open(f"{out}/neighbours.tsv").read())
