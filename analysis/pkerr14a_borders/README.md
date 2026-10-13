# pKerr14a (LT009732.1): is *acs* inside a T-region? (2026-10-13)

Question: the pKerr14a-type lineage groups with the *acs* of the Ti-type plasmid pKerr14a (CUX65485.1) and a cT-DNA fragment of *Parasponia*; pKerr14a is not in the 42-plasmid catalogue (draft assembly), so its T-region had never been delimited.

Method: `borders_pkerr.py` = the border scan of `analysis/refset/plasmids/build_plasmid_catalog.py` (25-mer TGGCAGGATATATTGTGGTGTAAAC, both strands, >= 18/25; also reported at 17 and 16) applied to the GenBank record `LT009732.gb` (NCBI efetch, 292,852 bp, 266 CDS); `extract.py` + DIAMOND blastp (env annot, --sensitive, 1e-5) of the 59 CDS in 60-130 kb against the catalogue proteins (`hits.tsv`).

Result (`borders_pkerr.out`):
- one border-like repeat at 82,583 bp (+), 20/25, core CAGGATATAT exact; *acs* 84,168-85,428 (+) starts 1.6 kb downstream;
- no second repeat >= 17/25 within 5-60 kb (the 16/25 hits are noise); T-region cannot be delimited, as for the 14 catalogue plasmids with a single border;
- neighbourhood: virD4 (67 kb) and a VirE2-family gene (70 kb) upstream; 14 kb of hypothetical proteins without catalogue homologs downstream of *acs*; then an opine-catabolic region (mas1/mas2 homologs at 33-39 %, hydantoinases, ABC transporters, 103-127 kb), traR, virC2/virC1/virD1/virG (128-142 kb). No iaa/ipt/6b homologs.

Conclusion used in the paper (Results 3.4, Discussion 4.2): "*acs* lies 1.6 kb downstream of a single border-like repeat; whether it belongs to a T-region is undetermined". The [check] about the position of *acs* on pKerr14a is closed.
