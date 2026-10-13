import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
#!/usr/bin/env python3
# per-locus RNA-seq read counts vs genome-wide background (reads per kb on registry loci vs all AUGUSTUS genes where available)
import csv, subprocess, glob, os, statistics
A = f'{GMO_ROOT}/analysis'
lt = [r for r in csv.DictReader(open(f'{A}/loci_pipeline/round2/loci_table.tsv'), delimiter='\t') if r['status'] in ('intact', 'pseudogene_like', 'fragment')]
def sh(c): return subprocess.run(c, shell=True, capture_output=True, text=True).stdout.strip()
out = open('expression_summary.tsv', 'w')
out.write('run\tgca\tlocus\tfam\tsubgroup\tstatus\tcontig\tstart\tend\tlen_bp\treads\treads_per_kb\tgenome_reads_per_kb_median_genes\tpercentile_among_genes\n')
for r in csv.reader(open('links.tsv'), delimiter='\t'):
    srr, sp, gca = r[0], r[1], r[2]
    bam = f'bam/{srr}.{gca}.bam'
    if not os.path.exists(bam + '.bai'): continue
    fna = glob.glob(f'{GMO_ROOT}/data/fungal_genomes/step2/ncbi_dataset/data/{gca}/*.fna')[0]
    # background: AUGUSTUS genes from alien_scan if present, else 2-kb windows on long contigs
    gff = glob.glob(f'{A}/alien_scan/{gca}/*.gff3') + glob.glob(f'{A}/alien_scan/{gca}/*/*.gff3')
    bg = []
    if gff:
        genes = [l.split('\t') for f in gff for l in open(f) if '\tgene\t' in l]
        genes = genes[::max(1, len(genes)//3000)]
        for g in genes:
            n = int(sh(f'samtools view -c -q 1 {bam} {g[0]}:{g[3]}-{g[4]}') or 0); bg.append(1000*n/(int(g[4])-int(g[3])+1))
        bgname = 'augustus_genes'
    else:
        fai = [l.split('\t') for l in open(fna + '.fai')] if os.path.exists(fna + '.fai') else []
        big = [f for f in fai if int(f[1]) >= 50000][:300]
        for f in big:
            for s in range(1, int(f[1]) - 3000, max(3000, int(f[1])//10)):
                n = int(sh(f'samtools view -c -q 1 {bam} {f[0]}:{s}-{s+2999}') or 0); bg.append(1000*n/3000)
        bgname = 'windows'
    bg.sort(); med = statistics.median(bg) if bg else float('nan')
    for L in lt:
        if L['gca'] != gca: continue
        cid = sh(f"grep -m1 -oP '^>\\K{L['contig']}\\S*' {fna}") or L['contig']
        s, e = sorted((int(L['contig_start']), int(L['contig_end'])))
        n = int(sh(f'samtools view -c -q 1 {bam} {cid}:{s}-{e}') or 0)
        rpk = 1000*n/(e-s+1)
        pct = 100*sum(1 for x in bg if x < rpk)/len(bg) if bg else float('nan')
        out.write('\t'.join(map(str, [srr, gca, L['locus_id'], L['fam'], L['subgroup'], L['status'], cid, s, e, e-s+1, n, f'{rpk:.2f}', f'{med:.2f}', f'{pct:.1f}'])) + '\n')
        print(srr, L['locus_id'], L['status'], 'reads', n, f'rpk {rpk:.2f} bg-median {med:.2f} pct {pct:.0f} ({bgname})')
out.close()
