#!/usr/bin/env python3
"""Per registry locus: contig depth vs genome, min depth inside locus, junction-spanning reads at locus start/end,
proper-pair fraction on contig vs genome, GC of contig vs assembly. Uses samtools (env annot)."""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import csv, subprocess, glob, os, statistics
from collections import defaultdict
D=f'{GMO_ROOT}'; LT=f'{D}/analysis/loci_pipeline/round2/loci_table.tsv'
loci=[r for r in csv.DictReader(open(LT),delimiter='\t') if r['status'] not in ('not_homolog','no_model','')]
targets={r[0]:r for r in csv.reader(open('links.tsv'),delimiter='\t')}
def sh(c): return subprocess.run(c,shell=True,capture_output=True,text=True).stdout
def gc(seq): s=seq.upper(); return (s.count('G')+s.count('C'))/max(1,len(s)-s.count('N'))
out=open('locus_coverage.tsv','w')
out.write('gca\tspecies\tlocus\tfam\tsubgroup\tstatus\tcontig\tcontig_len\tgenome_median_depth\tcontig_mean_depth\tcontig_ratio\tlocus_mean_depth\tlocus_min_depth\tlocus_ratio\tspan_left\tspan_right\tpp_frac_contig\tpp_frac_genome\tgc_contig\tgc_assembly\tmapq_contig\n')
for gca,row in targets.items():
    bam=f'bam/{gca}.bam'
    if not os.path.exists(f'depth/{gca}.coverage.tsv'): continue
    cov={r['#rname']:r for r in csv.DictReader(open(f'depth/{gca}.coverage.tsv'),delimiter='\t')}
    big=[float(r['meandepth']) for r in cov.values() if int(r['endpos'])>=10000]
    gmed=statistics.median(big) if big else float('nan')
    fs=open(f'bam/{gca}.flagstat').read()
    import re
    tot=int(re.search(r'(\d+) \+ \d+ in total',fs).group(1)); pp=int(re.search(r'(\d+) \+ \d+ properly paired',fs).group(1))
    ppg=pp/max(1,tot)
    # assembly GC
    fna=row[3]; seqs={}; cur=None
    for l in open(fna):
        if l.startswith('>'): cur=l[1:].split()[0]; seqs[cur]=[]
        else: seqs[cur].append(l.strip())
    seqs={k:''.join(v) for k,v in seqs.items()}
    tot_gc=gc(''.join(seqs.values()))
    for r in loci:
        if r['gca']!=gca: continue
        c=r['contig']; cid=[k for k in seqs if k.startswith(c)]
        cid=cid[0] if cid else c
        dep=f'depth/{gca}.{cid}.depth'
        if not os.path.exists(dep) or os.path.getsize(dep)==0:
            subprocess.run(f'samtools depth -a -r {cid} {bam} > {dep}',shell=True)
        d=[int(l.split('\t')[2]) for l in open(dep)]
        s,e=int(r['contig_start']),int(r['contig_end']); s,e=min(s,e),max(s,e)
        loc=d[s-1:e] if d else []
        cm=statistics.mean(d) if d else 0; lm=statistics.mean(loc) if loc else 0; lmin=min(loc) if loc else 0
        def span(pos):
            if pos<150 or pos>len(d)-150: return 'edge'
            return sh(f"samtools view -c -e 'pos<{pos-40} && endpos>{pos+40}' {bam} {cid}").strip()
        sl,sr=span(s),span(e)
        n_all=int(sh(f'samtools view -c {bam} {cid}').strip() or 0); n_pp=int(sh(f'samtools view -c -f 2 {bam} {cid}').strip() or 0)
        mq=sh(f"samtools view {bam} {cid} | awk '{{s+=$5;n++}} END{{if(n) printf \"%.1f\", s/n}}'")
        out.write('\t'.join(map(str,[gca,row[1],r['locus_id'],r['fam'],r['subgroup'],r['status'],cid,len(d),f'{gmed:.1f}',f'{cm:.1f}',f'{cm/gmed:.2f}' if gmed else 'nan',f'{lm:.1f}',lmin,f'{lm/gmed:.2f}' if gmed else 'nan',sl,sr,f'{n_pp/max(1,n_all):.2f}',f'{ppg:.2f}',f'{gc(seqs.get(cid,"")):.3f}',f'{tot_gc:.3f}',mq]))+'\n')
out.close(); print('locus_coverage.tsv written')
