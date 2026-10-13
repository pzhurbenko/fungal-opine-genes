import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
#!/usr/bin/env python3
import csv, json, glob, subprocess, os
from collections import defaultdict
A=f'{GMO_ROOT}/analysis'; D=f'{GMO_ROOT}/data/fungal_genomes'
cds={r['protein_id']:r for r in csv.DictReader(open(f'{A}/refset/plasmids/cds_table.tsv'),delimiter='\t')}
meta=list(csv.DictReader(open('island_genes_meta.tsv'),delimiter='\t'))
best={}
for q,s,pid,ln,cov,ev,bs in csv.reader(open('island_vs_plasmid.tsv'),delimiter='\t'):
    if q not in best or float(bs)>best[q][2]: best[q]=(s,float(pid),float(bs),float(cov),float(ev))
out=open('island_gene_order.tsv','w'); out.write('island_contig\tgca\tgene\tstart\tend\tstrand\texons\tkind\tbest_plasmid_protein\tplasmid\tplasmid_gene\tproduct\tin_TDNA\tplasmid_start\tplasmid_end\tpident\tqcov\tevalue\n')
by=defaultdict(list)
for m in meta:
    q=f"{m['contig']}|{m['gene']}"; b=best.get(q)
    row=dict(m)
    if b and b[4]<1e-5:
        pid=b[0].split('|')[0]; c=cds.get(pid,{})
        row.update({'best_plasmid_protein':pid,'plasmid':c.get('plasmid',b[0].split('|')[1]),'plasmid_gene':c.get('gene',''),'product':c.get('product',''),'in_TDNA':c.get('in_TDNA',''),'plasmid_start':c.get('start',''),'plasmid_end':c.get('end',''),'pident':b[1],'qcov':b[3],'evalue':b[4]})
    else: row.update({k:'' for k in ['best_plasmid_protein','plasmid','plasmid_gene','product','in_TDNA','plasmid_start','plasmid_end','pident','qcov','evalue']})
    by[m['contig']].append(row)
for ctg,rows in by.items():
    rows.sort(key=lambda r:int(r['start']))
    for r in rows: out.write('\t'.join(str(r.get(k,'')) for k in ['island','gca','gene','start','end','strand','exons','kind','best_plasmid_protein','plasmid','plasmid_gene','product','in_TDNA','plasmid_start','plasmid_end','pident','qcov','evalue'])+'\n')
out.close()
# collinearity: for each island, genes with T-DNA hits on the same plasmid — compare order
print('=== Island gene order (only genes with plasmid hits E<1e-5) ===')
for ctg,rows in by.items():
    hits=[r for r in rows if r['best_plasmid_protein']]
    print(f"\n## {ctg} ({rows[0]['gca']}): {len(rows)} genes, {len(hits)} with plasmid hits")
    for r in hits: print(f"  {r['gene']:22s} {r['start']:>9s}-{r['end']:<9s} {r['strand']} {r['kind'][:28]:28s} -> {r['best_plasmid_protein']:14s} {r['plasmid']:14s} {str(r['plasmid_gene'])[:14]:14s} {str(r['product'])[:38]:38s} TDNA={r['in_TDNA']:3s} pos={r['plasmid_start']}-{r['plasmid_end']} id={r['pident']} cov={r['qcov']}")
# border search
regions=json.load(open('regions.json'))
borders=set(r['seq'].upper() for r in csv.DictReader(open(f'{A}/refset/plasmids/borders.tsv'),delimiter='\t'))
cons='TGGCAGGATATATTGTGGTGTAAAC'  # canonical 25-bp RB consensus (pTi)
borders.add(cons)
def rc(s): return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def hd(a,b): return sum(1 for x,y in zip(a,b) if x!=y)
bout=open('border_like_hits.tsv','w'); bout.write('contig\tregion\tpos\tstrand\tmismatches\tseq\tbest_border\n')
print('\n=== T-DNA border-like 25-mers (<=5 mismatches to any catalogued border / consensus) in island regions ±3 kb ===')
for ctg,(gca,s,e) in regions.items():
    fna=glob.glob(f'{D}/step2/ncbi_dataset/data/{gca}/*.fna')+glob.glob(f'{D}/extra/ncbi_dataset/data/{gca}/*.fna')
    if not fna: print('no fna',gca); continue
    cid=subprocess.run(f"grep -m1 -oP '^>\\K{ctg}\\S*' {fna[0]}",shell=True,capture_output=True,text=True).stdout.strip() or ctg
    seq=subprocess.run(f'samtools faidx {fna[0]} {cid}:{s}-{e}',shell=True,capture_output=True,text=True).stdout
    seq=''.join(l.strip() for l in seq.split('\n')[1:]).upper()
    n=0; bestm=99
    for i in range(len(seq)-25):
        w=seq[i:i+25]
        for strand,ww in (('+',w),('-',rc(w))):
            m=min((hd(ww,b),b) for b in borders)
            bestm=min(bestm,m[0])
            if m[0]<=5:
                n+=1; bout.write('\t'.join(map(str,[ctg,f'{s}-{e}',s+i,strand,m[0],ww,m[1]]))+'\n')
    print(f"  {ctg}: region {s}-{e} ({len(seq)} bp): hits<=5mm: {n}; best mismatch count: {bestm}")
bout.close()
