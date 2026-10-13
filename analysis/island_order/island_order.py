#!/usr/bin/env python3
"""B1: gene order of Boletales islands vs T-DNA plasmids + T-DNA border-like 25-mers in island flanks.
Islands: bounded islands with >=1 bacterial-like gene or registry loci from key Boletales genomes."""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import csv, os, glob, subprocess, json, re
from collections import defaultdict
A=f'{GMO_ROOT}/analysis'; D=f'{GMO_ROOT}/data/fungal_genomes'
lt={r['locus_id']:r for r in csv.DictReader(open(f'{A}/loci_pipeline/round2/loci_table.tsv'),delimiter='\t')}
ISL=[('GCA_019364595.1','JAHRGV010000066'),('GCA_019364595.1','JAHRGV010000067'),('GCA_019364595.1','JAHRGV010000044'),('GCA_019364595.1','JAHRGV010000002'),
     ('GCA_038088715.1','JBBPJK010000045'),('GCA_038088715.1','JBBPJK010000239'),('GCA_041722695.1','JBBPWX010000022'),('GCA_041722695.1','JBBPWX010000076'),
     ('GCA_038089435.1','JBBPKX010000014'),('GCA_038089435.1','JBBPKX010000046'),('GCA_019096195.1','JAGTXK010000187'),('GCA_019096195.1','JAGTXK010000090'),
     ('GCA_038090635.1','JBBPNB010000203'),('GCA_021018615.1','JAEVFL010000001'),('GCA_003315995.1','QFEE01000702'),('GCA_038093215.1','JBBPRY010001778'),('GCA_038089235.1','JBBPKL010003317')]
isl=list(csv.DictReader(open(f'{A}/islands/islands.tsv'),delimiter='\t'))
# collect island genes: neighbours (any role) within [island_start-2kb, island_end+2kb] + registry targets
q=open('island_genes.faa','w'); meta=open('island_genes_meta.tsv','w'); meta.write('island\tgca\tcontig\tgene\tstart\tend\tstrand\texons\tkind\tlocus\n')
regions={}
for gca,ctg in ISL:
    rows=[r for r in isl if r['contig']==ctg]
    loci=[l for l in lt.values() if l['contig']==ctg and l['status'] not in ('not_homolog','no_model','')]
    if rows: s=min(int(r['island_start']) for r in rows); e=max(int(r['island_end']) for r in rows)
    else: s=min(int(l['contig_start']) for l in loci); e=max(int(l['contig_end']) for l in loci)
    s=max(1,s-3000); e=e+3000; regions[ctg]=(gca,s,e)
    seen=set()
    for l in loci:
        f=f"{A}/loci_pipeline/round2/loci/{l['locus_id']}/neighbors.faa"
        if not os.path.exists(f): continue
        cur=None
        for line in open(f):
            if line.startswith('>'):
                p=line[1:].strip().split('|'); c=p[2]; gs,ge=c[:-1].split('-'); gs,ge=int(gs),int(ge)
                cur=None
                if ge>=s and gs<=e and p[1] not in seen:
                    seen.add(p[1]); cur=f'{ctg}|{p[1]}'; q.write(f'>{cur}\n'); meta.write('\t'.join(map(str,[ctg,gca,ctg,p[1],gs,ge,c[-1],p[3].split("=")[1],'neighbor','']))+'\n')
            elif cur: q.write(line)
        tp=f"{A}/loci_pipeline/round2/loci/{l['locus_id']}/target_protein.faa"
        if os.path.exists(tp):
            seq=''.join(x.strip() for x in open(tp) if not x.startswith('>'))
            q.write(f">{ctg}|T_{l['locus_id']}\n{seq}\n"); meta.write('\t'.join(map(str,[ctg,gca,ctg,'T_'+l['locus_id'],l['contig_start'],l['contig_end'],l['strand'],l['n_exons'],'registry:'+l['fam']+':'+l['subgroup']+':'+l['status'],l['locus_id']]))+'\n')
q.close(); meta.close()
json.dump(regions,open('regions.json','w'))
print('queries written')
