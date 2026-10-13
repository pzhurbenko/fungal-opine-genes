import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
#!/usr/bin/env python3
import csv, os, glob
from Bio import SeqIO
D=f'{GMO_ROOT}'; LT=f'{D}/analysis/loci_pipeline/round2/loci_table.tsv'
SP={ 'Boletus_edulis':['GCA_038093215.1','GCA_003316165.1','GCA_054741165.1','GCA_015179015.1','GCA_014925795.1','GCA_022837065.1','GCA_038089135.1'],
 'Butyriboletus_roseoflavus':['GCA_019364595.1','GCA_003315995.1'],
 'Retiboletus_ornatipes':['GCA_003316065.1','GCA_038091215.1'],
 'Tylopilus_felleus':['GCA_038091295.1','GCA_054948845.1','GCA_038093275.1'],
 'Boletus_reticuloceps':['GCA_018397855.1','GCA_038093655.1'],
 'Austroboletus_festivus':['GCA_038091375.1','GCA_038088875.1'],
 'Boletellus_ananas':['GCA_038089235.1','GCA_038089895.1'],
 'Phylloporus_rhodoxanthus':['GCA_038091075.1','GCA_038091195.1'],
 'Boletus_variipes':['GCA_038088715.1','GCA_038088695.1','GCA_038092395.1'] }
def fna(g):
    for base in ['step2','extra']:
        f=glob.glob(f'{D}/data/fungal_genomes/{base}/ncbi_dataset/data/{g}/*.fna')
        if f: return f[0]
    raise SystemExit(f'no fna for {g}')
loci=[r for r in csv.DictReader(open(LT),delimiter='\t') if r['status'] not in ('not_homolog','no_model','')]
os.makedirs('contigs',exist_ok=True)
pairs=open('pairs.tsv','w'); genomes=open('genomes.tsv','w'); srcloci=open('source_loci.tsv','w')
srcloci.write('species\tgca\tlocus_id\tfam\tsubgroup\tstatus\tcontig\tcontig_len\tcontig_start\tcontig_end\n')
seen=set()
for sp,gcas in SP.items():
    for g in gcas:
        if g not in seen: genomes.write(f'{g}\t{sp}\t{fna(g)}\n'); seen.add(g)
    for src in gcas:
        L=[r for r in loci if r['gca']==src]
        if not L: continue
        recs=[]; done=set()
        for r in L:
            srcloci.write('\t'.join([sp,src,r['locus_id'],r['fam'],r['subgroup'],r['status'],r['contig'],r['contig_len'],r['contig_start'],r['contig_end']])+'\n')
            if r['contig'] in done: continue
            done.add(r['contig']); f=f"{D}/analysis/loci_pipeline/contigs/{r['contig']}.fa"
            if os.path.exists(f): recs+=list(SeqIO.parse(f,'fasta'))
            else: print('missing contig',r['contig'])
        SeqIO.write(recs,f'contigs/{sp}__{src}.fa','fasta')
        for tgt in gcas:
            if tgt!=src: pairs.write(f'{sp}\t{src}\t{tgt}\t{fna(tgt)}\n')
pairs.close(); genomes.close(); srcloci.close()
# all Boletales models of the 4 families, header = locus|subgroup|status
lt={r['locus_id']:r for r in csv.DictReader(open(LT),delimiter='\t')}
out=[]
for f in glob.glob(f'{D}/analysis/loci_pipeline/round2/proteins/*_models.faa'):
    for rec in SeqIO.parse(f,'fasta'):
        lid=rec.id.split('|')[0]; r=lt.get(lid)
        if r and r['order']=='Boletales' and r['status'] in ('intact','pseudogene_like'):
            rec.id=f"{lid}|{r['subgroup'] or 'NA'}|{r['status']}"; rec.description=''; out.append(rec)
SeqIO.write(out,'boletales_models.faa','fasta'); print('models',len(out),'pairs written')
