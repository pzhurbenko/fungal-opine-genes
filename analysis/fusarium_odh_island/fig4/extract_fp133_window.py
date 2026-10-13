"""Fp133 (GCA_046056315.1) genes in windows around the odh island copies on chr2 (CM100255.1) and chr3 (CM100256.1) -> fp133_window_genes.tsv, fp133_window.faa"""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import re
from Bio import SeqIO
from Bio.Seq import Seq
D=f'{GMO_ROOT}/data/fungal_genomes/fusarium/fus/ncbi_dataset/data/GCA_046056315.1/'
W={'CM100255.1':(9415000,9515000),'CM100256.1':(8370000,8470000)}
genome={r.id:r.seq for r in SeqIO.parse(D+'GCA_046056315.1_ASM4605631v1_genomic.fna','fasta') if r.id in W}
cds={}
for l in open(D+'genomic.gff'):
    if l.startswith('#'): continue
    p=l.rstrip('\n').split('\t')
    if len(p)<9 or p[2]!='CDS' or p[0] not in W: continue
    a=int(p[3]);b=int(p[4])
    if b<W[p[0]][0] or a>W[p[0]][1]: continue
    at=dict(x.split('=',1) for x in p[8].split(';'))
    k=at['locus_tag']; cds.setdefault(k,dict(c=p[0],s=p[6],ex=[],prot=at.get('protein_id'),prod=at.get('product','')))['ex'].append((a,b,int(p[7])))
out=open('fp133_window_genes.tsv','w'); out.write('contig\tlocus_tag\tprotein\tstart\tend\tstrand\texons\tproduct\n'); fa=open('fp133_window.faa','w')
for k,v in sorted(cds.items(),key=lambda kv:min(e[0] for e in kv[1]['ex'])):
    ex=sorted(v['ex']); s=''.join(str(genome[v['c']][a-1:b]) for a,b,_ in ex)
    sq=Seq(s); 
    if v['s']=='-': sq=sq.reverse_complement()
    pr=str(sq.translate()).rstrip('*')
    out.write(f"{v['c']}\t{k}\t{v['prot']}\t{ex[0][0]}\t{ex[-1][1]}\t{v['s']}\t{len(ex)}\t{v['prod']}\n"); fa.write(f'>{k}\n{pr}\n')
