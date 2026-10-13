import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import re,sys,collections
B=f'{GMO_ROOT}/data/fungal_genomes/extra2/x/ncbi_dataset/data'
J=[('GCA_058411305.1','KAO1061018.1','GCA_058411305.1_ASM5841130v1_genomic.fna'),('GCA_964275055.1','CAL9735043.1','GCA_964275055.1_MOSE0_genomic.fna')]
def fa(p):
    d={};n=None
    for l in open(p):
        if l[0]=='>': n=l[1:].split()[0]; d[n]=[]
        else: d[n].append(l.strip().upper())
    return {k:''.join(v) for k,v in d.items()}
def gc(s):
    a=sum(s.count(c) for c in 'ACGT'); return 100*(s.count('G')+s.count('C'))/a if a else float('nan')
for g,prot,fn in J:
    seq=fa(f'{B}/{g}/{fn}'); asm_gc=gc(''.join(seq.values()))
    mrna={};cds=collections.defaultdict(list);genes={}
    for l in open(f'{B}/{g}/genomic.gff'):
        if l[0]=='#': continue
        f=l.rstrip('\n').split('\t')
        if len(f)<9: continue
        at=dict(x.split('=',1) for x in f[8].split(';') if '=' in x)
        if f[2]=='CDS': cds[at['Parent']].append((int(f[3]),int(f[4]),f[6],at.get('protein_id'),at.get('product','')))
        elif f[2]=='gene': genes[at['ID']]=(f[0],int(f[3]),int(f[4]),f[6],at.get('locus_tag',''))
        elif f[2] in('mRNA',): mrna[at['ID']]=at['Parent']
    # transcript -> gene ; protein-coding genes with longest CDS set
    gl={}
    for tid,cl in cds.items():
        gid=mrna.get(tid,tid)
        if gid not in genes: # try parent of rna
            continue
        L=sum(b-a+1 for a,b,*_ in cl)
        if gid not in gl or L>gl[gid][0]: gl[gid]=(L,len(cl),cl[0][4],cl[0][3],tid)
    # fallback: RNA parents (non-mRNA types)
    print('===',g,prot,'genes',len(genes),'coding',len(gl),'asm_gc %.2f'%asm_gc)
    tgt=[gid for gid,v in gl.items() if v[3]==prot]
    if not tgt: print('NOT FOUND'); continue
    t=tgt[0]; c,s,e,st,lt=genes[t]; L,ne,prod,_,tid=gl[t]
    clen=len(seq[c]); print(f'target {lt} {c} len={clen} {s}-{e} {st} exons={ne} cds_nt={L} aa={L//3-1 if L%3==0 else "?"} prod={prod}')
    print(f'dist_left={s-1} dist_right={clen-e}')
    reg=seq[c][max(0,s-1-5000):e+5000]; print('GC gene+-5kb %.2f  gene only %.2f  window %d bp'%(gc(reg),gc(seq[c][s-1:e]),len(reg)))
    cg=sorted([(genes[x][1],genes[x][2],x) for x in gl if genes[x][0]==c])
    i=[k for k,v in enumerate(cg) if v[2]==t][0]
    for k in range(max(0,i-10),min(len(cg),i+11)):
        a,b,x=cg[k]; v=gl[x]; print(f'{k-i:+d}\t{genes[x][4]}\t{a}-{b}\t{genes[x][3]}\texons={v[1]}\taa={v[0]//3-1}\t{v[2]}\t{v[3]}')
    print(f'n_left_available={i} n_right_available={len(cg)-i-1}')
