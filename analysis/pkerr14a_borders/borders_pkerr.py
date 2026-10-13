from Bio import SeqIO
rec=SeqIO.read('LT009732.gb','genbank'); seq=str(rec.seq).upper(); L=len(seq)
BORDER='TGGCAGGATATATTGTGGTGTAAAC'
def rc(s): return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def find(seq,minmatch):
    hits=[]
    for strand,sq in (('+',seq),('-',rc(seq))):
        for st in range(len(sq)-25):
            w=sq[st:st+25]
            if w[3]!='C' and w[4]!='A': continue
            sc=sum(1 for a,b in zip(w,BORDER) if a==b)
            if sc>=minmatch:
                pos=st if strand=='+' else len(sq)-(st+25); hits.append((pos,strand,sc,w))
    hits.sort(); out=[]
    for h in hits:
        if out and abs(h[0]-out[-1][0])<10:
            if h[2]>out[-1][2]: out[-1]=h
        else: out.append(h)
    return out
cds=[]
for f in rec.features:
    if f.type=='CDS':
        pid=f.qualifiers.get('protein_id',[''])[0]; prod=f.qualifiers.get('product',[''])[0]; gene=f.qualifiers.get('gene',[''])[0]
        cds.append((int(f.location.start),int(f.location.end),f.location.strand,pid,gene,prod))
acs=[c for c in cds if c[3]=='CUX65485.1'][0]
print('length',L,'CDS',len(cds)); print('acs',acs[:3])
for mm in (18,17,16):
    h=find(seq,mm); print(f'\nborders >= {mm}/25: {len(h)}')
    for pos,strand,sc,w in h: print(f'  pos={pos:7} {strand} {sc}/25 {w} dist_to_acs={min(abs(pos-acs[0]),abs(pos-acs[1]))}')
print('\nvir/tra/rep/opine-related CDS:')
import re
for c in cds:
    if re.search(r'vir|tra[A-Z]|rep[A-C]|opine|nopaline|octopine|agrocinopine|mannopine|agropine|cucumopine|iaa|ipt|tryptophan|indole|cytokinin|T-DNA|6b|rol|border',c[5]+c[4],re.I):
        print(f'  {c[0]:7}-{c[1]:7} {"+" if c[2]==1 else "-"} {c[3]:12} {c[4]:6} {c[5]}')
print('\nneighbours of acs (±15 kb):')
for c in cds:
    if c[1]>acs[0]-15000 and c[0]<acs[1]+15000: print(f'  {c[0]:7}-{c[1]:7} {"+" if c[2]==1 else "-"} {c[3]:12} {c[4]:6} {c[5]}')
