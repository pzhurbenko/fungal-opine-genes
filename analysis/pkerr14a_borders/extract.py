from Bio import SeqIO
rec=SeqIO.read('LT009732.gb','genbank'); seq=str(rec.seq).upper()
with open('pkerr_60_130.faa','w') as o:
    for f in rec.features:
        if f.type=='CDS' and 60000<int(f.location.start)<130000:
            pid=f.qualifiers.get('protein_id',['?'])[0]; tr=f.qualifiers.get('translation',[''])[0]
            if tr: o.write(f'>{pid}|{int(f.location.start)}-{int(f.location.end)}|{f.location.strand}\n{tr}\n')
import re
for core in ('CAGGATATAT','GGATATATTG'):
    for strand,sq in (('+',seq),('-',str(rec.seq.reverse_complement()).upper())):
        for m in re.finditer(core,sq):
            pos=m.start() if strand=='+' else len(sq)-m.end()
            w=sq[max(0,m.start()-3):m.start()+22]
            print('core',core,strand,pos,w)
