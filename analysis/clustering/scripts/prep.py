import sys,re
inf,outf,mapf=sys.argv[1:4]
def reader(p):
    h=None;s=[]
    for line in open(p):
        line=line.rstrip('\n')
        if line.startswith('>'):
            if h is not None: yield h,''.join(s)
            h=line[1:];s=[]
        else: s.append(line.strip())
    if h is not None: yield h,''.join(s)
seen={}
n=0
with open(outf,'w') as o, open(mapf,'w') as m:
    m.write('clean_id\toriginal_header\n')
    for h,seq in reader(inf):
        cid=re.sub(r'_+','_',re.sub(r'[^A-Za-z0-9._-]+','_',h)).strip('_')
        if cid in seen: seen[cid]+=1; cid=f'{cid}__{seen[cid]}'
        else: seen[cid]=0
        o.write(f'>{cid}\n{seq}\n'); m.write(f'{cid}\t{h}\n'); n+=1
print(f'  prep: {n} seqs -> {outf}')
