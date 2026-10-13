import sys,os
clean,wd=sys.argv[1:3]
odir=os.path.join(wd,'of_input'); os.makedirs(odir,exist_ok=True)
recs=[];h=None;buf=[]
for line in open(clean):
    line=line.rstrip('\n')
    if line.startswith('>'):
        if h: recs.append((h,''.join(buf)))
        h=line[1:];buf=[]
    else: buf.append(line)
if h: recs.append((h,''.join(buf)))
m=open(os.path.join(wd,'of_species_map.tsv'),'w'); m.write('species_file\tgene_id\n')
for n,(h,s) in enumerate(recs):
    sp=f'sp{n:04d}'
    open(os.path.join(odir,f'{sp}.fa'),'w').write(f'>{h}\n{s}\n')
    m.write(f'{sp}\t{h}\n')
m.close()
print(f'  of_input: {len(recs)} single-seq files')
