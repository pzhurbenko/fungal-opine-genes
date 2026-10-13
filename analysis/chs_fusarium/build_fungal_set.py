import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import csv,sys
D=f'{GMO_ROOT}/data/fungal_proteomes/fungi_annotated/ncbi_dataset/data'
S=f'{GMO_ROOT}/analysis/refset/screen/seeds11/'
b=lambda g:g.split('_')[1].split('.')[0]
want=[] # (label,gca,protein)
seen=set()
for r in csv.DictReader(open(S+'seeds11_candidates.tsv'),delimiter='\t'):
    if r['seed']!='chs': continue
    if b(r['gca']) in seen: continue
    seen.add(b(r['gca']))
    lab='CHSlike_'+'_'.join(r['organism'].split()[:3]).replace('.','').replace('-','')+'_'+r['gca']+'_'+r['protein']
    want.append((lab,r['gca'],r['protein']))
for r in csv.DictReader(open(S+'seeds11_best.tsv'),delimiter='\t'):
    if r['seed']=='chs' and 'Dactylonectria estremocensis' in r['organism']: want.append(('CHSlike_Dactylonectria_estremocensis_'+r['protein'],r['gca'],r['protein']))
for r in csv.DictReader(open('ags_native_selected.tsv'),delimiter='\t'):
    if r['protein']=='KAH7120428.1': continue
    want.append(('AGSnat_'+r['genus']+'_'+r['class'].replace(' ','')+'_'+r['protein'],r['gca'],r['protein']))
need={}
for l,g,p in want: need.setdefault(g,{})[p]=l
out=open('fungal_part.faa','w')
for g,m in need.items():
    h=None
    for line in open(f'{D}/{g}/protein.faa'):
        if line[0]=='>':
            a=line[1:].split()[0]; h=m.get(a)
            if h: out.write('>'+h+'\n')
        elif h: out.write(line)
n=sum(1 for l in open('fungal_part.faa') if l[0]=='>'); print(n,len(want))
out.close()
