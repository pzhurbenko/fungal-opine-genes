#!/usr/bin/env python3
"""Download complete Ti/Ri plasmids (RefSeq preferred), parse CDS, detect T-DNA borders, HMM-classify all proteins.
Outputs: plasmids/<acc>.gb, plasmids/plasmids.tsv, plasmids/all_proteins.faa, plasmids/borders.tsv, plasmids/cds_table.tsv"""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import sys,re,csv,os,subprocess,glob
sys.path.insert(0,f'{GMO_ROOT}/pipeline/scripts'); from common import load_config, eutils, write_fasta, read_tsv
cfg=load_config(); os.chdir(f'{GMO_ROOT}/analysis/refset/plasmids')
cands=read_tsv('candidates.tsv')
# dedupe GenBank/RefSeq pairs: prefer NZ_/NC_ versions
by_base={}
for r in cands:
    base=re.sub(r'^NZ_','',r['acc']); pref=1 if r['acc'].startswith(('NZ_','NC_')) else 0
    if base not in by_base or pref>by_base[base][0]: by_base[base]=(pref,r)
plasmids=[r for _,r in by_base.values()]
extra=['AE007871.2','NC_002575.1','NC_010841.1','DQ058764.1','AB016260.1']
have={re.sub(r'^NZ_','',r['acc']).split('.')[0] for r in plasmids}
for a in extra:
    if a.split('.')[0] not in have: plasmids.append({'acc':a,'len':'','organism':'','title':'(canonical, added manually)'})
print('plasmids to process:',len(plasmids))
BORDER='TGGCAGGATATATTGTGGTGTAAAC'
def rc(s): return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def find_borders(seq,minmatch=18):
    """Full scan (no exact-core requirement): 25-mer with >= minmatch/25 identity to the border consensus, both strands.
    18/25 recovers TL+TR of octopine pTiAch5 (4 borders) and pRi1724 (2), while pTiC58 keeps exactly 2; pRi2659 has a divergent core."""
    hits=[]; B=BORDER
    for strand,sq in (('+',seq),('-',rc(seq))):
        L=len(sq)
        for st in range(L-25):
            w=sq[st:st+25]
            if w[3]!='C' and w[4]!='A': continue          # cheap prefilter on the 'CAGGATATAT' core start (allow 1 mismatch)
            sc=sum(1 for a,b in zip(w,B) if a==b)
            if sc>=minmatch:
                pos=st if strand=='+' else L-(st+25)
                hits.append((pos,strand,sc,w))
    # merge near-duplicate positions
    hits.sort(); out=[]
    for h in hits:
        if out and abs(h[0]-out[-1][0])<10:
            if h[2]>out[-1][2]: out[-1]=h
        else: out.append(h)
    return out
cds_rows=[]; border_rows=[]; prot=[]; plas_rows=[]
for r in plasmids:
    acc=r['acc']; gb=f'{acc}.gb'
    if not os.path.exists(gb):
        open(gb,'w').write(eutils(cfg,'efetch.fcgi',{'db':'nuccore','id':acc,'rettype':'gbwithparts','retmode':'text'}))
    txt=open(gb).read()
    if 'ORIGIN' not in txt: print('  no sequence:',acc); continue
    org=re.search(r'\n  ORGANISM  (.*)',txt); org=org.group(1).strip() if org else ''
    defn=re.search(r'DEFINITION  (.*?)\nACCESSION',txt,re.S); defn=' '.join(defn.group(1).split()) if defn else ''
    seq=''.join(re.findall(r'[acgtnACGTN]+',txt.split('ORIGIN')[-1])).upper()
    borders=find_borders(seq)
    # T-DNA intervals: consecutive borders 5-60 kb apart
    intervals=[]
    for i in range(len(borders)-1):
        a,b=borders[i][0],borders[i+1][0]
        if 5000<=b-a<=60000: intervals.append((a,b))
    virD2=[m.start() for m in re.finditer(r'VirD2|virD2',txt)]
    for pos,strand,sc,w in borders: border_rows.append({'plasmid':acc,'pos':pos,'strand':strand,'match':sc,'seq':w})
    # CDS features
    feats=re.split(r'\n     (?=CDS |gene |misc_feature|rRNA|tRNA|repeat_region|mobile_element|source)',txt.split('ORIGIN')[0])
    n=0
    for f in feats:
        if not f.startswith('CDS '): continue
        loc=f.split('\n')[0].split()[1]; nums=[int(x) for x in re.findall(r'\d+',loc)]
        if not nums: continue
        s,e=min(nums),max(nums); strand='-' if 'complement' in loc else '+'
        pid=re.search(r'/protein_id="([^"]+)"',f); gene=re.search(r'/gene="([^"]+)"',f); prod=re.search(r'/product="([^"]+)"',f)
        tr=re.search(r'/translation="([^"]+)"',f,re.S)
        if not tr: continue
        aa=re.sub(r'\s','',tr.group(1)); pid=pid.group(1) if pid else f'{acc}_{s}_{e}'
        intd=any(a<=s and e<=b for a,b in intervals)
        near=min((abs(s-p) for p,_,_,_ in borders),default=10**9)
        cds_rows.append({'plasmid':acc,'organism':org,'protein_id':pid,'gene':gene.group(1) if gene else '','product':' '.join(prod.group(1).split()) if prod else '',
                         'start':s,'end':e,'strand':strand,'aa_len':len(aa),'in_TDNA':'yes' if intd else 'no','dist_to_border':near if near<10**9 else ''})
        prot.append((f'{pid}|{acc}|{"TDNA" if intd else "nonTDNA"}',aa)); n+=1
    plas_rows.append({'plasmid':acc,'organism':org,'definition':defn[:100],'length':len(seq),'n_cds':n,'n_borders':len(borders),'n_TDNA_intervals':len(intervals),
                      'TDNA_intervals':';'.join(f'{a}-{b}' for a,b in intervals),'virD2_mentioned':'yes' if virD2 else 'no'})
    print(f'  {acc:18} {org[:32]:32} len={len(seq):7} cds={n:4} borders={len(borders):2} intervals={len(intervals)} {";".join(f"{a}-{b}" for a,b in intervals)[:60]}')
def wt(p,rows):
    with open(p,'w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0].keys()),delimiter='\t'); w.writeheader(); w.writerows(rows)
wt('plasmids.tsv',plas_rows); wt('borders.tsv',border_rows); wt('cds_table.tsv',cds_rows); write_fasta('all_proteins.faa',prot)
print('plasmids',len(plas_rows),'cds',len(cds_rows),'borders',len(border_rows))
