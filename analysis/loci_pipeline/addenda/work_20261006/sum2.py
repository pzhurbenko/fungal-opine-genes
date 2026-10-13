import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import re,os
from collections import defaultdict
D=f'{GMO_ROOT}/analysis/refset/screen/step2'
GAP=5000; FAMS=['ACS','CUS','SUS']
def ps(s):
    p=s.split('|')
    if p[0].startswith('FAM_'): return p[0][4:]
    return p[1]
for line in open(f'{D}/extra2.tsv'):
    g,fna,_=line.rstrip('\n').split('\t'); org=os.path.basename(fna).replace(g+'_','').replace('_genomic.fna','')
    od=f'{D}/out/{g}'
    coord={}
    for l in open(f'{od}/hit_orfs.faa'):
        m=re.match(r'>(\S+) source=(\S+) coords=(\d+)\.\.(\d+) length=(\d+) frame=(\d+)',l)
        if m: a,b=int(m[3]),int(m[4]); coord[m[1]]=(m[2],min(a,b),max(a,b),'+' if int(m[6])<=3 else '-')
    per=defaultdict(dict)
    for l in open(f'{od}/hmm.tbl'):
        if l[0]=='#': continue
        f=l.split(); per[f[0]][f[2]]=(float(f[5]),float(f[4]))
    conf=[];amb=0;outg=0
    for orf,d in per.items():
        for fam in FAMS:
            subs={k:v for k,v in d.items() if k.startswith(fam+'__') and not k.endswith('OUTGROUP')}
            og=d.get(f'{fam}__OUTGROUP',(0,99))[0]
            if not subs: continue
            bs=max(subs,key=lambda k:subs[k][0]); sc=subs[bs][0]; gap=sc-og
            if og>sc: outg+=1
            elif gap<20 or sc<30: amb+=1
            else:
                c=coord.get(orf,('?',0,0,'?')); conf.append(f'{fam}:{bs.split("__")[1]} score={sc:.0f} gap={gap:.0f} {c[0]}:{c[1]}-{c[2]}({c[3]}) {orf}')
    # diamond
    hs=defaultdict(list); clen={}
    for l in open(f'{od}/dmnd.tsv'):
        f=l.split('\t')
        if len(f)<13: continue
        q,s=f[0],f[1]; hs[(q,'+' if int(f[12])>0 else '-',s)].append((*sorted((int(f[6]),int(f[7]))),int(f[8]),int(f[9]),float(f[2]),float(f[11]),int(f[5])))
    chains=defaultdict(list)
    for (q,st,s),hl in hs.items():
        hl.sort(); cur=None
        for qs,qe,ss,se,pid,bits,slen in hl:
            if cur and qs-cur['end']<=GAP:
                cur['end']=max(cur['end'],qe);cur['bits']+=bits;cur['scov'].append((min(ss,se),max(ss,se)));cur['pid']=max(cur['pid'],pid)
            else:
                if cur: chains[(q,st)].append(cur)
                cur=dict(start=qs,end=qe,subj=s,bits=bits,scov=[(min(ss,se),max(ss,se))],pid=pid,slen=slen)
        chains[(q,st)].append(cur)
    strong=[]
    for (q,st),cl in chains.items():
        cl.sort(key=lambda c:c['start']); loci=[];cur=None
        for c in cl:
            if cur and c['start']<=cur['end']: cur['end']=max(cur['end'],c['end']);cur['m'].append(c)
            else:
                if cur: loci.append(cur)
                cur=dict(start=c['start'],end=c['end'],m=[c])
        loci.append(cur)
        for L in loci:
            b=max(L['m'],key=lambda m:m['bits']); mg=[]
            for a,e in sorted(b['scov']):
                if mg and a<=mg[-1][1]: mg[-1][1]=max(mg[-1][1],e)
                else: mg.append([a,e])
            sc=sum(e-a+1 for a,e in mg)/b['slen']; cat=ps(b['subj'])
            if b['pid']>=45 and sc>=0.5 and cat in('T1','T2','TDNA','opine_related'):
                strong.append(f'{q}:{L["start"]}-{L["end"]}({st}) {b["subj"][:90]} pid={b["pid"]:.0f} scov={sc:.2f}')
    print(f'=== {g} {org}\n  confident HMM calls ({len(conf)}); ambiguous={amb}; OUTGROUP={outg}')
    for c in conf: print('   ',c)
    print(f'  strong diamond loci ({len(strong)})')
    for s in strong: print('   ',s)
