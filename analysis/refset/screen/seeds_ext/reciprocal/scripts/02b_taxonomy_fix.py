"""Re-resolve unreliable organism names of the reverse search by the taxid of the hit protein itself (NCBI esummary).
Why: nr titles are older than the taxdump; 02_taxonomy.py falls back from the full name to the cleaned name/genus, and homonyms
then pick the wrong kingdom (e.g. 'Bosea sp. RAC05' -> plant genus Bosea, because the bacterial genus is now Allobosea).
Unreliable name = matched at a level other than 'exact', or the matched name has taxids in more than one kingdom.
Modes: --collect (bsub: list unreliable names + up to 5 hit accessions each) -> --fetch (login node: esummary/esearch, cached)
       -> --apply (bsub: rewrite organism_taxonomy.json, backup + taxonomy_fix_log.tsv)."""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..', '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import sys, re, json, time, collections, shutil, os, urllib.request, urllib.parse
R=f'{GMO_ROOT}/analysis/refset/screen/seeds_ext/reciprocal/'
BR=re.compile(r'\[([^\[\]]+)\]')
TRAIL=re.compile(r'(?:\s*\[[^\[\]]+\])+\s*$')  # the organism bracket(s) close an nr title; brackets followed by text are description fragments
NESTED=re.compile(r'\[\[([^\[\]]+)\]\s+')      # "[[Torrubiella] hemipterigena]" -> "[Torrubiella hemipterigena]"
K={4751:'Fungi',2:'Bacteria',33090:'Viridiplantae',2157:'Archaea',33208:'Metazoa',10239:'Viruses'}
def orgs_from_stitle(s):
    """Organism names of an nr title: only the trailing bracket(s) of each \\x01-separated title, so that
    'acyl-[acyl-carrier-protein] desaturase [Quercus suber]' yields ['Quercus suber'] and not 'acyl-carrier-protein'."""
    out=[]
    for seg in s.split('\x01'):
        m=TRAIL.search(NESTED.sub(r'[\1 ',seg))
        if m: out.extend(o.strip() for o in BR.findall(m.group(0)))
    return list(dict.fromkeys(out))
def clean(o):
    return re.sub(r'\s+(strain|str\.|isolate|sp\.|subsp\.|var\.|f\. sp\.|pv\.|cf\.).*$','',o).strip()
def candidates(o):  # same order as 02_taxonomy.py
    w=clean(o).split(); c=[('exact',o.lower()),('cleaned',clean(o).lower())]
    if len(w)>=2: c.append(('binomial',' '.join(w[:2]).lower()))
    if len(w)>=1: c.append(('genus',w[0].lower()))
    return c
def load_nodes():
    par={};rank={}
    for l in open(R+'taxdump/nodes.dmp'):
        f=l.split('\t|\t'); par[int(f[0])]=int(f[1]); rank[int(f[0])]=f[2]
    merged={int(l.split('\t|\t')[0]):int(l.split('\t|\t')[1].split('\t')[0]) for l in open(R+'taxdump/merged.dmp')}
    return par,rank,merged
def lineage(t,par):
    out=[];seen=0
    while t in par and t!=par[t] and seen<100: out.append(t); t=par[t]; seen+=1
    return out
def kingdom(t,par):
    L=lineage(t,par)
    for a in L:
        if a in K: return K[a]
    return 'Eukaryota_other' if 2759 in L else 'NA'
def collect():
    tax=json.load(open(R+'organism_taxonomy.json'))
    want={}
    for o,d in tax.items():
        for lv,x in candidates(o):
            if lv==d.get('level'): want[o]=x; break
    names=set(want.values()); n2t=collections.defaultdict(set)
    for l in open(R+'taxdump/names.dmp'):
        f=l.rstrip('\n').split('\t|\t')
        if f[1].lower() in names: n2t[f[1].lower()].add(int(f[0]))
    par,_,_=load_nodes()
    unsafe={}
    for o,d in tax.items():
        lv=d.get('level')
        if lv!='exact': unsafe[o]=dict(old=d,reason=f'level={lv}',accessions=[],n_hit_lines=0); continue
        if len(o.split())==1: unsafe[o]=dict(old=d,reason='genus-only name',accessions=[],n_hit_lines=0); continue  # e.g. '[Bosea]': plant genus now, bacterial genus in older nr titles
        ts=n2t.get(want.get(o,''),set()); ks={kingdom(t,par) for t in ts}
        if len(ts)>1 and len(ks)>1: unsafe[o]=dict(old=d,reason='homonym:'+'/'.join(sorted(ks)),accessions=[],n_hit_lines=0)
    seen=set()
    for l in open(R+'diamond_nr.tsv'):
        f=l.rstrip('\n').split('\t',8)
        for o in orgs_from_stitle(f[8]):
            seen.add(o)
            if o in unsafe:
                u=unsafe[o]; u['n_hit_lines']+=1
                if f[1] not in u['accessions'] and len(u['accessions'])<5: u['accessions'].append(f[1])
    qnames={l.split('\t')[3] for l in open(R+'queries.tsv')}
    junk=[o for o in unsafe if o not in seen and o not in qnames]  # bracket fragments of descriptions, never a trailing organism bracket
    for o in junk: unsafe.pop(o)
    json.dump(junk,open(R+'nonorganism_brackets.json','w'),indent=0)
    json.dump(unsafe,open(R+'unsafe_names.json','w'),indent=0)
    print(len(junk),'non-organism bracket names dropped (nonorganism_brackets.json), e.g.',junk[:8])
    print(len(unsafe),'unreliable names;',collections.Counter(u['reason'].split(':')[0] for u in unsafe.values()),
          '; hit lines:',sum(u['n_hit_lines'] for u in unsafe.values()),'; names without accession:',sum(1 for u in unsafe.values() if not u['accessions']))
def eutils(tool,params,tries=3):
    url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/%s.fcgi?'%tool+urllib.parse.urlencode(params)
    for i in range(tries):
        try:
            with urllib.request.urlopen(url,timeout=60) as h: return json.load(h)
        except Exception as e:
            print('retry',i,tool,str(e)[:80]); time.sleep(2+2*i)
    return None
def fetch():
    unsafe=json.load(open(R+'unsafe_names.json'))
    cf=R+'eutils_cache.json'; cache=json.load(open(cf)) if os.path.exists(cf) else {}
    accs=sorted({a for u in unsafe.values() for a in u['accessions'] if a not in cache})
    print(len(accs),'accessions to fetch')
    for i in range(0,len(accs),150):
        b=accs[i:i+150]; j=eutils('esummary',dict(db='protein',id=','.join(b),retmode='json'))
        got={}
        if j and 'result' in j:
            for uid in j['result'].get('uids',[]):
                d=j['result'][uid]; got[d.get('caption','')]=dict(taxid=d.get('taxid'),organism=d.get('organism'))
        for a in b:
            cache[a]=got.get(a.split('.')[0],dict(taxid=None,organism=None,note='not returned'))
        json.dump(cache,open(cf,'w'),indent=0); time.sleep(0.4)
    nt=R+'name_taxid.json'; names=json.load(open(nt)) if os.path.exists(nt) else {}
    for o,u in unsafe.items():
        if any(cache.get(a,{}).get('taxid') for a in u['accessions']) or o in names: continue
        names[o]=None
        for term in dict.fromkeys([o,clean(o)]):  # full name, then the cleaned name; never the bare genus (homonyms)
            if len(term.split())<2: continue
            j=eutils('esearch',dict(db='taxonomy',term=term,retmode='json')); ids=(j or {}).get('esearchresult',{}).get('idlist',[]); time.sleep(0.4)
            if len(ids)==1: names[o]=int(ids[0]); break
    json.dump(names,open(nt,'w'),indent=0)
    print('unresolved by accession:',sum(1 for a in accs if not cache[a].get('taxid')),'; names via esearch:',sum(1 for v in names.values() if v),'of',len(names))
def apply():
    unsafe=json.load(open(R+'unsafe_names.json')); cache=json.load(open(R+'eutils_cache.json'))
    nt=R+'name_taxid.json'; names=json.load(open(nt)) if os.path.exists(nt) else {}
    par,rank,merged=load_nodes(); need=set(); new={}
    for o,u in unsafe.items():
        ts=[merged.get(int(cache[a]['taxid']),int(cache[a]['taxid'])) for a in u['accessions'] if cache.get(a,{}).get('taxid')]
        lv='accession'
        if not ts and names.get(o): ts=[merged.get(names[o],names[o])]; lv='esearch'
        if not ts and u['reason']=='genus-only name' and u['old'].get('level')=='exact':  # e.g. 'bacterium': exact taxdump match, no accession evidence against it
            old=u['old']; new[o]=dict(kingdom=old['kingdom'],cls=old['cls'],order=old['order'],genus=old['genus'],level='exact',taxid=None,flag='kept_exact',organism_ncbi=''); continue
        if not ts: new[o]=dict(kingdom='NA',cls='',order='',genus=clean(o).split()[0] if o else '',level='unresolved',taxid=None,flag='unresolved',organism_ncbi=''); continue
        ks=collections.Counter(kingdom(t,par) for t in ts); k,nk=ks.most_common(1)[0]
        t=collections.Counter(t for t in ts if kingdom(t,par)==k).most_common(1)[0][0]
        d={}
        for a in lineage(t,par):
            if rank.get(a) in('genus','order','class') and rank[a] not in d: d[rank[a]]=a; need.add(a)
        org=next((cache[a]['organism'] for a in u['accessions'] if str(cache.get(a,{}).get('taxid'))==str(t)),'')
        new[o]=dict(kingdom=k,cls=d.get('class'),order=d.get('order'),genus=d.get('genus'),level=lv,taxid=t,organism_ncbi=org,
                    flag='conflict:'+'/'.join(f'{a}={b}' for a,b in ks.items()) if len(ks)>1 else '')
    sci={}
    for l in open(R+'taxdump/names.dmp'):
        f=l.rstrip('\n').split('\t|\t')
        if 'scientific' in f[3] and int(f[0]) in need: sci[int(f[0])]=f[1]
    bak=R+'organism_taxonomy.json.bak_20261013'  # the output of 02_taxonomy.py; made once, re-runs start from it
    if not os.path.exists(bak): shutil.copy(R+'organism_taxonomy.json',bak)
    tax=json.load(open(bak))
    trans=collections.Counter(); lines=collections.Counter(); changed=[]
    with open(R+'taxonomy_fix_log.tsv','w') as f:
        f.write('\t'.join(['name','reason','old_level','old_kingdom','old_genus','new_level','new_kingdom','new_cls','new_order','new_genus','taxid','organism_ncbi','n_accessions','n_hit_lines','flag'])+'\n')
        for o,n in new.items():
            u=unsafe[o]; old=u['old']
            for k2 in('cls','order','genus'):
                if isinstance(n[k2],int): n[k2]=sci.get(n[k2],'')
                n[k2]=n[k2] or ''
            if not n['genus']: n['genus']=clean(o).split()[0] if o else ''
            tax[o]=dict(kingdom=n['kingdom'],cls=n['cls'],order=n['order'],genus=n['genus'],level=n['level'])
            trans[(old['kingdom'],n['kingdom'])]+=1; lines[(old['kingdom'],n['kingdom'])]+=u['n_hit_lines']
            if old['kingdom']!=n['kingdom']: changed.append((o,old['kingdom'],n['kingdom'],u['n_hit_lines']))
            f.write('\t'.join(str(x) for x in [o,u['reason'],old['level'],old['kingdom'],old['genus'],n['level'],n['kingdom'],n['cls'],n['order'],n['genus'],n['taxid'] or '',n['organism_ncbi'] or '',len(u['accessions']),u['n_hit_lines'],n['flag']])+'\n')
    json.dump(tax,open(R+'organism_taxonomy.json','w'))
    print('kingdom transitions (names / hit lines):')
    for k2,v in sorted(trans.items(),key=lambda kv:-kv[1]): print(' ',k2[0],'->',k2[1],v,'/',lines[k2])
    print('names with kingdom change:',len(changed))
    for c in sorted(changed,key=lambda c:-c[3]): print(' ',*c,sep='\t')
    print('unresolved:',sum(1 for n in new.values() if n['level']=='unresolved'),'; conflicts:',sum(1 for n in new.values() if n['flag'].startswith('conflict')))
    for o in('Bosea sp. RAC05','Bacillus','Bosea sp. (in: a-proteobacteria)'):
        if o in tax: print(' check',o,tax[o])
if __name__=='__main__':
    {'--collect':collect,'--fetch':fetch,'--apply':apply}[sys.argv[1]]()
