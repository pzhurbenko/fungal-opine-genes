import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import csv,collections
from Bio import Phylo
TR=f'{GMO_ROOT}/analysis/report_A/trees'
SIX=['JBBPKL010003317','JBBPKL010018021','JBBPNB010000203','JBBPNB010000771','JBBPOG010012756','JBBPRY010001778']
def load(fam):
    rows=list(csv.DictReader(open(f'{TR}/{fam}_leaf_labels.tsv'),delimiter='\t'))
    L={}
    for x in rows:L.setdefault(x['label'],x)
    t=Phylo.read(f'{TR}/{fam}_final_MAFFT_IQTREE.nwk','newick')
    og=[l.name for l in t.get_terminals() if l.name.startswith('OG__')]
    t.root_with_outgroup(og)
    return t,L
def grp(L,n): return L[n]['group'] if n in L else 'OG'
def leaves(t,L,groups,excl=''):
    return [l.name for l in t.get_terminals() if grp(L,l.name) in groups and not(excl and excl in l.name)]
def clade(t,names): return {x.name for x in t.common_ancestor(names).get_terminals()}
def sup(c):
    s=c.name if c.name else (str(c.confidence) if c.confidence is not None else '')
    return s
def categories(fam,t,L):
    """returns dict leaf->category key"""
    cat={}
    T=t.get_terminals()
    if fam=='SUS':
        c=clade(t,leaves(t,L,{'REF_sus'}))
        for l in T:
            g=grp(L,l.name)
            if l.name in c: k='sus'
            elif g=='REF_ocd_rhizobiales_chromosomal': k='chromOCD'
            elif g=='REF_ocd_other_bacteria': k='bactOCD'
            elif g=='FUNGI_sus': k='fungOCD'
            elif g=='NATIVE': k='native'
            else: k='OG'
            cat[l.name]=k
    elif fam=='CUS':
        ct=clade(t,leaves(t,L,{'REF_cus'})); ac=clade(t,leaves(t,L,{'FUNGI_cus_actinobacteria'})); mi=clade(t,leaves(t,L,{'FUNGI_mis'}))
        for l in T:
            n=l.name;g=grp(L,n)
            if n in ct:k='cus'
            elif n in ac:k='actino'
            elif n in mi:k='mis'
            elif g=='OG':k='OG'
            elif g.startswith('FUNGI'):k='fungother'
            else:k='bact'
            cat[n]=k
    else:
        td=clade(t,leaves(t,L,{'REF_acs'},excl='PON40589')); ac=clade(t,leaves(t,L,{'REF_gdpd_meso'}))
        fung=[l.name for l in T if grp(L,l.name).startswith('FUNGI')]
        six=[n for n in fung if any(s in n for s in SIX)]
        kerr=six+leaves(t,L,{'REF_kerr14'})+[l.name for l in T if 'PON40589' in l.name]
        sx=clade(t,kerr)
        nat=[n for n in t.get_nonterminals() if n.count_terminals()==94 and all(grp(L,x.name)=='FUNGI_acs' for x in n.get_terminals())]
        natc={x.name for x in nat[0].get_terminals()}
        for l in T:
            n=l.name;g=grp(L,n)
            if n in td:k='tdna'
            elif n in ac:k='acc'
            elif n in sx:k='kerr'
            elif n in natc:k='natfung'
            elif g=='OG':k='OG'
            elif g.startswith('FUNGI'):k='fungunres'
            elif g=='NATIVE':k='native'
            else:k='bact'
            cat[n]=k
    return cat
def collapse(node,cat):
    """returns nested: ('leafset',cat,[names]) or ('node',support,[children])"""
    names=[x.name for x in node.get_terminals()]
    cs={cat[n] for n in names}
    if len(cs)==1: return ('tri',cs.pop(),names,node)
    return ('node',sup(node),[collapse(c,cat) for c in node.clades],node)
