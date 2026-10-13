#!/usr/bin/env python3
"""Figure 3 (redesign 2026-10-13): phylograms (real branch lengths) of the six families; focal structure open, rest collapsed to grey triangles.
Run: python make_fig3.py
Rooting / categories as in paper/tables/build_table1.py (OG__ leaves, SEED_ags, SEED_ocs/vis). Previous version: make_fig3.py.bak_20261013x."""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import os
import sys, os, re, collections, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from common import load, categories, grp, sup, collapse

plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
                     'font.size': 6.0, 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
FS = 6.0
W = 17.0  # cm
CM = 1 / 2.54
RED, BLUE, GREEN, PURP = '#D7301F', '#0072B2', '#009E73', '#7B3294'
LRED, LBLUE, LGREEN, LPURP = '#F4B6AE', '#A9D0E8', '#A6DCCB', '#CDB7DD'
DG, MG, LG, NAT = '#6E6E6E', '#A0A0A0', '#D9D9D9', '#E6E6E6'
ORANGE, LORANGE = '#E69F00', '#F6D9A3'   # chs-like
PINK, LPINK = '#CC79A7', '#EBC3DA'       # nos-ODH Boletales
BROWN, LBROWN = '#A0522D', '#DDB89F'     # nos-ODH Rhypophila
INDIGO, LINDIGO = '#332288', '#B8B0DC'   # odh island
TEAL, LTEAL = '#1B7F8C', '#A8D5DB'       # A. larrymoorei plasmid ODH (donor)
INFO = {}  # numbers for report
FUNG_DOM = {'native', 'fungunres', 'fungOCD', 'fungother', 'agsnat'}  # fungal grey blocks -> light grey; bacterial/reference -> mid grey
MERGED = []  # report of merged ladder blocks

def esc(s): return s
# ---------- rich text: *italic*, **bold** segments ----------
def rich(ax, x, y, text, fs=FS, color='k', bold=False, va='center'):
    """draw text with *italic* parts starting at data x; returns end x (data)"""
    fig = ax.figure; r = fig.canvas.get_renderer()
    parts = re.split(r'(\*[^*]+\*)', text)
    cx = x
    for p in parts:
        if not p: continue
        it = p.startswith('*') and p.endswith('*')
        s = p[1:-1] if it else p
        t = ax.text(cx, y, s, fontsize=fs, color=color, fontweight='bold' if bold else 'normal',
                    fontstyle='italic' if it else 'normal', ha='left', va=va)
        bb = t.get_window_extent(r)
        inv = ax.transData.inverted()
        x1 = inv.transform((bb.x1, bb.y0))[0]
        # keep trailing spaces: matplotlib extent ignores them -> add space width manually
        cx = x1
        if s.startswith(' '): pass
    return cx

# ---------- family definitions ----------
NAMES = {
 'SUS': {'chromOCD': 'Rhizobiales chromosomal OCD', 'bactOCD': 'other bacterial OCD', 'fungOCD': 'other fungal OCD-like',
         'native': 'native fungal OCD', 'OG': 'outgroup'},
 'CUS': {'bact': 'other bacterial CUS', 'fungother': 'native fungal gene', 'OG': 'outgroup'},
 'ACS': {'bact': 'bacterial GDPD', 'fungunres': 'unresolved fungal GDPD', 'native': 'other native fungal GDPD', 'OG': 'outgroup'},
 'CHS': {'bact': 'other bacterial M24 aminopeptidases', 'agsnat': 'native fungal Ags homologs', 'OG': 'outgroup (Ags seed)'},
 'NOS': {'noslike': 'rhizobial nos-like (non-T-DNA)', 'bact': 'other bacterial ODH', 'native': 'native fungal ODH', 'OG': 'outgroup (ocs/vis seeds)'},
 'ODH': {'bact': 'other bacterial ODH', 'native': 'native fungal ODH', 'OG': 'outgroup (ocs/vis seeds)'},
}
KEY = {'SUS': {'sus'}, 'CUS': {'cus', 'actino', 'mis'}, 'ACS': {'tdna', 'acc', 'kerr', 'natfung'},
       'CHS': {'fus', 'dac', 'tichs'}, 'NOS': {'nosbol', 'nosrhy', 'tdnanos', 'isl'}, 'ODH': {'isl', 'larry', 'nossub'}}
COL = {'sus': (RED, LRED), 'cus': (BLUE, LBLUE), 'acc': (GREEN, LGREEN), 'kerr': (PURP, LPURP), 'tdna': (DG, MG), 'actino': (DG, MG),
       'mis': ('#8C8C8C', '#BDBDBD'), 'natfung': ('#8C8C8C', LG),
       'fus': (ORANGE, LORANGE), 'dac': (ORANGE, LORANGE), 'tichs': (DG, MG),
       'nosbol': (PINK, LPINK), 'nosrhy': (BROWN, LBROWN), 'tdnanos': (DG, MG), 'nossub': (DG, MG),
       'isl': (INDIGO, LINDIGO), 'larry': (TEAL, LTEAL)}

def n_plant(names): return sum('Plant' in n for n in names)

def key_label(fam, k, names, L):
    g = collections.Counter(grp(L, n) for n in names)
    nf = sum(v for a, v in g.items() if a.startswith('FUNGI')); nref = len(names) - nf
    if fam in GEN:
        nf = sum(is_fungal_leaf(n) for n in names); nref = len(names) - nf
    npl = n_plant([n for n in names if not grp(L, n).startswith('FUNGI')])
    INFO[(fam, k)] = dict(n=len(names), fungal=nf, nonfungal=nref, plant=npl, groups=dict(g))
    if fam == 'SUS':
        return ['T-DNA *sus* clade', f'{nref - npl} Ti/Ri plasmid + {npl} plant cT-DNA references; {nf} fungal (*sus*-like)']
    if fam == 'CUS':
        if k == 'cus':
            nrc = g['REF_cus']; plc = n_plant([n for n in names if grp(L, n) == 'REF_cus']); oth = nref - nrc
            return ['T-DNA *cus* clade', f'{nrc - plc} Ti/Ri plasmid + {plc} plant cT-DNA references (+{oth} other bacterial); {nf} fungal (*cus*-like)',
                    'fungal sequences interleaved with references']
        if k == 'actino':
            return ['Actinobacterial-type *cus* clade', f'{nf} fungal sequences only (no bacterial or plant references)',
                    'contains the *Beauveria bassiana* gene *Fcs1* (Li et al. [61])']
        if k == 'mis':
            nm = g['FUNGI_mis'] + g['FUNGI_cusrel181']
            return ['*mis* lineage', f'{nm} fungal + {nref} reference sequences (Ti/Ri plasmid, plant cT-DNA, bacterial)']
    if fam == 'ACS':
        if k == 'tdna': return ['T-DNA *acs* clade', f'{g["REF_acs"]} Ti/Ri plasmid + plant cT-DNA references + {g["BACT_acc"]} bacterial; no fungal']
        if k == 'acc': return ['Acc-type GDPD', f'{nf} fungal (Boletales) + {nref} bacterial acc-locus references']
        if k == 'kerr':
            return ['pKerr14a-type *acs*-like', f'{nf} fungal (Boletales) + pKerr14a *acs* + Parasponia cT-DNA fragment']
        if k == 'natfung': return ['Native fungal GDPD (main clade)', f'{nf} fungal sequences only (no bacterial or plant references)']

    if fam == 'CHS':
        if k == 'fus': return ['*Fusarium* *chs*-like', f'{nf} sequences ({len({n.split("_")[2] for n in names})} species)']
        if k == 'dac': return ['*Dactylonectria estremocensis* *chs*-like', '1 sequence (KAH7120428.1)']
        if k == 'tichs': return ['Ti-plasmid *chs* clade', f'*chs* seed (pTi14-818c) + {len(names) - 1} other agrobacterial *chs*']
    if fam == 'NOS':
        if k == 'nosbol': return ['*nos*-subfamily ODH, Boletales', f'{nf} loci (*Boletus*, *Butyriboletus*)']
        if k == 'nosrhy': return ['*Rhypophila* *nos* homolog (not assigned)', f'{nf} sequences (2 species)']
        if k == 'tdnanos':
            nt = sum(n.startswith('NOS_TDNA') for n in names); nl = len(names) - nt
            return ['T-DNA *nos* clade', f'{nt} Ti-plasmid *nos* references + {nl} rhizobial *nos*-like']
        if k == 'isl': return ['*odh*-island ODH (*F. poae*) + *A. larrymoorei* plasmid ODH', 'not in the *nos* subfamily']
    if fam == 'ODH':
        if k == 'isl': return ['*odh*-island ODH (*Fusarium*)', f'{nf} sequences (*F. poae* ×2, *F. langsethiae*)']
        if k == 'larry': return ['*A. larrymoorei* plasmid ODH', 'WP_309770759 (donor locus)']
        if k == 'nossub':
            return ['*nos* subfamily (T-DNA type)', 'T-DNA *nos* seeds (2), *Rhypophila* (2), *Butyriboletus* (1), WP_309770762']

def mixed_label(fam, names, cat, L):
    c = collections.Counter(cat[n] for n in names)
    nm = NAMES[fam]
    items = [f'{nm.get(k, k)} ({v})' for k, v in c.most_common()]
    return items


# ---------- generic schematic trees (chs / nos-subfamily / odh island) ----------
from Bio import Phylo
A = f'{GMO_ROOT}/analysis'
GEN = {
 'CHS': dict(tree=f'{A}/chs_fusarium/chs_tree.iqtree.treefile', og=lambda n: n.startswith('SEED_ags')),
 'NOS': dict(tree=f'{A}/fusarium_odh_island/nos_check/nos_tree.iqtree.treefile', og=lambda n: n.startswith('SEED_ocs') or n.startswith('SEED_vis')),
 'ODH': dict(tree=f'{A}/fusarium_odh_island/tree/odh_tree.iqtree.treefile', og=lambda n: n.startswith('SEED_ocs') or n.startswith('SEED_vis')),
}
def is_fungal_leaf(n):
    return n.startswith(('CHSlike_', 'AGSnat_', 'BOL_', 'FUNG_', 'NAT_', 'ISL_'))
def _clade_of(t, names): return {x.name for x in t.common_ancestor(names).get_terminals()}
def categories_generic(fam, t):
    T = [l.name for l in t.get_terminals()]; cat = {}
    if fam == 'CHS':
        tichs = _clade_of(t, [n for n in T if n.startswith('SEED_chs')] + ['BACT_WP_084631803.1'])
        for n in T:
            if n in tichs: k = 'tichs'
            elif n.startswith('CHSlike_Fusarium'): k = 'fus'
            elif n.startswith('CHSlike_Dactylonectria'): k = 'dac'
            elif n.startswith('AGSnat_'): k = 'agsnat'
            elif n.startswith('SEED_ags'): k = 'OG'
            else: k = 'bact'
            cat[n] = k
    elif fam == 'NOS':
        tdna = _clade_of(t, [n for n in T if n.startswith('NOS_TDNA')])
        isl = _clade_of(t, [n for n in T if n.startswith('ISL_') or 'larrymoorei' in n])
        for n in T:
            if n in tdna: k = 'tdnanos'
            elif n.startswith('BOL_'): k = 'nosbol'
            elif n.startswith('FUNG_Rhypophila'): k = 'nosrhy'
            elif n in isl: k = 'isl'
            elif n.startswith('NOSLIKE_'): k = 'noslike'
            elif n.startswith('NAT_'): k = 'native'
            elif n.startswith('SEED_'): k = 'OG'
            else: k = 'bact'
            cat[n] = k
    elif fam == 'ODH':
        nossub = _clade_of(t, [n for n in T if n.startswith('SEED_nos') or n.startswith('FUNG_') or n.startswith('BACT2_')])
        for n in T:
            if n.startswith('ISL_'): k = 'isl'
            elif n == 'BACT_Agrobacterium_larrymoorei_WP_309770759': k = 'larry'
            elif n in nossub: k = 'nossub'
            elif n.startswith('NAT_'): k = 'native'
            elif n.startswith('SEED_'): k = 'OG'
            else: k = 'bact'
            cat[n] = k
    return cat
def load_generic(fam):
    t = Phylo.read(GEN[fam]['tree'], 'newick')
    og = [l.name for l in t.get_terminals() if GEN[fam]['og'](l.name)]
    t.root_with_outgroup(og)
    INFO[(fam, 'n_leaves')] = len(t.get_terminals())
    return t, {}

# =====================================================================================
# Phylogram drawing (redesign 2026-10-13): real branch lengths; focal structure open, the rest = grey triangles
# =====================================================================================
COL.update({'accref': (DG, '#DCEFE8'), 'mis': ('#5B7FA6', '#C9D6E6'), 'actino': ('#3A6EA5', '#B7CCE3')})
GREYNAME = {'OG': 'outgroup', 'chromOCD': 'Rhizobiales chromosomal OCD', 'bactOCD': 'other bacterial OCD',
            'fungOCD': 'other fungal OCD-like', 'native': 'native fungal OCD'}
GREYNAME_FAM = {
 'SUS': GREYNAME,
 'CUS': {'bact': 'other bacterial CUS', 'fungother': 'native fungal CUS', 'OG': 'outgroup'},
 'ACS': {'bact': 'bacterial GDPD', 'fungunres': 'unresolved fungal GDPD', 'native': 'native fungal GDPD', 'natfung': 'native fungal GDPD', 'OG': 'outgroup'},
 'CHS': {'bact': 'other bacterial M24', 'agsnat': 'native fungal Ags homologs', 'OG': 'outgroup'},
 'NOS': {'noslike': 'rhizobial nos-like', 'bact': 'other bacterial ODH', 'native': 'native fungal ODH', 'OG': 'outgroup'},
 'ODH': {'bact': 'other bacterial ODH', 'native': 'native fungal ODH', 'OG': 'outgroup'},
}
FOCAL_LABEL = {'sus': '*sus* lineage, fungal', 'cus': '*cus* lineage, fungal', 'tdna': 'T-DNA references', 'actino': 'Actinobacterial-type *cus*, fungal',
               'mis': '*mis* lineage', 'acc': 'Acc-type GDPD, fungal', 'accref': 'Acc-locus bacterial references', 'kerr': 'pKerr14a-type, fungal',
               'fus': '*Fusarium chs*-like', 'dac': '*Dactylonectria chs*-like', 'tichs': 'Ti-plasmid *chs* references',
               'nosbol': '*nos*-subfamily ODH, Boletales', 'nosrhy': '*Rhypophila nos* homolog', 'isl': '*odh*-island ODH',
               'larry': '*A. larrymoorei* plasmid ODH'}
FOCAL_LABEL_FAM = {('NOS', 'tdna'): 'T-DNA *nos* + rhizobial *nos*-like', ('NOS', 'tdnanos'): 'T-DNA *nos* + rhizobial *nos*-like',
                   ('ODH', 'tdna'): 'T-DNA *nos* seeds / *nos*-type', ('CHS', 'tdna'): 'Ti-plasmid *chs* references',
                   ('SUS', 'tdna'): 'T-DNA *sus* references'}
FOCAL_SETS = {'SUS': {'sus'}, 'CUS': {'cus', 'actino', 'mis'}, 'ACS': {'tdna', 'acc', 'kerr'}, 'CHS': {'fus', 'dac', 'tichs'},
              'NOS': {'tdnanos', 'nosbol', 'nosrhy', 'isl'}, 'ODH': {'isl', 'larry', 'nossub'}}
TRI_CAP = 4.0
ROWH = 0.2     # cm per tip row
MAX_H = 23.9     # cm, total figure height limit (240 mm minus margin)
TREE_W = 6.4     # cm, width of the tree area
X0 = 1.0
XL = X0 + TREE_W + 0.2   # label column
XB = XL + 3.15           # bracket column
SPECIES_SKIP = {'Plant', 'Agro', 'T1', 'T2', 'T2x', 'MULTISPECIES', 'Rhizobiaceae', 'Gillisella'}

def short_name(n):
    plant = 'Plant' in n.split('__')[-1].split('_')[:1] or n.startswith(('REF_T1_sus__Plant', 'REF_T1_acs__Plant', 'REF_T2_acs__Plant')) or '__Plant_' in n
    tail = n.split('__')[-1] if '__' in n else re.sub(r'^(BACT2|BACT|SEED|ISL|FUNG|BOL|NOS_TDNA|NOSLIKE|CHSlike|AGSnat|NAT)_', '', n)
    toks = tail.split('_')
    nm = None
    for a, b in zip(toks, toks[1:]):
        if a in SPECIES_SKIP: continue
        if re.fullmatch(r'[A-Z][a-z]+', a) and re.fullmatch(r'[a-z]+\.?', b) and b not in ('protein', 'family', 'synthase', 'isolate', 'chromosome', 'extra', 'hypothetical'):
            nm = f'{a} {b}'; break
    if nm is None:
        for a in toks:
            if re.fullmatch(r'[A-Z][a-z]{3,}', a) and a not in SPECIES_SKIP:
                nm = a + ' sp.'; break
    if nm is None:
        m_ = re.search(r'(WP_\d{6,}(?:[._]\d)?)', tail)
        if m_: nm = m_.group(1).replace('_1', '.1') if m_.group(1).endswith('_1') else m_.group(1)
    if nm is None:
        nm = re.sub(r'^(Plant|Agro)_', '', tail)[:24].replace('_', ' ')
    nm = nm.replace('.', '') if False else nm
    if nm.endswith(' aff.'):
        mm = re.search(r'__[A-Z]+_([A-Z0-9]{6,16})_L\d', n)
        if mm: nm += f' ({mm.group(1)[:6]})'
    if plant: nm += ' (plant)'
    if n.startswith('BACT2_'): nm += ' (nos-type)'
    if len(nm) > 25:
        p = nm.split(' ')
        nm = f'{p[0][0]}. {" ".join(p[1:])}'
        if len(nm) > 25: nm = nm[:24] + '…'
    return nm

def style_of(fam, n, c, L):
    """colour/style key of a leaf, or None if it is not part of the focal structure"""
    if c not in FOCAL_SETS[fam]: return None
    fung = n.startswith(('FUNGI', 'CHSlike_', 'AGSnat_', 'BOL_', 'FUNG_', 'NAT_', 'ISL_'))
    if fam == 'CUS':
        if c in ('actino', 'mis'): return c
        return 'cus' if fung else ('tdna' if grp(L, n) == 'REF_cus' else None)
    if c in ('sus', 'cus', 'acc', 'kerr'):
        return c if fung else ('accref' if c == 'acc' else 'tdna')
    if c == 'nossub': return 'nosrhy' if 'Rhypophila' in n else ('nosbol' if n.startswith('BOL_') else 'tdna')
    if 'larrymoorei' in n and c == 'isl': return 'larry'
    if c == 'tdnanos': return 'tdna'
    return c

def build_all(fam):
    if fam in GEN:
        t, L = load_generic(fam); cat = categories_generic(fam, t)
    else:
        t, L = load(fam); cat = categories(fam, t, L)
        if fam == 'CUS':
            t.prune('FUNGI_cus_intact__CUS_JBBPKX010002805_L1__Boletus_austroedulis')
    cat = dict(cat); sk = {}
    for n, c in cat.items():
        sk[n] = style_of(fam, n, c, L)
        if sk[n] is None and fam == 'CUS' and c == 'cus': cat[n] = 'bact'
    X = {id(t.root): 0.0}
    def setx(nd):
        for c in nd.clades:
            X[id(c)] = X[id(nd)] + max(c.branch_length or 0.0, 0.0); setx(c)
    setx(t.root)
    leafs = {}
    def lf(nd):
        leafs[id(nd)] = [nd.name] if nd.is_terminal() else [x for c in nd.clades for x in lf(c)]
        return leafs[id(nd)]
    lf(t.root)
    maxx = {}
    def mx(nd):
        maxx[id(nd)] = X[id(nd)] if nd.is_terminal() else max(mx(c) for c in nd.clades)
        return maxx[id(nd)]
    mx(t.root)
    return t, L, cat, sk, X, leafs, maxx

class D: pass
def make_display(fam, thr):
    t, L, cat, sk, X, leafs, maxx = build_all(fam)
    FS_ = FOCAL_SETS[fam]
    def grey_label(names):
        c = collections.Counter(cat[x] for x in names).most_common(1)[0][0]
        return f'{GREYNAME_FAM[fam].get(c, c)} ({len(names)})', c
    def rec(nd, ingroup=None):
        d = D(); names = leafs[id(nd)]; ks = {sk[x] for x in names}
        d.clade = nd; d.x = X[id(nd)]; d.n = len(names); d.group = None; d.kids = []; d.names = names
        d.depth = maxx[id(nd)]
        if ks == {None} and ingroup is None:
            d.kind = 'tri'; d.style = None; d.label, d.dom = grey_label(names); d.rows = min(TRI_CAP, max(0.85, len(names) / 30))
            return d
        if len(ks) == 1 and ingroup is None:
            s = next(iter(ks))
            if len(names) > thr:
                d.kind = 'tri'; d.style = s; d.rows = min(TRI_CAP, max(0.85, len(names) / 30)); d.label = None; return d
            d.group = s
            ingroup = s
        if nd.is_terminal():
            d.kind = 'tip'; d.style = sk[nd.name]; d.rows = 1.0; d.label = short_name(nd.name); return d
        d.kind = 'node'; d.kids = [rec(c, ingroup) for c in nd.clades]
        d.style = None; d.rows = 0
        return d
    root = rec(t.root)
    # y layout
    y = [0.0]
    def lay(d):
        if d.kids:
            for k in d.kids: lay(k)
            d.y = (d.kids[0].y + d.kids[-1].y) / 2
        else:
            d.y = y[0] + d.rows / 2; y[0] += d.rows + 0.0
    lay(root)
    # supports: focal stems, cat-level focal clades, and their parents
    S = {}
    def marks(d, parent_cat=False, inside=False):
        if d.kind == 'tip' or inside: return False
        cs = {cat[x] for x in d.names}
        catlev = len(cs) == 1 and cs <= FS_
        f = (d.kind == 'tri' and d.style is not None) or d.group is not None or (catlev and not parent_cat)
        any_kid = False
        for k in d.kids:
            if marks(k, catlev or parent_cat, d.group is not None): any_kid = True
        if (f or any_kid) and d is not root: S[id(d.clade)] = 0 if (catlev and not parent_cat) else (1 if f else 2)
        return f
    marks(root)
    return t, root, y[0], S, X, maxx

def count_items(d, c):
    if d.kind == 'tip': c['tips'] += 1
    elif d.kind == 'tri': c['tri_focal' if d.style else 'tri_grey'] += 1
    for k in d.kids: count_items(k, c)

def draw_panel(ax, fam, root, rows, S, hp, letter, title):
    maxdepth = [0.0]
    def md(d):
        maxdepth[0] = max(maxdepth[0], d.depth if d.kind == 'tri' else d.x)
        for k in d.kids: md(k)
    md(root)
    sc = (TREE_W - 0.15) / maxdepth[0]
    top = hp - 0.4
    PX = lambda v: X0 + v * sc
    PY = lambda r: top - r * ROWH
    lc = '#333333'
    brackets = []; cand = []
    ax.text(0.05, hp - 0.08, letter, fontsize=9, fontweight='bold', va='top')
    rich(ax, 0.75, hp - 0.1, title, fs=FS + 1, bold=True, va='top')
    def stylec(s): return COL[s][0] if s else '#808080'
    def rec(d, px):
        x = PX(d.x); y = PY(d.y)
        ax.plot([px, x], [y, y], color=lc, lw=0.7, solid_capstyle='butt')
        if d.kind == 'node':
            ys = [PY(k.y) for k in d.kids]
            ax.plot([x, x], [min(ys), max(ys)], color=lc, lw=0.7, solid_capstyle='butt')
            for k in d.kids: rec(k, x)
            if d.group is not None:
                brackets.append((PY(d.kids[0].y if False else min_y(d)), PY(max_y(d)), d))
        elif d.kind == 'tri':
            w = max((d.depth - d.x) * sc, 0.2); h = d.rows * ROWH * 0.94
            if d.style:
                fc, ec = COL[d.style][1], COL[d.style][0]
            else:
                fc = 'white' if d.dom == 'OG' else (LG if d.dom in FUNG_DOM else MG); ec = '#808080'
            ax.add_patch(Polygon([[x, y], [x + w, y + h / 2], [x + w, y - h / 2]], closed=True, fc=fc, ec=ec, lw=0.7, zorder=2))
            lab = d.label if d.style is None else f'{FOCAL_LABEL_FAM.get((fam, d.style), FOCAL_LABEL[d.style])} ({d.n})'
            rich(ax, XL, y, lab, color=(COL[d.style][0] if d.style else 'k'), bold=bool(d.style))
            if d.style and d.style in ('tdna',): pass
            ax.plot([x + w + 0.03, XL - 0.05], [y, y], color='#BBBBBB', lw=0.4, ls=(0, (1, 2)), zorder=1)
        else:
            col = stylec(d.style)
            ax.plot([x], [y], marker='o', ms=2.4, color=col, mew=0, zorder=3)
            ax.text(XL, y, d.label, fontsize=FS, va='center', ha='left', fontstyle='italic', color='#222222' if d.style == 'tdna' else col)
            ax.plot([x + 0.06, XL - 0.05], [y, y], color='#BBBBBB', lw=0.4, ls=(0, (1, 2)), zorder=1)
        # supports (collected; placed after the tree with collision avoidance, focal nodes first)
        if id(d.clade) in S and d.clade.name and d.kind != 'tip':
            cand.append((S[id(d.clade)], x, y, d.clade.name))
    def min_y(d): return min([d.y] if d.kind != 'node' else [min_y(k) for k in d.kids]) if d.kind != 'tri' else d.y - d.rows / 2
    def max_y(d): return max([d.y] if d.kind != 'node' else [max_y(k) for k in d.kids]) if d.kind != 'tri' else d.y + d.rows / 2
    # root stem
    ax.plot([X0 - 0.25, PX(root.x)], [PY(root.y)] * 2, color=lc, lw=0.7)
    rec(root, X0 - 0.25)
    placed = []
    for pr, x, y, txt in sorted(cand, key=lambda c: c[0]):
        below = y + 0.2 > hp - 0.36
        bx = (x - 1.15, x - 0.02, y - 0.42, y + 0.22) if below else (x - 1.15, x - 0.02, y - 0.22, y + 0.42)
        if any(not (bx[1] < p[0] or bx[0] > p[1] or bx[3] < p[2] or bx[2] > p[3]) for p in placed): continue
        placed.append(bx)
        ax.text(x - 0.04, y - 0.015 if below else y + 0.015, txt, fontsize=6, ha='right', va='top' if below else 'bottom', color='#222222',
                bbox=dict(facecolor='white', edgecolor='none', pad=0.2, alpha=0.85), zorder=5)
    for ya, yb, d in brackets:
        s = d.group; col = COL[s][0]
        ax.plot([XB - 0.05, XB, XB, XB - 0.05], [ya + 0.09, ya + 0.09, yb - 0.09, yb - 0.09], color=col, lw=1.0, solid_capstyle='butt')
        lab = f'{FOCAL_LABEL_FAM.get((fam, s), FOCAL_LABEL[s])} ({d.n})'
        rich(ax, XB + 0.12, (ya + yb) / 2, lab, color=col, bold=True)
    # scale bar
    L_ = next(v for v in (0.02, 0.05, 0.1, 0.2, 0.5, 1.0) if v * sc >= 0.9 or v == 1.0)
    # choose the largest nice length <= 40% of the tree width
    cand = [v for v in (0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0) if v * sc <= 0.4 * TREE_W]
    L_ = cand[-1] if cand else 0.01
    ysb = hp - 0.2
    ax.plot([XB + 2.0, XB + 2.0 + L_ * sc], [ysb, ysb], color='k', lw=0.9, solid_capstyle='butt')
    ax.text(XB + 2.0 + L_ * sc + 0.1, ysb, f'{L_:g} subst./site', fontsize=FS, va='center', ha='left')

def render(specs, base, rowh, thr_list=(8,), tri_cap=4.0):
    global ROWH, TRI_CAP
    ROWH = rowh; TRI_CAP = tri_cap
    fig = plt.figure(figsize=(W * CM, 22.0 * CM)); fig.canvas.draw()
    for thr in thr_list:
        ps = [make_display(f, thr) for f, _, _ in specs]
        hts = [p[2] * ROWH + PAD for p in ps]
        total = sum(hts) + 0.1
        if total <= MAX_H: break
    print(base, 'tip threshold used', thr, 'heights', [round(h, 2) for h in hts], 'total', round(total, 2))
    fig.set_size_inches(W * CM, total * CM)
    cur = total
    for (fam, letter, title), (t, root, rows, S, X, maxx), hp in zip(specs, ps, hts):
        ax = fig.add_axes([0, (cur - hp) / total, 1, hp / total]); ax.set_xlim(0, W); ax.set_ylim(0, hp); ax.axis('off')
        draw_panel(ax, fam, root, rows, S, hp, letter, title)
        c = collections.Counter(); count_items(root, c)
        print('PANEL', fam, 'leaves', len(t.get_terminals()), 'open tips', c['tips'], 'focal triangles', c['tri_focal'], 'grey triangles', c['tri_grey'], 'rows', round(rows, 1))
        cur -= hp
    for ext, kw in (('png', dict(dpi=int(os.environ.get("FIG_DPI", "300")))), ('svg', {}), ('pdf', {})):
        fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), f'{base}.{ext}'), **kw)
    plt.close(fig)
    print(base, 'fig size cm', W, round(total, 2))

PAD = 0.45
def main():
    global TRI_CAP, ROWH
    specs = (('SUS', 'a', 'OCD-fold family (*Sus*)'), ('CUS', 'b', 'Cucumopine synthase family (*Cus*)'),
             ('CHS', 'c', 'M24 aminopeptidase family (*Chs*, *Ags*)'), ('ACS', 'd', 'GDPD family'),
             ('NOS', 'e', 'Opine dehydrogenase family, *nos* subfamily'), ('ODH', 'f', 'Opine dehydrogenase family, *odh* island'))
    for cap in (4.0, 3.0):
        for thr in (4, 3, 2, 1):
            TRI_CAP = cap
            rows = sum(make_display(f, thr)[2] for f, _, _ in specs)
            for rowh in (0.24, 0.23, 0.22):
                tot = rows * rowh + PAD * 6 + 0.1
                print('try cap', cap, 'thr', thr, 'rowh', rowh, 'rows', round(rows, 1), 'total', round(tot, 2))
                if tot <= MAX_H:
                    render(specs, 'Fig3_trees', rowh, (thr,), cap); return
    print('NO FIT'); render(specs, 'Fig3_trees', 0.22, (1,), 3.0)

if __name__ == '__main__':
    main()
