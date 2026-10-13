#!/usr/bin/env python3
"""Figure 3: schematic (collapsed) trees of the six families (sus, cus, chs, GDPD, nos subfamily, odh island).
The AU-test panel was moved to make_figS_AU.py (Supplementary Figure; data from Table S5) on 2026-10-12 (decision 6.3b).
Run: python make_fig3.py
Clade definitions are repeated from paper/tables/build_table1.py (rooted on OG__ leaves)."""
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
PINK, LPINK = '#ff7f0e', '#FFD0A8'       # nos-ODH Boletales (orange, as in Figures 4 and 5)
BROWN, LBROWN = '#A0522D', '#DDB89F'     # nos-ODH Rhypophila
INDIGO, LINDIGO = '#332288', '#B8B0DC'   # odh island
TEAL, LTEAL = '#1B7F8C', '#A8D5DB'       # A. larrymoorei plasmid ODH (source)
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
 'SUS': {'chromOCD': 'Rhizobiales chromosomal OCD', 'bactOCD': 'other bacterial OCD', 'fungOCD': 'native fungal OCD',
         'native': 'native fungal OCD', 'OG': 'outgroup'},
 'CUS': {'bact': 'other bacterial Cus-like', 'fungother': 'native fungal Cus-like', 'OG': 'outgroup'},
 'ACS': {'bact': 'other bacterial GDPD', 'fungunres': 'native fungal GDPD', 'native': 'native fungal GDPD', 'OG': 'outgroup'},
 'CHS': {'bact': 'other bacterial M24 aminopeptidases', 'agsnat': 'native fungal Ags-like aminopeptidases', 'OG': 'outgroup: T-DNA Ags'},
 'NOS': {'noslike': 'rhizobial nos-like, outside the T-DNA', 'bact': 'other bacterial ODH', 'native': 'native fungal ODH', 'OG': 'outgroup: T-DNA Ocs and Vis'},
 'ODH': {'bact': 'other bacterial ODH', 'native': 'native fungal ODH', 'OG': 'outgroup: T-DNA Ocs and Vis'},
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
            return ['Actinobacterial-type *cus* clade', f'{nf} fungal sequences; their closest relatives in nr are actinobacterial and are not in this tree',
                    'contains the *Beauveria bassiana* gene *Fcs1* (Li et al.)']
        if k == 'mis':
            nm = g['FUNGI_mis'] + g['FUNGI_cusrel181']
            return ['*mis* group', f'{nm} fungal + {nref} reference sequences (Ti/Ri plasmid, plant cT-DNA, bacterial)']
    if fam == 'ACS':
        if k == 'tdna': return ['T-DNA *acs* clade', f'{g["REF_acs"]} Ti/Ri plasmid + plant cT-DNA references + {g["BACT_acc"]} other bacterial; no fungal']
        if k == 'acc': return ['Acc-type GDPD', f'{nf} fungal (Boletales) + {nref} GDPD of rhizobial acc-like loci']
        if k == 'kerr':
            return ['pKerr14a-type *acs*-like', f'{nf} fungal (Boletales) + pKerr14a *acs* + Parasponia cT-DNA fragment']
        if k == 'natfung': return ['Native fungal GDPD (main clade)', f'{nf} fungal sequences only (no bacterial or plant references)']

    if fam == 'CHS':
        if k == 'fus': return ['*Fusarium* *chs*-like', f'{nf} sequences ({len({n.split("_")[2] for n in names})} species)']
        if k == 'dac': return ['*Dactylonectria estremocensis* *chs*-like', '1 sequence (KAH7120428.1)']
        if k == 'tichs': return ['Ti-plasmid *chs* clade', 'chrysopine synthases of Ti plasmids, incl. pTi14-818c']
    if fam == 'NOS':
        if k == 'nosbol': return ['*nos*-subfamily ODH, Boletales', f'{nf} loci (*Boletus*, *Butyriboletus*)']
        if k == 'nosrhy': return ['*Rhypophila* *nos*-like genes (origin not resolved)', f'{nf} sequences (2 species)']
        if k == 'tdnanos':
            nt = sum(n.startswith('NOS_TDNA') for n in names); nl = len(names) - nt
            return ['T-DNA *nos* clade', f'{nt} Ti-plasmid *nos* references + {nl} rhizobial *nos*-like']
        if k == 'isl': return ['*odh*-island ODH (*F. poae*) + *A. larrymoorei* plasmid ODH', 'not in the *nos* subfamily']
    if fam == 'ODH':
        if k == 'isl': return ['*odh*-island ODH (*Fusarium*)', f'{nf} sequences (*F. poae* ×2, *F. langsethiae*)']
        if k == 'larry': return ['*A. larrymoorei* plasmid ODH', 'WP_309770759 (source locus)']
        if k == 'nossub':
            return ['*nos* subfamily', 'T-DNA Nos (2), *Rhypophila* (2), *Butyriboletus* (1), *A. larrymoorei* *nos*-type gene']

def mixed_label(fam, names, cat, L):
    c = collections.Counter(cat[n] for n in names)
    nm = NAMES[fam]
    agg = {}   # merge categories that share a display name (e.g. screen-derived and reference native fungal genes)
    for k, v in c.most_common():
        agg[nm.get(k, k)] = agg.get(nm.get(k, k), 0) + v
    items = [f'{name} ({v})' for name, v in sorted(agg.items(), key=lambda x: -x[1])]
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

# ---------- build display tree ----------
class N: pass
def build(fam):
    if fam in GEN:
        t, L = load_generic(fam); cat = categories_generic(fam, t)
    else:
        t, L = load(fam); cat = categories(fam, t, L)
        if fam == 'CUS':   # author decision 2026-10-13: the excluded native fungal leaf is omitted from the schematic (kept in Fig. S3)
            t.prune('FUNGI_cus_intact__CUS_JBBPKX010002805_L1__Boletus_austroedulis')
    KEYS = KEY[fam]
    def haskey(node): return any(cat[x.name] in KEYS for x in node.get_terminals())
    def rec(node, depth_sup):
        # collapsing rule (2026-10-13): a collapsed block is exactly one monophyletic clade with a single category;
        # monophyletic clade without focal leaves = one block (may mix non-focal categories); no merging of siblings / grades
        names = [x.name for x in node.get_terminals()]
        n = N(); n.sup = sup(node) if not node.is_terminal() else ''; n.names = names; n.kids = []
        cs = {cat[x] for x in names}
        if len(cs) == 1 and next(iter(cs)) in KEYS:
            c0 = next(iter(cs)); n.kind = 'key'; n.cat = c0; n.label = key_label(fam, c0, names, L); return n
        if not (cs & KEYS):   # maximal monophyletic clade without any focal leaf -> one grey block labelled by composition
            c = collections.Counter(cat[x] for x in names)
            n.kind = 'mixed'; n.dom = c.most_common(1)[0][0]; n.cat = 'OG' if set(c) == {'OG'} else 'mixed'
            n.label = mixed_label(fam, names, cat, L); return n
        n.kind = 'node'; n.kids = [rec(c, depth_sup) for c in node.clades]
        return n
    root = rec(t.root, 0)
    return root, t, L, cat

def height(n): return 0 if not n.kids else 1 + max(height(k) for k in n.kids)
def tri_h(n):
    k = len(n.names); return 0.13 + 0.22 * math.log10(max(k, 1))
def nlines(n): return len(n.label)  # approx, wrapping handled in label()

def wrap_label(n, maxc):
    """return list of lines for display: title (bold) + wrapped details"""
    import textwrap
    if n.kind == 'key':
        out = [('b', n.label[0])]
        for s in n.label[1:]:
            for ln in textwrap.wrap(s, maxc): out.append(('n', ln))
        return out
    s = '; '.join(n.label)
    return [('n', ln) for ln in textwrap.wrap(s, maxc)]

LH = 0.215  # line height cm
def merge_ladders(root, fam, cat, L):
    """Consecutive (vertical order) mixed tips with the same category set that sit on one ladder
    (the MRCA of two neighbours has one of them as a direct child) are drawn as one block."""
    tips = []; parent = {}
    def collect(n):
        if not n.kids: tips.append(n)
        for k in n.kids: parent[id(k)] = n; collect(k)
    collect(root)
    def base(t):
        if t.kind != 'mixed': return None
        b = frozenset(cat[x] for x in t.names)
        return b if len(b) == 1 else None   # grades only for single-category non-focal blocks (not outgroup)
    def anc(t):
        a = []
        while id(t) in parent: t = parent[id(t)]; a.append(t)
        return a
    def ladder(t1, t2):
        a1 = anc(t1); a2 = set(map(id, anc(t2)))
        m = next(a for a in a1 if id(a) in a2)
        return t1 in m.kids or t2 in m.kids
    groups = []; cur = [tips[0]]
    for t in tips[1:]:
        if base(t) is not None and base(t) == base(cur[-1]) and ladder(cur[-1], t): cur.append(t)
        else: groups.append(cur); cur = [t]
    groups.append(cur)
    for g in groups:
        g[0].grp = g
        for t in g[1:]: t.skip = True
        if len(g) > 1:
            names = [x for t in g for x in t.names]; b = len(g)
            c = collections.Counter(cat[x] for x in names); nm = NAMES[fam]
            k0 = next(iter(c)); items = [f'{nm.get(k0, k0)} ({len(names)} in {b} branches)']
            g[0].label = items
            MERGED.append((fam, '; '.join(items), [len(t.names) for t in g]))
    return groups

def layout(root, maxc, groups):
    y = 0.0
    for g in groups:
        for t in g: t.h = tri_h(t)
        n = g[0]; n.lines = wrap_label(n, maxc)
        STUB = 0.5 * LH
        rh = (max(len(g) * STUB, len(n.lines) * LH) if len(g) > 1 else max(sum(t.h for t in g), len(n.lines) * LH)) + 0.025
        if len(g) > 1:
            y0 = y + (rh - len(g) * STUB) / 2
            for i, t in enumerate(g):
                t.y = y0 + STUB * (i + 0.5); t.h = STUB
        else: n.y = y + rh / 2
        y += rh
    def setx(n):
        if n.kids:
            for k in n.kids: setx(k)
            n.y = (n.kids[0].y + n.kids[-1].y) / 2
    setx(root)
    return y

def draw_tree(ax, root, fam, ytop, total, x0=1.3, xw=3.6, tw=0.8, xl=5.9, nested=None):
    D = max(height(root), 1); step = xw / D
    def X(n): return x0 + (D - height(n)) * step
    def Y(n): return ytop - n.y
    xtip = x0 + D * step
    lc = '#333333'
    def rec(n, parent_x):
        x, y = X(n), Y(n)
        if n.kids:
            ys = [Y(k) for k in n.kids]
            ax.plot([x, x], [min(ys), max(ys)], color=lc, lw=0.8, solid_capstyle='butt')
            for k in n.kids:
                ax.plot([x, X(k)], [Y(k), Y(k)], color=lc, lw=0.8, solid_capstyle='butt'); rec(k, x)
        else:
            col = (COL[n.cat][1] if n.kind == 'key' else (LG if getattr(n, 'dom', '') in FUNG_DOM else MG))
            edge = COL[n.cat][0] if n.kind == 'key' else '#808080'
            if n.kind == 'mixed' and n.cat == 'OG': col = 'white'
            h = n.h
            if getattr(n, 'skip', False): return
            g = getattr(n, 'grp', [n])
            if len(g) > 1:  # paraphyletic grade: short stubs + bracket on the label side, one label
                for t in g:
                    ax.plot([x, x + 0.35 * tw], [Y(t), Y(t)], color=lc, lw=0.8, solid_capstyle='butt')
                    ax.plot([x + 0.35 * tw], [Y(t)], marker='o', ms=1.8, color=DG, mew=0, zorder=3)
                xb = xl - 0.12; y1, y2 = Y(g[0]) + 0.05, Y(g[-1]) - 0.05
                ax.plot([xb - 0.05, xb, xb, xb - 0.05], [y1, y1, y2, y2], color='#555555', lw=0.7, solid_capstyle='butt')
                y = (Y(g[0]) + Y(g[-1])) / 2
            else:
                ax.add_patch(Polygon([[x, y], [x + tw, y + h / 2], [x + tw, y - h / 2]], closed=True, fc=col, ec=edge, lw=0.7))
            if nested and n.kind == 'key' and n.cat == nested[0]:
                f = nested[1]; hh = h * f
                ax.add_patch(Polygon([[x, y], [x + tw * f, y + hh / 2], [x + tw * f, y - hh / 2]], closed=True, fc=COL[n.cat][0], ec='none', alpha=0.9))
            lx = xl
            yy = y + (len(n.lines) - 1) * LH / 2
            for st, ln in n.lines:
                rich(ax, lx, yy, ln, bold=(st == 'b'), color=(COL[n.cat][0] if (n.kind == 'key' and st == 'b') else 'k'))
                yy -= LH
    rec(root, x0 - 0.3)
    ax.plot([x0 - 0.3, x0], [Y(root), Y(root)], color=lc, lw=0.8)
    # supports: parents of key clades and key stems
    def sups(n):
        for k in n.kids:
            if k.kind == 'key' and k.sup:
                last = (k is n.kids[-1])
                ax.text(X(k) - 0.05, Y(k) + (-0.04 if last else 0.04), k.sup, fontsize=4.9, ha='right', va='top' if last else 'bottom', color=COL[k.cat][0], fontweight='bold', bbox=dict(facecolor='white', edgecolor='none', pad=0.3), zorder=5)
            sups(k)
        if n is not root and n.kids and any((k.kind == 'key') for k in n.kids) and n.sup and not any(len(getattr(k, 'grp', [k])) > 1 for k in n.kids):
            ax.text(X(n) - 0.06, Y(n) + 0.03, n.sup, fontsize=4.9, ha='right', va='bottom', color='#222222', bbox=dict(facecolor='white', edgecolor='none', pad=0.3), zorder=5)
    sups(root)

def panel(fam, letter, title, fig, ytop_frac, maxc=62):
    pass

# AU-test data/panel: see make_figS_AU.py (reads supplementary/TableS5_*.tsv)

# layout knobs (Fig. 3 must fit one page: width 17.5 cm, height <= 24 cm)
MAXC_WIDE = 90    # characters per label line for the 4 narrow-tree panels
MAXC_MAIN = 100   # for sus/cus panels
PANEL_PAD = 0.38

def main():
    fig = plt.figure(figsize=(W * CM, 22.0 * CM))
    fig.canvas.draw()
    panels = []
    for fam, letter, title, nested in (('SUS', 'a', 'OCD-fold family (*Sus*)', None), ('CUS', 'b', 'Cucumopine synthase family (*Cus*)', None),
                                       ('CHS', 'c', 'M24 aminopeptidase family (*Chs*, *Ags*)', None), ('ACS', 'd', 'GDPD family', None),
                                       ('NOS', 'e', 'Opine dehydrogenase family, *nos* subfamily', None),
                                       ('ODH', 'f', 'Opine dehydrogenase family, *odh* island', None)):
        root, t, L, cat = build(fam)
        tips_ = []
        def _c(n):
            if not n.kids: tips_.append(n)
            for k in n.kids: _c(k)
        _c(root)
        groups = merge_ladders(root, fam, cat, L)   # grades only for same-category non-focal ladders
        tot = layout(root, MAXC_WIDE if fam in ('ACS', 'CHS', 'NOS', 'ODH') else MAXC_MAIN, groups)
        panels.append((fam, letter, title, nested, root, tot, len(groups)))
        INFO[(fam, 'height_cm')] = tot
        print('PANEL', fam, 'rows(blocks)', len(groups), 'height', round(tot, 2))
    plt.close(fig)
    hts = [p[5] + PANEL_PAD for p in panels]
    total = sum(hts) + 0.15
    print('single-figure height', round(total, 2), 'cm')
    parts = [('Fig3_trees', panels)]
    for base, ps in parts:
        tt = sum(p[5] + PANEL_PAD for p in ps) + 0.15 + 0.4
        fig = plt.figure(figsize=(W * CM, tt * CM)); fig.canvas.draw()
        cur = tt
        for fam, letter, title, nested, root, tot, _n in ps:
            hp = tot + PANEL_PAD
            ax = fig.add_axes([0, (cur - hp) / tt, 1, hp / tt]); ax.set_xlim(0, W); ax.set_ylim(0, hp); ax.axis('off')
            ax.text(0.05, hp - 0.1, letter, fontsize=9, fontweight='bold', va='top')
            rich(ax, 0.75, hp - 0.1, title, fs=FS + 1, bold=True, va='top')
            draw_tree(ax, root, fam, hp - PANEL_PAD, tot, nested=nested, **({'xw': 4.6, 'xl': 7.0} if fam in ('ACS', 'CHS', 'NOS', 'ODH') else {'xw': 3.4, 'xl': 5.95}))
            cur -= hp
        fig.text(0.0035, 0.21 / tt, 'Triangle: one monophyletic clade. Bracketed stubs: consecutive branches of one category that do not form a clade (a paraphyletic grade).', fontsize=FS, va='center', ha='left')
        for ext, kw in (('png', dict(dpi=int(os.environ.get("FIG_DPI", "300")))), ('svg', {}), ('pdf', {})):
            fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), f'{base}.{ext}'), **kw)
        plt.close(fig)
        print(base, 'fig size cm', W, round(tt, 2))

if __name__ == '__main__':
    main()
