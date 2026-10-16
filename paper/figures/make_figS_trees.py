#!/usr/bin/env python3
"""Supplementary Figures: full maximum-likelihood trees of the six families drawn in Fig. 3 (all leaves, SH-aLRT/UFBoot
supports >= 50 at the nodes, leaf colours = lineage categories of Fig. 3, same rooting as Fig. 3).
Sources (same treefiles as make_fig3.py): trees_final/<F>/<F>.mafft.iqtree.treefile (F = SUS, CUS, ACS),
chs_fusarium/chs_tree.iqtree.treefile, fusarium_odh_island/nos_check/nos_tree.iqtree.treefile,
fusarium_odh_island/tree/odh_tree.iqtree.treefile.
Output: supplementary/figures/FigureS{2..7}_tree_{sus,cus,chs,gdpd,nos,odh}.{png,svg,pdf} (Supplementary Figures S2-S7)
Run: python make_figS_trees.py [family ...]"""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import sys, os, re, csv, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import make_fig3 as F3
from common import load, categories, grp

HERE = os.path.dirname(os.path.abspath(__file__))
OUTD = os.path.join(os.path.dirname(HERE), 'supplementary', 'figures')
A = F3.A
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
                     'pdf.fonttype': 42, 'svg.fonttype': 'none'})
CM = 1 / 2.54; W = 17.5
GREY, DGREY, NAT, OGC = '#7A7A7A', '#444444', '#8C6D2F', '#A5A5A5'
# colours: key lineage categories as in Fig. 3; other categories grey/brown
def colour(cat):
    if cat in F3.COL: return F3.COL[cat][0]
    return {'native': NAT, 'natfung': NAT, 'agsnat': NAT, 'OG': OGC, 'fungother': '#B8A57C', 'fungunres': '#B8A57C'}.get(cat, DGREY)
LEG = {  # legend text per category per family (order of appearance)
 'sus': 'T-DNA sus clade (fungal sus-like + Ti/Ri plasmid and plant cT-DNA references)', 'chromOCD': 'Rhizobiales chromosomal OCD',
 'bactOCD': 'other bacterial OCD-fold proteins', 'fungOCD': 'other fungal OCD-like', 'native': 'native fungal homologs', 'OG': 'outgroup',
 'cus': 'T-DNA cus clade (fungal cus-like + references)', 'actino': 'actinobacterial-type cus clade (fungal)', 'mis': 'mis group (fungal + references)',
 'bact': 'other bacterial homologs', 'fungother': 'native fungal gene, excluded',
 'tdna': 'T-DNA acs clade (references only)', 'acc': 'acc-type GDPD (fungal + acc-locus references)', 'kerr': 'pKerr14a-type acs-like (fungal + pKerr14a acs + Parasponia)',
 'natfung': 'native fungal GDPD (main clade)', 'fungunres': 'unresolved fungal GDPD',
 'fus': 'Fusarium chs-like', 'dac': 'Dactylonectria chs-like', 'tichs': 'Ti-plasmid chs clade', 'agsnat': 'native fungal Ags homologs',
 'nosbol': 'nos-subfamily ODH, Boletales', 'nosrhy': 'Rhypophila nos homolog (not assigned)', 'tdnanos': 'T-DNA nos clade', 'noslike': 'rhizobial nos-like (non-T-DNA)',
 'isl': 'odh-island ODH (Fusarium)', 'larry': 'A. larrymoorei plasmid ODH', 'nossub': 'nos subfamily (T-DNA type)'}
TITLE = {'SUS': 'OCD-fold family (Sus): full maximum-likelihood tree', 'CUS': 'Cucumopine synthase family (Cus): full maximum-likelihood tree',
         'ACS': 'GDPD family (Acs): full maximum-likelihood tree', 'CHS': 'M24 aminopeptidase family (Chs, Ags): full maximum-likelihood tree',
         'NOS': 'Opine dehydrogenase family, nos subfamily: full maximum-likelihood tree', 'ODH': 'Opine dehydrogenase family, odh island: full maximum-likelihood tree'}
OUTNAME = {'SUS': 'S2_tree_sus', 'CUS': 'S3_tree_cus', 'ACS': 'S5_tree_gdpd', 'CHS': 'S4_tree_chs', 'NOS': 'S6_tree_nos', 'ODH': 'S7_tree_odh'}  # numbering 2026-10-12q
ROW = {'SUS': 0.26, 'CUS': 0.135, 'ACS': 0.19, 'CHS': 0.34, 'NOS': 0.32, 'ODH': 0.34}  # cm per leaf
IQ = {'SUS': f'{A}/trees_final/SUS/SUS.mafft.iqtree.iqtree', 'CUS': f'{A}/trees_final/CUS/CUS.mafft.iqtree.iqtree', 'ACS': f'{A}/trees_final/ACS/ACS.mafft.iqtree.iqtree',
      'CHS': f'{A}/chs_fusarium/chs_tree.iqtree.iqtree', 'NOS': f'{A}/fusarium_odh_island/nos_check/nos_tree.iqtree.iqtree', 'ODH': f'{A}/fusarium_odh_island/tree/odh_tree.iqtree.iqtree'}

def model_of(fam):
    m = n = ''
    for l in open(IQ[fam]):
        if l.startswith('Model of substitution:'): m = l.split(':', 1)[1].strip()
        if l.startswith('Input data:'): n = l.split(':', 1)[1].strip()
    return m, n

# ---------- organism lookup for bare accessions (chs / nos / odh bacterial leaves) ----------
def org_map():
    M = {}
    p = f'{A}/chs_fusarium/bact_selected.tsv'
    for r in csv.DictReader(open(p), delimiter='\t'): M[r['acc']] = r['organism']
    for p in (f'{A}/fusarium_odh_island/nos_check/bact_selected.tsv', f'{A}/fusarium_odh_island/tree/bact_selected.tsv'):
        if os.path.exists(p):
            for r in csv.DictReader(open(p), delimiter='\t'):
                if 's' in r and 'org' in r: M[r['s']] = r['org']
    return M
ORG = org_map()
def sp(s): return s.replace('_', ' ').strip()

def label(fam, n, L):
    """readable leaf label"""
    if fam in ('SUS', 'CUS', 'ACS'):
        g = grp(L, n); parts = n.split('__')
        if g.startswith('FUNGI') or g.startswith('NATIVE'):
            if len(parts) >= 3:
                status = 'pseudogene-like' if 'pseudogene' in parts[0] else 'intact'
                tag = 'native ' if g.startswith('NATIVE') else ''
                return f'{tag}{sp(parts[2])}  {parts[1]}  [{status}]'
        tail = sp(parts[-1]) if len(parts) > 1 else sp(n)
        tail = re.sub(r'^(Agro|Plant) ', lambda m: 'plant cT-DNA: ' if m.group(1) == 'Plant' else 'plasmid: ', tail)
        tail = re.sub(r'(WP \d+) (\d)\b', r'\1.\2', tail)  # WP_x_1 -> WP_x.1
        tail = tail.replace('WP ', 'WP_')
        tag = {'OG': 'outgroup: '}.get(g, '')
        if g.startswith('REF_T2x'): tag = 'bacterial: '
        elif g.startswith('REF_T2'): tag = 'reference: '
        elif g.startswith('REF_T1'): tag = 'reference (T-DNA/cT-DNA): '
        return (tag + tail)[:70]
    if fam == 'CHS':
        if n.startswith('CHSlike_'):
            m = re.match(r'CHSlike_(.+?)_(GCA_\d+\.\d)_(\S+)$', n)
            return f'{sp(m.group(1))}  {m.group(2)}  {m.group(3)}' if m else sp(n[8:])
        if n.startswith('AGSnat_'):
            m = re.match(r'AGSnat_([A-Za-z]+)_([A-Za-z]*)_(\S+)$', n)
            return f'{m.group(1)} ({m.group(2)})  {m.group(3)}' if (m and m.group(2)) else (f'{m.group(1)}  {m.group(3)}' if m else sp(n[7:]))
        if n.startswith('BACT_'):
            acc = n[5:]; return f'{ORG.get(acc, "")}  {acc}'.strip()
        if n.startswith('SEED_ags'): return f'ags seed  {n[9:]}'
        if n.startswith('SEED_chs'): return f'chs seed (pTi14-818c)  {n[9:]}'
        return sp(n)
    # NOS / ODH: PREFIX_Organism_..._ACC
    m = re.match(r'^(ISL|FUNG|BOL|NOS_TDNA|NOSLIKE|NAT|BACT2|BACT|SEED)_(.+?)_([A-Z]{2,3}_?\d+(?:\.\d)?|[A-Z]\d[A-Z0-9]{3,}\d|P\d{5})$', n)
    if m:
        pre, org, acc = m.groups()
        tag = {'SEED': 'seed: ', 'NOS_TDNA': 'T-DNA nos: ', 'NOSLIKE': 'rhizobial nos-like: ', 'NAT': 'native: ', 'BACT': '', 'BACT2': '', 'ISL': '', 'FUNG': '', 'BOL': ''}[pre]
        org = sp(org)
        org = re.sub(r'\b(\d{9}\.\d)\b', r'GCA_\1', org)  # BOL assembly ids
        if pre == 'BACT' and org in ('Rhypophila',): pass
        return f'{tag}{org}  {acc}'
    return sp(n)

# ---------- tree layout ----------
def parse_sup(c):
    s = c.name if c.name else (str(c.confidence) if c.confidence is not None else '')
    if not s: return None, None, ''
    if '/' in s:
        a, b = s.split('/', 1)
        try: return float(a), float(b), s
        except ValueError: return None, None, s
    try: return None, float(s), s
    except ValueError: return None, None, s

def draw(fam):
    if fam in F3.GEN:
        t, L = F3.load_generic(fam); cat = F3.categories_generic(fam, t)
    else:
        t, L = load(fam); cat = categories(fam, t, L)
    t.ladderize(reverse=True)
    leaves = t.get_terminals(); n = len(leaves)
    row = ROW[fam]; fs = max(3.0, min(7.0, row * 23))  # pt
    Y = {}
    for i, l in enumerate(leaves): Y[l] = i
    for cl in t.get_nonterminals(order='postorder'): Y[cl] = sum(Y[c] for c in cl.clades) / len(cl.clades)
    X = t.depths()  # cumulative branch length
    maxd = max(X.values()) or 1.0
    model, nsites = model_of(fam)
    tw = 6.0 if fam != 'CUS' else 5.5  # tree width cm
    xscale = tw / maxd
    top, bottom = 0.9, 2.1
    H = top + n * row + bottom
    fig = plt.figure(figsize=(W * CM, H * CM))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis('off')
    x0 = 0.4
    def px(c): return x0 + X[c] * xscale
    def py(c): return H - top - (Y[c] + 0.5) * row
    lw = 0.45 if n > 200 else 0.6
    for cl in t.get_nonterminals():
        x, y = px(cl), py(cl)
        ys = [py(c) for c in cl.clades]
        ax.plot([x, x], [min(ys), max(ys)], color='#333333', lw=lw, solid_capstyle='butt')
        for c in cl.clades: ax.plot([x, px(c)], [py(c), py(c)], color='#333333', lw=lw, solid_capstyle='butt')
        if cl is not t.root:
            sh, uf, s = parse_sup(cl)
            if s and ((uf is not None and uf >= 50) or (sh is not None and sh >= 50 and uf is None)):
                ax.text(x - 0.03, y + row * 0.08, s, fontsize=fs * 0.8, ha='right', va='bottom', color='#222222')
    ax.plot([x0 - 0.25, x0], [py(t.root), py(t.root)], color='#333333', lw=lw)
    counts = collections.Counter()
    for l in leaves:
        c = cat[l.name]; counts[c] += 1
        ax.text(px(l) + 0.08, py(l), label(fam, l.name, L) + (' [native fungal gene, excluded]' if l.name == 'FUNGI_cus_intact__CUS_JBBPKX010002805_L1__Boletus_austroedulis' else ''), fontsize=fs, va='center', color=colour(c),
                fontweight='bold' if c in F3.KEY[fam] else 'normal')
    # title
    ax.text(0.3, H - 0.15, TITLE[fam], fontsize=8, fontweight='bold', va='top')
    ax.text(0.3, H - 0.5, f'{n} sequences; {nsites}; model {model}; node labels SH-aLRT (%) / UFBoot (%), shown when UFBoot ≥ 50; '
            f'rooted as in Fig. 3', fontsize=6, va='top', color='#333333')
    # scale bar
    sb = 0.5 if maxd > 2 else 0.2
    ax.plot([x0, x0 + sb * xscale], [1.55, 1.55], color='k', lw=0.8)
    ax.text(x0 + sb * xscale + 0.1, 1.55, f'{sb} substitutions/site', fontsize=6, ha='left', va='center')
    # legend
    order = [k for k in list(F3.KEY[fam]) + list(F3.NAMES[fam].keys()) + list(counts) if k in counts]
    seen = []; hs = []
    for k in order:
        if k in seen: continue
        seen.append(k)
        hs.append(Line2D([], [], marker='s', ls='', color=colour(k), ms=5, label=f'{LEG.get(k, k)} ({counts[k]})'))
    ax.legend(handles=hs, loc='lower left', bbox_to_anchor=(x0 / W + 0.0, 0.0), ncol=2 if len(hs) > 5 else 1, fontsize=6, frameon=False,
              handletextpad=0.4, labelspacing=0.25, columnspacing=1.2, borderaxespad=0.2)
    # (legend sits in the bottom margin below the scale bar; adjust bottom if it overlaps)
    os.makedirs(OUTD, exist_ok=True)
    out = os.path.join(OUTD, f'Figure{OUTNAME[fam]}')
    for ext, kw in (('png', dict(dpi=300 if n > 200 else 200)), ('svg', {}), ('pdf', {})):
        fig.savefig(f'{out}.{ext}', **kw)
    plt.close(fig)
    print(f'{fam}: {n} leaves, {dict(counts)}, size cm {W} x {H:.1f}, font {fs:.1f} pt, model {model}')

if __name__ == '__main__':
    fams = sys.argv[1:] or ['SUS', 'CUS', 'ACS', 'CHS', 'NOS', 'ODH']
    for f in fams: draw(f)
