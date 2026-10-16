#!/usr/bin/env python3
"""Build TableS1_loci_catalog.tsv, Table1_lineages.tsv, numbers_for_text.txt (+ lineage_check.tsv).
Run: python3 build_table1.py   (needs Biopython; no heavy compute)

Lineage rules (trees: analysis/trees_final/<FAM>/<FAM>.mafft.iqtree.treefile, rooted on the OG__ leaves):
 SUS: 'sus-like (T-DNA type)' = MRCA of all REF_sus leaves.  Other fungal leaves -> 'OCD-like, outside T-DNA sus clade'.
 CUS: 'cus-like (T-DNA type)' = MRCA(REF_cus); 'actinobacterial-type' = MRCA(FUNGI_cus_actinobacteria leaves)
      (pure fungal clade); 'mis lineage' = MRCA(FUNGI_mis leaves) (also contains the HMM-cusrel181 fungal leaves);
      remaining fungal leaves -> 'other bacterial-type CUS'.
 ACS: 'T-DNA acs clade' = MRCA(REF_acs without Parasponia PON40589) (expected: no fungi);
      'acc-type GDPD' = MRCA(REF_gdpd_meso); 'pKerr14a-type acs-like' = the six Boletales leaves
      (contig list from analysis/trees_final/au_lineage_acs.py) + 13 non-tree loci listed in analysis/funnel/carriers_tree_based.tsv;
      'native fungal GDPD' = all remaining FUNGI_acs leaves (94-leaf pure clade + 25 elsewhere).
 Duplicate loci (leaf_labels row whose label is the tree leaf of a representative) inherit the representative's clade.
"""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import csv, json, os, re, collections, sys
from Bio import Phylo

G = f'{GMO_ROOT}'
OUT = f'{G}/paper/tables'
TR = f'{G}/analysis/trees_final'
LT = f'{G}/analysis/loci_pipeline/round2/loci_table.tsv'
ADD = f'{G}/analysis/loci_pipeline/addenda/loci_added_20261006.tsv'
CLU = f'{G}/analysis/boletales_mash_517/clusters.tsv'
CAR = f'{G}/analysis/funnel/carriers_tree_based.tsv'
CSUM = f'{G}/analysis/funnel/carriers_summary.tsv'
REG = [f'{G}/analysis/loci_pipeline/registry/loci.tsv', f'{G}/analysis/loci_pipeline/registry/loci_step2_added.tsv']
BLAST = f'{G}/data/tsv/merged/all_blast_stats.tsv'

def rd(p): return list(csv.DictReader(open(p), delimiter='\t'))

NOTPLACED = {  # HMM subgroup -> name used when locus is not on the tree
 'acs': 'acs-HMM group', 'gdpd_meso': 'acc-type GDPD', 'kerr14': 'pKerr14a-type acs-like', 'sus': 'sus-like (T-DNA type)',
 'cus': 'cus-like (T-DNA type)', 'mis': 'mis lineage', 'cus_actinobacteria': 'actinobacterial-type', 'cusrel181': 'cus-related small ORF (cusrel181)',
 'rcorf12': 'rcorf12', 'cus_betaproteobacteria': 'other bacterial-type CUS (beta-proteobacteria)', 'cus_other_bacteria': 'other bacterial-type CUS',
 'ocd_other_bacteria': 'OCD-like bacterial', 'OUTGROUP': 'outgroup-HMM', '': 'no HMM subgroup'}
NOTHOM = 'not a homolog'
QUAR = {'CUS_JBBPKX010002805_L1': 'quarantined: native fungal gene'}  # analysis/curated/quarantine.tsv (2026-09-07)
ORDER_LIN = ['sus-like (T-DNA type)', 'cus-like (T-DNA type)', 'acc-type GDPD', 'pKerr14a-type acs-like']
MAIN = [('SUS','sus-like (T-DNA type)','T-DNA'),('CUS','cus-like (T-DNA type)','T-DNA'),('CUS','mis lineage','T-DNA, direction unresolved'),
 ('CUS','actinobacterial-type','bacterial, non-rhizobial'),('ACS','acc-type GDPD','acc-like catabolic locus'),('ACS','pKerr14a-type acs-like','Ti-type plasmid pKerr14a'),
 ('ACS','native fungal GDPD (main clade)','fungal'),('ACS','other fungal GDPD (position unresolved)','fungal')]

# ---------------- tree-based clade assignment ----------------
def load_tree(fam):
    rows = rd(f'{TR}/{fam}/{fam}_labels.tsv')
    t = Phylo.read(f'{TR}/{fam}/{fam}.mafft.iqtree.treefile', 'newick')
    L = {}
    for x in rows: L.setdefault(x['label'], x)
    og = [l.name for l in t.get_terminals() if l.name.startswith('OG__')]
    t.root_with_outgroup(og)
    return rows, t, L

def grp(L, n): return L[n]['group'] if n in L else 'OG'
def leaves(t, L, groups, excl=''):
    return [l.name for l in t.get_terminals() if grp(L, l.name) in groups and not (excl and excl in l.name)]
def clade(t, names): return {x.name for x in t.common_ancestor(names).get_terminals()}

def assign_tree(fam):
    """returns {leaf_label: (lineage, note)} for fungal leaves + info dict"""
    rows, t, L = load_tree(fam)
    fung = [l.name for l in t.get_terminals() if grp(L, l.name).startswith('FUNGI')]
    res, info = {}, {}
    if fam == 'SUS':
        c = clade(t, leaves(t, L, {'REF_sus'}))
        for n in fung: res[n] = ('sus-like (T-DNA type)', '') if n in c else ('OCD-like, outside T-DNA sus clade', 'tree: basal, sister only to outgroup leaves')
    elif fam == 'CUS':
        ct = clade(t, leaves(t, L, {'REF_cus'}))
        ac = clade(t, leaves(t, L, {'FUNGI_cus_actinobacteria'}))
        mi = clade(t, leaves(t, L, {'FUNGI_mis'}))
        for n in fung:
            if n in ct: res[n] = ('cus-like (T-DNA type)', '')
            elif n in ac: res[n] = ('actinobacterial-type', '')
            elif n in mi: res[n] = ('mis lineage', 'HMM cusrel181; inside mis region' if grp(L, n) == 'FUNGI_cusrel181' else '')
            else: res[n] = ('other bacterial-type CUS', 'tree: inside mixed bacterial-ref clade, outside T-DNA cus clade')
    else:
        SIX = ['JBBPKL010003317', 'JBBPKL010018021', 'JBBPNB010000203', 'JBBPNB010000771', 'JBBPOG010012756', 'JBBPRY010001778']
        six = [n for n in fung if any(s in n for s in SIX)]
        td = clade(t, leaves(t, L, {'REF_acs'}, excl='PON40589'))
        ac = clade(t, leaves(t, L, {'REF_gdpd_meso'}))
        sx = clade(t, six)
        info['six_leaves'] = len(six)
        nat = [n for n in t.get_nonterminals() if n.count_terminals() == 94 and all(grp(L, x.name) == 'FUNGI_acs' for x in n.get_terminals())]
        natc = {x.name for x in nat[0].get_terminals()} if nat else set()
        for n in fung:
            if n in td: res[n] = ('T-DNA acs clade', '')
            elif n in ac: res[n] = ('acc-type GDPD', '')
            elif n in sx: res[n] = ('pKerr14a-type acs-like', '')
            elif n in natc: res[n] = ('native fungal GDPD (main clade)', 'core 94-leaf pure fungal clade')
            else: res[n] = ('other fungal GDPD (position unresolved)', 'outside 94-leaf clade; not in T-DNA/acc/six clades')
    return rows, t, L, res

# ---------------- taxonomy ----------------
def taxmaps():
    byloc, bygen = {}, {}
    for p in REG:
        for r in rd(p):
            byloc[r['locus_id']] = (r.get('division', ''), r.get('family_tax', ''))
            g = r.get('genus') or r['organism'].split()[0]
            if r.get('division') or r.get('family_tax'): bygen.setdefault(g, (r.get('division', ''), r.get('family_tax', '')))
    bg2 = {}
    for r in csv.DictReader(open(BLAST), delimiter='\t'):
        if r.get('Genus'): bg2.setdefault(r['Genus'], (r.get('Division', ''), r.get('Family', '')))
    return byloc, bygen, bg2

def binomial(org):
    t = org.replace('_', ' ').split()
    if not t: return ''
    if t[0] in ('aff.', 'cf.') and len(t) > 2: return ' '.join(t[:3])
    if len(t) > 2 and t[1] in ('sp.', 'aff.', 'cf.'): return ' '.join(t[:3])
    return ' '.join(t[:2])

def main():
    lt = rd(LT); add = rd(ADD)
    lt = [r for r in lt if r['fam'] != 'ORF358']
    for r in add: r['_added'] = 1
    # carriers table
    car = {r['locus_id']: r['lineage'] for r in rd(CAR)}
    # clusters
    gca2sp, gca2cl = {}, {}
    cl_names = collections.defaultdict(set)
    for r in rd(CLU):
        cl_names[r['cluster_name']].add(r['cluster_id'])
        for g in r['genomes'].split(';'): gca2sp[g] = r['cluster_name']; gca2cl[g] = r['cluster_id']
    dupname = {n for n, c in cl_names.items() if len(c) > 1}  # same accepted name, different genomic clusters
    byloc, bygen, bg2 = taxmaps()
    # tree assignment
    pl = {}  # locus_id -> (lineage, note)
    for fam in ('ACS', 'CUS', 'SUS'):
        rows, t, L, res = assign_tree(fam)
        for x in rows:
            if x['group'].startswith('FUNGI') and x['label'] in res:
                lin, note = res[x['label']]
                if x['duplicate_of'] and x['label'] != x['locus'] and not x['label'].split('__')[1] == x['locus']:
                    note = (note + '; ' if note else '') + 'duplicate of tree leaf'
                pl[x['locus']] = (lin, note, x['label'])
    out = []
    for r in lt + add:
        fam, lid, st = r['fam'], r['locus_id'], r['status']
        note = ''; leaf = ''
        if lid in QUAR: lin, basis, note = NOTHOM, 'hmm_only', QUAR[lid]
        elif r.get('_added'): lin, basis, note = 'cus-like (T-DNA type)', 'similarity', 'addendum 20261006; not on tree (HMM cus + identical/similar to tree FUNGI_cus)'
        elif st in ('not_homolog', 'not_run'): lin, basis = NOTHOM, 'hmm_only'
        elif lid in pl: (lin, note, leaf), basis = pl[lid], 'tree'
        elif fam == 'ACS' and car.get(lid, '').startswith('ACS kerr14'): lin, basis, note = 'pKerr14a-type acs-like', 'similarity', '>=60% identity to the six tree-placed members'
        else: lin, basis = NOTPLACED.get(r['subgroup'], r['subgroup']) + ' (not placed)', 'hmm_only'
        gca = r['gca']; org = r['organism']
        sp = gca2sp.get(gca) if r['order'] == 'Boletales' else None
        scl = gca2cl.get(gca, '') if sp else ''
        sp_src = 'mash_cluster' if sp else 'binomial'
        sp = sp or binomial(org)
        if sp == 'Kazachstania servazzii': sp = 'Monosporozyma servazzii'; note = (note + '; ' if note else '') + 'syn. Kazachstania servazzii'
        skey = scl if scl else sp
        genus = org.replace('_', ' ').split()[0]
        div, fa = byloc.get(lid, ('', ''))
        if not div: div = bygen.get(genus, ('', ''))[0] or bg2.get(genus, ('', ''))[0]
        if not fa: fa = bygen.get(genus, ('', ''))[1] or bg2.get(genus, ('', ''))[1]
        if scl and sp in dupname: sp = f'{sp} (cluster {scl})'
        out.append(dict(species_key=skey, species_cluster=scl, locus_id=lid, family=fam, lineage_final=lin, basis=basis, hmm_subgroup=r['subgroup'], status=st, organism=org,
            species=sp, division=div, **{'class': r['class']}, order=r['order'], family_tax=fa, gca=gca, contig=r['contig'], contig_len=r['contig_len'],
            start=r['contig_start'] or r['hsp_start'] if 'hsp_start' in r else r['contig_start'], end=r['contig_end'] or r.get('hsp_end', ''), strand=r['strand'],
            n_exons=r['n_exons'], protein_len=r['protein_len'], context_call=r['context_call'] if 'context_call' in r else '',
            identity_T1T2=r.get('identity_T1T2', ''), was_deleted=r['was_deleted'] if 'was_deleted' in r else '', note=note, species_source=sp_src, leaf=leaf))
    cols = ['locus_id', 'family', 'lineage_final', 'basis', 'hmm_subgroup', 'status', 'organism', 'species', 'division', 'class', 'order', 'family_tax',
            'gca', 'species_cluster', 'contig', 'contig_len', 'start', 'end', 'strand', 'n_exons', 'protein_len', 'context_call', 'identity_T1T2', 'was_deleted', 'note', 'species_source']
    with open(f'{OUT}/TableS1_loci_catalog.tsv', 'w', newline='') as f:
        w = csv.DictWriter(f, cols, delimiter='\t', extrasaction='ignore'); w.writeheader(); w.writerows(out)
    summarize(out)
    main_table(out, cl_names)

def top3(rs):
    c = collections.Counter(r['order'] or 'NA' for r in rs)
    return '; '.join(f'{k} ({v})' for k, v in sorted(c.items(), key=lambda kv: (-kv[1], kv[0]))[:3])

def stats(rs):
    s = collections.Counter(r['status'] for r in rs)
    bol = [r for r in rs if r['order'] == 'Boletales']
    return dict(n_loci=len(rs), intact=s['intact'], pseudogene_like=s['pseudogene_like'], fragment_remnant=s['fragment'] + s['remnant'],
        ambiguous=s['ambiguous'], other_status=len(rs) - s['intact'] - s['pseudogene_like'] - s['fragment'] - s['remnant'] - s['ambiguous'],
        n_genomes=len({r['gca'] for r in rs}), n_species=len({r['species_key'] for r in rs}),
        n_divisions=len({r['division'] for r in rs if r['division']}), n_classes=len({r['class'] for r in rs if r['class']}),
        n_orders=len({r['order'] for r in rs if r['order']}), n_families=len({r['family_tax'] for r in rs if r['family_tax']}),
        top3_orders=top3(rs), loci_Boletales=len(bol), species_Boletales=len({r['species_key'] for r in bol}))

def summarize(out):
    cols = ['block', 'lineage_final', 'family', 'n_loci', 'intact', 'pseudogene_like', 'fragment_remnant', 'ambiguous', 'other_status', 'n_genomes', 'n_species',
            'n_divisions', 'n_classes', 'n_orders', 'n_families', 'top3_orders', 'loci_Boletales', 'species_Boletales']
    lines = []
    for block, sel in (('all_loci', lambda r: True), ('tree_placed_only', lambda r: r['basis'] == 'tree')):
        rs0 = [r for r in out if sel(r)]
        lins = sorted({(r['lineage_final'], r['family']) for r in rs0}, key=lambda k: (k[1], k[0]))
        for lin, fam in lins:
            lines.append(dict(block=block, lineage_final=lin, family=fam, **stats([r for r in rs0 if r['lineage_final'] == lin and r['family'] == fam])))
        lines.append(dict(block=block, lineage_final='TOTAL', family='ACS+CUS+SUS', **stats(rs0)))
    with open(f'{OUT}/Table1_lineages.tsv', 'w', newline='') as f:
        w = csv.DictWriter(f, cols, delimiter='\t'); w.writeheader(); w.writerows(lines)
    # ---- checks and numbers
    chk = []
    allrows = [l for l in lines if l['block'] == 'all_loci' and l['lineage_final'] != 'TOTAL']
    tot = [l for l in lines if l['block'] == 'all_loci' and l['lineage_final'] == 'TOTAL'][0]
    chk.append(f"sum over lineages = {sum(l['n_loci'] for l in allrows)}; catalog rows = {len(out)}; TOTAL row = {tot['n_loci']}")
    car = rd(CAR); csum = {r['lineage']: r for r in rd(CSUM)}
    m = {'SUS T-DNA clade': 'sus-like (T-DNA type)', 'CUS cus clade': 'cus-like (T-DNA type)', 'ACS acc-locus (gdpd_meso) clade': 'acc-type GDPD',
         'ACS kerr14/acc-like ("six") lineage': 'pKerr14a-type acs-like'}
    for old, new in m.items():
        mine = {r['locus_id'] for r in out if r['lineage_final'] == new}
        theirs = {r['locus_id'] for r in car if r['lineage'] == old}
        s = csum[old]
        mg = {r['gca'] for r in out if r['lineage_final'] == new}
        chk.append(f"{new}: mine={len(mine)} (genomes {len(mg)}) vs carriers={len(theirs)} (summary n_loci {s['n_loci']}, genomes {s['n_genomes']}); only_mine={sorted(mine-theirs)}; only_carriers={sorted(theirs-mine)}")
    open(f'{OUT}/selfcheck.txt', 'w').write('\n'.join(chk) + '\n')
    # numbers_for_text
    N = []
    hom = [r for r in out if r['lineage_final'] != NOTHOM]
    N.append(f"Loci in catalog without ORF358: {len(out)} (round2 ACS+CUS+SUS = {len(out)-2}, + 2 addenda); homologs (not not_homolog/not_run): {len(hom)}; not a homolog: {len(out)-len(hom)}")
    N.append(f"Homologs placed on tree (basis=tree): {sum(1 for r in out if r['basis']=='tree')}; hmm_only homologs: {sum(1 for r in hom if r['basis']=='hmm_only')}; similarity: {sum(1 for r in out if r['basis']=='similarity')}")
    N.append('\nPer lineage (all loci): loci / intact / genomes / species')
    for l in allrows: N.append(f"  {l['family']} | {l['lineage_final']}: {l['n_loci']} / {l['intact']} / {l['n_genomes']} / {l['n_species']}")
    N.append('\nFour paper lineages, tree-placed + listed non-tree (as in funnel/carriers_tree_based.tsv): outside Boletales')
    for lin in ORDER_LIN:
        rs = [r for r in out if r['lineage_final'] == lin]
        nb = [r for r in rs if r['order'] != 'Boletales']
        N.append(f"  {lin}: {len(rs)} loci, {sum(1 for r in rs if r['status']=='intact')} intact, {len({r['gca'] for r in rs})} genomes, {len({r['species_key'] for r in rs})} species; "
                 f"outside Boletales: {len(nb)} loci, {len({r['species_key'] for r in nb})} species")
        c = collections.Counter((r['species'], r['order'], r['status']) for r in nb)
        for (sp, o, s), v in sorted(c.items()): N.append(f"      {sp} [{o}] {s} x{v}")
    nofam = sum(1 for r in out if not r['family_tax']); nodiv = sum(1 for r in out if not r['division'])
    N.append(f"\nTaxonomy: loci without family_tax: {nofam}; without division: {nodiv}")
    N.append(f"Species via mash clusters (Boletales): {sum(1 for r in out if r['species_source']=='mash_cluster')}; Boletales loci without cluster: {sum(1 for r in out if r['order']=='Boletales' and r['species_source']!='mash_cluster')}")
    open(f'{OUT}/numbers_for_text.txt', 'w').write('\n'.join(N) + '\n')

def main_table(out, cl_names_g):
    cols = ['family','lineage','origin','n_loci','n_nonredundant_tree_leaves','n_intact','n_pseudogene','n_fragments','n_genomes','n_species','n_orders','orders','n_species_Boletales','additional_fragments_profile_only']
    rows = []
    for fam, lin, origin in MAIN:
        rs = [r for r in out if r['lineage_final'] == lin and r['family'] == fam and r['basis'] in ('tree', 'similarity')]
        oc = collections.defaultdict(set)
        for r in rs: oc[r['order'] or 'NA'].add(r['species_key'])
        ords = '; '.join(f'{k} ({len(v)})' for k, v in sorted(oc.items(), key=lambda kv: (-len(kv[1]), kv[0])))
        extra = sum(1 for r in out if r['lineage_final'] == lin + ' (not placed)' and r['family'] == fam)
        rows.append(dict(family=fam, lineage=lin, origin=origin, n_loci=len(rs), n_nonredundant_tree_leaves=len({r['leaf'] for r in rs if r['leaf']}),
            n_intact=sum(1 for r in rs if r['status'] == 'intact'), n_pseudogene=sum(1 for r in rs if r['status'] == 'pseudogene_like'),
            n_fragments=sum(1 for r in rs if r['status'] in ('fragment', 'remnant')),
            n_genomes=len({r['gca'] for r in rs}), n_species=len({r['species_key'] for r in rs}), n_orders=len(oc), orders=ords,
            n_species_Boletales=len({r['species_key'] for r in rs if r['order'] == 'Boletales'}), additional_fragments_profile_only=extra))
    with open(f'{OUT}/Table1_main.tsv', 'w', newline='') as f:
        w = csv.DictWriter(f, cols, delimiter='\t'); w.writeheader(); w.writerows(rows)
    with open(f'{OUT}/selfcheck.txt', 'a') as f:
        for r in rows:
            ok = r['n_intact'] + r['n_pseudogene'] + r['n_fragments'] == r['n_loci']
            f.write(f"Table1_main {r['lineage']}: intact {r['n_intact']} + pseudogene {r['n_pseudogene']} + fragments {r['n_fragments']} = {r['n_intact']+r['n_pseudogene']+r['n_fragments']} vs n_loci {r['n_loci']}: {'OK' if ok else 'MISMATCH'}\n")
        dn = {n: sorted(c) for n, c in cl_names_g.items() if len(c) > 1}
        f.write(f"species key = mash cluster id (Boletales) / species name (others); cluster names shared by several clusters: {dn}\n")
    hom = [r for r in out if r['lineage_final'] != NOTHOM]
    hp = sum(1 for r in hom if r['basis'] == 'hmm_only')
    open(f'{OUT}/Table1_main_footnotes.txt', 'w').write(
        f"Total loci in catalog (ACS+CUS+SUS, ORF358 excluded): {len(out)}\nHomologs: {len(hom)}\nNot a homolog (incl. 1 quarantined): {len(out)-len(hom)}\n"
        f"Homologs assigned by profile (HMM) only, not placed on tree: {hp}\nHomologs placed on tree: {sum(1 for r in out if r['basis']=='tree')}; by similarity: {sum(1 for r in out if r['basis']=='similarity')}\n"
        "species counted by mash cluster id (Boletales) or name (others); n_fragments = fragment+remnant; n_loci = basis tree + similarity (duplicates of tree leaves included); n_nonredundant_tree_leaves = distinct tree leaves.\n")

if __name__ == '__main__':
    main()
