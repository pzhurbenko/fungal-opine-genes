"""Summary of step-2 screen over 517 Boletales assemblies for ALL T-DNA / opine / backbone subjects. Reads existing out/<GCA>/ only;
does not touch shared step2_*.tsv. Chaining/locus rules replicate aggregate_step2.py / summarize_ba.py (GAP 5000, overlapping chains merge into loci,
locus best subject = max bits). strong: pident>=45 & subject cov>=0.5; background: pident>=30 & cov>=0.5."""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import re, os, csv, statistics as st
from collections import defaultdict, Counter
D = f'{GMO_ROOT}/analysis/refset/screen/step2'; OUT = f'{D}/boletales517_tdna_summary'
GEN = f'{GMO_ROOT}/analysis/boletales_mash_517/genomes.tsv'
S1 = f'{GMO_ROOT}/paper/tables/TableS1_loci_catalog.tsv'
QF = f'{GMO_ROOT}/analysis/refset/step2/queries.faa'
CL = f'{GMO_ROOT}/analysis/refset/plasmids/tdna_clusters.tsv'
GAP = 5000; BACK = ('vir', 'tra', 'rep_partition', 'mobile'); CONTAM_PID = 60.0
CARR = {'sus-like (T-DNA type)', 'cus-like (T-DNA type)', 'acc-type GDPD', 'pKerr14a-type acs-like'}
def sv(c): return re.sub(r'\.\d+$', '', c)
def num(g): return re.sub(r'^GC[AF]_', '', g).split('.')[0]

def ps(s):
    p = s.split('|')
    if p[0].startswith('FAM_'): return dict(cat='FAM_' + p[0][4:], lin=p[1], acc=p[2] if len(p) > 2 else '', plasmid='', gene='', prod='')
    return dict(cat=p[1], lin=p[4] if len(p) > 4 and p[4] != '-' else '', acc=p[0], plasmid=p[2] if len(p) > 2 else '',
                gene=p[4] if len(p) > 4 and p[4] != '-' else '', prod=p[5] if len(p) > 5 else '')
def gkey(s):
    p = ps(s)
    if p['cat'].startswith('FAM_'): return f"{p['cat']}|{p['lin']}"
    return p['gene'] or p['prod'] or p['acc']
def scov_of(iv, slen):
    mg = []
    for a, e in sorted(iv):
        if mg and a <= mg[-1][1]: mg[-1][1] = max(mg[-1][1], e)
        else: mg.append([a, e])
    return sum(e - a + 1 for a, e in mg) / slen
def w(name, rows, cols=None):
    with open(f'{OUT}/{name}', 'w') as f:
        if not rows: return
        ks = cols or list(rows[0].keys()); f.write('\t'.join(ks) + '\n')
        for r in rows: f.write('\t'.join(str(r.get(k, '')) for k in ks) + '\n')

genomes = [l.rstrip('\n').split('\t') for l in open(GEN)]
carr_g = set(); carr_c = set(); carr_loci = defaultdict(list)
for r in csv.DictReader(open(S1), delimiter='\t'):
    if r['lineage_final'] in CARR:
        carr_g.add(num(r['gca'])); carr_c.add((num(r['gca']), sv(r['contig']))); carr_loci[(num(r['gca']), sv(r['contig']))].append(r['lineage_final'])
qdb = {}
for l in open(QF):
    if l[0] == '>': h = l[1:].strip(); qdb[h] = ps(h)
notdone = [g[0] for g in genomes if not os.path.exists(f'{D}/out/{g[0]}/done')]

strong_rows = []; best_subj = defaultdict(dict)   # subj -> gca -> dict(best strong/bg pid)
cat_stat = defaultdict(lambda: defaultdict(dict)); key_stat = defaultdict(lambda: defaultdict(dict))
chain_best = defaultdict(dict)   # subj -> gca -> (pid, cov) best chain with cov>=0.5
dh_rows = []; contam = []; hmm_stat = defaultdict(dict); genome_rows = []
def upd(dct, g, pid, cov):
    o = dct.get(g)
    if o is None: o = dct[g] = dict(n=0, bg=0, st=0, max=0.0, maxbg=0.0, maxall=0.0)
    o['n'] += 1; o['maxall'] = max(o['maxall'], pid)
    if cov >= 0.5 and pid >= 30: o['bg'] += 1; o['maxbg'] = max(o['maxbg'], pid)
    if cov >= 0.5 and pid >= 45: o['st'] += 1
DH = re.compile(r'dehydrogenase', re.I)
def is_dh(s):
    p = ps(s)
    return p['cat'] in ('TDNA', 'opine_related') and (DH.search(p['prod']) or p['gene'] in ('nos', 'ocs', 'noxA', 'noxB', 'odh', 'ons'))
for gca, org, genus, fna in genomes:
    od = f'{D}/out/{gca}'
    if gca in notdone: continue
    hs = defaultdict(list); clen = {}
    for l in open(f'{od}/dmnd.tsv'):
        f = l.rstrip('\n').split('\t')
        if len(f) < 13: continue
        q = f[0]; clen[q] = int(f[4]); qs, qe = sorted((int(f[6]), int(f[7])))
        hs[(q, '+' if int(f[12]) > 0 else '-', f[1])].append((qs, qe, int(f[8]), int(f[9]), float(f[2]), float(f[11]), int(f[5])))
    chains = defaultdict(list)
    for (q, s, sub), hl in hs.items():
        hl.sort(); cur = None
        for qs, qe, ss, se, pid, bits, slen in hl:
            if cur and qs - cur['end'] <= GAP:
                cur['end'] = max(cur['end'], qe); cur['bits'] += bits; cur['scov'].append((min(ss, se), max(ss, se))); cur['pid'] = max(cur['pid'], pid)
            else:
                if cur: chains[(q, s)].append(cur)
                cur = dict(start=qs, end=qe, subj=sub, bits=bits, scov=[(min(ss, se), max(ss, se))], pid=pid, slen=slen)
        chains[(q, s)].append(cur)
    cb = Counter(); cb_strict = Counter(); cb_best = {}
    gloc = gstrong = 0
    for (q, s), cl in chains.items():
        cl.sort(key=lambda c: c['start']); loci = []; cur = None
        for c in cl:
            c['sc'] = scov_of(c['scov'], c['slen'])
            if c['sc'] >= 0.5:
                o = chain_best[c['subj']].get(gca)
                if o is None or c['pid'] > o[0]: chain_best[c['subj']][gca] = (c['pid'], c['sc'])
            if is_dh(c['subj']) and c['pid'] >= 35 and c['sc'] >= 0.5:
                dh_rows.append(dict(gca=gca, organism=org, contig=q, contig_len=clen[q], start=c['start'], end=c['end'], strand=s, subject=c['subj'],
                                    pident=f"{c['pid']:.1f}", subj_cov=f"{c['sc']:.2f}", genome_is_carrier=int(num(gca) in carr_g),
                                    contig_has_carrier_locus=int((num(gca), sv(q)) in carr_c), carrier_lineages_on_contig=';'.join(sorted(set(carr_loci.get((num(gca), sv(q)), []))))))
            if cur and c['start'] <= cur['end']: cur['end'] = max(cur['end'], c['end']); cur['m'].append(c)
            else:
                if cur: loci.append(cur)
                cur = dict(start=c['start'], end=c['end'], m=[c])
        loci.append(cur)
        for L in loci:
            b = max(L['m'], key=lambda x: x['bits']); p = ps(b['subj']); gloc += 1
            if p['cat'] in BACK and b['pid'] >= CONTAM_PID:
                cb[q] += 1
                if b['sc'] >= 0.5: cb_strict[q] += 1
            upd(best_subj[b['subj']], gca, b['pid'], b['sc']); upd(cat_stat[p['cat']], gca, b['pid'], b['sc']); upd(key_stat[(p['cat'], gkey(b['subj']))], gca, b['pid'], b['sc'])
            if b['pid'] >= 45 and b['sc'] >= 0.5 and not p['cat'].startswith('FAM_'):
                gstrong += 1
                strong_rows.append(dict(gca=gca, organism=org, contig=q, contig_len=clen[q], start=L['start'], end=L['end'], strand=s, subject=b['subj'],
                                        gene_product=gkey(b['subj']), plasmid=p['plasmid'], category=p['cat'], pident=f"{b['pid']:.1f}", subj_cov=f"{b['sc']:.2f}"))
    for q, n in cb.items():
        contam.append(dict(gca=gca, organism=org, contig=q, contig_len=clen[q], n_backbone_loci_ge60=n, n_backbone_loci_ge60_cov05=cb_strict.get(q, 0),
                           genome_is_carrier=int(num(gca) in carr_g), contig_has_carrier_locus=int((num(gca), sv(q)) in carr_c)))
    for l in open(f'{od}/hmm.tbl'):
        if l[0] == '#': continue
        f = l.split(); e = float(f[4]); sc = float(f[5])
        o = hmm_stat[f[2]].get(gca)
        if o is None or sc > o[0]: hmm_stat[f[2]][gca] = (sc, e)
    genome_rows.append(dict(gca=gca, organism=org, n_loci=gloc, n_strong_nonfamily=gstrong, n_contam_contigs_agg=len(cb), n_contam_contigs_cov05=len(cb_strict), is_carrier=int(num(gca) in carr_g)))
ng = len(genomes) - len(notdone)
def stats(dct):
    ms = [o['maxall'] for o in dct.values()]
    return dict(n_asm_strong=sum(1 for o in dct.values() if o['st']), n_loci_strong=sum(o['st'] for o in dct.values()),
                n_asm_bg30=sum(1 for o in dct.values() if o['bg']), n_loci_bg30=sum(o['bg'] for o in dct.values()),
                n_asm_any_locus=len(dct), median_best_pident=f'{st.median(ms):.1f}' if ms else '', max_pident=f'{max(ms):.1f}' if ms else '')
# --- per category / per gene key / per subject ---
cat_rows = [dict(category=c, **stats(d)) for c, d in sorted(cat_stat.items())]
key_rows = [dict(category=c, gene_or_product=k, **stats(d)) for (c, k), d in sorted(key_stat.items())]
sub_rows = []
for h, p in qdb.items():
    d = best_subj.get(h, {})
    sub_rows.append(dict(subject=h, category=p['cat'], plasmid=p['plasmid'], gene_or_product=gkey(h), **stats(d)))
w('category_summary.tsv', cat_rows); w('gene_summary.tsv', key_rows); w('subject_summary.tsv', sub_rows)
# --- HMM lineages ---
cl = {r['rep']: r for r in csv.DictReader(open(CL), delimiter='\t')}
MAN = {'iaaH_iaaH': 'AAF77122.1', 'ipt_ipt': 'AAF77124.1', 'mas1__mas1_': 'AAF77131.1', 'nos_D-nopaline_dehydroge': 'WP_010974826.1', 'ocs_ocs': 'AAF77127.1',
       'ons_ons': 'AAF77125.1', 'rolB_C__protein': 'AAK90962.1', 'gene_': 'AAF77128.1'}
hrows = []
for name, d in sorted(hmm_stat.items()):
    if '__' in name and name.split('__')[0] in ('ACS', 'CUS', 'SUS', 'ORF358', 'OG') or name.startswith('OG|'): typ = 'family_subgroup'
    else: typ = 'tdna_lineage'
    m = re.match(r'^([A-Z]+_?\d+\.\d+)', name); rep = m.group(1) if m else MAN.get(name, '')
    c = cl.get(rep, {})
    ok = {g: v for g, v in d.items() if v[1] <= 1e-5}
    hrows.append(dict(hmm=name, type=typ, cluster_rep=rep, cluster_genes=c.get('genes', ''), cluster_products=c.get('products', ''),
                      n_asm_E_le_1e5=len(ok), best_score=f"{max((v[0] for v in d.values()), default=0):.1f}",
                      best_score_asm=max(d, key=lambda g: d[g][0]) if d else '', rep_mapping='manual(unverified)' if name in MAN else ''))
w('hmm_profile_summary.tsv', hrows)
# --- strong loci ---
strong_rows.sort(key=lambda r: (r['category'], r['gene_product'], r['gca'], r['contig'], int(r['start'])))
w('strong_loci.tsv', strong_rows)
grp = defaultdict(list)
for r in strong_rows: grp[(r['category'], r['gene_product'], r['plasmid'])].append(r)
w('strong_loci_by_gene.tsv', [dict(category=k[0], gene_or_product=k[1], plasmid=k[2], n_loci=len(v), n_asm=len({r['gca'] for r in v}),
    max_pident=max(float(r['pident']) for r in v), organisms='; '.join(sorted({f"{r['organism']} ({r['gca']})" for r in v}))) for k, v in sorted(grp.items())])
# --- opine dehydrogenases ---
dh_rows.sort(key=lambda r: (r['gca'], r['contig'], int(r['start']))); w('opine_dehydrogenase_loci.tsv', dh_rows)
dg = defaultdict(list)
for r in dh_rows: dg[r['gca']].append(r)
w('opine_dehydrogenase_genomes.tsv', [dict(gca=g, organism=v[0]['organism'], n_chains=len(v), max_pident=max(float(r['pident']) for r in v), genome_is_carrier=v[0]['genome_is_carrier'],
    n_on_carrier_contig=sum(r['contig_has_carrier_locus'] for r in v), subjects=';'.join(sorted({r['subject'].split('|')[0] + ':' + gkey(r['subject']) for r in v}))) for g, v in sorted(dg.items())])
contam.sort(key=lambda r: (r['gca'], -r['contig_len'])); w('contamination_contigs.tsv', contam)
w('genome_summary.tsv', genome_rows)
# --- classic genes (manual mapping of gene -> query accessions; see report) ---
def qs_by(pred): return [h for h, p in qdb.items() if pred(h, p)]
fam_acs = qs_by(lambda h, p: p['cat'] == 'FAM_ACS')
CL_MAP = {'rolB (AAK90962, gene=rolB)': ['AAK90962.1'], 'rolC? (AAK90962 product "C_protein"; same query as rolB)': ['AAK90962.1'],
          '5 (AAK90960 "5_protein")': ['AAK90960.1'], '6a (AAF77125 "ons"; cluster product 6a_protein)': ['AAF77125.1'],
          '6b/tml (WP_174089376 cluster tml/tiorf186; AAZ50400 Tml)': ['WP_174089376.1', 'AAZ50400.1'],
          'plast D (WP_174168441)': ['WP_174168441.1'], 'plast E (WP_172690589)': ['WP_172690589.1'],
          'RolB_family_protein (all queries)': qs_by(lambda h, p: p['prod'] == 'RolB_family_protein'),
          'orf13/orf14': [], 'iaaM (WP_010974824 FAD-dependent oxidoreductase, cluster iaaM/tms1)': ['WP_010974824.1'],
          'iaaH': qs_by(lambda h, p: p['gene'] == 'iaaH'), 'ipt': qs_by(lambda h, p: p['gene'] == 'ipt' or p['prod'].lower() == 'ipt'),
          'mas1': qs_by(lambda h, p: p['gene'].startswith('mas1') or p['prod'] == 'Mannopine_synthase'), 'mas2 (WP_262512305 SIS domain; cluster mas2)': ['WP_262512305.1'],
          'ags (WP_080864006 M24 peptidase; cluster ags)': ['WP_080864006.1'], 'ocs': qs_by(lambda h, p: p['gene'] == 'ocs'), 'nos': qs_by(lambda h, p: p['gene'] == 'nos'),
          'acs agrocinopine synthase (WP_010974816 acsx, cluster)': ['WP_010974816.1'], 'acs family queries FAM_ACS T1/T2': fam_acs}
crow = []
for lab, ids in CL_MAP.items():
    full = [h for h in qdb if h.split('|')[0] in ids or h in ids] if not ids or '|' not in ids[0] else ids
    asm_st = set(); mx = 0.0; mxloc = 0.0; asm_loc = set(); asm_chain45 = set(); mxchain = 0.0
    for h in full:
        for g, o in best_subj.get(h, {}).items():
            mxloc = max(mxloc, o['maxall'])
            if o['st']: asm_st.add(g)
        for g, (pid, cov) in chain_best.get(h, {}).items():
            mxchain = max(mxchain, pid)
            if pid >= 45: asm_chain45.add(g)
    crow.append(dict(gene=lab, n_queries=len(full), queries=';'.join(x.split('|')[0] + (':' + x.split('|')[2] if x.startswith('FAM_') else '') for x in full)[:200],
                     n_asm_strong_locus_best_subject=len(asm_st), max_pident_best_locus=f'{mxloc:.1f}' if full else 'n/a',
                     n_asm_chain_pid45_cov05=len(asm_chain45), max_chain_pident_cov05=f'{mxchain:.1f}' if full else 'n/a'))
w('classic_genes.tsv', crow)
print('assemblies', len(genomes), 'done', ng, 'notdone', notdone)
print('strong loci', len(strong_rows), 'asm', len({r['gca'] for r in strong_rows}), 'dh chains', len(dh_rows), 'dh asm', len(dg), 'contam contigs', len(contam), 'asm', len({r['gca'] for r in contam}))
