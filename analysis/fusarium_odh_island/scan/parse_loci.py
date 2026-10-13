#!/usr/bin/env python3
"""Parse miniprot runs (mp_out/<GCA>/*.gff) for the F. poae ODH island (3 seeds: ODH KAL3595992.1 + neighbours KAL3595993.1, KAL3595994.1)
-> loci_all.tsv (every hit cluster per seed), loci.tsv (island-like: identity >= ID_MIN), genomes.tsv, species_summary.tsv, neighbors.tsv.
Adapted from chs_fusarium/stage2/parse_loci.py. Usage: parse_loci.py <genomes_files.tsv: GCA fna gff>"""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import sys, os, re, csv, collections
S = f'{GMO_ROOT}/analysis/fusarium_odh_island/scan'
CHS = f'{GMO_ROOT}/analysis/chs_fusarium/stage2'
SEEDS = {'KAL3595992.1': 353, 'KAL3595993.1': 383, 'KAL3595994.1': 259}
SNAME = {'KAL3595992.1': 'odh', 'KAL3595993.1': 'pepP', 'KAL3595994.1': 'g3'}
ID_MIN = 0.60   # island-like; native fungal OCD-family homologs are <= ~0.40
_src = open(f'{CHS}/parse_loci.py').read()
exec(_src[:_src.index("H = ['gca'")])  # FSSC, MOB, gl, species(), complex_of(), parse_mp(), read_contig(), gc(), parse_gff() from the chs script
S = f'{GMO_ROOT}/analysis/fusarium_odh_island/scan'; SEEDS = {'KAL3595992.1': 353, 'KAL3595993.1': 383, 'KAL3595994.1': 259}  # restore after exec

H = ['gca','organism','species','complex','level','seed','gene','locus_id','contig','start','end','strand','contig_len','dist_left','dist_right','identity','coverage',
     'cls','status','qfrom','qto','seed_len','exons','frameshifts','stops','missing','edge','gene_gc','assembly_gc',
     'annot_gene','annot_prod','annot_protein','annot_exons','annot_pseudo','annot_overlap_bp','annot_aa','stop_pos_in_seed','note']
files = [l.rstrip('\n').split('\t') for l in open(sys.argv[1])]
fo = open(f'{S}/loci_all.tsv', 'w'); fo.write('\t'.join(H) + '\n')
fn = open(f'{S}/neighbors.tsv', 'w'); fn.write('gca\tlocus_id\trel\tgene\tprotein\tproduct\texons\tdist_bp\tmobile\n')
gen = []; nmiss = 0
for gca, fna, gff in files:
    d = f'{S}/mp_out/{gca}'
    if not os.path.exists(d + '/done'): nmiss += 1; continue
    org = gl[gca]['organism']; sp, _, _ = species(org); cx = complex_of(org)
    ctg = {}; tot = 0; gcs = 0
    for l in open(d + '/contigs.tsv'):
        n, ln, g = l.split('\t')[:3]; ln = int(ln); g = float(g) if g.strip() not in ('', 'nan') else 0; ctg[n] = ln; tot += ln; gcs += ln * g
    agc = gcs / tot if tot else float('nan')
    rows = []
    for sd in SEEDS:
        hits = parse_mp(f'{d}/{sd}.gff', sd)
        hits.sort(key=lambda h: (h['contig'], h['start']))
        clusters = []
        for h in hits:
            for c in clusters:
                if c['contig'] == h['contig'] and c['strand'] == h['strand'] and h['start'] <= c['end'] and h['end'] >= c['start']:
                    c['hits'].append(h); c['start'] = min(c['start'], h['start']); c['end'] = max(c['end'], h['end']); break
            else: clusters.append(dict(contig=h['contig'], strand=h['strand'], start=h['start'], end=h['end'], hits=[h]))
        for k, c in enumerate(clusters, 1):
            prim = max(c['hits'], key=lambda h: h['score'])
            cls = 'island_like' if prim['ident'] >= ID_MIN and prim['cov'] >= .25 else 'other'
            L = ctg.get(c['contig'], 0); dl, dr = prim['start'] - 1, L - prim['end']; ql = SEEDS[sd]
            if prim['fs'] or prim['stop']: status = 'pseudogene'
            elif prim['cov'] >= 0.9: status = 'intact'
            else: status = 'truncated'
            nmis = prim['qfrom'] - 1; cmis = ql - prim['qto']; miss = []; edge = ''
            if status == 'truncated':
                plus = prim['strand'] == '+'
                for lab, n, side in (('N', nmis, dl if plus else dr), ('C', cmis, dr if plus else dl)):
                    if n > 5: miss.append(f'{lab}:{n}aa')
                    if n > 5 and side < 3 * n + 30: edge += f'{lab}-end blocked by contig edge ({side}nt); '
            if cls == 'other': status = 'NA'
            rows.append(dict(gca=gca, organism=org, species=sp, complex=cx, level=gl[gca]['level'], seed=sd, gene=SNAME[sd], locus_id=f'{gca}_{SNAME[sd]}{k}', contig=c['contig'],
                 start=prim['start'], end=prim['end'], strand=prim['strand'], contig_len=L, dist_left=dl, dist_right=dr, identity=prim['ident'], coverage=round(prim['cov'], 3),
                 cls=cls, status=status, qfrom=prim['qfrom'], qto=prim['qto'], seed_len=ql, exons=len(prim['cds']), frameshifts=prim['fs'], stops=prim['stop'],
                 missing=','.join(miss), edge=edge.strip('; '), gene_gc='', assembly_gc=round(agc, 2), annot_gene='', annot_prod='', annot_protein='', annot_exons='', annot_pseudo='',
                 annot_overlap_bp='', annot_aa='', stop_pos_in_seed=','.join(map(str, prim['stoppos'])), note='', _cds=prim['cds']))
    sel = [r for r in rows if r['cls'] == 'island_like']
    if sel:
        bc = parse_gff(gff) if gff and os.path.exists(gff) else {}
        for r in sel:
            seq = read_contig(fna, r['contig'])
            cds = ''.join(seq[s - 1:e] for s, e in sorted(r['_cds']))
            r['gene_gc'] = round(gc(cds), 2)
            gs = bc.get(r['contig'], [])
            ov = [(min(G['end'], r['end']) - max(G['start'], r['start']) + 1, G) for G in gs if G['start'] <= r['end'] and G['end'] >= r['start']]
            if ov:
                o, G = max(ov, key=lambda x: x[0]); r.update(annot_gene=G['name'], annot_prod=G['prod'], annot_protein=G['prot'], annot_exons=G['ex'], annot_pseudo=int(G['pseudo']), annot_overlap_bp=o, annot_aa=G.get('aa', ''))
            if r['status'] == 'pseudogene':
                r['note'] = (f"in-frame stop(s) at seed aa {r['stop_pos_in_seed'] or 'NA'}" if r['stops'] else '') + ('; ' if r['stops'] and r['frameshifts'] else '') + (f"{r['frameshifts']} frameshift(s) vs seed" if r['frameshifts'] else '')
            ids = {G['id'] for _, G in ov}
            left = [G for G in gs if G['end'] < r['start'] and G['id'] not in ids][-10:][::-1]
            right = [G for G in gs if G['start'] > r['end'] and G['id'] not in ids][:10]
            for sgn, lst in ((-1, left), (1, right)):
                for i, G in enumerate(lst, 1):
                    dist = r['start'] - G['end'] if sgn < 0 else G['start'] - r['end']
                    fn.write('\t'.join(map(str, [gca, r['locus_id'], sgn * i, G['name'], G['prot'], G['prod'], G['ex'], dist, int(bool(MOB.search(G['prod'] + ' ' + G['btype'])))])) + '\n')
    for r in rows: fo.write('\t'.join(str(r[h]) for h in H) + '\n')
    g = dict(gca=gca, organism=org, species=sp, complex=cx, level=gl[gca]['level'], n_contigs=len(ctg), assembly_len=tot, assembly_gc=round(agc, 2))
    for sd, nm in SNAME.items():
        st = [r['status'] for r in sel if r['seed'] == sd]
        g[nm] = 'absent' if not st else ('intact' if 'intact' in st else 'truncated' if 'truncated' in st else 'pseudogene')
        g[nm + '_n'] = len(st)
    gen.append(g)
fo.close(); fn.close()
with open(f'{S}/loci.tsv', 'w') as f:
    f.write('\t'.join(H) + '\n')
    for l in open(f'{S}/loci_all.tsv'):
        if '\tisland_like\t' in l: f.write(l)
with open(f'{S}/genomes.tsv', 'w') as f:
    ks = list(gen[0]); f.write('\t'.join(ks) + '\n')
    for g in gen: f.write('\t'.join(str(g[k]) for k in ks) + '\n')
cnt = collections.defaultdict(collections.Counter)
for g in gen: cnt[g['species']]['total'] += 1; cnt[g['species']]['odh_' + g['odh']] += 1; cnt[g['species']]['pepP_' + g['pepP']] += 1; cnt[g['species']]['g3_' + g['g3']] += 1
with open(f'{S}/species_summary.tsv', 'w') as f:
    cols = ['odh_intact','odh_truncated','odh_pseudogene','odh_absent','pepP_intact','pepP_truncated','pepP_pseudogene','pepP_absent','g3_intact','g3_truncated','g3_pseudogene','g3_absent']
    f.write('species\tcomplex\tgenomes\t' + '\t'.join(cols) + '\n')
    cxm = {g['species']: g['complex'] for g in gen}
    for s in sorted(cnt, key=lambda s: -cnt[s]['total']):
        c = cnt[s]; f.write(f"{s}\t{cxm[s]}\t{c['total']}\t" + '\t'.join(str(c[k]) for k in cols) + '\n')
print('genomes parsed', len(gen), 'missing mp', nmiss)
