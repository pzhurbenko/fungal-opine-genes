#!/usr/bin/env python3
"""Island composition around every homologous locus of round 2 (paper plan 6.2).
Uses the round-2 outputs only: neighbours (all AUGUSTUS genes in the context, up to 250 kb) and their diamond-nr hits.
Gene classification (hits of the SAME GENUS are excluded — annotated genomes hit their own proteins first):
  top kingdom of remaining hits; 'bacterial_like' = top kingdom Bacteria OR >=3/10 bacterial hits OR gene overlaps a registry locus.
Island = genes between the nearest 'fungal anchor' on each side of the target (anchor = top kingdom Fungi, not bacterial_like, >=2 exons).
Outputs: islands.tsv (one row per locus), islands_genes.tsv (genes inside islands), summary on stdout."""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import os, re, csv, json, glob, sys
from collections import defaultdict, Counter
OUT = f'{GMO_ROOT}/analysis/loci_pipeline'; RD = f'{OUT}/round2'; HERE = os.path.dirname(os.path.abspath(__file__))
HOM = {'intact', 'pseudogene_like', 'fragment', 'remnant'}
def tsv(p): return list(csv.DictReader(open(p), delimiter='\t'))
T = tsv(f'{RD}/loci_table.tsv'); loci = [r for r in T if r['status'] in HOM]
cache = json.load(open(f'{OUT}/registry/taxonomy_cache.json'))
def kingdom(org):
    o = re.sub(r'\s+(strain|str\.|isolate|sp\.).*$', '', org).strip(); return cache.get(o, 'other')
reg = defaultdict(list)
for r in T: reg[r['contig']].append((r['locus_id'], r['fam'], int(r['hsp_start']), int(r['hsp_end']), r['status'], r['subgroup']))
asm_gc = {}
for r in tsv(f'{OUT}/registry/assemblies.tsv'):
    if r['gc_percent']: asm_gc[r['wgs_prefix']] = float(r['gc_percent'])
wgs = lambda c: (re.match(r'^([A-Z]{6}\d{2}|[A-Z]{4}\d{2})', c) or [None, None])[1] if re.match(r'^([A-Z]{6}\d{2}|[A-Z]{4}\d{2})', c) else None
# ---- diamond hits per query (only neighbours/targets of homolog loci)
need = {r['locus_id'] for r in loci}
hits = defaultdict(list)
for line in open(f'{RD}/diamond_nr.tsv'):
    f = line.rstrip('\n').split('\t')
    if len(f) < 9 or f[0].split('|')[0] not in need: continue
    org = re.findall(r'\[([^\]]+)\]', f[8]); org = org[-1] if org else ''
    hits[f[0]].append((float(f[7]), float(f[2]), org, f[8][:100]))
sys.stderr.write(f'loci: {len(loci)}; queries with hits: {len(hits)}\n')
def gc_of(seq): 
    s = seq.upper(); n = sum(s.count(c) for c in 'ACGT'); return 100.0 * (s.count('G') + s.count('C')) / n if n else float('nan')
seqcache = {}
def contig_seq(contig):
    if contig not in seqcache:
        p = f'{OUT}/contigs/{contig}.fa'; h, s = open(p).read().split('\n', 1); s = s.replace('\n', '')
        seqcache[contig] = (int(re.search(r'seq_start=(\d+)', h).group(1)), s)
    return seqcache[contig]
rows = []; grows = []
for L in loci:
    lid = L['locus_id']; d = f'{RD}/loci/{lid}'
    if not os.path.exists(f'{d}/neighbors.faa'): continue
    genus = L['organism'].split()[0].strip('[]')
    genes = []
    for line in open(f'{d}/neighbors.faa'):
        if not line.startswith('>'): continue
        m = re.match(r'>(\S+?)\|(g\d+)\|(\d+)-(\d+)([+-])\|exons=(\d+)', line)
        if m: genes.append({'q': line[1:].split()[0], 'id': m.group(2), 's': int(m.group(3)), 'e': int(m.group(4)), 'strand': m.group(5), 'exons': int(m.group(6))})
    ts, te = int(L['contig_start'] or L['hsp_start']), int(L['contig_end'] or L['hsp_end'])
    others = [(l2, f2, a, b, st, sg) for l2, f2, a, b, st, sg in reg[L['contig']] if l2 != lid]
    for g in genes:
        hs = sorted(hits.get(g['q'], []), key=lambda x: -x[0])[:10]
        ex = [h for h in hs if not h[2].split()[0].strip('[]') == genus] if hs else []
        ks = [kingdom(h[2]) if h[2] else 'other' for h in ex]
        g['n_hits'] = len(hs); g['top_kingdom'] = ks[0] if ks else ('self_only' if hs else 'no_hit')
        g['n_bact'] = sum(1 for k in ks if k == 'Bacteria'); g['top_org'] = ex[0][2] if ex else (hs[0][2] if hs else '')
        g['top_pident'] = ex[0][1] if ex else ''; g['top_title'] = re.sub(r'^\S+\s', '', ex[0][3]) if ex else ''
        ov = [l2 for l2, f2, a, b, st, sg in reg[L['contig']] if g['s'] <= b and g['e'] >= a]
        g['registry_locus'] = ','.join(ov)
        g['is_target'] = g['s'] <= te and g['e'] >= ts
        g['bacterial_like'] = g['top_kingdom'] == 'Bacteria' or g['n_bact'] >= 3 or bool(ov)
        g['anchor'] = (not g['bacterial_like']) and g['top_kingdom'] == 'Fungi' and g['exons'] >= 2
    genes.sort(key=lambda g: g['s'])
    # island boundaries: nearest anchors left/right of the target
    left = [g for g in genes if g['e'] < ts]; right = [g for g in genes if g['s'] > te]
    la = next((g for g in reversed(left) if g['anchor']), None); ra = next((g for g in right if g['anchor']), None)
    seq0, seq = contig_seq(L['contig'])
    isl_s = la['e'] + 1 if la else seq0; isl_e = ra['s'] - 1 if ra else seq0 + len(seq) - 1
    inside = [g for g in genes if g['s'] >= isl_s and g['e'] <= isl_e]
    gc_isl = gc_of(seq[isl_s - seq0: isl_e - seq0 + 1]); gc_ctx = gc_of(seq)
    c = Counter('target' if g['is_target'] else 'bacterial_like' if g['bacterial_like'] else g['top_kingdom'] for g in inside)
    fam_inside = sorted({f2 for l2, f2, a, b, st, sg in others if a >= isl_s and b <= isl_e})
    rows.append({'locus_id': lid, 'fam': L['fam'], 'subgroup': L['subgroup'], 'status': L['status'], 'organism': L['organism'], 'class': L['class'],
                 'contig': L['contig'], 'contig_len': L['contig_len'], 'island_start': isl_s, 'island_end': isl_e, 'island_len': isl_e - isl_s + 1,
                 'bounded_left': 'yes' if la else 'no', 'bounded_right': 'yes' if ra else 'no',
                 'n_genes_inside': len(inside), 'n_bacterial_like': c['bacterial_like'], 'n_no_hit': c['no_hit'] + c['self_only'],
                 'n_fungal_singleexon': sum(1 for g in inside if not g['is_target'] and not g['bacterial_like'] and g['top_kingdom'] == 'Fungi'),
                 'n_other': sum(v for k, v in c.items() if k not in ('target', 'bacterial_like', 'no_hit', 'self_only', 'Fungi')),
                 'other_fams_inside': ','.join(fam_inside), 'n_registry_loci_inside': sum(1 for l2, f2, a, b, st, sg in others if a >= isl_s and b <= isl_e),
                 'gc_island': f'{gc_isl:.1f}', 'gc_context': f'{gc_ctx:.1f}', 'gc_assembly': asm_gc.get(wgs(L['contig']) or '', ''),
                 'context_call': L['context_call']})
    for g in inside:
        grows.append({'locus_id': lid, 'organism': L['organism'], 'contig': L['contig'], 'gene': g['id'], 'start': g['s'], 'end': g['e'], 'strand': g['strand'],
                      'exons': g['exons'], 'role': 'TARGET' if g['is_target'] else ('bacterial_like' if g['bacterial_like'] else g['top_kingdom']),
                      'registry_locus': g['registry_locus'], 'n_hits': g['n_hits'], 'n_bact_of_top10': g['n_bact'], 'top_org': g['top_org'],
                      'top_pident': g['top_pident'], 'top_title': g['top_title']})
def w(p, rs):
    with open(p, 'w') as f:
        wr = csv.DictWriter(f, fieldnames=list(rs[0].keys()), delimiter='\t'); wr.writeheader(); wr.writerows(rs)
w(f'{HERE}/islands.tsv', rows); w(f'{HERE}/islands_genes.tsv', grows)
print('loci:', len(rows))
print('\n== островки, ограниченные грибными якорями с обеих сторон: длина и число не-грибных генов внутри ==')
import statistics as st
for fam in ['ACS', 'CUS', 'SUS', 'ORF358']:
    rs = [r for r in rows if r['fam'] == fam and r['bounded_left'] == 'yes' and r['bounded_right'] == 'yes']
    if not rs: continue
    L = [r['island_len'] for r in rs]; nb = [r['n_bacterial_like'] for r in rs]
    print(f"{fam:7s} n={len(rs):4d}  длина медиана {st.median(L)/1e3:.1f} кб (q90 {sorted(L)[int(len(L)*.9)]/1e3:.0f} кб); только целевой ген внутри: {sum(1 for r in rs if r['n_genes_inside']==1)}; "
          f">=1 доп. бактериальноподобный: {sum(1 for x in nb if x>=1)}; >=3: {sum(1 for x in nb if x>=3)}; др. семейство внутри: {sum(1 for r in rs if r['other_fams_inside'])}")
print('\n== крупнейшие островки с >=2 бактериальноподобными генами (не считая цели), по классам ==')
big = sorted([r for r in rows if r['n_bacterial_like'] >= 2], key=lambda r: -r['n_bacterial_like'])
for r in big[:40]:
    print(f"  {r['organism'][:28]:28s} {r['class'][:15]:15s} {r['locus_id']:32s} {r['island_len']/1e3:6.1f} кб  генов {r['n_genes_inside']:3d} бакт {r['n_bacterial_like']:2d} nohit {r['n_no_hit']:2d} fungi1ex {r['n_fungal_singleexon']:2d}  GC {r['gc_island']}/{r['gc_context']}/{r['gc_assembly']}  fams {r['other_fams_inside']} bounded {r['bounded_left']}/{r['bounded_right']}")
