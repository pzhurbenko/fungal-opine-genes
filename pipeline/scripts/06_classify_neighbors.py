#!/usr/bin/env python3
"""Stage 06 — classify neighbouring genes (and targets) by their best nr hits.

  06_classify_neighbors.py <round> collect   -> round<N>/diamond_query.faa (all neighbors.faa + target proteins, deduplicated)
  06_classify_neighbors.py <round> parse     -> round<N>/neighbors_class.tsv, round<N>/locus_context.tsv
(the diamond run itself is submitted to LSF by run_pipeline.sh between the two calls)

The cluster's diamond nr has NO taxonomy map, so the organism is parsed from the hit title "[Organism]"
and mapped to a kingdom via NCBI Taxonomy (esearch/efetch, cached in out/registry/taxonomy_cache.json).
Per neighbour: top-hit kingdom, majority kingdom over top hits, best e-value, TE flag (keyword regex on titles).
Per locus: counts of fungal / bacterial / plant / TE / no-hit neighbours inside the adaptive window.
"""
import os, sys, re, json, glob
from collections import Counter, defaultdict
from common import load_config, round_dir, read_fasta, write_fasta, read_tsv, write_tsv, eutils, log

cfg = load_config(); P = cfg['_paths']; rnd, mode = sys.argv[1], sys.argv[2]; RD = round_dir(cfg, rnd)
Q = os.path.join(RD, 'diamond_query.faa'); OUT = os.path.join(RD, 'diamond_nr.tsv')

if mode == 'collect':
    recs = {}
    for f in glob.glob(os.path.join(RD, 'loci', '*', 'neighbors.faa')) + glob.glob(os.path.join(RD, 'loci', '*', 'target_protein.faa')):
        for h, s in read_fasta(f):
            if s: recs[h.split()[0]] = s.replace('*', '')
    write_fasta(Q, sorted(recs.items())); log(f'diamond query: {len(recs)} proteins -> {Q}'); sys.exit()

# ---------------------------------------------------------------- parse
TE = re.compile(cfg['validate']['te_keywords'], re.I)
cache_p = os.path.join(P['registry'], 'taxonomy_cache.json')
cache = json.load(open(cache_p)) if os.path.exists(cache_p) else {}
KINGDOM_RANKS = {'Fungi': 'Fungi', 'Bacteria': 'Bacteria', 'Viridiplantae': 'Plants', 'Metazoa': 'Metazoa', 'Archaea': 'Archaea', 'Viruses': 'Viruses'}

def lineage_kingdom(org):
    org = re.sub(r'\s+(strain|str\.|isolate|sp\.).*$', '', org).strip()
    if org in cache: return cache[org]
    k = 'other'
    try:
        ids = eutils(cfg, 'esearch.fcgi', {'db': 'taxonomy', 'term': org, 'retmode': 'json'}, as_json=True)['esearchresult']['idlist']
        if not ids and ' ' in org:  # retry with genus only
            ids = eutils(cfg, 'esearch.fcgi', {'db': 'taxonomy', 'term': org.split()[0], 'retmode': 'json'}, as_json=True)['esearchresult']['idlist']
        if ids:
            xml = eutils(cfg, 'efetch.fcgi', {'db': 'taxonomy', 'id': ids[0], 'retmode': 'xml'})
            lin = re.search(r'<Lineage>(.*?)</Lineage>', xml, re.S); lin = lin.group(1) if lin else ''
            for key, val in KINGDOM_RANKS.items():
                if key in lin.split('; '): k = val; break
            if k == 'other' and 'Eukaryota' in lin: k = 'other_eukaryote'
    except Exception as ex:
        log('taxonomy lookup failed', org, ex)
    cache[org] = k; return k

hits = defaultdict(list)
if os.path.exists(OUT):
    for line in open(OUT):
        f = line.rstrip('\n').split('\t')
        if len(f) < 9: continue
        q, s, pid, alen, qlen, slen, ev, bits, title = f[:9]
        org = re.findall(r'\[([^\]]+)\]', title); org = org[-1] if org else ''
        hits[q].append({'sseqid': s, 'pident': float(pid), 'evalue': float(ev), 'bits': float(bits), 'title': title, 'org': org})
orgs = sorted({h['org'] for hs in hits.values() for h in hs if h['org']})
log(f'unique organisms to classify: {len(orgs)} (cached {sum(1 for o in orgs if o in cache)})')
for i, o in enumerate(orgs):
    lineage_kingdom(o)
    if (i + 1) % 200 == 0: json.dump(cache, open(cache_p, 'w')); log(f'  taxonomy {i+1}/{len(orgs)}')
json.dump(cache, open(cache_p, 'w'))

rows = []
for h, s in read_fasta(Q):
    q = h.split()[0]; hs = sorted(hits.get(q, []), key=lambda x: -x['bits'])[:cfg['diamond']['max_target_seqs']]
    locus = q.split('|')[0]; kind = 'target' if '|miniprot|' in q else 'neighbor'
    if not hs:
        rows.append({'query': q, 'locus_id': locus, 'kind': kind, 'n_hits': 0, 'top_kingdom': 'no_hit', 'majority_kingdom': 'no_hit', 'top_pident': '', 'top_evalue': '', 'top_org': '', 'top_title': '', 'is_TE': 'no'}); continue
    ks = [lineage_kingdom(x['org']) if x['org'] else 'other' for x in hs]
    top = hs[0]
    rows.append({'query': q, 'locus_id': locus, 'kind': kind, 'n_hits': len(hs), 'top_kingdom': ks[0], 'majority_kingdom': Counter(ks).most_common(1)[0][0],
                 'top_pident': top['pident'], 'top_evalue': f"{top['evalue']:.1e}", 'top_org': top['org'], 'top_title': top['title'][:120],
                 'is_TE': 'yes' if any(TE.search(x['title']) for x in hs[:3]) else 'no'})
write_tsv(os.path.join(RD, 'neighbors_class.tsv'), rows)

# ---- per-locus context summary (neighbours inside the adaptive window only)
# Neighbours that overlap ANOTHER registered locus on the same contig (e.g. the CUS next to an ACS in a T-DNA-like
# cluster) are T-DNA candidates themselves, not evidence of bacterial contamination -> counted separately as nb_tdna_loci.
reg_loci = defaultdict(list)
for r in read_tsv(os.path.join(P['registry'], 'loci.tsv')):
    reg_loci[r['contig']].append((r['locus_id'], int(r['hsp_start']), int(r['hsp_end'])))
ctx = []
for jp in glob.glob(os.path.join(RD, 'loci', '*', 'locus.json')):
    J = json.load(open(jp)); lid = J['locus_id']; ad = J['adaptive_neighbors']
    hs_s, hs_e = J['hsp_local']; cs = J['context_local'][0]; s0 = J['seq_start']
    lo = hs_s - ad.get('left', {}).get('distance_bp', 0) + s0 - 1; hi = hs_e + ad.get('right', {}).get('distance_bp', 0) + s0 - 1
    others = [(a, b) for l2, a, b in reg_loci[J['contig']] if l2 != lid]
    nb = [r for r in rows if r['locus_id'] == lid and r['kind'] == 'neighbor']
    inwin = []; tdna = 0
    for r in nb:
        m = re.search(r'\|(\d+)-(\d+)[+-]\|', r['query'])
        if not (m and int(m.group(2)) >= lo and int(m.group(1)) <= hi): continue
        gs, ge = int(m.group(1)), int(m.group(2))
        if any(gs <= b and ge >= a for a, b in others): tdna += 1; continue
        inwin.append(r)
    c = Counter(r['majority_kingdom'] for r in inwin)
    ctx.append({'locus_id': lid, 'neighbors_in_window': len(inwin), 'nb_fungi': c.get('Fungi', 0), 'nb_bacteria': c.get('Bacteria', 0),
                'nb_plants': c.get('Plants', 0), 'nb_other': sum(v for k, v in c.items() if k not in ('Fungi', 'Bacteria', 'Plants', 'no_hit')),
                'nb_no_hit': c.get('no_hit', 0), 'nb_TE': sum(1 for r in inwin if r['is_TE'] == 'yes'), 'nb_tdna_loci': tdna,
                'window_left_bp': ad.get('left', {}).get('distance_bp', ''), 'window_right_bp': ad.get('right', {}).get('distance_bp', ''),
                'target_top_org': next((r['top_org'] for r in rows if r['locus_id'] == lid and r['kind'] == 'target'), ''),
                'target_top_kingdom': next((r['top_kingdom'] for r in rows if r['locus_id'] == lid and r['kind'] == 'target'), '')})
write_tsv(os.path.join(RD, 'locus_context.tsv'), ctx)
log('done:', len(rows), 'proteins classified;', len(ctx), 'loci summarised')
