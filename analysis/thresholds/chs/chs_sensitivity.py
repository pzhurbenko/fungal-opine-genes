#!/usr/bin/env python3
"""Sensitivity of chs-locus thresholds (analysis/chs_fusarium/stage2/parse_loci.py:109-110).
Re-parses ALL miniprot hits (mp_out/*) with the same clustering as parse_loci.py, keeps ags_like_or_other rows
(which loci_all.tsv does not contain), and evaluates thresholds on id_chs and id_EEU (F. vanettenii seed)."""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import csv, collections, os
S = f'{GMO_ROOT}/analysis/chs_fusarium/stage2'
O = f'{GMO_ROOT}/analysis/thresholds/chs'
SEEDS = {'chs_A0A2Z2PGA9': 407, 'EEU41297.1': 405, 'KNB14377.1': 260}
gl = {r['gca']: r for r in csv.DictReader(open(f'{S}/genomes_list.tsv'), delimiter='\t')}
def parse_mp(path, seed):   # identical to parse_loci.py
    out = []; cur = None
    for line in open(path):
        if line[0] == '#': continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9: continue
        if f[2] == 'mRNA':
            a = dict(kv.split('=', 1) for kv in f[8].split(';') if '=' in kv); t = a['Target'].split()
            cur = dict(seed=seed, contig=f[0], start=int(f[3]), end=int(f[4]), score=int(f[5]), strand=f[6], ident=float(a['Identity']),
                       qfrom=int(t[1]), qto=int(t[2]), fs=int(a.get('Frameshift', 0)), stop=int(a.get('StopCodon', 0)))
            cur['cov'] = (cur['qto'] - cur['qfrom'] + 1) / SEEDS[seed]; out.append(cur)
    return out
old = {r['locus_id']: r for r in csv.DictReader(open(f'{S}/loci_all.tsv'), delimiter='\t')}
rows = []
for gca, fna, gff in (l.rstrip('\n').split('\t') for l in open(f'{S}/genomes_files.tsv')):
    d = f'{S}/mp_out/{gca}'
    hits = []
    for sd in SEEDS: hits += parse_mp(f'{d}/{sd}.gff', sd)
    hits.sort(key=lambda h: (h['contig'], h['start']))
    cl = []
    for h in hits:
        for c in cl:
            if c['contig'] == h['contig'] and c['strand'] == h['strand'] and h['start'] <= c['end'] and h['end'] >= c['start']:
                c['hits'].append(h); c['start'] = min(c['start'], h['start']); c['end'] = max(c['end'], h['end']); break
        else: cl.append(dict(contig=h['contig'], strand=h['strand'], start=h['start'], end=h['end'], hits=[h]))
    for k, c in enumerate(cl, 1):
        best = {}
        for h in c['hits']:
            if h['seed'] not in best or h['score'] > best[h['seed']]['score']: best[h['seed']] = h
        g = lambda s, key: best[s][key] if s in best else 0
        ic, ie = g('chs_A0A2Z2PGA9', 'ident'), g('EEU41297.1', 'ident')
        mx = max(g('chs_A0A2Z2PGA9', 'cov'), g('EEU41297.1', 'cov'), g('KNB14377.1', 'cov'))
        cls = 'chs_like' if (ic >= .45 or ie >= .70) and mx >= .25 else 'borderline' if (ic >= .35 or ie >= .55) and mx >= .25 else 'ags_like_or_other'
        prim = best.get('EEU41297.1') or max(best.values(), key=lambda h: h['score'])
        status = 'pseudogene' if (prim['fs'] or prim['stop']) else 'intact' if prim['cov'] >= .9 else 'truncated'
        if cls == 'ags_like_or_other': status = 'NA'
        lid = f'{gca}_L{k}'
        rows.append(dict(locus_id=lid, gca=gca, organism=gl[gca]['organism'], contig=c['contig'], start=c['start'], end=c['end'],
                         pid_to_chs_seed=ic, cov_chs=round(g('chs_A0A2Z2PGA9', 'cov'), 3), pid_to_vanettenii_EEU41297=ie, cov_EEU=round(g('EEU41297.1', 'cov'), 3),
                         pid_to_KNB14377=g('KNB14377.1', 'ident'), max_cov=round(mx, 3), cls=cls, status=status,
                         in_loci_all=int(lid in old and old[lid]['cls'] == cls)))
H = list(rows[0])
with open(f'{O}/chs_loci_miniprot_all.tsv', 'w') as f:
    f.write('\t'.join(H) + '\n')
    for r in sorted(rows, key=lambda r: -r['pid_to_chs_seed']): f.write('\t'.join(str(r[h]) for h in H) + '\n')
print('clusters', len(rows), collections.Counter(r['cls'] for r in rows), 'chs_like consistent w/ loci_all:', sum(r['in_loci_all'] for r in rows if r['cls']=='chs_like'), len(old))
