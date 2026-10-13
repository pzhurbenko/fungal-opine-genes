#!/usr/bin/env python3
"""Plant cT-DNA vs fungal rule-B proteins: identity to bacterial T-DNA seeds of the same gene.
Inputs: queries.faa (headers gene|id|group|origin), bact_seeds.faa (seeds_ext bacterial), Ti/Ri plasmid proteome.
Outputs: pairwise_hits.tsv (per query best same-gene seed + canonical seed11 + best plasmid protein), summary.tsv."""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import subprocess, csv, statistics, collections, os
D = f'{GMO_ROOT}/analysis/ctdna_plants_identity'
DM = 'diamond'
PL = f'{GMO_ROOT}/analysis/refset/plasmids/all_proteins.faa'
os.chdir(D)
FMT = 'qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore'.split()

def run(db_fa, out):
    subprocess.run([DM, 'makedb', '--in', db_fa, '-d', out + '.db', '--quiet'], check=True)
    subprocess.run([DM, 'blastp', '-q', 'queries.faa', '-d', out + '.db', '-o', out + '.m8', '--very-sensitive',
                    '-k', '500', '--max-hsps', '1', '-e', '1e-3', '-p', '4', '--quiet', '--outfmt', '6'] + FMT, check=True)
    rows = collections.defaultdict(list)
    for l in open(out + '.m8'):
        r = dict(zip(FMT, l.rstrip('\n').split('\t')))
        for k in ('pident', 'bitscore'): r[k] = float(r[k])
        for k in ('length', 'qlen', 'slen', 'qstart', 'qend', 'sstart', 'send'): r[k] = int(r[k])
        r['qcov'] = (r['qend'] - r['qstart'] + 1) / r['qlen']; r['scov'] = (r['send'] - r['sstart'] + 1) / r['slen']
        rows[r['qseqid']].append(r)
    return rows

norm = lambda g: 'acs' if g.startswith('acs') else g
seedhits = run('bact_seeds.faa', 'vs_seeds')
plhits = run(PL, 'vs_plasmids')
with open('plant_refs.faa', 'w') as o:
    keep = False
    for l in open('queries.faa'):
        if l.startswith('>'): keep = '|plant_cTDNA|' in l
        if keep: o.write(l)
planthits = run('plant_refs.faa', 'vs_plants')
EXCL = {'KHG03718.1': 'native SDR homolog (31%, partial), not cT-DNA', 'XP_044455008.1': 'native SDR homolog (31%), not cT-DNA',
        'KAB1203964.1': 'native ocs/vis-family homolog (26%), not cT-DNA',
        'GMD07831.1': 'only 32% of query aligned to mas1 (chimeric/partial model)'}
seed11 = {l.split('\t')[0]: l.split('\t')[1] for l in open(f'{GMO_ROOT}/analysis/refset/screen/seeds11/seeds11.tsv') if not l.startswith('gene')}
queries = [l[1:].strip() for l in open('queries.faa') if l.startswith('>')]
seqs = {}; _k = None
for l in open('queries.faa'):
    if l.startswith('>'): _k = l[1:].split('|')[1]; seqs[_k] = ''
    else: seqs[_k] += l.strip()

out = []
for q in queries:
    gene, qid, group, origin = q.split('|', 3)
    hs = sorted(seedhits.get(q, []), key=lambda r: -r['bitscore'])
    top = hs[0] if hs else None
    if gene == 'auto':
        gene = norm(top['sseqid'].split('|')[0]) if top else 'none'
    same = [r for r in hs if norm(r['sseqid'].split('|')[0]) == gene]
    b = same[0] if same else None
    c = next((r for r in hs if r['sseqid'].split('|')[1] == seed11.get(gene)), None)
    p = sorted(plhits.get(q, []), key=lambda r: -r['bitscore'])
    p = p[0] if p else None
    ph = [r for r in sorted(planthits.get(q, []), key=lambda r: -r['bitscore']) if r['sseqid'] != q and r['sseqid'].split('|')[0] in (gene, 'auto')] if group == 'fungal_ruleB' else []
    ph = ph[0] if ph else None
    out.append(dict(group=group, gene=gene, query=qid, status=('exclude: ' + EXCL[qid]) if qid in EXCL else 'include', origin=origin, qlen=(b or top or {}).get('qlen', ''),
        best_seed=b['sseqid'].split('|')[1] if b else '', pident=round(b['pident'], 1) if b else '',
        qcov=round(b['qcov'], 3) if b else '', scov=round(b['scov'], 3) if b else '', bits=b['bitscore'] if b else '',
        evalue=b['evalue'] if b else '',
        seed11=seed11.get(gene, ''), pident_seed11=round(c['pident'], 1) if c else '', qcov_seed11=round(c['qcov'], 3) if c else '',
        top_any_seed=top['sseqid'].split('|')[0] + '|' + top['sseqid'].split('|')[1] if top else '',
        top_any_pident=round(top['pident'], 1) if top else '',
        best_plasmid=p['sseqid'].split('|')[0] + '|' + '|'.join(p['sseqid'].split('|')[2:4]) if p else '',
        pident_plasmid=round(p['pident'], 1) if p else '', qcov_plasmid=round(p['qcov'], 3) if p else '',
        best_plant_cTDNA=ph['sseqid'].split('|')[1] if ph else '', pident_plant_cTDNA=round(ph['pident'], 1) if ph else '',
        qcov_plant_cTDNA=round(ph['qcov'], 3) if ph else ''))
# second pass: fungal best hit to included plant cT-DNA proteins of the SAME resolved gene
pg = {r['query']: r['gene'] for r in out if r['group'] == 'plant_cTDNA' and r['status'] == 'include'}
for r in out:
    if r['group'] != 'fungal_ruleB': continue
    q = next(x for x in queries if x.split('|')[1] == r['query'])
    hs = [h for h in sorted(planthits.get(q, []), key=lambda h: -h['bitscore']) if pg.get(h['sseqid'].split('|')[1]) == r['gene']]
    h = hs[0] if hs else None
    r['best_plant_cTDNA'] = h['sseqid'].split('|')[1] if h else ''
    r['pident_plant_cTDNA'] = round(h['pident'], 1) if h else ''
    r['qcov_plant_cTDNA'] = round(h['qcov'], 3) if h else ''
with open('pairwise_hits.tsv', 'w') as f:
    w = csv.DictWriter(f, fieldnames=list(out[0]), delimiter='\t'); w.writeheader(); w.writerows(out)

# summary: by gene x group; 'all' and 'qcov>=0.5' subsets
def stats(v):
    return (len(v), round(statistics.median(v), 1), min(v), max(v)) if v else (0, '', '', '')
summ = []
for gene in sorted({r['gene'] for r in out}):
    for group in ('plant_cTDNA', 'fungal_ruleB'):
        rs = [r for r in out if r['gene'] == gene and r['group'] == group and r['pident'] != '' and r['status'] == 'include']
        if not rs: continue
        a = stats([r['pident'] for r in rs]); f = stats([r['pident'] for r in rs if r['scov'] >= 0.7])
        seen = set(); ru = [r for r in rs if not (seqs[r['query']] in seen or seen.add(seqs[r['query']]))]; u = stats([r['pident'] for r in ru])
        s = stats([r['pident_seed11'] for r in rs if r['pident_seed11'] != ''])
        summ.append(dict(gene=gene, group=group, n=a[0], median_pident=a[1], min=a[2], max=a[3],
            median_qcov=round(statistics.median([r['qcov'] for r in rs]), 2),
            n_fulllen_scov70=f[0], median_pident_fulllen=f[1], min_fulllen=f[2], max_fulllen=f[3],
            n_unique_seq=u[0], median_pident_unique=u[1],
            median_pident_seed11=s[1], range_seed11=f'{s[2]}-{s[3]}' if s[0] else '',
            median_pident_to_plant_cTDNA=(round(statistics.median([r['pident_plant_cTDNA'] for r in rs if r['pident_plant_cTDNA'] != '']), 1) if any(r['pident_plant_cTDNA'] != '' for r in rs) else ''),
            queries=','.join(r['query'] for r in rs)))
with open('summary.tsv', 'w') as f:
    w = csv.DictWriter(f, fieldnames=list(summ[0]), delimiter='\t'); w.writeheader(); w.writerows(summ)
print(open('summary.tsv').read())
