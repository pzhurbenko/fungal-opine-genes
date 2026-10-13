import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
# genes following the rel -2 ortholog (direction away from rel -3 ortholog) in the 4 FSSC targets
import csv, re, collections
G = f'{GMO_ROOT}/data/fungal_genomes/fusarium/fus/ncbi_dataset/data'
MOB = re.compile(r'transpos|integrase|retro|helicase|reverse transcriptase', re.I)
fl = list(csv.DictReader(open('step6_vanettenii_flanks.tsv'), delimiter='\t'))
out = open('step6_vanettenii_next_genes.tsv', 'w'); out.write('target\tcontig\torient\tidx\tgene\tstart\tend\tproduct\texons\tdist_from_rel-2_ortholog\tmobile\n')
for T in sorted({r['target'] for r in fl}):
    p = {int(r['rel']): r for r in fl if r['target'] == T and r['mutual'] == 'mutual'}
    a, b = p[-2], p[-3]; c = a['contig']; sgn = 1 if int(a['start']) > int(b['start']) else -1
    m = {}
    for l in open(f'{G}/{T}/genomic.gff'):
        if l[0] == '#': continue
        f = l.split('\t')
        if len(f) < 9 or f[0] != c: continue
        at = dict(kv.split('=', 1) for kv in f[8].strip().split(';') if '=' in kv)
        if f[2] in ('gene', 'pseudogene'): m[at['ID']] = dict(name=at.get('Name', at['ID']), s=int(f[3]), e=int(f[4]))
    # simpler second pass for product/exons
    prod = {}; ex = collections.Counter(); par = {}
    for l in open(f'{G}/{T}/genomic.gff'):
        f = l.split('\t')
        if len(f) < 9 or f[0] != c: continue
        at = dict(kv.split('=', 1) for kv in f[8].strip().split(';') if '=' in kv)
        if f[2] in ('mRNA', 'transcript'): par[at['ID']] = at['Parent']
        elif f[2] == 'CDS':
            g = par.get(at.get('Parent'), at.get('Parent')); ex[g] += 1; prod[g] = at.get('product', '')
    L = sorted(m.items(), key=lambda kv: kv[1]['s'])
    ref = int(a['end']) if sgn > 0 else int(a['start'])
    nxt = [x for x in L if (x[1]['s'] > ref if sgn > 0 else x[1]['e'] < ref)]
    if sgn < 0: nxt = nxt[::-1]
    for i, (gid, g) in enumerate(nxt[:5], 1):
        d = g['s'] - ref if sgn > 0 else ref - g['e']
        out.write('\t'.join(map(str, [T, c, '+' if sgn > 0 else '-', i, g['name'], g['s'], g['e'], prod.get(gid, ''), ex[gid], d, int(bool(MOB.search(prod.get(gid, ''))))])) + '\n')
