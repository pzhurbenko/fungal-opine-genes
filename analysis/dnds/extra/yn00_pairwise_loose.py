import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
#!/usr/bin/env python3
# yn00 drops every codon column with a gap in ANY sequence, so yn00 on the <=20%-gap alignment equals the strict one.
# Here: yn00 per pair on the <=20%-gap alignment with pairwise deletion (only columns gapped in that pair removed).
import sys, os, re, subprocess, itertools, tempfile
s = sys.argv[1]; d = f'{GMO_ROOT}/analysis/dnds/extra/{s}/loose'
seqs = []; name = None
for l in open(f'{d}/{s}.codon.fna'):
    l = l.strip()
    if l.startswith('>'): name = l[1:]; seqs.append([name, ''])
    else: seqs[-1][1] += l
out = open(f'{d}/{s}.yn00_pairwise.tsv', 'w'); out.write('seq1\tseq2\tcodons\tomega\tdN\tdS\n')
with tempfile.TemporaryDirectory() as tmp:
    for (a, x), (b, y) in itertools.combinations(seqs, 2):
        keep = [i for i in range(0, len(x), 3) if '-' not in x[i:i+3] + y[i:i+3] and 'N' not in x[i:i+3] + y[i:i+3]]
        X = ''.join(x[i:i+3] for i in keep); Y = ''.join(y[i:i+3] for i in keep)
        open(f'{tmp}/p.phy', 'w').write(f' 2 {len(X)}\nA  {X}\nB  {Y}\n')
        open(f'{tmp}/p.ctl', 'w').write('seqfile = p.phy\noutfile = p.out\nverbose = 0\nicode = 0\nweighting = 0\ncommonf3x4 = 0\n')
        subprocess.run(['yn00', 'p.ctl'], cwd=tmp, capture_output=True)
        m = re.search(r'^\s*2\s+1\s+[\d.]+\s+[\d.]+\s+[-\d.]+\s+[-\d.]+\s+([-\d.]+)\s+([-\d.]+)\s+\+-\s+[\d.]+\s+([-\d.]+)', open(f'{tmp}/p.out').read(), re.M)
        out.write(f'{a}\t{b}\t{len(keep)}\t' + ('\t'.join(m.groups()) if m else '\t\t') + '\n')
