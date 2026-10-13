#!/usr/bin/env python3
import sys
fam, model = sys.argv[1], sys.argv[2]
M = {'M0': (0, 0, '0'), 'M7M8': (0, 0, '7 8'), 'M8a': (0, 0, '8')}
extra = '' if model != 'M8a' else 'fix_omega = 1\nomega = 1\n'
nsm = {'M0': '0', 'M7M8': '7 8', 'M8a': '8'}[model]
open(f'{fam}.{model}.ctl', 'w').write(f'''seqfile = {fam}.codon.phy
treefile = {fam}.tree.nwk
outfile = {fam}.{model}.out
noisy = 0
verbose = 0
runmode = 0
seqtype = 1
CodonFreq = 2
model = 0
NSsites = {nsm}
icode = 0
fix_kappa = 0
kappa = 2
{extra if extra else "fix_omega = 0"}
{"" if extra else "omega = 0.4"}
cleandata = 0
''')
