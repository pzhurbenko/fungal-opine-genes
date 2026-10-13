#!/usr/bin/env python3
"""Threshold tables. Inputs: chs_identity.tsv (miniprot loci, from chs_sensitivity.py) and
diamond_hits_to_chs_seed.raw.txt (= grep 'chs|A0A2Z2PGA9' refset/screen/seeds_ext/out/*.tsv; all fungal proteins
hitting the chs seed in the annotated-proteome DIAMOND screen, --very-sensitive -e 1e-3)."""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import csv, collections
O = f'{GMO_ROOT}/analysis/thresholds/chs'
org = {r['gca']: r['organism'] for r in csv.DictReader(open(f'{GMO_ROOT}/analysis/refset/screen/seeds_ext/seeds_ext_best.tsv'), delimiter='\t')}
FUS = ('Fusarium', 'Neocosmospora', 'Dactylonectria')
loci = list(csv.DictReader(open(f'{O}/chs_loci_miniprot_all.tsv'), delimiter='\t'))
for r in loci:
    for k in ('pid_to_chs_seed', 'pid_to_vanettenii_EEU41297', 'max_cov'): r[k] = float(r[k])
dia = []
for l in open(f'{O}/diamond_hits_to_chs_seed.raw.txt'):
    fn, rest = l.split(':', 1); f = rest.rstrip('\n').split('\t'); g = fn.split('/')[-1][:-4]; o = org.get(g, '?')
    dia.append(dict(source='diamond_protein', id=f[0], gca=g, organism=o, pid_to_chs_seed=float(f[2]) / 100, cov_seed=round((int(f[9]) - int(f[8]) + 1) / int(f[5]), 3),
                    aln_len=f[3], qlen=f[4], bits=f[11], cls='chs_lineage(Fusarium/Neocosmospora/Dactylonectria)' if o.split()[0] in FUS else 'native_fungal_homolog(ags-like/other)'))
# combined identity table (overwrite chs_identity.tsv with both sources; keep miniprot columns)
H = ['source', 'id', 'gca', 'organism', 'pid_to_chs_seed', 'cov_seed', 'pid_to_vanettenii_EEU41297', 'max_cov', 'cls', 'status', 'note']
out = []
for r in loci:
    out.append(dict(source='miniprot_locus', id=r['locus_id'], gca=r['gca'], organism=r['organism'], pid_to_chs_seed=r['pid_to_chs_seed'], cov_seed=r['cov_chs'],
                    pid_to_vanettenii_EEU41297=r['pid_to_vanettenii_EEU41297'], max_cov=r['max_cov'], cls=r['cls'], status=r['status'],
                    note='no chs-seed alignment reported by miniprot (id=0); chs_like via id_EEU>=0.70' if r['pid_to_chs_seed'] == 0 else ''))
for d in dia: out.append(dict(d, pid_to_vanettenii_EEU41297='', max_cov='', status='', note=f"aln {d['aln_len']} aa, qlen {d['qlen']}, bits {d['bits']}"))
with open(f'{O}/chs_identity.tsv', 'w') as f:
    f.write('\t'.join(H) + '\n')
    for r in sorted(out, key=lambda r: (r['source'] != 'miniprot_locus', -r['pid_to_chs_seed'])): f.write('\t'.join(str(r.get(h, '')) for h in H) + '\n')

def chs_like(r, tc=.45, te=.70): return (r['pid_to_chs_seed'] >= tc or r['pid_to_vanettenii_EEU41297'] >= te) and r['max_cov'] >= .25
base = {r['locus_id'] for r in loci if chs_like(r)}
with open(f'{O}/chs_thresholds.tsv', 'w') as f:
    f.write('threshold_chs\tmp_loci_chs_like(rule as implemented, te=0.70)\tmp_genomes_chs_like\tmp_loci_lost_vs_0.45\tmp_loci_with_id_chs>=T(alone)\t'
            'dmnd_proteins>=T(all)\tdmnd_chs_lineage>=T\tdmnd_native>=T\tdmnd_genomes>=T\tdmnd_native_names>=T\n')
    for T in (.30, .35, .40, .45, .50, .55, .60):
        s = [r for r in loci if chs_like(r, tc=T)]; lost = sorted(base - {r['locus_id'] for r in s})
        d = [x for x in dia if x['pid_to_chs_seed'] >= T]; nat = [x for x in d if x['cls'].startswith('native')]
        f.write(f"{int(T*100)}\t{len(s)}\t{len({r['gca'] for r in s})}\t{','.join(lost) or '-'}\t{sum(r['pid_to_chs_seed'] >= T for r in loci)}\t"
                f"{len(d)}\t{len(d)-len(nat)}\t{len(nat)}\t{len({x['gca'] for x in d})}\t{','.join(sorted({x['organism'] for x in nat})) or '-'}\n")
# vanettenii
dep = [r for r in loci if r['pid_to_chs_seed'] < .45 and r['cls'] == 'chs_like']  # loci that are chs_like only thanks to id_EEU>=0.70
with open(f'{O}/vanettenii_identity.tsv', 'w') as f:
    f.write('locus_id\tgca\torganism\tpid_to_vanettenii_EEU41297\tcov_EEU\tpid_to_chs_seed\tstatus\tdepends_on_70pct_rule\n')
    for r in sorted(loci, key=lambda r: -r['pid_to_vanettenii_EEU41297']):
        f.write(f"{r['locus_id']}\t{r['gca']}\t{r['organism']}\t{r['pid_to_vanettenii_EEU41297']}\t{r['cov_EEU']}\t{r['pid_to_chs_seed']}\t{r['status']}\t{int(r in dep)}\n")
with open(f'{O}/vanettenii_thresholds.tsv', 'w') as f:
    f.write('threshold_vanettenii\tmp_loci_chs_like\tdependent_loci_kept(of %d)\tdependent_kept_by_status\tloci_lost_vs_0.70\n' % len(dep))
    for T in (.55, .60, .65, .70, .75, .80, .85, .90):
        s = [r for r in loci if chs_like(r, te=T)]; ids = {r['locus_id'] for r in s}; lost = sorted(base - ids)
        kept = [r for r in dep if r['locus_id'] in ids]
        f.write(f"{int(T*100)}\t{len(s)}\t{len(kept)}\t{dict(collections.Counter(r['status'] for r in kept))}\t{','.join(lost) or '-'}\n")
print('dep', len(dep), collections.Counter(r['status'] for r in dep), sorted({r['pid_to_vanettenii_EEU41297'] for r in dep}))
print('all ie', sorted(r['pid_to_vanettenii_EEU41297'] for r in loci))
