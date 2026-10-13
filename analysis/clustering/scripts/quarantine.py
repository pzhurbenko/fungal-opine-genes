#!/usr/bin/env python3
"""Move flagged sequences out of clean.fasta into a quarantine, write curated inputs.

Input : analysis/clustering/<FAM>/clean.fasta (sanitised, deduplicated)
Output: analysis/curated/<FAM>.fasta          — clean.fasta minus quarantined ids
        analysis/quarantine/<FAM>_quarantine.fasta — removed sequences
        analysis/quarantine/quarantine.tsv     — ledger (append-only): who/why/when
"""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import os, sys, datetime
BASE = f'{GMO_ROOT}/analysis'
FAMS = ['ACS', 'CUS', 'ORF358', 'SUS']

# (family, substring uniquely identifying the clean_id, reason)
QUARANTINE = [
    ('SUS', 'JBBPRU010007955',
     'False-positive tblastn hit, not a SUS homolog. tblastn E=0.035, id 33%, query cover 32%, '
     'caught only by plant query SUS_Camellia (no Agro query); singleton in user Tree_clade and in mmseqs -c0.6; '
     'Pfam sensitive scan: 13 WD40/beta-propeller models over aa 1-118 (best E 1e-7), no OCD/Shikimate_DH. '
     'Best agro hit in rep_seqs = TIR-domain protein (unrelated). 122 aa WD40 repeat fragment.'),
    ('SUS', 'JAKRWS020000013',
     'Native fungal OCD/mu-crystallin-family protein, not T-DNA sus. tblastn E=0.0065, id 26%, cover 44% (plant query SUS_Camellia only). '
     'Best agro hit = WP_151986804 chromosomal delta1-pyrroline-2-carboxylate reductase (Rhizobium sp.), NOT the pRi T-DNA sus '
     '(WP_269710210 etc.) that all other SUS reps hit; fragment 41% id to that reductase vs 26% to T-DNA sus. '
     'Pfam sensitive: OCD_Mu_crystall + Shikimate_DH + F420_oxidored + Sacchrp_dh_NADP (architecture absent in any SUS). '
     'miniprot: none of 23 Agro/Plant SUS guides align to the 19.5 kb window; AUGUSTUS (ustilago) predicts no gene. Confirmed by user 2026-09-05.'),
]

def read_fasta(p):
    recs, h, s = [], None, []
    for line in open(p):
        line = line.rstrip('\n')
        if line.startswith('>'):
            if h is not None: recs.append((h, ''.join(s)))
            h, s = line[1:], []
        else:
            s.append(line)
    if h is not None: recs.append((h, ''.join(s)))
    return recs

os.makedirs(f'{BASE}/curated', exist_ok=True)
os.makedirs(f'{BASE}/quarantine', exist_ok=True)
ledger = f'{BASE}/quarantine/quarantine.tsv'
new_ledger = not os.path.exists(ledger)
today = datetime.date.today().isoformat()
already = set()
if not new_ledger:  # idempotent: never re-append an id that is already in the ledger
    for line in open(ledger):
        p = line.rstrip('\n').split('\t')
        if len(p) >= 3 and p[0] != 'date': already.add(p[2])

with open(ledger, 'a') as L:
    if new_ledger:
        L.write('date\tfamily\tclean_id\tlength\treason\n')
    for fam in FAMS:
        recs = read_fasta(f'{BASE}/clustering/{fam}/clean.fasta')
        rules = [(k, r) for f, k, r in QUARANTINE if f == fam]
        keep, drop = [], []
        for h, s in recs:
            hit = [r for k, r in rules if k in h]
            (drop if hit else keep).append((h, s, hit[0] if hit else None))
        # every rule must match exactly one record
        for k, _ in rules:
            n = sum(1 for h, _, _ in drop if k in h)
            assert n == 1, f'{fam}: rule {k} matched {n} records'
        with open(f'{BASE}/curated/{fam}.fasta', 'w') as out:
            for h, s, _ in keep: out.write(f'>{h}\n{s}\n')
        if drop:
            with open(f'{BASE}/quarantine/{fam}_quarantine.fasta', 'w') as q:
                for h, s, r in drop:
                    q.write(f'>{h} QUARANTINED {today}: {r}\n{s}\n')
                    if h not in already:
                        L.write(f'{today}\t{fam}\t{h}\t{len(s)}\t{r}\n')
        print(f'{fam}: in={len(recs)} kept={len(keep)} quarantined={len(drop)}')
