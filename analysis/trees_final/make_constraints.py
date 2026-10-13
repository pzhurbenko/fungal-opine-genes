#!/usr/bin/env python3
"""Constraint trees for the AU test: H_tdna = the target fungal lineage is monophyletic WITH the T-DNA references;
H_alt = the same fungal lineage is monophyletic with the alternative donor / native paralogs instead.
Rejecting H_alt (and not H_tdna) is the statistical statement 'these fungal genes sit inside the T-DNA clade'.
usage: make_constraints.py <FAM> <alignment.fasta>   -> <FAM>.H_tdna.constraint.nwk, <FAM>.H_alt.constraint.nwk"""
import sys, os, csv
FAM, ALN = sys.argv[1], sys.argv[2]
HERE=os.path.dirname(os.path.abspath(ALN)) or '.'
names=[l[1:].split()[0] for l in open(ALN) if l.startswith('>')]
lab={r['label']:r for r in csv.DictReader(open(f'{HERE}/{FAM}_labels.tsv'), delimiter='\t')}
TDNA={'ACS':['acs'],'CUS':['cus','mis','rcorf12','cusrel181'],'SUS':['sus'],'ORF358':['orf358']}[FAM]
TARGET={'ACS':['acs','gdpd_meso'],'CUS':['cus','mis'],'SUS':['sus'],'ORF358':['orf358']}[FAM]
def sub(n): return lab.get(n,{}).get('subgroup','')
fungi=[n for n in names if n.startswith('FUNGI_') and sub(n) in TARGET]
refs_tdna=[n for n in names if n.startswith('REF_') and any(f'_{t}__' in n for t in TDNA)]
alt=[n for n in names if n.startswith('OG__') or n.startswith('NATIVE_') or n.startswith('BACT_') or
     (n.startswith('REF_') and not any(f'_{t}__' in n for t in TDNA))]
rest=[n for n in names if n not in set(fungi)|set(refs_tdna)|set(alt)]
def grp(x): return '(' + ','.join(x) + ')'
out={}
if fungi and refs_tdna: out['H_tdna']=f"({grp(fungi+refs_tdna)},{','.join(alt+rest)});"
if fungi and alt:       out['H_alt'] =f"({grp(fungi+alt)},{','.join(refs_tdna+rest)});"
for k,v in out.items():
    open(f'{HERE}/{FAM}.{k}.constraint.nwk','w').write(v+'\n')
print(f'{FAM}: fungi {len(fungi)}, T-DNA refs {len(refs_tdna)}, alt {len(alt)}, rest {len(rest)} -> {list(out)}', file=sys.stderr)
