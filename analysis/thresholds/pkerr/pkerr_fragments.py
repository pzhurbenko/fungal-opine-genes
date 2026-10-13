#!/usr/bin/env python3
"""Sensitivity of the '>=60% identity to the six tree-placed pKerr14a-type members' rule
(paper/tables/build_table1.py:145 reads the 'ACS kerr14' label from analysis/report_A/tables/carriers_tree_based.tsv).
Recomputed: DIAMOND 2.2.6 blastp --very-sensitive -e 1e-3 --max-target-seqs 6 --max-hsps 1, all 332 ACS proteins
(loci_pipeline/round2/proteins/ACS_models.faa + ACS_remnants_fragments.faa) vs the six -> acs_vs_six.diamond.tsv."""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import csv, collections
P = f'{GMO_ROOT}/analysis/thresholds/pkerr'
SIX = ['ACS_JBBPKL010003317_L1','ACS_JBBPNB010000771_L1','ACS_JBBPKL010018021_L1','ACS_JBBPOG010012756_L1','ACS_JBBPNB010000203_L1','ACS_JBBPRY010001778_L1']
cat = {r['locus_id']: r for r in csv.DictReader(open(f'{GMO_ROOT}/paper/tables/TableS1_loci_catalog.tsv'), delimiter='\t')}
hits = collections.defaultdict(list)
for l in open(f'{P}/acs_vs_six.diamond.tsv'):
    f = l.rstrip('\n').split('\t')
    if f[0] == f[1]: continue
    hits[f[0]].append(dict(s=f[1], pid=float(f[2]), alen=int(f[3]), qlen=int(f[4]), qcov=round((int(f[7])-int(f[6])+1)/int(f[4]), 3), bits=float(f[11])))
rows = []
for q, hs in hits.items():
    hs.sort(key=lambda h: (-h['bits'], -h['pid']))   # nearest = DIAMOND top hit (highest bitscore); max pid over six reported separately
    c = cat.get(q, {}); b1 = hs[0]; b2 = hs[1] if len(hs) > 1 else None
    hb = max(hs, key=lambda h: h['pid'])
    rows.append(dict(locus_id=q, is_six=int(q in SIX), lineage_final=c.get('lineage_final', '?'), basis=c.get('basis', '?'), status=c.get('status', '?'), organism=c.get('organism', '?'),
                     qlen=b1['qlen'], nearest_intact=b1['s'], pid1=b1['pid'], alen1=b1['alen'], qcov1=b1['qcov'],
                     second_intact=b2['s'] if b2 else '', pid2=b2['pid'] if b2 else '', bits1=b1['bits'], second_qcov=b2['qcov'] if b2 else '', bits2=b2['bits'] if b2 else '', maxpid_intact=hb['s'], maxpid=hb['pid'], maxpid_alen=hb['alen'], n_six_hit=len(hs)))
rows.sort(key=lambda r: -r['pid1'])
H = list(rows[0])
with open(f'{P}/fragment_identity.tsv', 'w') as f:
    f.write('\t'.join(H) + '\n')
    for r in rows: f.write('\t'.join(str(r[h]) for h in H) + '\n')
orig = {k for k, c in cat.items() if c['basis'] == 'similarity' and c['lineage_final'] == 'pKerr14a-type acs-like'}
nonsix = [r for r in rows if not r['is_six']]
with open(f'{P}/fragment_thresholds.tsv', 'w') as f:
    f.write('threshold\tn_assigned(non-six, any basis)\tn_assigned_not_tree_placed\tof_original13_kept\toriginal13_lost\tnew_not_in_original13(locus:lineage_final:pid)\tnearest_intact_same_as_at_60\tALT_maxpid_def:n_assigned\tALT_maxpid_def:original13_lost\tALT_maxpid_def:new\n')
    for T in (40, 45, 50, 55, 60, 65, 70, 75, 80):
        a = [r for r in nonsix if r['pid1'] >= T]; b = [r for r in nonsix if r['maxpid'] >= T]; nt = [r for r in a if r['basis'] != 'tree']
        ids = {r['locus_id'] for r in a}
        new = [f"{r['locus_id']}:{r['lineage_final']}:{r['pid1']}" for r in a if r['locus_id'] not in orig]
        f.write(f"{T}\t{len(a)}\t{len(nt)}\t{len(orig & ids)}\t{','.join(sorted(orig - ids)) or '-'}\t{';'.join(new) or '-'}\tyes(nearest fixed per fragment)\t{len(b)}\t{','.join(sorted(orig - {r['locus_id'] for r in b})) or '-'}\t{';'.join(r['locus_id']+':'+str(r['maxpid']) for r in b if r['locus_id'] not in orig) or '-'}\n")
print('orig13 present in diamond:', len(orig & set(hits)), 'of', len(orig))
for r in rows:
    if not r['is_six'] and (r['locus_id'] in orig or r['pid1'] >= 40):
        print(r['locus_id'], r['lineage_final'][:22], r['basis'], r['status'], r['qlen'], r['nearest_intact'][4:], r['pid1'], r['alen1'], r['second_intact'][4:], r['pid2'], r['qcov1'], 'maxpid', r['maxpid_intact'][4:], r['maxpid'], r['maxpid_alen'])
nb = sorted(r['pid1'] for r in nonsix if r['locus_id'] not in orig)
print('non-assigned max pid1 (top5):', nb[-5:], 'n', len(nb))
