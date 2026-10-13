#!/usr/bin/env python3
"""vir/tra/rep(/mobile) loci at >=60% identity on contigs of the Boletales assemblies listed in contamination_contigs.tsv (audit A2, 2026-10-13).
Same chaining/locus rules as summarize_517.py (GAP 5000; locus best subject = max bits); lists each such locus with its best subject.
Output: vir_tra_rep_hits.tsv (one row per locus)."""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import csv, re
from collections import defaultdict
D = f'{GMO_ROOT}/analysis/refset/screen/step2'; OUT = f'{D}/boletales517_tdna_summary'
GAP = 5000; BACK = ('vir', 'tra', 'rep_partition', 'mobile')
def scov_of(iv, slen):
    mg = []
    for a, e in sorted(iv):
        if mg and a <= mg[-1][1]: mg[-1][1] = max(mg[-1][1], e)
        else: mg.append([a, e])
    return sum(e - a + 1 for a, e in mg) / slen
rows = []
for cc in csv.DictReader(open(f'{OUT}/contamination_contigs.tsv'), delimiter='\t'):
    gca, tgt = cc['gca'], cc['contig']
    hs = defaultdict(list)
    for l in open(f'{D}/out/{gca}/dmnd.tsv'):
        f = l.rstrip('\n').split('\t')
        if len(f) < 13 or f[0] != tgt: continue
        qs, qe = sorted((int(f[6]), int(f[7])))
        hs[('+' if int(f[12]) > 0 else '-', f[1])].append((qs, qe, int(f[8]), int(f[9]), float(f[2]), float(f[11]), int(f[5]), int(f[4])))
    chains = defaultdict(list)
    for (s, sub), hl in hs.items():
        hl.sort(); cur = None
        for qs, qe, ss, se, pid, bits, slen, clen in hl:
            if cur and qs - cur['end'] <= GAP:
                cur['end'] = max(cur['end'], qe); cur['bits'] += bits; cur['scov'].append((min(ss, se), max(ss, se))); cur['pid'] = max(cur['pid'], pid)
            else:
                if cur: chains[s].append(cur)
                cur = dict(start=qs, end=qe, subj=sub, bits=bits, scov=[(min(ss, se), max(ss, se))], pid=pid, slen=slen)
        chains[s].append(cur)
    for s, cl in chains.items():
        cl.sort(key=lambda c: c['start']); loci = []; cur = None
        for c in cl:
            c['sc'] = scov_of(c['scov'], c['slen'])
            if cur and c['start'] <= cur['end']: cur['end'] = max(cur['end'], c['end']); cur['m'].append(c)
            else:
                if cur: loci.append(cur)
                cur = dict(start=c['start'], end=c['end'], m=[c])
        loci.append(cur)
        for L in loci:
            b = max(L['m'], key=lambda x: x['bits']); p = b['subj'].split('|')
            if p[1] in BACK and b['pid'] >= 60:
                rows.append(dict(gca=gca, organism=cc['organism'], contig=tgt, contig_len=cc['contig_len'], start=L['start'], end=L['end'], strand=s,
                                 subject=b['subj'], category=p[1], pident=f"{b['pid']:.1f}", subj_cov=f"{b['sc']:.2f}",
                                 genome_is_carrier=cc['genome_is_carrier'], contig_has_carrier_locus=cc['contig_has_carrier_locus']))
rows.sort(key=lambda r: (r['gca'], r['start']))
with open(f'{OUT}/vir_tra_rep_hits.tsv', 'w') as f:
    ks = list(rows[0].keys()); f.write('\t'.join(ks) + '\n')
    for r in rows: f.write('\t'.join(str(r[k]) for k in ks) + '\n')
print(len(rows))
