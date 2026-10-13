#!/usr/bin/env python3
"""Summarise backbone signal per contig and the fraction flagged at 40/50/60/70/80 %.

Backbone categories = vir, tra, rep_partition, mobile (as aggregate_step2.py, line 129).
(1) blastp (prodigal proteins vs 7,105-protein labelled catalogue): each protein takes its best hit (max bitscore);
    a protein is a backbone hit when that best hit is a backbone protein. best_backbone_pid = max identity over such
    proteins; n_backbone_hits_ge60 = number of such proteins with identity >= 60.
    best_any_bb_pid = max identity of ANY backbone subject among the top-25 hits (lenient variant).
(2) bx_* : exact re-implementation of the step-2 rule on DIAMOND blastx vs the 389 representatives:
    HSPs chained per (contig,strand,subject) with gap <= 5 kb (chain pid = max HSP pid), overlapping chains merged into
    loci, locus best member = max summed bits; contig flagged if a locus' best member is backbone with pid >= threshold.
Outputs: backbone_contigs.tsv, backbone_thresholds.tsv, backbone_scan_1pct.tsv (in analysis/thresholds/backbone/)
"""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import csv
from collections import defaultdict

B = f'{GMO_ROOT}/analysis/thresholds/backbone'; W = f'{B}/work'
BB = {'vir', 'tra', 'rep_partition', 'mobile'}
GAP = 5000
THR = [40, 50, 60, 70, 80]


def cat_of(s):
    p = s.split('|')
    return p[0][4:] if p[0].startswith('FAM_') else (p[1] if len(p) > 1 else '')


def label_of(s):
    p = s.split('|')
    return f"{p[1]}:{p[4] if len(p) > 4 and p[4] != '-' else (p[5][:40] if len(p) > 5 else '')}:{p[3] if len(p) > 3 else ''}"


man = list(csv.DictReader(open(f'{W}/manifest.tsv'), delimiter='\t'))
genes = {r['contig_key']: r for r in csv.DictReader(open(f'{W}/genes.tsv'), delimiter='\t')}

# ---- (1) blastp ----
best = {}; any_bb = defaultdict(float); any_bb_lab = {}
for l in open(f'{W}/blastp_vs_catalog.tsv'):
    f = l.rstrip('\n').split('\t'); q, s, pid, bits = f[0], f[1], float(f[2]), float(f[11])
    if q not in best or bits > best[q][2]: best[q] = (s, pid, bits)
    ck = q.rsplit('__', 1)[0]
    if cat_of(s) in BB and pid > any_bb[ck]: any_bb[ck] = pid; any_bb_lab[ck] = label_of(s)
per = defaultdict(lambda: dict(n_hit=0, bb=[]))
for q, (s, pid, bits) in best.items():
    ck = q.rsplit('__', 1)[0]; per[ck]['n_hit'] += 1
    if cat_of(s) in BB: per[ck]['bb'].append((pid, label_of(s)))

# ---- (2) blastx, step-2 locus logic ----
hsps = defaultdict(list)
for l in open(f'{W}/blastx_vs_q389.tsv'):
    f = l.rstrip('\n').split('\t')
    q, s = f[0], f[1]; pid = float(f[2]); qs, qe = sorted((int(f[6]), int(f[7]))); bits = float(f[11])
    hsps[(q, '+' if int(f[12]) > 0 else '-', s)].append((qs, qe, pid, bits))
chains = defaultdict(list)
for (q, st, s), hl in hsps.items():
    hl.sort(); cur = None
    for qs, qe, pid, bits in hl:
        if cur and qs - cur['end'] <= GAP:
            cur['end'] = max(cur['end'], qe); cur['bits'] += bits; cur['pid'] = max(cur['pid'], pid)
        else:
            if cur: chains[(q, st)].append(cur)
            cur = dict(start=qs, end=qe, subj=s, bits=bits, pid=pid)
    chains[(q, st)].append(cur)
bx = defaultdict(list)   # contig -> [(pid,label)] of loci whose best member is backbone
for (q, st), cl in chains.items():
    cl.sort(key=lambda c: c['start']); loci = []; cur = None
    for c in cl:
        if cur and c['start'] <= cur['end']: cur['end'] = max(cur['end'], c['end']); cur['m'].append(c)
        else:
            if cur: loci.append(cur)
            cur = dict(end=c['end'], m=[c])
    loci.append(cur)
    for L in loci:
        b = max(L['m'], key=lambda m: m['bits'])
        if cat_of(b['subj']) in BB: bx[q].append((b['pid'], label_of(b['subj'])))

rows = []
for m in man:
    ck = f"{m['sample']}__{m['contig']}"; g = genes[ck]; L = int(m['length']); ng = int(g['n_genes'])
    bb = per[ck]['bb']; top = max(bb) if bb else (0.0, '')
    xb = bx.get(ck, []); xtop = max(xb) if xb else (0.0, '')
    r = dict(sample=m['sample'], contig=m['contig'], source=m['source'], lineage=m['lineage'], gca=m['gca'],
             organism=m['organism'], length=L, GC=m['GC'], assembly_gc=m['assembly_gc'], n_genes=ng,
             genes_per_kb=f'{1000 * ng / L:.3f}', coding_frac=f"{int(g['coding_bp']) / L:.3f}",
             frac_prot_catalog_hit=f"{per[ck]['n_hit'] / ng:.3f}" if ng else '0',
             best_backbone_pid=f'{top[0]:.1f}', best_backbone_label=top[1],
             n_backbone_hits_ge60=sum(p >= 60 for p, _ in bb), n_backbone_hits_any=len(bb),
             best_any_bb_pid=f'{any_bb.get(ck, 0):.1f}', best_any_bb_label=any_bb_lab.get(ck, ''),
             bx_best_backbone_pid=f'{xtop[0]:.1f}', bx_best_backbone_label=xtop[1],
             bx_n_backbone_loci_ge60=sum(p >= 60 for p, _ in xb), annot_genes_on_contig=m['annot_genes_on_contig'])
    rows.append(r)
cols = list(rows[0])
with open(f'{B}/backbone_contigs.tsv', 'w') as fo:
    w = csv.DictWriter(fo, cols, delimiter='\t'); w.writeheader(); w.writerows(rows)

# unique contigs per sample for threshold fractions
uniq = {}
for r in rows: uniq[(r['sample'], r['contig'])] = r


def frac(sample, key, t):
    xs = [float(r[key]) for (s, c), r in uniq.items() if s == sample]
    return sum(x >= t for x in xs), len(xs)


with open(f'{B}/backbone_thresholds.tsv', 'w') as fo:
    fo.write('threshold\tfrac_bacterial_flagged\tfrac_fungal_flagged\tn_bacterial\tn_fungal\t'
             'bx_frac_bacterial_flagged\tbx_frac_fungal_flagged\tanyhit_frac_bacterial_flagged\tanyhit_frac_fungal_flagged\n')
    for t in THR:
        a, na = frac('bacterial', 'best_backbone_pid', t); f_, nf = frac('fungal', 'best_backbone_pid', t)
        xa, _ = frac('bacterial', 'bx_best_backbone_pid', t); xf, _ = frac('fungal', 'bx_best_backbone_pid', t)
        ya, _ = frac('bacterial', 'best_any_bb_pid', t); yf, _ = frac('fungal', 'best_any_bb_pid', t)
        fo.write(f'{t}\t{a / na:.3f}\t{f_ / nf:.3f}\t{na}\t{nf}\t{xa / na:.3f}\t{xf / nf:.3f}\t{ya / na:.3f}\t{yf / nf:.3f}\n')
with open(f'{B}/backbone_scan_1pct.tsv', 'w') as fo:
    fo.write('threshold\tn_bact_flagged\tn_fung_flagged\tbx_n_bact_flagged\tbx_n_fung_flagged\n')
    for t in range(20, 101):
        fo.write(f"{t}\t{frac('bacterial', 'best_backbone_pid', t)[0]}\t{frac('fungal', 'best_backbone_pid', t)[0]}\t"
                 f"{frac('bacterial', 'bx_best_backbone_pid', t)[0]}\t{frac('fungal', 'bx_best_backbone_pid', t)[0]}\n")
print('contigs', len(uniq))
