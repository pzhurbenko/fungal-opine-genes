#!/usr/bin/env python3
"""Aggregate step-2 screening over whole fungal assemblies (diamond blastx + six-frame hmmsearch).

Per genome:
  * diamond HSPs -> chained per (contig, strand, subject) with gap <= GAP -> overlapping chains merged into LOCI;
    each locus keeps the best subject (bits), its category (TDNA/opine_related/vir/tra/rep_partition/T1/T2 = our families)
    and lineage, pident, subject coverage, and the set of categories hit.
  * hmmsearch on six-frame ORFs -> per ORF competitive family call (same rule as step 3: OUTGROUP if og>sub; ambiguous if
    gap<20 or score<30) + best T-DNA lineage HMM; ORF coords mapped back to the contig and attached to overlapping loci.
  * loci are compared with the round-1 registry (registry/loci.tsv) -> known / new.
Outputs (in step2/): step2_loci.tsv, step2_family_calls.tsv, step2_contamination.tsv, step2_genome_summary.tsv, step2_new_loci.tsv
"""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import os, glob, csv, re, sys
from collections import defaultdict, Counter
D = f'{GMO_ROOT}/analysis/refset/screen/step2'
REG = f'{GMO_ROOT}/analysis/loci_pipeline/registry'
GAP = 5000          # max nt gap when chaining HSPs of the same subject (introns / frameshift breaks)
CONTAM_PID = 60.0   # plasmid backbone hit >= this identity -> bacterial contamination signal
CAND_PID, CAND_COV = 45.0, 0.5
SHORT_CONTIG = 20000  # bp   # diamond-only candidate threshold (same as step-3 outlier rule)
FAMS = ['ACS', 'CUS', 'SUS', 'ORF358']
FAM_OF_T = {}       # T1/T2 query id -> family (FAM_ACS|T1|...)

def strip_ver(c): return re.sub(r'\.\d+$', '', c)

# ---------- registry ----------
reg = defaultdict(list)        # contig(no version) -> [(start,end,locus_id,fam)]
asm_org = {}; wgs_tax = {}; genus_tax = {}; asm_tax = {}
for r in csv.DictReader(open(f'{REG}/loci.tsv'), delimiter='\t'):
    reg[r['contig']].append((int(r['hsp_start']), int(r['hsp_end']), r['locus_id'], r['fam']))
    wgs_tax[r['genome'][:6]] = (r['division'], r['class'], r['order'])
    if r['genus'] and r['class']: genus_tax[r['genus']] = (r['division'], r['class'], r['order'])
for r in csv.DictReader(open(f'{REG}/assemblies.tsv'), delimiter='\t'):
    asm_org[r['gca']] = r['organism']
    genus = r['organism'].split()[0].strip('[]') if r['organism'] else ''
    asm_tax[r['gca']] = wgs_tax.get(r['wgs_prefix'][:6]) or genus_tax.get(genus, ('', '', ''))

def overlaps(a, b, c, d): return min(b, d) - max(a, c) >= 0

# ---------- helpers ----------
def parse_subject(s):
    p = s.split('|')
    if p[0].startswith('FAM_'):                       # FAM_ACS|T1|Agro_xxx
        return dict(cat=p[0][4:], lineage=p[1], acc=p[2] if len(p) > 2 else '')
    lineage = p[4] if len(p) > 4 and p[4] != '-' else (p[5][:30] if len(p) > 5 else '')
    return dict(cat=p[1], lineage=lineage, acc=p[0])

loci_rows, fam_rows, contam_rows, genome_rows, new_rows = [], [], [], [], []
for done in sorted(glob.glob(f'{D}/out/*/done')):
    g = os.path.basename(os.path.dirname(done)); od = os.path.dirname(done)
    # ---- diamond: chain HSPs per (contig, strand, subject) ----
    hsps = defaultdict(list); contig_len = {}
    for l in open(f'{od}/dmnd.tsv'):
        f = l.rstrip('\n').split('\t')
        if len(f) < 13: continue
        q, s = f[0], f[1]; pid = float(f[2]); alen = int(f[3]); slen = int(f[5]); contig_len[q] = int(f[4])
        qs, qe = sorted((int(f[6]), int(f[7]))); ss, se = int(f[8]), int(f[9]); ev = float(f[10]); bits = float(f[11])
        strand = '+' if int(f[12]) > 0 else '-'
        hsps[(q, strand, s)].append((qs, qe, ss, se, pid, bits, ev, slen))
    chains = defaultdict(list)   # (contig,strand) -> list of chains
    for (q, strand, s), hl in hsps.items():
        hl.sort()
        cur = None
        for qs, qe, ss, se, pid, bits, ev, slen in hl:
            if cur and qs - cur['end'] <= GAP:
                cur['end'] = max(cur['end'], qe); cur['bits'] += bits; cur['n'] += 1
                cur['scov'].append((min(ss, se), max(ss, se))); cur['pid'] = max(cur['pid'], pid); cur['ev'] = min(cur['ev'], ev)
            else:
                if cur: chains[(q, strand)].append(cur)
                cur = dict(start=qs, end=qe, subj=s, bits=bits, n=1, scov=[(min(ss, se), max(ss, se))], pid=pid, ev=ev, slen=slen)
        chains[(q, strand)].append(cur)
    # ---- merge overlapping chains (any subject) into loci ----
    loci = []
    for (q, strand), cl in chains.items():
        cl.sort(key=lambda c: c['start'])
        cur = None
        for c in cl:
            if cur and c['start'] <= cur['end']:
                cur['end'] = max(cur['end'], c['end']); cur['members'].append(c)
            else:
                if cur: loci.append(cur)
                cur = dict(contig=q, strand=strand, start=c['start'], end=c['end'], members=[c])
        loci.append(cur)
    # ---- HMM on six-frame ORFs ----
    orf_coord = {}
    for l in open(f'{od}/hit_orfs.faa'):
        if l.startswith('>'):
            m = re.match(r'>(\S+) source=(\S+) coords=(\d+)\.\.(\d+) length=(\d+) frame=(\d+)', l)
            if m:
                a, b = int(m.group(3)), int(m.group(4))
                orf_coord[m.group(1)] = (m.group(2), min(a, b), max(a, b), '+' if int(m.group(6)) <= 3 else '-')
    per = defaultdict(dict)
    for l in open(f'{od}/hmm.tbl'):
        if l.startswith('#'): continue
        f = l.split(); per[f[0]][f[2]] = (float(f[5]), float(f[4]))
    orf_calls = {}   # orf -> dict(fam_call, fam_sub, fam_score, gap, tdna_best, tdna_score)
    for orf, d in per.items():
        rec = {}
        for fam in FAMS:
            subs = {k: v for k, v in d.items() if k.startswith(fam + '__') and not k.endswith('OUTGROUP')}
            og = d.get(f'{fam}__OUTGROUP', (0, 99))[0]
            if not subs:
                if og > 0: rec.setdefault('outgroup_only', []).append(f'{fam}:{og:.0f}')
                continue
            bs = max(subs, key=lambda k: subs[k][0]); sc = subs[bs][0]; gap = sc - og
            call = 'OUTGROUP' if og > sc else ('ambiguous' if (gap < 20 or sc < 30) else bs.split('__')[1])
            rec.setdefault('fam', []).append((fam, call, sc, gap))
        td = {k: v for k, v in d.items() if '__' not in k}
        if td:
            bt = max(td, key=lambda k: td[k][0]); rec['tdna'] = (bt, td[bt][0])
        orf_calls[orf] = rec
    # ---- attach ORF calls to loci; HMM-only loci ----
    orf_by_contig = defaultdict(list)
    for orf, (c, a, b, st) in orf_coord.items(): orf_by_contig[c].append((a, b, st, orf))
    used = set()
    for L in loci:
        L['orfs'] = [o for a, b, st, o in orf_by_contig.get(L['contig'], []) if overlaps(a, b, L['start'], L['end'])]
        used.update(L['orfs'])
    for orf, (c, a, b, st) in orf_coord.items():
        if orf in used: continue
        rec = orf_calls.get(orf, {})
        if rec.get('fam') or rec.get('tdna'):
            loci.append(dict(contig=c, strand=st, start=a, end=b, members=[], orfs=[orf]))
    # ---- contamination signals per contig (before emitting) ----
    contig_backbone = Counter(); contig_hi = Counter()
    for L in loci:
        if not L['members']: continue
        b = max(L['members'], key=lambda m: m['bits']); cat = parse_subject(b['subj'])['cat']
        if cat in ('vir', 'tra', 'rep_partition', 'mobile') and b['pid'] >= CONTAM_PID: contig_backbone[L['contig']] += 1
        if b['pid'] >= CONTAM_PID: contig_hi[L['contig']] += 1
    genome_contam = len(contig_backbone)
    # ---- emit loci ----
    gsum = Counter(); gsum['gca'] = g
    for L in loci:
        L['members'].sort(key=lambda m: -m['bits'])
        best = L['members'][0] if L['members'] else None
        cats = Counter(parse_subject(m['subj'])['cat'] for m in L['members'])
        # subject coverage of best chain
        scov = 0.0
        if best:
            iv = sorted(best['scov']); merged = []
            for a, b in iv:
                if merged and a <= merged[-1][1]: merged[-1][1] = max(merged[-1][1], b)
                else: merged.append([a, b])
            scov = sum(b - a + 1 for a, b in merged) / best['slen']
        # HMM summary
        fam_calls = []; td_best = ('', 0.0); og_only = []
        for o in L['orfs']:
            rec = orf_calls.get(o, {})
            fam_calls += [(fam, call, sc, gap, o) for fam, call, sc, gap in rec.get('fam', [])]
            og_only += rec.get('outgroup_only', [])
            if rec.get('tdna') and rec['tdna'][1] > td_best[1]: td_best = rec['tdna']
        fam_calls.sort(key=lambda x: -x[2])
        fc = fam_calls[0] if fam_calls else None
        bsub = parse_subject(best['subj']) if best else dict(cat='', lineage='', acc='')
        tax = asm_tax.get(g, ('', '', ''))
        row = dict(gca=g, organism=asm_org.get(g, ''), division=tax[0], tax_class=tax[1], tax_order=tax[2],
                   contig=L['contig'], contig_nover=strip_ver(L['contig']), contig_len=contig_len.get(L['contig'], ''), strand=L['strand'],
                   start=L['start'], end=L['end'], span=L['end'] - L['start'] + 1,
                   best_subject=best['subj'] if best else '', best_cat=bsub['cat'], best_lineage=bsub['lineage'],
                   best_pident=f"{best['pid']:.1f}" if best else '', best_bits=f"{best['bits']:.0f}" if best else '',
                   best_evalue=f"{best['ev']:.1e}" if best else '', subj_cov=f'{scov:.2f}', n_hsp=best['n'] if best else 0,
                   n_subjects=len(L['members']), cats=','.join(f'{k}:{v}' for k, v in cats.most_common()),
                   hmm_family=fc[0] if fc else '', hmm_call=fc[1] if fc else '', hmm_score=f'{fc[2]:.0f}' if fc else '',
                   hmm_gap=f'{fc[3]:.0f}' if fc else '', hmm_outgroup_only=';'.join(sorted(set(og_only))),
                   hmm_tdna_lineage=td_best[0][:40], hmm_tdna_score=f'{td_best[1]:.0f}' if td_best[0] else '', n_orfs=len(L['orfs']))
        # contamination flags: contig carries plasmid-backbone genes >=60% / >=60%-hits of any category; genome-level count
        row['contig_backbone_hits'] = contig_backbone.get(L['contig'], 0); row['contig_hi_id_hits'] = contig_hi.get(L['contig'], 0)
        row['genome_contam_contigs'] = genome_contam
        # short_contig: bacterial-like (>=60% to a plasmid protein) hit on a contig < SHORT_CONTIG bp — typical for contaminating fragments
        clen = contig_len.get(L['contig'], 10**9)
        row['contam_flag'] = ('contig' if contig_backbone.get(L['contig'], 0) > 0 else 'genome' if genome_contam >= 3
                              else 'short_contig' if (best and best['pid'] >= CONTAM_PID and clen < SHORT_CONTIG) else 'no')
        # registry match
        known = [lid for a, b, lid, fam in reg.get(row['contig_nover'], []) if overlaps(a, b, L['start'], L['end'])]
        row['registry_loci'] = ','.join(known); row['is_new'] = 'no' if known else 'yes'
        loci_rows.append(row)
        for fam, call, sc, gap, o in fam_calls:
            if call != 'OUTGROUP':
                fam_rows.append(dict(gca=g, organism=row['organism'], contig=L['contig'], start=L['start'], end=L['end'], strand=L['strand'],
                                     family=fam, call=call, score=f'{sc:.0f}', gap=f'{gap:.0f}', orf=o,
                                     best_subject_cat=row['best_cat'], best_pident=row['best_pident'], is_new=row['is_new'],
                                     registry_loci=row['registry_loci'], contam_flag=row['contam_flag'], tax_class=row['tax_class'],
                                     contig_len=row['contig_len']))
        # candidate = confident competitive family call, OR a T-DNA/opine/family subject hit that is clearly above the native
        # background (identity >= CAND_PID over >= CAND_COV of the subject); lineage-wise outlier analysis is done later across genomes
        strong = best is not None and best['pid'] >= CAND_PID and scov >= CAND_COV and bsub['cat'] in ('T1', 'T2', 'TDNA', 'opine_related')
        row['candidate'] = 'yes' if (row['hmm_call'] not in ('', 'OUTGROUP', 'ambiguous') or strong) else 'no'
        gsum['n_loci'] += 1; gsum['n_' + (bsub['cat'] or 'hmm_only')] += 1
        if row['candidate'] == 'yes': gsum['n_candidates'] += 1
        if row['is_new'] == 'yes' and row['candidate'] == 'yes': gsum['n_new_candidates'] += 1
    for c, n in contig_backbone.items():
        contam_rows.append(dict(gca=g, organism=asm_org.get(g, ''), contig=c, contig_len=contig_len.get(c, ''), n_backbone_hits_ge60=n,
                                n_hi_id_hits=contig_hi.get(c, 0)))
    gsum['n_contam_contigs'] = len(contig_backbone); gsum['n_hi_id_contigs'] = len(contig_hi); gsum['organism'] = asm_org.get(g, '')
    genome_rows.append(dict(gsum))

def write(name, rows, cols=None):
    if not rows: open(f'{D}/{name}', 'w').close(); return
    cols = cols or list(rows[0].keys())
    with open(f'{D}/{name}', 'w') as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter='\t', extrasaction='ignore'); w.writeheader(); w.writerows(rows)
write('step2_loci.tsv', loci_rows)
write('step2_family_calls.tsv', fam_rows)
write('step2_contamination.tsv', contam_rows)
gcols = ['gca', 'organism', 'n_loci', 'n_candidates', 'n_new_candidates', 'n_contam_contigs', 'n_hi_id_contigs'] + sorted({k for r in genome_rows for k in r if k.startswith('n_') and k not in ('n_loci', 'n_candidates', 'n_new_candidates', 'n_contam_contigs', 'n_hi_id_contigs')})
write('step2_genome_summary.tsv', genome_rows, gcols)
write('step2_new_loci.tsv', [r for r in loci_rows if r['is_new'] == 'yes' and r['candidate'] == 'yes'])
print('genomes:', len(genome_rows), 'loci:', len(loci_rows), 'family calls (non-OUTGROUP):', len(fam_rows),
      'candidates:', sum(1 for r in loci_rows if r['candidate'] == 'yes'),
      'new candidates:', sum(1 for r in loci_rows if r['is_new'] == 'yes' and r['candidate'] == 'yes'))
