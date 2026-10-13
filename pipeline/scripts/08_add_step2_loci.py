#!/usr/bin/env python3
"""Stage 08 — add NEW loci found by the step-2 whole-assembly screen (diamond blastx + six-frame hmmsearch) to the registry,
so that round 2 annotates them together with the round-1 loci.

Input : analysis/refset/screen/step2/step2_family_calls.tsv (one row per ORF call; is_new=yes -> not overlapping a registry locus)
        analysis/refset/screen/step2/genomes.tsv (gca -> local genome fasta)
Output: registry/loci.tsv (+ appended rows; backup registry/loci.round1.tsv), registry/contigs.tsv, contigs/<acc>.fa + contigs/contigs.tsv
        (local extraction from the downloaded assemblies, same full/window rule as stage 02), registry/loci_step2_added.tsv
Locus id: <FAM>_<contig>_S2L<i>. Taxonomy (division/class/order) via genus -> registry map. Calls incl. 'ambiguous' are added:
round-2 validation with the v1 reference set (competitive OUTGROUPs) makes the final call.
"""
import os, re, sys, shutil, subprocess
from collections import defaultdict
from common import load_config, read_tsv, write_tsv, read_fasta, wgs_prefix, log

cfg = load_config(); P = cfg['_paths']
S2 = os.path.join(cfg['project_dir'], 'analysis/refset/screen/step2')
FULLMAX = cfg['fetch']['contig_full_max']; FLANK = cfg['fetch']['flank_context']; FLANK_MODEL = cfg['annotate']['flank_model']
GAP = cfg['registry']['merge_gap']

loci = read_tsv(os.path.join(P['registry'], 'loci.tsv'))
if not os.path.exists(os.path.join(P['registry'], 'loci.round1.tsv')):
    shutil.copy(os.path.join(P['registry'], 'loci.tsv'), os.path.join(P['registry'], 'loci.round1.tsv'))
    shutil.copy(os.path.join(P['registry'], 'contigs.tsv'), os.path.join(P['registry'], 'contigs.round1.tsv'))
    shutil.copy(os.path.join(P['contigs'], 'contigs.tsv'), os.path.join(P['contigs'], 'contigs.round1.tsv'))
existing_ids = {l['locus_id'] for l in loci}
if any('_S2L' in i for i in existing_ids):
    log('registry already contains step-2 loci; nothing to do'); sys.exit(0)
genus_tax = {}
for l in loci:
    if l['genus'] and l['class']: genus_tax.setdefault(l['genus'], (l['division'], l['class'], l['order'], l['family_tax']))
asm = read_tsv(os.path.join(P['registry'], 'assemblies.tsv')); gca_org = {a['gca']: a['organism'] for a in asm}
genomes = {r['gca']: r['fna'] if 'fna' in r else r[list(r)[1]] for r in read_tsv(os.path.join(S2, 'genomes.tsv'))} if False else {}
for line in open(os.path.join(S2, 'genomes.tsv')):
    f = line.rstrip('\n').split('\t'); genomes[f[0]] = f[1]

# ---- collect new loci: dedupe ORF calls per (gca, contig, start, end, family), then merge same-family loci closer than GAP
calls = defaultdict(list)
for r in read_tsv(os.path.join(S2, 'step2_family_calls.tsv')):
    if r['is_new'] != 'yes': continue
    calls[(r['gca'], r['contig'], r['family'])].append(r)
new = []; by_contig_fam = defaultdict(list)
for (gca, contig_v, fam), rs in calls.items():
    regions = defaultdict(list)
    for r in rs: regions[(int(r['start']), int(r['end']))].append(r)
    items = []
    for (s, e), rr in sorted(regions.items()):
        best = max(rr, key=lambda r: float(r['score']))
        items.append({'start': s, 'end': e, 'best': best, 'n': len(rr)})
    merged = []
    for it in items:
        if merged and it['start'] - merged[-1]['end'] <= GAP:
            m = merged[-1]; m['end'] = max(m['end'], it['end']); m['n'] += it['n']
            if float(it['best']['score']) > float(m['best']['score']): m['best'] = it['best']
        else: merged.append(dict(it))
    for m in merged: by_contig_fam[(gca, contig_v, fam)].append(m)

contig_fams = defaultdict(set)
for l in loci: contig_fams[l['contig']].add(l['fam'])
rows = []
for (gca, contig_v, fam), ms in sorted(by_contig_fam.items()):
    contig = re.sub(r'\.\d+$', '', contig_v); org = gca_org.get(gca) or ms[0]['best']['organism']
    genus = org.split()[0].strip('[]') if org else ''
    div, cls, order, ftax = genus_tax.get(genus, ('', '', '', ''))
    wp = wgs_prefix(contig)
    genome = (wp[:-2] + '0' * (len(contig) - len(wp) + 2)) if wp else gca
    for i, m in enumerate(ms, 1):
        b = m['best']
        rows.append({'locus_id': f'{fam}_{contig}_S2L{i}', 'fam': fam, 'contig': contig, 'genome': genome, 'organism': org, 'genus': genus,
                     'division': div, 'class': cls, 'order': order, 'family_tax': ftax,
                     'hsp_start': m['start'], 'hsp_end': m['end'], 'n_hsp': m['n'],
                     'best_identity': f"{float(b['best_pident'] or 0) / 100:.2f}", 'best_evalue': '',
                     'queries': f"step2:{b['call']}:score{b['score']}:gap{b['gap']}", 'strategies': 'step2_blastx_hmm',
                     'was_deleted': 'no', 'reasons': '', 'tree_clade': f"S2_{b['call']}", 'is_rep_contig': 'no',
                     'other_fams_on_contig': '', 'gca': gca, 'contig_v': contig_v, 'contam_flag': b['contam_flag']})
        contig_fams[contig].add(fam)
log(f'new loci from step 2: {len(rows)} on {len({r["contig"] for r in rows})} contigs, genomes {len({r["gca"] for r in rows})}')

# ---- contigs: extract locally (full <= 1 Mb, else window around ALL loci of the contig), update contigs/contigs.tsv
cmeta_p = os.path.join(P['contigs'], 'contigs.tsv'); cmeta = {r['contig']: r for r in read_tsv(cmeta_p)}
need = defaultdict(set)   # gca -> contig_v
span = defaultdict(lambda: [10**12, 0])
for l in loci: 
    s = span[l['contig']]; s[0] = min(s[0], int(l['hsp_start'])); s[1] = max(s[1], int(l['hsp_end']))
for r in rows:
    s = span[r['contig']]; s[0] = min(s[0], r['hsp_start']); s[1] = max(s[1], r['hsp_end'])
for r in rows:
    c = cmeta.get(r['contig'])
    if c and (c['fetched'] == 'full' or (int(c['seq_start']) <= span[r['contig']][0] - FLANK_MODEL and int(c['seq_end']) >= span[r['contig']][1] + FLANK_MODEL)):
        continue
    need[r['gca']].add(r['contig_v'])
log(f'contigs to extract/re-extract locally: {sum(len(v) for v in need.values())} from {len(need)} genomes')
n_done = 0
for gca, cvs in need.items():
    ids = os.path.join(P['contigs'], f'_ids_{gca}.txt'); open(ids, 'w').write('\n'.join(sorted(cvs)) + '\n')
    tmp = os.path.join(P['contigs'], f'_tmp_{gca}.fa')
    subprocess.run(['seqkit', 'grep', '-f', ids, genomes[gca], '-o', tmp], check=True, stderr=subprocess.DEVNULL)
    for h, seq in read_fasta(tmp):
        cv = h.split()[0]; contig = re.sub(r'\.\d+$', '', cv); L = len(seq); org = gca_org.get(gca, '')
        hmin, hmax = span[contig]
        if L <= FULLMAX: s, e, mode = 1, L, 'full'
        else: s, e, mode = max(1, hmin - FLANK), min(L, hmax + FLANK), 'window'
        path = os.path.join(P['contigs'], f'{contig}.fa')
        with open(path, 'w') as f:
            f.write(f'>{contig} {cv} {org} seq_start={s} seq_end={e} contig_len={L}\n')
            sub = seq[s-1:e]
            for i in range(0, len(sub), 80): f.write(sub[i:i+80] + '\n')
        cmeta[contig] = {'contig': contig, 'acc_ver': cv, 'length': L, 'organism': org, 'fetched': mode, 'seq_start': s, 'seq_end': e,
                         'path': path, 'wgs_prefix': wgs_prefix(contig) or ''}
        n_done += 1
    os.remove(ids); os.remove(tmp)
write_tsv(cmeta_p, list(cmeta.values()))
log(f'contigs written: {n_done}; total in contigs.tsv: {len(cmeta)}')
missing = [r['contig'] for r in rows if r['contig'] not in cmeta]
if missing: log(f'WARNING: {len(missing)} new loci have no contig sequence (first: {missing[:5]})')

# ---- registry: append loci, recompute other_fams_on_contig, rebuild registry/contigs.tsv
cols = list(loci[0].keys()) + ['gca']
for l in loci: l.setdefault('gca', '')
all_loci = loci + [{k: r.get(k, '') for k in cols} for r in rows]
for l in all_loci: l['other_fams_on_contig'] = ','.join(sorted(contig_fams[l['contig']] - {l['fam']}))
write_tsv(os.path.join(P['registry'], 'loci.tsv'), all_loci, cols)
write_tsv(os.path.join(P['registry'], 'loci_step2_added.tsv'), rows)
by_c = defaultdict(list)
for l in all_loci: by_c[l['contig']].append(l)
write_tsv(os.path.join(P['registry'], 'contigs.tsv'),
          [{'contig': c, 'fams': ','.join(sorted({l['fam'] for l in ls})), 'n_loci': len(ls), 'organism': ls[0]['organism'],
            'hsp_min': min(int(l['hsp_start']) for l in ls), 'hsp_max': max(int(l['hsp_end']) for l in ls)} for c, ls in sorted(by_c.items())])
from collections import Counter
log('registry now:', len(all_loci), 'loci;', dict(Counter(l['fam'] for l in all_loci)), '; added by call:', dict(Counter(r['tree_clade'] for r in rows)))
