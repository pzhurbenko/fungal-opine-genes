#!/usr/bin/env python
"""GC/GC3/RSCU of nos-subfamily ODH-like loci of Boletales (Table S12b, 6 loci; miniprot models with KAG8213183.1)
vs own-genome background, using exactly the procedure of composition.py (its code is executed up to the per-locus
block with OUT redirected to a scratch dir, so the original outputs are not touched; bolpool background is identical).
Output: analysis/composition/per_locus_nosODH.tsv, bg_sources_nosODH.tsv, background_genes_nosODH.tsv.gz,
cds_loci_nosODH.fna
"""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import re, gzip, sys
SRC = f'{GMO_ROOT}/analysis/composition/scripts/composition.py'
WORK = f'{GMO_ROOT}/analysis/composition/nosodh_work'
code = open(SRC).read().split('# ---------------- per locus')[0]
code = code.replace("OUT = f'{ROOT}/analysis/composition'", f"OUT = '{WORK}/orig_rerun'; DL_DIR = f'{{ROOT}}/analysis/composition/cds_dl'")
code = code.replace("dl = f'{OUT}/cds_dl'", "dl = DL_DIR")
import os; os.makedirs(f'{WORK}/orig_rerun', exist_ok=True)
exec(code)                                  # defines bg, S, S_pool, pool_genes, eval_locus, prep, ... (globals)

FINAL = f'{GMO_ROOT}/analysis/composition'
# ---- nos-ODH loci: miniprot best model per contig
contigs = dict(read_fasta(f'{WORK}/contigs.fa'))
mrna, cds = {}, {}
for l in open(f'{WORK}/mp.gff'):
    if l[0] == '#': continue
    p = l.rstrip('\n').split('\t')
    if len(p) < 9: continue
    a = dict(kv.split('=', 1) for kv in p[8].split(';') if '=' in kv)
    if p[2] == 'mRNA':
        if p[0] in mrna: continue                         # keep best rank per contig
        mrna[p[0]] = dict(id=a['ID'], s=int(p[3]), e=int(p[4]), st=p[6], ident=a['Identity'],
                          fs=a.get('Frameshift', '0'), stop=a.get('StopCodon', '0'))
    elif p[2] == 'CDS' and p[0] in mrna and a['Parent'] == mrna[p[0]]['id']:
        cds.setdefault(p[0], []).append((int(p[3]), int(p[4]), p[6], int(p[7])))

meta = {}   # contig -> organism (Table S12b)
for x in csv.DictReader(open(f'{GMO_ROOT}/paper/supplementary/TableS12b_Opine_dehydrogenase_like_loci_in_Boletales_assembl.tsv'), delimiter='\t'):
    meta[(x['assembly'], x['contig'])] = x['organism']

loci = []
for key, m in mrna.items():
    g, c = key.split('__')
    seq = contigs[key].upper()
    fs = sorted(cds[key], key=lambda f: f[0], reverse=(m['st'] == '-'))
    nt, third = [], []
    for s, e, st, ph in fs:     # phase-aware per segment (handles miniprot frameshift splits in pseudogenes)
        seg = seq[s - 1:e]; seg = rc(seg) if st == '-' else seg
        seg = seg[ph:]; seg = seg[:len(seg) - len(seg) % 3]
        nt.append(seg)
    full = ''.join(nt)
    pseudo = m['fs'] != '0' or m['stop'] != '0'
    lid = f'nosODH_{c.split(".")[0]}_L1'
    reg[(g, re.sub(r'\.\d+$', '', c))].append((m['s'], m['e']))       # exclude locus itself from background
    loci.append(dict(lid=lid, g=g, c=c, s=m['s'], e=m['e'], st=m['st'], full=full, pseudo=pseudo, m=m,
                     org=meta.get((g, c), '')))

# ---- background for these genomes (same rules as composition.py)
new_g = sorted({x['g'] for x in loci})
bgn, srcn = {}, {}
for g in new_g:
    if os.path.exists(f'{dl}/{g}/ncbi_dataset/data/{g}/cds_from_genomic.fna'):
        bgn[g], s = bg_from_datasets(g); srcn[g] = ('ncbi_annotation_CDS', s)
    else:
        bgn[g], s = bg_from_augustus(g); srcn[g] = ('augustus_neighbors', s)
    print(g, srcn[g][0], len(bgn[g]), flush=True)
Sn = {}
for g in new_g:
    if len(bgn[g]) >= 5:
        Sn[g] = prep(bgn[g]); Sn[g]['d'] = bg_distances(Sn[g]['cnt']); Sn[g]['rs'] = bg_rscu_vec(Sn[g]['tot'][None, :])

fmt = lambda v: ('%.4f' % v) if isinstance(v, (float, np.floating)) else str(v)
cols = ['locus_id', 'family', 'lineage', 'organism', 'order', 'gca', 'cds_len_nt', 'cds_source', 'gc', 'gc3',
        'bg_source', 'n_bg_genes', 'bg_small_flag', 'bg_med_gc', 'bg_q1_gc', 'bg_q3_gc', 'bg_med_gc3', 'bg_q1_gc3', 'bg_q3_gc3',
        'delta_gc3', 'gc3_pctile', 'gc3_z', 'gc_pctile', 'codon_dist_1minusr_rscu', 'codon_dist_pctile', 'n_len_matched_bg',
        'bolpool_n_genes', 'bolpool_med_gc3', 'bolpool_delta_gc3', 'bolpool_gc3_pctile', 'bolpool_gc3_z',
        'bolpool_codon_dist_pctile', 'contig', 'start', 'end', 'strand', 'status', 'miniprot_identity']
rows = []
for x in loci:
    g = x['g']; full = x['full']
    # codons with stop -> drop stops (pseudogenes) for codon-usage; GC/GC3 on all aligned codons
    cod = [full[i:i + 3] for i in range(0, len(full), 3)]
    sense = ''.join(c for c in cod if c in CODE and CODE[c] != '*')
    gc, gc3 = gc_gc3(full)
    r = dict(locus_id=x['lid'], family='opine dehydrogenase (nos subfamily)', lineage='nos-subfamily ODH (Boletales)',
             organism=x['org'], order='Boletales', gca=g, cds_len_nt=len(full),
             cds_source='miniprot(KAG8213183.1) CDS, phase-aware' + (' (pseudogene: stop/frameshift)' if x['pseudo'] else ''),
             gc=gc, gc3=gc3, bg_source=srcn[g][0], n_bg_genes=len(bgn[g]), contig=x['c'], start=x['s'], end=x['e'],
             strand=x['st'], status='pseudogene' if x['pseudo'] else 'intact', miniprot_identity=x['m']['ident'])
    if g in Sn:
        e = eval_locus(sense, Sn[g])
        # GC/GC3 percentiles must use the full aligned sequence -> recompute on gc/gc3 directly
        st = Sn[g]; n = len(st['gc3'])
        pct3 = 100 * (np.sum(st['gc3'] < gc3) + 0.5 * np.sum(st['gc3'] == gc3)) / n
        pctgc = 100 * (np.sum(st['gc'] < gc) + 0.5 * np.sum(st['gc'] == gc)) / n
        r.update(bg_small_flag='bg_small(<30)' if len(bgn[g]) < MIN_BG else 'ok', bg_med_gc=e['med_gc'], bg_q1_gc=e['q1_gc'],
                 bg_q3_gc=e['q3_gc'], bg_med_gc3=e['med_gc3'], bg_q1_gc3=e['q1_gc3'], bg_q3_gc3=e['q3_gc3'],
                 delta_gc3=gc3 - e['med_gc3'], gc3_pctile=pct3, gc3_z=(gc3 - st['gc3'].mean()) / st['gc3'].std(ddof=1),
                 gc_pctile=pctgc, codon_dist_1minusr_rscu=e['codon_dist'], codon_dist_pctile=e['codon_dist_pct'],
                 n_len_matched_bg=e['n_len_matched'])
    else:
        r['bg_small_flag'] = 'bg_unavailable(<5)'
    e = eval_locus(sense, S_pool); n = len(S_pool['gc3'])
    r.update(bolpool_n_genes=len(pool_genes), bolpool_med_gc3=e['med_gc3'], bolpool_delta_gc3=gc3 - e['med_gc3'],
             bolpool_gc3_pctile=100 * (np.sum(S_pool['gc3'] < gc3) + 0.5 * np.sum(S_pool['gc3'] == gc3)) / n,
             bolpool_gc3_z=(gc3 - S_pool['gc3'].mean()) / S_pool['gc3'].std(ddof=1), bolpool_codon_dist_pctile=e['codon_dist_pct'])
    rows.append(r)
rows.sort(key=lambda r: (r['gca'], r['contig']))
with open(f'{FINAL}/per_locus_nosODH.tsv', 'w') as f:
    f.write('\t'.join(cols) + '\n')
    for r in rows: f.write('\t'.join(fmt(r.get(c, '')) for c in cols) + '\n')
with open(f'{FINAL}/cds_loci_nosODH.fna', 'w') as f:
    for x in loci: f.write(f">{x['lid']} {x['g']} {x['c']}:{x['s']}-{x['e']}{x['st']}\n{x['full']}\n")
with open(f'{FINAL}/bg_sources_nosODH.tsv', 'w') as f:
    f.write('gca\tbg_source\tn_bg_genes\tnotes\n')
    for g in new_g: f.write(f'{g}\t{srcn[g][0]}\t{len(bgn[g])}\t{srcn[g][1]}\n')
with gzip.open(f'{FINAL}/background_genes_nosODH.tsv.gz', 'wt') as f:
    f.write('gca\tgene\tcds_len_nt\tgc\tgc3\n')
    for g in new_g:
        for gid, c in bgn[g]:
            a, b = gc_gc3(c); f.write(f'{g}\t{gid}\t{len(c)}\t{a:.4f}\t{b:.4f}\n')
print('DONE')
