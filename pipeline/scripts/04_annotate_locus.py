#!/usr/bin/env python3
"""Stage 04 — annotate ONE locus (run as LSF job array element or directly).

  usage: 04_annotate_locus.py <round> <locus_id | index>   (index = 1-based line in round<N>/loci_list.txt)

Steps
  1. model window  = HSP ± annotate.flank_model (15 kb) on the downloaded contig sequence
     context       = whole downloaded sequence if <= annotate.context_max, else HSP ± flank_context
  2. miniprot (guides T1+T2+T3 of the family) on the model window  -> gene structure, frameshift/stop flags,
     protein (##STA), CDS nucleotide (from CDS lines), best model = highest miniprot score overlapping the HSP
  3. AUGUSTUS ab initio on the context, species chosen by taxonomic class, hints = miniprot CDS (CDSpart)
     -> independent model of the target + all neighbouring genes (exon counts, proteins)
  4. neighbours: adaptive window — expand 15 -> 50 -> 100 kb -> edge until >= N genes with >= 2 exons per side
  5. GC profile in 10-kb bins over the context; GC of the model window
Outputs: round<N>/loci/<locus_id>/{model_window.fa, context.fa, guides.faa->symlink, miniprot.gff, hints.gff, augustus.gff3,
         target_protein.faa, target_cds.fna, neighbors.faa, locus.json}
"""
import os, sys, re, json, subprocess, shutil
from common import load_config, round_dir, read_tsv, read_fasta, write_fasta, revcomp, gc, log

cfg = load_config(); P = cfg['_paths']; A = cfg['annotate']
rnd = sys.argv[1]; key = sys.argv[2]
RD = round_dir(cfg, rnd)
if key.isdigit():
    key = open(os.path.join(RD, 'loci_list.txt')).read().split()[int(key) - 1]
loci = {r['locus_id']: r for r in read_tsv(os.path.join(P['registry'], 'loci.tsv'))}
L = loci[key]
cmeta = {r['contig']: r for r in read_tsv(os.path.join(P['contigs'], 'contigs.tsv'))}
C = cmeta[L['contig']]
out = os.path.join(RD, 'loci', key); os.makedirs(out, exist_ok=True)

hdr, seq = read_fasta(C['path'])[0]
seq_start = int(C['seq_start'])                    # contig coordinate of seq[0]
def to_local(x): return x - seq_start + 1          # contig -> 1-based local
hs, he = to_local(int(L['hsp_start'])), to_local(int(L['hsp_end']))
n = len(seq)
ms, me = max(1, hs - A['flank_model']), min(n, he + A['flank_model'])
if n <= A['context_max']: cs, ce = 1, n
else: cs, ce = max(1, hs - cfg['fetch']['flank_context']), min(n, he + cfg['fetch']['flank_context'])
mw = seq[ms-1:me]; ctx = seq[cs-1:ce]
write_fasta(os.path.join(out, 'model_window.fa'), [(f'{key}_mw local_start={ms}', mw)])
write_fasta(os.path.join(out, 'context.fa'), [(f'{key}_ctx local_start={cs}', ctx)])

# ---------------------------------------------------------------- miniprot
guides = os.path.join(P['guides'], f"{L['fam']}_guides.faa")
mp_gff = os.path.join(out, 'miniprot.gff')
subprocess.run(['miniprot', '-t', str(A['miniprot']['threads']), '-G', str(A['miniprot']['max_intron']), '--gff', '--trans',
                os.path.join(out, 'model_window.fa'), guides], stdout=open(mp_gff, 'w'), stderr=open(os.path.join(out, 'miniprot.log'), 'w'), check=True)

def parse_miniprot(path):
    models = []; cur = None; sta = ''
    for line in open(path):
        if line.startswith('##STA'): sta = line.split('\t')[1].strip(); continue
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9: continue
        if f[2] == 'mRNA':
            a = dict(kv.split('=', 1) for kv in f[8].split(';') if '=' in kv); t = a['Target'].split()
            cur = {'id': a['ID'], 'start': int(f[3]), 'end': int(f[4]), 'strand': f[6], 'score': int(f[5]),
                   'guide': t[0], 'guide_tier': t[0].split('|')[0], 'guide_from': int(t[1]), 'guide_to': int(t[2]),
                   'identity': float(a.get('Identity', 0)), 'positive': float(a.get('Positive', 0)),
                   'frameshift': int(a.get('Frameshift', 0)), 'stop_codons': int(a.get('StopCodon', 0)),
                   'cds': [], 'protein': sta}
            sta = ''; models.append(cur)
        elif f[2] == 'CDS' and cur is not None:
            cur['cds'].append((int(f[3]), int(f[4])))
    return models

mp = parse_miniprot(mp_gff)
# HSP in model-window coordinates
hw_s, hw_e = hs - ms + 1, he - ms + 1
overl = [m for m in mp if m['start'] <= hw_e + 300 and m['end'] >= hw_s - 300]
best = max(overl, key=lambda m: m['score']) if overl else None
best_by_tier = {}
for m in overl:
    if m['guide_tier'] not in best_by_tier or m['score'] > best_by_tier[m['guide_tier']]['score']:
        best_by_tier[m['guide_tier']] = m

def cds_seq(m, window):
    parts = [window[s-1:e] for s, e in sorted(m['cds'])]
    s = ''.join(parts)
    return revcomp(s) if m['strand'] == '-' else s

model_source = 'miniprot' if best else None
if best:
    write_fasta(os.path.join(out, 'target_protein.faa'),
                [(f"{key}|{L['fam']}|miniprot|{best['guide_tier']}|exons={len(best['cds'])}|id={best['identity']:.2f}|fs={best['frameshift']}|stop={best['stop_codons']}|guide={best['guide']}", best['protein'])])
    write_fasta(os.path.join(out, 'target_cds.fna'), [(f"{key}|{L['fam']}|miniprot|{best['guide_tier']}|strand={best['strand']}", cds_seq(best, mw))])

# ---------------------------------------------------------------- AUGUSTUS
def species_for(L):
    return A['augustus']['by_class'].get(L['class']) or A['augustus']['by_division'].get(L['division']) or A['augustus']['default_species']
sp = species_for(L)
hints = os.path.join(out, 'hints.gff')
with open(hints, 'w') as h:
    for m in mp:
        for s, e in m['cds']:   # model-window -> context coordinates
            h.write(f"{key}_ctx\tminiprot\tCDSpart\t{s + ms - cs}\t{e + ms - cs}\t.\t{m['strand']}\t.\tsrc=P;pri=4;grp={m['id']}\n")
cfgdir = subprocess.run(['bash', '-c', 'echo $AUGUSTUS_CONFIG_PATH'], capture_output=True, text=True).stdout.strip() or os.path.join(os.path.dirname(shutil.which('augustus')), '..', 'config')
extr = os.path.join(cfgdir, 'extrinsic', 'extrinsic.M.RM.E.W.P.cfg')
cmd = ['augustus', f'--species={sp}', '--strand=both', '--gff3=on', '--protein=on', '--codingseq=on', '--UTR=off',
       '--softmasking=0', os.path.join(out, 'context.fa')]
if os.path.getsize(hints) > 0 and os.path.exists(extr):
    cmd[1:1] = [f'--hintsfile={hints}', f'--extrinsicCfgFile={extr}']
aug_gff = os.path.join(out, 'augustus.gff3')
subprocess.run(cmd, stdout=open(aug_gff, 'w'), stderr=open(os.path.join(out, 'augustus.log'), 'w'))

def parse_augustus(path):
    genes = {}; prot = {}; cur = None; buf = None
    for line in open(path):
        if line.startswith('#'):
            m = re.match(r'# protein sequence = \[(.*)', line)
            if m:
                s = m.group(1).rstrip()
                if s.endswith(']'): prot[cur] = s[:-1]
                else: buf = [s]
            elif buf is not None:
                s = line[2:].strip()
                if s.endswith(']'): buf.append(s[:-1]); prot[cur] = ''.join(buf); buf = None
                else: buf.append(s)
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9: continue
        if f[2] == 'gene':
            gid = re.search(r'ID=([^;]+)', f[8]).group(1); cur = gid
            genes[gid] = {'id': gid, 'start': int(f[3]), 'end': int(f[4]), 'strand': f[6], 'exons': 0}
        elif f[2] == 'transcript': cur = re.search(r'ID=([^;]+)', f[8]).group(1).rsplit('.t', 1)[0]
        elif f[2] == 'CDS':
            gid = re.search(r'Parent=([^;]+)', f[8]).group(1).rsplit('.t', 1)[0]; genes[gid]['exons'] += 1
    for g in genes: genes[g]['protein'] = prot.get(g, '')
    return list(genes.values())

ag = parse_augustus(aug_gff) if os.path.exists(aug_gff) else []
hc_s, hc_e = hs - cs + 1, he - cs + 1                       # HSP in context coordinates
target_ag = [g for g in ag if g['start'] <= hc_e and g['end'] >= hc_s]
neighbors = [g for g in ag if not (g['start'] <= hc_e and g['end'] >= hc_s)]
# AUGUSTUS target gene(s): always exported (independent model); used as THE model when miniprot found nothing
write_fasta(os.path.join(out, 'augustus_target.faa'),
            [(f"{key}|{L['fam']}|augustus|{g['id']}|exons={g['exons']}|{g['start'] + cs - 1 + seq_start - 1}-{g['end'] + cs - 1 + seq_start - 1}{g['strand']}", g['protein']) for g in target_ag if g['protein']])
if not best and target_ag:
    g = max(target_ag, key=lambda g: len(g['protein']))
    model_source = 'augustus'
    write_fasta(os.path.join(out, 'target_protein.faa'),
                [(f"{key}|{L['fam']}|augustus|NA|exons={g['exons']}|id=0.00|fs=0|stop=0|guide={g['id']}", g['protein'])])
    # CDS from AUGUSTUS gff3 (context coordinates)
    cds = []
    for line in open(aug_gff):
        f = line.rstrip('\n').split('\t')
        if len(f) >= 9 and f[2] == 'CDS' and re.search(r'Parent=' + re.escape(g['id']) + r'\.t', f[8]):
            cds.append(ctx[int(f[3]) - 1:int(f[4])])
    s = ''.join(cds); s = revcomp(s) if g['strand'] == '-' else s
    write_fasta(os.path.join(out, 'target_cds.fna'), [(f"{key}|{L['fam']}|augustus|NA|strand={g['strand']}", s)])
if model_source is None:
    # Third fallback (a): the user's original tblastn HSP translation(s) for this contig/region (data/Fasta/<FAM>_tblastn.fas,
    # header 'desc|ACC|start:end|len'); '*' = in-frame stops are kept as information (counted, then replaced by X for HMM).
    tfas = os.path.join(cfg['project_dir'], 'data', 'Fasta', f"{L['fam']}_tblastn.fas")
    hsp_recs = []
    if os.path.exists(tfas):
        for h, s in read_fasta(tfas):
            parts = h.split('|')
            if len(parts) >= 3 and parts[1].split('.')[0] == L['contig']:
                a, b = (int(x) for x in parts[2].split(':')); a, b = min(a, b), max(a, b)
                if a <= int(L['hsp_end']) + 5000 and b >= int(L['hsp_start']) - 5000: hsp_recs.append((a, b, s))
    if hsp_recs:
        hsp_recs.sort(); prot = ''.join(s for _, _, s in hsp_recs); nstop = prot.count('*')
        model_source = 'hsp_concat'
        write_fasta(os.path.join(out, 'target_protein.faa'),
                    [(f"{key}|{L['fam']}|hsp_concat|NA|exons=1|id=0.00|fs=0|stop={nstop}|guide=tblastn_hsp_x{len(hsp_recs)}", prot.replace('*', 'X'))])
if model_source is None:
    # Third fallback (b): translate the HSP region (6 frames), keep the longest stop-free segment.
    # Such 'remnants' (typically 80-130 aa) are evidence of past presence/decay; validated by HMM in stage 05 -> status 'remnant'.
    CODON = {a+b+c: aa for a, b, c, aa in zip(*[iter(''.join(x*16 for x in 'TCAG') + ''.join((x*4 for x in 'TCAG')) * 4 + 'TCAG' * 16)]*1, 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')} if False else None
    import itertools
    bases = 'TCAG'; aas = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    CODON = {b1 + b2 + b3: aas[i] for i, (b1, b2, b3) in enumerate(itertools.product(bases, repeat=3))}
    def translate(s):
        return ''.join(CODON.get(s[i:i+3].upper(), 'X') for i in range(0, len(s) - 2, 3))
    region = seq[max(0, hs - 1 - 30):he + 30]
    best_seg, best_frame = '', ''
    for strand, sq in (('+', region), ('-', revcomp(region))):
        for f in range(3):
            for seg in translate(sq[f:]).split('*'):
                if len(seg) > len(best_seg): best_seg, best_frame = seg, f'{strand}{f}'
    if len(best_seg) >= 30:
        model_source = 'hsp_translation'
        write_fasta(os.path.join(out, 'target_protein.faa'),
                    [(f"{key}|{L['fam']}|hsp_translation|NA|exons=1|id=0.00|fs=0|stop=0|guide=frame{best_frame}", best_seg)])
write_fasta(os.path.join(out, 'neighbors.faa'),
            [(f"{key}|{g['id']}|{g['start'] + cs - 1 + seq_start - 1}-{g['end'] + cs - 1 + seq_start - 1}{g['strand']}|exons={g['exons']}", g['protein']) for g in neighbors if g['protein']])

# ---------------------------------------------------------------- adaptive neighbour window
N = A['neighbors']['min_multiexon_each_side']
def count_side(side, dist):
    lo, hi = (hc_s - dist, hc_s) if side == 'left' else (hc_e, hc_e + dist)
    sel = [g for g in neighbors if g['end'] >= lo and g['start'] <= hi]
    return len(sel), sum(1 for g in sel if g['exons'] >= 2)
adaptive = {}
for side in ['left', 'right']:
    edge = hc_s - 1 if side == 'left' else (ce - cs + 1) - hc_e
    for step in [k * 1000 for k in A['neighbors']['steps_kb']] + [edge]:
        d = min(step, edge); tot, multi = count_side(side, d)
        if multi >= N or d >= edge:
            adaptive[side] = {'distance_bp': d, 'genes': tot, 'multiexon_genes': multi, 'reached_edge': d >= edge}; break

# ---------------------------------------------------------------- GC
B = A['neighbors']['gc_bin']
gc_bins = [{'local_start': cs + i, 'gc': round(gc(ctx[i:i+B]), 1)} for i in range(0, len(ctx), B) if len(ctx[i:i+B]) >= B // 2]

info = {
    'locus_id': key, 'fam': L['fam'], 'contig': L['contig'], 'organism': L['organism'], 'class': L['class'], 'model_source': model_source,
    'contig_len': int(C['length']), 'fetched': C['fetched'], 'seq_start': seq_start,
    'hsp_local': [hs, he], 'model_window_local': [ms, me], 'context_local': [cs, ce], 'augustus_species': sp,
    'miniprot_n_models_overlapping': len(overl),
    'miniprot_best': ({k: v for k, v in best.items() if k not in ('protein',)} | {'protein_len': len(best['protein']), 'n_exons': len(best['cds']),
                       'contig_start': best['start'] + ms - 1 + seq_start - 1, 'contig_end': best['end'] + ms - 1 + seq_start - 1}) if best else None,
    'miniprot_best_by_tier': {t: {'score': m['score'], 'identity': m['identity'], 'guide': m['guide'], 'n_exons': len(m['cds']), 'protein_len': len(m['protein']),
                                  'frameshift': m['frameshift'], 'stop_codons': m['stop_codons']} for t, m in best_by_tier.items()},
    'augustus_target': [{'id': g['id'], 'exons': g['exons'], 'protein_len': len(g['protein']), 'strand': g['strand'],
                         'contig_start': g['start'] + cs - 1 + seq_start - 1, 'contig_end': g['end'] + cs - 1 + seq_start - 1} for g in target_ag],
    'augustus_n_genes_context': len(ag), 'n_neighbors': len(neighbors),
    'neighbors_exons': sorted(g['exons'] for g in neighbors), 'adaptive_neighbors': adaptive,
    'gc_model_window': round(gc(mw), 1), 'gc_context': round(gc(ctx), 1), 'gc_bins': gc_bins,
}
json.dump(info, open(os.path.join(out, 'locus.json'), 'w'), indent=1)
print(f"{key}\tbest={'none' if not best else f'{best['guide_tier']} id={best['identity']:.2f} ex={len(best['cds'])} fs={best['frameshift']} stop={best['stop_codons']} len={len(best['protein'])}'}\taug_target={len(target_ag)}\tneighbors={len(neighbors)}\tsp={sp}")
