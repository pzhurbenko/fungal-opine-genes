#!/usr/bin/env python3
"""Stage 05 — validate every target model DIRECTLY against the T-DNA references (independent of the guide that built it).

  1. per family: MAFFT (T1+T2) -> hmmbuild -> family HMM   (round<N>/hmm/<FAM>.hmm)
  2. hmmsearch of all target proteins (round<N>/loci/*/target_protein.faa) + NCBI protein subjects (registry) against the HMM
  3. status per locus:
       no_model          miniprot produced nothing overlapping the HSP
       not_homolog       HMM E-value > validate.hmm_evalue AND miniprot identity < 0.25
       fragment          protein shorter than validate.fragment_frac x median T1 length
       pseudogene_like   frameshift > 0 or internal stop codons > 0
       intact            everything else
Outputs: round<N>/validation.tsv ; round<N>/hmm/ ; round<N>/protein_subjects_hmm.tsv
"""
import os, sys, json, glob, subprocess, statistics
from common import load_config, round_dir, read_tsv, write_tsv, read_fasta, write_fasta, eutils, log

cfg = load_config(); P = cfg['_paths']; V = cfg['validate']
rnd = sys.argv[1]; RD = round_dir(cfg, rnd); HD = os.path.join(RD, 'hmm'); os.makedirs(HD, exist_ok=True)

# ---- HMMs: one per SUBGROUP (guides/<FAM>__<sub>_seed.faa, written by stage 03); family-level = min over subgroups
import glob as _glob
median_len = {}; subgroups = {}
for fam in cfg['families']:
    t1 = read_fasta(os.path.join(P['guides'], f'{fam}_T1.faa')); t2 = read_fasta(os.path.join(P['guides'], f'{fam}_T2.faa'))
    median_len[fam] = statistics.median(len(s) for _, s in t1)
    seeds = sorted(_glob.glob(os.path.join(P['guides'], f'{fam}__*_seed.faa')))
    subgroups[fam] = []
    for seed in seeds:
        sg = os.path.basename(seed)[len(fam) + 2:-len('_seed.faa')]
        if sg == 'unassigned': continue
        recs = read_fasta(seed)
        if len(recs) < 1: continue
        aln = os.path.join(HD, f'{fam}__{sg}.aln'); hmm = os.path.join(HD, f'{fam}__{sg}.hmm')
        if not os.path.exists(hmm):
            if len(recs) == 1: write_fasta(aln, recs)          # single sequence -> hmmbuild --singlemx
            else: subprocess.run(['mafft', '--auto', '--quiet', '--thread', '4', seed], stdout=open(aln, 'w'), check=True)
            cmd = ['hmmbuild', '--amino', '-n', f'{fam}__{sg}'] + (['--singlemx'] if len(recs) == 1 else []) + [hmm, aln]
            subprocess.run(cmd, stdout=subprocess.DEVNULL, check=True)
        subgroups[fam].append(sg)
    og_seed = os.path.join(P['guides'], f'{fam}_outgroup_seed.faa')
    if os.path.exists(og_seed) and os.path.getsize(og_seed) > 0:
        recs = read_fasta(og_seed); aln = os.path.join(HD, f'{fam}__OUTGROUP.aln'); hmm = os.path.join(HD, f'{fam}__OUTGROUP.hmm')
        if not os.path.exists(hmm):
            if len(recs) == 1: write_fasta(aln, recs)
            else: subprocess.run(['mafft', '--auto', '--quiet', '--thread', '4', og_seed], stdout=open(aln, 'w'), check=True)
            subprocess.run(['hmmbuild', '--amino', '-n', f'{fam}__OUTGROUP'] + (['--singlemx'] if len(recs) == 1 else []) + [hmm, aln], stdout=subprocess.DEVNULL, check=True)
        subgroups[fam].append('OUTGROUP')
    log(f'{fam}: subgroup HMMs {subgroups[fam]} ; T1={len(t1)} T2={len(t2)}, median T1 length {median_len[fam]}')

# ---- collect targets
targets = []
for f in sorted(glob.glob(os.path.join(RD, 'loci', '*', 'target_protein.faa'))):
    targets += read_fasta(f)
allt = os.path.join(RD, 'all_targets.faa'); write_fasta(allt, [(h, s.replace('*', 'X')) for h, s in targets])

def hmmsearch_one(hmm, fasta, out):
    if os.path.getsize(fasta) == 0: return {}
    subprocess.run(['hmmsearch', '--tblout', out, '--noali', '-E', '10', hmm, fasta], stdout=subprocess.DEVNULL, check=True)
    res = {}
    for line in open(out):
        if line.startswith('#'): continue
        f = line.split()
        res[f[0]] = {'evalue': float(f[4]), 'score': float(f[5])}
    return res

def hmmsearch(fam, fasta, out):
    """Search all subgroup HMMs + OUTGROUP of a family. Per query: best NON-outgroup subgroup ('score','evalue','subgroup'),
    'og_score' (OUTGROUP, 0 if no hit) and 'gap' = score - og_score. 'subgroup' is set to OUTGROUP when the outgroup wins."""
    per = {}
    for sg in subgroups[fam]:
        r = hmmsearch_one(os.path.join(HD, f'{fam}__{sg}.hmm'), fasta, out.replace('.tbl', f'__{sg}.tbl'))
        for q, v in r.items(): per.setdefault(q, {})[sg] = v
    res = {}
    for q, d in per.items():
        og = d.get('OUTGROUP', {'score': 0.0})['score']
        subs = {k: v for k, v in d.items() if k != 'OUTGROUP'}
        if subs:
            bs = max(subs, key=lambda k: subs[k]['score']); b = subs[bs]
            res[q] = {'score': b['score'], 'evalue': b['evalue'], 'best_sub': bs, 'og_score': og, 'gap': b['score'] - og,
                      'subgroup': 'OUTGROUP' if og > b['score'] else bs}
        else:
            res[q] = {'score': 0.0, 'evalue': 99.0, 'best_sub': '', 'og_score': og, 'gap': -og, 'subgroup': 'OUTGROUP'}
    return res

hmm_hits = {fam: hmmsearch(fam, allt, os.path.join(HD, f'{fam}_targets.tbl')) for fam in cfg['families']}

# ---- status per locus
loci = read_tsv(os.path.join(P['registry'], 'loci.tsv'))
rows = []
for L in loci:
    lid = L['locus_id']; fam = L['fam']; jp = os.path.join(RD, 'loci', lid, 'locus.json')
    if not os.path.exists(jp): continue
    J = json.load(open(jp)); b = J['miniprot_best']; src = J.get('model_source')
    r = {'locus_id': lid, 'fam': fam, 'organism': L['organism'], 'contig': L['contig'], 'contig_len': J['contig_len'],
         'augustus_species': J['augustus_species'], 'was_deleted': L['was_deleted'], 'model_source': src or ''}
    if not src:
        r.update(status='no_model'); rows.append(r); continue
    if not b and src == 'augustus':   # AUGUSTUS fallback model
        g = max(J['augustus_target'], key=lambda g: g['protein_len'])
        b = {'guide_tier': 'NA', 'guide': g['id'], 'identity': 0.0, 'n_exons': g['exons'], 'protein_len': g['protein_len'], 'frameshift': 0, 'stop_codons': 0,
             'contig_start': g['contig_start'], 'contig_end': g['contig_end'], 'strand': g['strand']}
    elif not b:   # hsp_concat / hsp_translation remnant
        tp = read_fasta(os.path.join(RD, 'loci', lid, 'target_protein.faa'))
        plen = len(tp[0][1]) if tp else 0
        import re as _re
        nst = _re.search(r'\|stop=(\d+)\|', tp[0][0]) if tp else None
        b = {'guide_tier': 'NA', 'guide': tp[0][0].split('guide=')[-1] if tp else '', 'identity': 0.0, 'n_exons': '', 'protein_len': plen, 'frameshift': '', 'stop_codons': int(nst.group(1)) if nst else '',
             'contig_start': L['hsp_start'], 'contig_end': L['hsp_end'], 'strand': ''}
    hdr = [h for h, _ in targets if h.startswith(lid + '|')]
    hh = hmm_hits[fam].get(hdr[0]) if hdr else None
    own_E = hh['evalue'] if hh else 99.0
    other = min((hmm_hits[f2].get(hdr[0], {'evalue': 99.0})['evalue'] for f2 in cfg['families'] if f2 != fam), default=99.0) if hdr else 99.0
    bt = J.get('miniprot_best_by_tier', {}); t12 = [bt[k] for k in ('T1', 'T2') if k in bt]
    best12 = max(t12, key=lambda m: m['identity']) if t12 else None
    r.update(guide_tier=b['guide_tier'], guide=b['guide'], identity=f"{b['identity']:.2f}", n_exons=b['n_exons'], protein_len=b['protein_len'],
             frameshift=b['frameshift'], stop_codons=b['stop_codons'], contig_start=b['contig_start'], contig_end=b['contig_end'], strand=b['strand'],
             identity_T1T2=f"{best12['identity']:.2f}" if best12 else '', guide_T1T2=best12['guide'][:60] if best12 else '',
             hmm_evalue=f'{own_E:.1e}', hmm_score=hh['score'] if hh else '', hmm_other_family_evalue=f'{other:.1e}',
             subgroup=hh['subgroup'] if hh else '', best_subgroup=hh['best_sub'] if hh else '',
             og_score=f"{hh['og_score']:.1f}" if hh else '', hmm_gap_bits=f"{hh['gap']:.1f}" if hh else '',
             len_frac_of_T1_median=f"{b['protein_len'] / median_len[fam]:.2f}",
             augustus_target_exons=';'.join(str(g['exons']) for g in J['augustus_target']), augustus_target_len=';'.join(str(g['protein_len']) for g in J['augustus_target']),
             n_neighbors=J['n_neighbors'], neighbors_multiexon=sum(1 for e in J['neighbors_exons'] if e >= 2),
             adaptive_left_bp=J['adaptive_neighbors'].get('left', {}).get('distance_bp', ''), adaptive_right_bp=J['adaptive_neighbors'].get('right', {}).get('distance_bp', ''),
             gc_model_window=J['gc_model_window'], gc_context=J['gc_context'])
    GAP = V.get('ambiguous_gap_bits', 20); MINS = V.get('min_score_bits', 30)
    if src in ('hsp_translation', 'hsp_concat'):                          # remnant: no gene model, only translated HSP
        st = 'remnant' if (hh and hh['evalue'] <= V['hmm_evalue'] and hh['gap'] > 0) else 'not_homolog'
    elif hh and hh['gap'] <= -GAP: st = 'not_homolog'                    # decoy HMM (native/chromosomal paralog) clearly wins
    elif hh and (abs(hh['gap']) < GAP or hh['score'] < MINS): st = 'ambiguous'   # nobody wins clearly / both weak -> tree decides
    elif own_E > V['hmm_evalue'] and b['identity'] < 0.25: st = 'not_homolog'
    elif own_E > V['hmm_evalue'] and src == 'augustus': st = 'not_homolog'  # AUGUSTUS-only model must pass the HMM
    elif not hh and src == 'augustus': st = 'not_homolog'
    elif b['protein_len'] < V['fragment_frac'] * median_len[fam]: st = 'fragment'
    elif (b['frameshift'] or 0) > 0 or (b['stop_codons'] or 0) > 0: st = 'pseudogene_like'
    else: st = 'intact'
    r['status'] = st; rows.append(r)
write_tsv(os.path.join(RD, 'validation.tsv'), rows)
from collections import Counter
log('status:', dict(Counter((r['fam'], r['status']) for r in rows)))

# ---- NCBI protein subjects (already gene models): fetch once, hmmsearch, mark T3 eligibility
ps = read_tsv(os.path.join(P['registry'], 'protein_subjects.tsv'))
pf = os.path.join(P['registry'], 'protein_subjects.faa')
if ps and not os.path.exists(pf):
    accs = sorted(set(r['protein_acc'] for r in ps)); recs = []
    for i in range(0, len(accs), 150):
        txt = eutils(cfg, 'efetch.fcgi', {'db': 'protein', 'id': ','.join(accs[i:i+150]), 'rettype': 'fasta', 'retmode': 'text'})
        tmp = pf + '.part'; open(tmp, 'w').write(txt); recs += read_fasta(tmp)
    write_fasta(pf, recs); log(f'protein subjects fetched: {len(recs)}/{len(accs)}')
if os.path.exists(pf):
    fam_of = {}
    for r in ps: fam_of.setdefault(r['protein_acc'].split('.')[0], set()).add(r['fam'])
    out = []
    for fam in cfg['families']:
        hits = hmmsearch(fam, pf, os.path.join(HD, f'{fam}_protein_subjects.tbl'))
        for h, s in read_fasta(pf):
            acc = h.split()[0]
            if fam in fam_of.get(acc.split('.')[0], set()):
                e = hits.get(acc, {'evalue': 99.0, 'subgroup': '', 'score': 0.0, 'og_score': 0.0, 'gap': -99.0}); out.append({'fam': fam, 'protein_acc': acc, 'length': len(s), 'hmm_evalue': f"{e['evalue']:.1e}", 'subgroup': e.get('subgroup', ''),
                                                                 'og_score': f"{e.get('og_score', 0):.1f}", 'gap': f"{e.get('gap', 0):.1f}",
                                                                 'hmm_pass': 'yes' if (e['evalue'] <= V['hmm_evalue'] and e.get('gap', -99) >= V.get('ambiguous_gap_bits', 20) and e.get('score', 0) >= V.get('min_score_bits', 30)) else ('ambiguous' if e.get('subgroup') and abs(e.get('gap', -99)) < V.get('ambiguous_gap_bits', 20) else 'no'), 'title': h.split(' ', 1)[1] if ' ' in h else ''})
    write_tsv(os.path.join(RD, 'protein_subjects_hmm.tsv'), out)
    log('protein subjects HMM pass:', dict(Counter((r['fam'], r['hmm_pass']) for r in out)))
