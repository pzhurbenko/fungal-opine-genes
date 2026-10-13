#!/usr/bin/env python3
"""Stage 07 — final per-locus table, per-family full-length protein/CDS sets, T3 guides for the next round.

Outputs (round<N>/):
  loci_table.tsv          registry + validation + context, one row per locus, with final 'status' and 'context_call'
  proteins/<FAM>_models.faa   full-length target proteins (status intact|pseudogene_like), header = locus_id|status|tier|...
  proteins/<FAM>_models.fna   matching CDS nucleotide sequences (for GC3 / codon usage / dN/dS)
  ../guides/<FAM>_T3.faa      validated fungal models -> guides for round N+1 (intact, HMM pass, identity >= tier3_min_identity)
context_call: fungal_context | mixed | bacterial_context | no_neighbors | short_contig
"""
import os, sys, glob, json
from collections import Counter
from common import load_config, round_dir, read_tsv, write_tsv, read_fasta, write_fasta, log

cfg = load_config(); P = cfg['_paths']; rnd = sys.argv[1]; RD = round_dir(cfg, rnd)
reg = {r['locus_id']: r for r in read_tsv(os.path.join(P['registry'], 'loci.tsv'))}
val = {r['locus_id']: r for r in read_tsv(os.path.join(RD, 'validation.tsv'))}
ctxp = os.path.join(RD, 'locus_context.tsv'); ctx = {r['locus_id']: r for r in read_tsv(ctxp)} if os.path.exists(ctxp) else {}
asm_p = os.path.join(P['registry'], 'assemblies.tsv'); asm = {r['wgs_prefix']: r for r in read_tsv(asm_p)} if os.path.exists(asm_p) else {}
asm_by_gca = {r['gca']: r for r in asm.values() if r.get('gca')}   # step-2 loci on non-WGS (chromosome) accessions carry 'gca' in loci.tsv
from common import wgs_prefix

rows = []
for lid, L in reg.items():
    v = val.get(lid, {}); c = ctx.get(lid, {}); a = asm.get(wgs_prefix(L['contig']) or '', {}) or asm_by_gca.get(L.get('gca', ''), {})
    r = dict(L); r.update({k: v.get(k, '') for k in v if k not in r}); r.update({k: c.get(k, '') for k in c if k not in r})
    r['gca'] = a.get('gca', ''); r['assembly_annotated'] = a.get('annotated', ''); r['assembly_gc'] = a.get('gc_percent', '')
    r['status'] = v.get('status', 'not_run')   # intact | pseudogene_like | fragment | remnant | ambiguous | not_homolog | no_model
    r['origin'] = (cfg['guides'].get('origin', {}) or {}).get(v.get('subgroup', ''), '')   # T-DNA | catabolic_locus | bacterial_chromosomal | ...
    nf, nb, nn = int(c.get('nb_fungi', 0) or 0), int(c.get('nb_bacteria', 0) or 0), int(c.get('neighbors_in_window', 0) or 0)
    if not c: call = ''
    elif int(v.get('contig_len', 0) or 0) < 10000: call = 'short_contig'
    elif nn == 0: call = 'no_neighbors'
    elif nf >= 2 and nb == 0: call = 'fungal_context'
    elif nb >= 2 and nf == 0: call = 'bacterial_context'
    else: call = 'mixed'
    r['context_call'] = call
    try: r['gc_delta_vs_assembly'] = f"{float(v['gc_context']) - float(a['gc_percent']):.1f}" if v.get('gc_context') and a.get('gc_percent') else ''
    except Exception: r['gc_delta_vs_assembly'] = ''
    rows.append(r)
cols = ['locus_id', 'fam', 'subgroup', 'origin', 'status', 'context_call', 'model_source', 'organism', 'class', 'order', 'genome', 'gca', 'assembly_annotated', 'contig', 'contig_len',
        'hsp_start', 'hsp_end', 'contig_start', 'contig_end', 'strand', 'n_exons', 'protein_len', 'len_frac_of_T1_median', 'frameshift', 'stop_codons',
        'guide_tier', 'identity', 'guide', 'identity_T1T2', 'guide_T1T2', 'hmm_evalue', 'hmm_score', 'og_score', 'hmm_gap_bits', 'best_subgroup', 'hmm_other_family_evalue', 'augustus_species', 'augustus_target_exons', 'augustus_target_len',
        'n_neighbors', 'neighbors_multiexon', 'neighbors_in_window', 'nb_fungi', 'nb_bacteria', 'nb_plants', 'nb_TE', 'nb_no_hit', 'nb_tdna_loci', 'window_left_bp', 'window_right_bp',
        'gc_model_window', 'gc_context', 'assembly_gc', 'gc_delta_vs_assembly', 'target_top_org', 'target_top_kingdom',
        'best_identity', 'best_evalue', 'n_hsp', 'queries', 'was_deleted', 'reasons', 'tree_clade', 'is_rep_contig', 'other_fams_on_contig']
for r in rows:
    for k in cols: r.setdefault(k, '')
write_tsv(os.path.join(RD, 'loci_table.tsv'), rows, cols)
log('status:', dict(Counter((r['fam'], r['status']) for r in rows)))
log('context:', dict(Counter(r['context_call'] for r in rows if r['context_call'])))

# ---- per-family model sets + T3 guides
os.makedirs(os.path.join(RD, 'proteins'), exist_ok=True)
MIN3 = cfg['guides']['tier3_min_identity']; HE = cfg['validate']['hmm_evalue']
for fam in cfg['families']:
    prot, cds, t3 = [], [], []
    for r in rows:
        if r['fam'] != fam or r['status'] not in ('intact', 'pseudogene_like'): continue
        d = os.path.join(RD, 'loci', r['locus_id'])
        pf, cf = os.path.join(d, 'target_protein.faa'), os.path.join(d, 'target_cds.fna')
        if not os.path.exists(pf): continue
        h, s = read_fasta(pf)[0]; hdr = f"{r['locus_id']}|{r['subgroup']}|{r['status']}|{r['guide_tier']}|{r['organism'].replace(' ', '_')}|len={len(s)}|ex={r['n_exons']}|id={r['identity']}"
        prot.append((hdr, s))
        if os.path.exists(cf): cds.append((hdr, read_fasta(cf)[0][1]))
        if r['status'] == 'intact' and float(r['hmm_evalue'] or 99) <= HE and float(r['identity_T1T2'] or r['identity'] or 0) >= MIN3 and r['context_call'] not in ('bacterial_context',):
            t3.append((f"T3|{r['locus_id']}|{r['organism'].replace(' ', '_')}", s.replace('*', '')))
    write_fasta(os.path.join(RD, 'proteins', f'{fam}_models.faa'), prot); write_fasta(os.path.join(RD, 'proteins', f'{fam}_models.fna'), cds)
    rem = []
    for r in rows:
        if r['fam'] == fam and r['status'] in ('remnant', 'fragment'):
            pf = os.path.join(RD, 'loci', r['locus_id'], 'target_protein.faa')
            if os.path.exists(pf):
                h, s = read_fasta(pf)[0]; rem.append((f"{r['locus_id']}|{r['subgroup']}|{r['status']}|{r['model_source']}|{r['organism'].replace(' ', '_')}|len={len(s)}", s))
    write_fasta(os.path.join(RD, 'proteins', f'{fam}_remnants_fragments.faa'), rem)
    # T3 also from NCBI protein subjects that pass the HMM
    psh = os.path.join(RD, 'protein_subjects_hmm.tsv'); pf = os.path.join(P['registry'], 'protein_subjects.faa')
    if os.path.exists(psh) and os.path.exists(pf):
        ok = {r['protein_acc'] for r in read_tsv(psh) if r['fam'] == fam and r['hmm_pass'] == 'yes'}
        t3 += [(f"T3|{h.split()[0]}|ncbi_protein", s.replace('*', '')) for h, s in read_fasta(pf) if h.split()[0] in ok]
    write_fasta(os.path.join(P['guides'], f'{fam}_T3.faa'), t3)
    log(f'{fam}: models {len(prot)} (cds {len(cds)}), T3 guides for next round {len(t3)}')
