#!/usr/bin/env python3
"""Stage 03 — guide proteins for miniprot, in three tiers (provenance is recorded per model later).

  T1  Agro_/Plant_ references from analysis/curated/<FAM>.fasta (the original query set, after quarantine)
  T2  Agrobacterium/Rhizobium best hits listed in rep_seqs (Subj_Accession_agro), fetched from NCBI protein,
      kept only if identity >= guides.tier2_min_identity to at least one T1 protein
      -> removes chromosomal paralogs of the same fold (OCD, reductases, TIR ...). Review table is written.
  T3  fungal full-length models validated in a previous round (written by stage 07), if present.

Header convention (miniprot reports Target=<header>):  T1|<original header> ; T2|<acc>|<organism> ; T3|<locus_id>
Outputs (out/guides/): <FAM>_T1.faa <FAM>_T2.faa <FAM>_T2_review.tsv <FAM>_guides.faa (T1+T2+T3 concatenated)
"""
import os, re, glob, sys
from common import load_config, read_fasta, write_fasta, read_tsv, write_tsv, eutils, log
from Bio import Align
from Bio.Align import substitution_matrices

cfg = load_config(); P = cfg['_paths']; os.makedirs(P['guides'], exist_ok=True)
MIN_ID = cfg['guides']['tier2_min_identity']
al = Align.PairwiseAligner(); al.substitution_matrix = substitution_matrices.load('BLOSUM62')
al.open_gap_score = -10; al.extend_gap_score = -0.5; al.mode = 'local'

def pid(a, b):
    """Identity normalised by the SHORTER sequence (identical residues / min length).
    Plain local identity over aligned columns is inflated for a short conserved core shared by all
    proteins of a fold (OCD vs T-DNA sus gave 0.57 that way; normalised = 0.35)."""
    x = al.align(a, b)[0]; cols = ident = 0
    for p, q in zip(x[0], x[1]):
        if p != '-' and q != '-': cols += 1; ident += (p == q)
    return ident / min(len(a), len(b)), cols

quarantined = set()
if os.path.exists(P['quarantine_tsv']):
    quarantined = {r['clean_id'] for r in read_tsv(P['quarantine_tsv'])}

for fam in cfg['families']:
    # ---- T1
    excl = cfg['guides'].get('exclude_T1', {}).get(fam, []) or []
    t1_all = [(h, s.replace('*', '')) for h, s in read_fasta(os.path.join(P['curated_dir'], f'{fam}.fasta'))
              if h.startswith(('Agro_', 'Plant_')) and h not in quarantined]
    t1 = [(f'T1|{h}', s) for h, s in t1_all if not any(x in h for x in excl)]
    extra = cfg['guides'].get('extra_T1', {}).get(fam, []) or []
    if extra:   # additional T1 references by NCBI accession (e.g. cT-DNA fragments found by context inspection)
        txt = eutils(cfg, 'efetch.fcgi', {'db': 'protein', 'id': ','.join(extra), 'rettype': 'fasta', 'retmode': 'text'})
        tmp = os.path.join(P['guides'], f'{fam}_extraT1_raw.faa'); open(tmp, 'w').write(txt)
        for h, s in read_fasta(tmp):
            org = re.search(r'\[(.*?)\]', h); org = org.group(1).replace(' ', '_') if org else 'unknown'
            kind = 'Agro' if re.search(r'Rhizob|Agrobact|Ensifer|Sinorhiz|Mesorhiz|Neorhiz', org) else 'Plant'
            t1.append((f'T1|{kind}_{h.split()[0].replace(".", "_")}_{org}_extra', s.replace('*', '')))
        log(f'{fam}: extra T1 added: {len(extra)}')
    t1_ex = [(f'T1excluded|{h}', s) for h, s in t1_all if any(x in h for x in excl)]
    write_fasta(os.path.join(P['guides'], f'{fam}_T1.faa'), t1)
    write_fasta(os.path.join(P['guides'], f'{fam}_T1_excluded.faa'), t1_ex)   # outgroups for trees, not guides
    if t1_ex: log(f'{fam}: excluded from T1 (config guides.exclude_T1): {len(t1_ex)}')
    # ---- T2 candidates from rep_seqs
    rp = glob.glob(os.path.join(P['rep_seqs_dir'], f'{fam}_rep_seqs.tsv')) + glob.glob(os.path.join(P['rep_seqs_dir'], f'{fam.lower()}_rep_seqs.tsv'))
    cands = {}
    for r in read_tsv(rp[0]) if rp else []:
        a = r.get('Subj_Accession_agro', '').strip(); n = r.get('Subj_Name_agro', '').strip()
        if a: cands.setdefault(a, {'acc': a, 'name': n, 'n_reps': 0}); cands[a]['n_reps'] += 1
    have = {re.sub(r'\.\d+$', '', m.group(1)) for h, _ in t1 for m in [re.search(r'(?:Agro|Plant)_([A-Z]{2,3}_?\d{5,9}(?:\.\d)?)', h)] if m}
    new = [a for a in cands if re.sub(r'\.\d+$', '', a) not in have]
    t2, review = [], []
    if new:
        txt = eutils(cfg, 'efetch.fcgi', {'db': 'protein', 'id': ','.join(new), 'rettype': 'fasta', 'retmode': 'text'})
        tmp = os.path.join(P['guides'], f'{fam}_T2_raw.faa'); open(tmp, 'w').write(txt)
        for h, s in read_fasta(tmp):
            acc = h.split()[0]; base = re.sub(r'\.\d+$', '', acc)
            best = max(((pid(s, t)[0], th) for th, t in t1), default=(0, ''))
            org = re.search(r'\[(.*?)\]', h); org = org.group(1) if org else ''
            c = cands.get(base) or cands.get(acc) or {'name': h, 'n_reps': 0}
            keep = best[0] >= MIN_ID
            review.append({'fam': fam, 'acc': acc, 'organism': org, 'ncbi_title': h.split(' ', 1)[1] if ' ' in h else '',
                           'n_reps_best_hit': c['n_reps'], 'max_identity_to_T1': f'{best[0]:.2f}', 'closest_T1': best[1][:60],
                           'length': len(s), 'kept_T2': 'yes' if keep else 'no'})
            if keep: t2.append((f'T2|{acc}|{org.replace(" ", "_")}', s.replace('*', '')))
    write_fasta(os.path.join(P['guides'], f'{fam}_T2.faa'), t2)
    write_tsv(os.path.join(P['guides'], f'{fam}_T2_review.tsv'), sorted(review, key=lambda r: -float(r['max_identity_to_T1'])))
    # ---- subgroups_extra: additional lineages given as NCBI accessions (e.g. chromosomal OCD of Rhizobiales, CUS-like genes of
    #      Actinobacteria) -> fetched, added to guides as tier T2x and to their OWN subgroup HMM (competitive with the T-DNA lineages)
    extra_sg = cfg['guides'].get('subgroups_extra', {}).get(fam, {}) or {}
    t2x = []; t2x_assign = {}
    for sgname, accs in extra_sg.items():
        if not accs: continue
        txt = eutils(cfg, 'efetch.fcgi', {'db': 'protein', 'id': ','.join(accs), 'rettype': 'fasta', 'retmode': 'text'})
        tmp = os.path.join(P['guides'], f'{fam}_{sgname}_raw.faa'); open(tmp, 'w').write(txt)
        for h, s in read_fasta(tmp):
            org = re.search(r'\[(.*?)\]', h); org = org.group(1).replace(' ', '_') if org else ''
            hh = f'T2x|{h.split()[0]}|{org}'; t2x.append((hh, s.replace('*', ''))); t2x_assign[hh] = sgname
        log(f'{fam}: subgroups_extra {sgname}: {sum(1 for v in t2x_assign.values() if v == sgname)}')
    # ---- subgroups: T1 by config substrings, T2 by closest T1; per-subgroup seed fasta for the HMMs (stage 05)
    SG = cfg['guides'].get('subgroups', {}).get(fam) or {fam.lower(): None}
    def sub_of_t1(h):
        for sg, keys in SG.items():
            if keys is None or any(k in h for k in keys): return sg
        return 'unassigned'
    assign = {h: sub_of_t1(h) for h, _ in t1}
    t1d = dict(t1)
    for h, s in t2:
        best = max(((pid(s, t1d[h1])[0], h1) for h1 in t1d), default=(0, None))
        assign[h] = assign[best[1]] if best[1] and best[0] >= MIN_ID else 'unassigned'
    assign.update(t2x_assign)
    sg_rows = [{'fam': fam, 'subgroup': sg, 'tier': h.split('|')[0], 'header': h[:100]} for h, sg in assign.items()]
    write_tsv(os.path.join(P['guides'], f'{fam}_subgroups.tsv'), sorted(sg_rows, key=lambda r: (r['subgroup'], r['tier'])))
    for old_seed in glob.glob(os.path.join(P['guides'], f'{fam}__*_seed.faa')): os.remove(old_seed)   # stale subgroups
    for sg in set(assign.values()):
        recs = [(h, s) for h, s in t1 + t2 + t2x if assign[h] == sg]
        write_fasta(os.path.join(P['guides'], f'{fam}__{sg}_seed.faa'), recs)
    from collections import Counter
    log(f'{fam}: subgroups ' + ', '.join(f'{k}={v}' for k, v in sorted(Counter(assign.values()).items())))
    # ---- outgroup seed (decoy HMM): excluded T1 refs by substring + NCBI accessions
    og_keys = cfg['guides'].get('outgroups', {}).get(fam, []) or []
    og = [(f'OG|{h.split("|",1)[1]}', s) for h, s in t1_ex if any(k in h for k in og_keys)]
    accs = [k for k in og_keys if re.match(r'^[A-Z]{2,3}_?\d{5,9}\.\d+$', k)]
    if accs:
        txt = eutils(cfg, 'efetch.fcgi', {'db': 'protein', 'id': ','.join(accs), 'rettype': 'fasta', 'retmode': 'text'})
        tmp = os.path.join(P['guides'], f'{fam}_outgroup_raw.faa'); open(tmp, 'w').write(txt)
        og += [(f'OG|{h.split()[0]}|{(re.search(r"\[(.*?)\]", h) or [None, ""])[1].replace(" ", "_")}', s.replace('*', '')) for h, s in read_fasta(tmp)]
    write_fasta(os.path.join(P['guides'], f'{fam}_outgroup_seed.faa'), og)
    log(f'{fam}: outgroup seed {len(og)}')
    # ---- T3 (from previous round, optional)
    t3p = os.path.join(P['guides'], f'{fam}_T3.faa'); t3 = read_fasta(t3p) if os.path.exists(t3p) else []
    write_fasta(os.path.join(P['guides'], f'{fam}_guides.faa'), t1 + t2 + t2x + t3)
    log(f'{fam}: T1={len(t1)} T2={len(t2)}/{len(review)} T2x={len(t2x)} T3={len(t3)} -> {fam}_guides.faa')
