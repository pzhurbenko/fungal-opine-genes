#!/usr/bin/env python3
"""Stage 02 — download contig sequences + assembly metadata for every contig in the registry.

Rules (config.fetch):
  * contig length <= contig_full_max  -> whole contig
  * else                              -> window [min(hsp)-flank_context, max(hsp)+flank_context]
Resumable: contigs already present in contigs.tsv are skipped. Rate-limited via common.eutils.

Gotcha (2026-09): freshly released WGS records are NOT in the esearch/esummary index for weeks,
but efetch by 'ACC.1' works — we fall back to efetch and derive length from the sequence.

Outputs (out/contigs/): <acc>.fa ; contigs.tsv (acc, acc_ver, length, organism, fetched=full|window, seq_start, seq_end, path)
         (out/registry/) assemblies.tsv (wgs_prefix, gca, organism, status, annotated, refseq, submitter, gc_percent, total_len)
"""
import os, sys, json, argparse
from common import load_config, read_tsv, write_tsv, eutils, wgs_prefix, log

ap = argparse.ArgumentParser(); ap.add_argument('--only', help='file with contig accessions to restrict to (testing)')
ap.add_argument('--skip-assemblies', action='store_true'); args = ap.parse_args()

cfg = load_config(); P = cfg['_paths']; os.makedirs(P['contigs'], exist_ok=True)
FULLMAX = cfg['fetch']['contig_full_max']; FLANK = cfg['fetch']['flank_context']
reg = read_tsv(os.path.join(P['registry'], 'contigs.tsv'))
if args.only:
    keep = set(l.strip() for l in open(args.only) if l.strip()); reg = [r for r in reg if r['contig'] in keep]
meta_path = os.path.join(P['contigs'], 'contigs.tsv')
done = {r['contig']: r for r in read_tsv(meta_path)} if os.path.exists(meta_path) else {}
log(f'contigs to fetch: {len(reg)}, already done: {sum(1 for r in reg if r["contig"] in done)}')

def summary(acc):
    s = eutils(cfg, 'esummary.fcgi', {'db': 'nuccore', 'id': acc, 'retmode': 'json'}, as_json=True)['result']
    if not s.get('uids'): return None
    d = s[s['uids'][0]]
    return {'acc_ver': d['accessionversion'], 'length': int(d['slen']), 'organism': d.get('organism', '')}

def fetch_seq(acc_ver, start=None, stop=None):
    p = {'db': 'nuccore', 'id': acc_ver, 'rettype': 'fasta', 'retmode': 'text'}
    if start: p.update(seq_start=start, seq_stop=stop)
    txt = eutils(cfg, 'efetch.fcgi', p)
    if not txt.startswith('>'): raise RuntimeError(f'efetch returned no FASTA for {acc_ver}')
    hdr, seq = txt.split('\n', 1)
    return hdr[1:], seq.replace('\n', '')

n_new = 0
for r in reg:
    acc = r['contig']
    if acc in done: continue
    try:
        info = summary(acc)
        if info is None:                      # not indexed yet -> efetch whole record directly
            hdr, seq = fetch_seq(acc + '.1')
            info = {'acc_ver': hdr.split()[0], 'length': len(seq), 'organism': ' '.join(hdr.split()[1:3])}
            full = (hdr, seq)
        else:
            full = None
        L = info['length']; hmin, hmax = int(r['hsp_min']), int(r['hsp_max'])
        if L <= FULLMAX:
            if full is None: full = fetch_seq(info['acc_ver'])
            hdr, seq = full; s, e, mode = 1, L, 'full'
        else:
            s, e = max(1, hmin - FLANK), min(L, hmax + FLANK); mode = 'window'
            if full is not None: hdr, seq = full[0], full[1][s-1:e]
            else: hdr, seq = fetch_seq(info['acc_ver'], s, e)
        path = os.path.join(P['contigs'], f'{acc}.fa')
        with open(path, 'w') as f:
            f.write(f'>{acc} {info["acc_ver"]} {info["organism"]} seq_start={s} seq_end={e} contig_len={L}\n')
            for i in range(0, len(seq), 80): f.write(seq[i:i+80] + '\n')
        done[acc] = {'contig': acc, 'acc_ver': info['acc_ver'], 'length': L, 'organism': info['organism'],
                     'fetched': mode, 'seq_start': s, 'seq_end': e, 'path': path, 'wgs_prefix': wgs_prefix(acc)}
        n_new += 1
        if n_new % 25 == 0:
            write_tsv(meta_path, list(done.values())); log(f'  fetched {n_new} ...')
    except Exception as ex:
        log(f'  FAILED {acc}: {ex}')
write_tsv(meta_path, list(done.values()))
log(f'contigs fetched now: {n_new}; total in table: {len(done)}')

# ---------------------------------------------------------------- assemblies (annotation availability, GC)
if not args.skip_assemblies:
    asm_path = os.path.join(P['registry'], 'assemblies.tsv')
    asm = {r['wgs_prefix']: r for r in read_tsv(asm_path)} if os.path.exists(asm_path) else {}
    prefixes = sorted(set(d['wgs_prefix'] for d in done.values() if d.get('wgs_prefix')) - set(asm))
    log(f'assemblies to resolve: {len(prefixes)}')
    for i, p in enumerate(prefixes):
        try:
            ids = eutils(cfg, 'esearch.fcgi', {'db': 'assembly', 'term': f'{p}[WGS]', 'retmode': 'json'}, as_json=True)['esearchresult']['idlist']
            if not ids: asm[p] = {'wgs_prefix': p, 'gca': '', 'organism': '', 'status': 'UNRESOLVED', 'annotated': '', 'refseq': '', 'submitter': '', 'gc_percent': '', 'total_len': ''}; continue
            d = eutils(cfg, 'esummary.fcgi', {'db': 'assembly', 'id': ids[0], 'retmode': 'json'}, as_json=True)['result'][ids[0]]
            gcp = ''
            try:   # GC% and total length live in the Datasets API report
                import urllib.request
                rep = json.load(urllib.request.urlopen(f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{d['assemblyaccession']}/dataset_report", timeout=60))
                st = rep['reports'][0]['assembly_stats']; gcp = st.get('gc_percent', ''); tl = st.get('total_sequence_length', '')
            except Exception: tl = ''
            asm[p] = {'wgs_prefix': p, 'gca': d.get('assemblyaccession', ''), 'organism': d.get('speciesname', ''),
                      'status': d.get('assemblystatus', ''), 'annotated': 'yes' if 'has_annotation' in d.get('propertylist', []) else 'no',
                      'refseq': 'yes' if d.get('rsuid') else 'no', 'submitter': d.get('submitterorganization', ''),
                      'gc_percent': gcp, 'total_len': tl}
        except Exception as ex:
            log(f'  assembly FAILED {p}: {ex}')
        if (i + 1) % 25 == 0: write_tsv(asm_path, list(asm.values())); log(f'  assemblies {i+1}/{len(prefixes)}')
    write_tsv(asm_path, list(asm.values()))
    log(f'assemblies: {len(asm)}, annotated: {sum(1 for a in asm.values() if a["annotated"]=="yes")}')
