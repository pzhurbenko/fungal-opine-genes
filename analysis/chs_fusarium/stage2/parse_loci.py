#!/usr/bin/env python3
"""Parse miniprot runs (mp_out/<GCA>/*.gff) -> loci_all.tsv, loci.tsv (chs-like + borderline), neighbors.tsv, genomes.tsv, species_summary.tsv.
Usage: parse_loci.py <genomes_files.tsv: GCA fna gff> """
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import sys, os, re, csv, collections
S = f'{GMO_ROOT}/analysis/chs_fusarium/stage2'
SEEDS = {'chs_A0A2Z2PGA9': 407, 'EEU41297.1': 405, 'KNB14377.1': 260}
FSSC = {'solani','vanettenii','phalaenopsidis','nematophilum','solani-melongenae','petroliphilum','keratoplasticum','falciforme','pisi','striatum',
        'suttonianum','lichenicola','metavorans','cuneirostrum','ambrosium','euwallaceae','kuroshium','oligoseptatum','mori','haematococcum','piperis','paranaense','longipes','illudens','cucurbitae','robiniae','rhododendri','epiphyta','geoasiaticum','tonkinense','brasiliensis','croci','gamsii','plagianthi','noneumartii','catenisporum','liriodendri','bataticola','falciforme','vasinfecta','phaseoli','martii'}
MOB = re.compile(r'transpos|integrase|retro|helicase|reverse transcriptase', re.I)
gl = {r['gca']: r for r in csv.DictReader(open(f'{S}/genomes_list.tsv'), delimiter='\t')}

def species(org):
    t = org.replace('[', '').replace(']', '').split()
    return ' '.join(t[:2]) if len(t) > 1 else org, t[0], (t[1] if len(t) > 1 else '')
def complex_of(org):
    sp, g, e = species(org)
    if g == 'Neocosmospora' or (g == 'Fusarium' and e in FSSC): return 'FSSC/Neocosmospora'
    if sp == 'Fusarium oxysporum': return 'F. oxysporum'
    if g == 'Dactylonectria': return 'Dactylonectria'
    return 'other Fusarium'

def parse_mp(path, seed):
    out = []; cur = None; pend = ''
    for line in open(path):
        if line.startswith('##STA'): pend = line.split('\t')[1].strip(); continue
        if line[0] == '#': continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9: continue
        if f[2] == 'mRNA':
            a = dict(kv.split('=', 1) for kv in f[8].split(';') if '=' in kv); t = a['Target'].split()
            cur = dict(seed=seed, contig=f[0], start=int(f[3]), end=int(f[4]), score=int(f[5]), strand=f[6], ident=float(a['Identity']),
                       qfrom=int(t[1]), qto=int(t[2]), fs=int(a.get('Frameshift', 0)), stop=int(a.get('StopCodon', 0)), cds=[], stoppos=[i for i, ch in enumerate(pend[:-1]) if ch == '*'])
            cur['stoppos'] = [cur['qfrom'] + i for i in cur['stoppos']]; pend = ''
            cur['cov'] = (cur['qto'] - cur['qfrom'] + 1) / SEEDS[seed]; out.append(cur)
        elif f[2] == 'CDS' and cur: cur['cds'].append((int(f[3]), int(f[4])))
    return out

def read_contig(fna, name):
    seq = []; on = False
    for l in open(fna):
        if l[0] == '>': 
            if on: break
            on = l[1:].split()[0] == name
        elif on: seq.append(l.strip())
    return ''.join(seq).upper()
def gc(s):
    s = [c for c in s if c in 'ACGT']; return 100 * sum(c in 'GC' for c in s) / len(s) if s else float('nan')
comp = str.maketrans('ACGT', 'TGCA')

def parse_gff(path):
    genes = {}; mrna = {}; cdsn = collections.Counter(); cdl = collections.Counter(); prod = {}; pid = {}
    for l in open(path):
        if l[0] == '#': continue
        f = l.rstrip('\n').split('\t')
        if len(f) < 9: continue
        a = dict(kv.split('=', 1) for kv in f[8].split(';') if '=' in kv)
        if f[2] in ('gene', 'pseudogene'):
            genes[a['ID']] = dict(contig=f[0], start=int(f[3]), end=int(f[4]), strand=f[6], name=a.get('Name', a['ID']), pseudo=(f[2] == 'pseudogene' or a.get('pseudo') == 'true'), prod='', ex=0, prot='', btype=a.get('gene_biotype', ''))
        elif f[2] in ('mRNA', 'transcript', 'tRNA', 'rRNA', 'ncRNA', 'lnc_RNA', 'snoRNA', 'snRNA'):
            mrna[a['ID']] = a.get('Parent'); 
            if a.get('Parent') in genes and not genes[a['Parent']]['prod']: genes[a['Parent']]['prod'] = a.get('product', f[2])
        elif f[2] == 'CDS':
            p = a.get('Parent'); g = mrna.get(p, p)
            if g in genes:
                G = genes[g]; 
                if not G['prot']: G['prot'] = a.get('Name', '')
                if not G['prod'] or G['prod'] in ('mRNA',): G['prod'] = a.get('product', G['prod'])
                if 'product' in a: G['prod'] = a['product']
                cdsn[(g, p)] += 1; cdl[(g, p)] += int(f[4]) - int(f[3]) + 1
    for (g, p), n in cdsn.items(): genes[g]['ex'] = max(genes[g]['ex'], n); genes[g]['aa'] = max(genes[g].get('aa', 0), cdl[(g, p)] // 3 - 1)
    bc = collections.defaultdict(list)
    for gid, G in genes.items(): G['id'] = gid; bc[G['contig']].append(G)
    for c in bc: bc[c].sort(key=lambda g: (g['start'], g['end']))
    return bc

H = ['gca','organism','species','complex','level','locus_id','contig','start','end','strand','contig_len','dist_left','dist_right','id_chs','cov_chs','id_EEU','cov_EEU','id_KNB','cov_KNB',
     'cls','status','primary_seed','prim_qfrom','prim_qto','prim_len','exons','frameshifts','stops','missing','edge','gene_gc','assembly_gc',
     'annot_gene','annot_prod','annot_protein','annot_exons','annot_pseudo','annot_overlap_bp','annot_aa','seed_fs_stop','stop_pos_in_primary_seed','note']
files = [l.rstrip('\n').split('\t') for l in open(sys.argv[1])]
fo = open(f'{S}/loci_all.tsv', 'w'); fo.write('\t'.join(H) + '\n')
fn = open(f'{S}/neighbors.tsv', 'w'); fn.write('gca\tlocus_id\trel\tgene\tprotein\tproduct\texons\tdist_bp\tmobile\n')
gen = []; nmiss = 0
for gca, fna, gff in files:
    d = f'{S}/mp_out/{gca}'
    if not os.path.exists(d + '/done'): nmiss += 1; continue
    org = gl[gca]['organism']; sp, _, _ = species(org); cx = complex_of(org)
    ctg = {}; tot = 0; gcs = 0
    for l in open(d + '/contigs.tsv'):
        n, ln, g = l.split('\t')[:3]; ln = int(ln); g = float(g) if g.strip() not in ('', 'nan') else 0; ctg[n] = ln; tot += ln; gcs += ln * g
    agc = gcs / tot if tot else float('nan')
    hits = []
    for sd in SEEDS: hits += parse_mp(f'{d}/{sd}.gff', sd)
    hits.sort(key=lambda h: (h['contig'], h['start']))
    clusters = []
    for h in hits:  # merge overlapping hits on same contig and strand
        for c in clusters:
            if c['contig'] == h['contig'] and c['strand'] == h['strand'] and h['start'] <= c['end'] and h['end'] >= c['start']:
                c['hits'].append(h); c['start'] = min(c['start'], h['start']); c['end'] = max(c['end'], h['end']); break
        else: clusters.append(dict(contig=h['contig'], strand=h['strand'], start=h['start'], end=h['end'], hits=[h]))
    rows = []
    for k, c in enumerate(clusters, 1):
        best = {}
        for h in c['hits']:
            if h['seed'] not in best or h['score'] > best[h['seed']]['score']: best[h['seed']] = h
        g = lambda s, key: best[s][key] if s in best else 0
        ic, ie = g('chs_A0A2Z2PGA9', 'ident'), g('EEU41297.1', 'ident')
        cc, ce = g('chs_A0A2Z2PGA9', 'cov'), g('EEU41297.1', 'cov')
        mx = max(cc, ce, g('KNB14377.1', 'cov'))
        if (ic >= .45 or ie >= .70) and mx >= .25: cls = 'chs_like'
        elif (ic >= .35 or ie >= .55) and mx >= .25: cls = 'borderline'
        else: cls = 'ags_like_or_other'
        prim = best.get('EEU41297.1') or max(best.values(), key=lambda h: h['score'])
        L = ctg.get(c['contig'], 0); dl, dr = prim['start'] - 1, L - prim['end']
        pc = prim['cov']; ql = SEEDS[prim['seed']]
        if prim['fs'] or prim['stop']: status = 'pseudogene'
        elif pc >= 0.9: status = 'intact'
        else: status = 'truncated'
        nmis = prim['qfrom'] - 1; cmis = ql - prim['qto']; miss = []; edge = ''
        if status == 'truncated':
            plus = prim['strand'] == '+'
            for lab, n, side in (('N', nmis, dl if plus else dr), ('C', cmis, dr if plus else dl)):
                if n > 5: miss.append(f'{lab}:{n}aa'); 
                if n > 5 and side < 3 * n + 30: edge += f'{lab}-end blocked by contig edge ({side}nt); '
        if cls == 'ags_like_or_other': status = 'NA'
        r = dict(gca=gca, organism=org, species=sp, complex=cx, level=gl[gca]['level'], locus_id=f'{gca}_L{k}', contig=c['contig'], start=prim['start'], end=prim['end'], strand=prim['strand'],
                 contig_len=L, dist_left=dl, dist_right=dr, id_chs=ic, cov_chs=round(cc, 3), id_EEU=ie, cov_EEU=round(ce, 3), id_KNB=g('KNB14377.1', 'ident'), cov_KNB=round(g('KNB14377.1', 'cov'), 3),
                 cls=cls, status=status, primary_seed=prim['seed'], prim_qfrom=prim['qfrom'], prim_qto=prim['qto'], prim_len=ql, exons=len(prim['cds']), frameshifts=prim['fs'], stops=prim['stop'],
                 missing=','.join(miss), edge=edge.strip('; '), gene_gc='', assembly_gc=round(agc, 2), annot_gene='', annot_prod='', annot_protein='', annot_exons='', annot_pseudo='', annot_overlap_bp='', annot_aa='', seed_fs_stop=' '.join(f"{k[:3]}:{best[k]['fs']}/{best[k]['stop']}" for k in best), stop_pos_in_primary_seed=','.join(map(str, prim['stoppos'])), note='', _cds=prim['cds'])
        rows.append(r)
    sel = [r for r in rows if r['cls'] != 'ags_like_or_other']
    if sel:
        bc = parse_gff(gff) if gff and os.path.exists(gff) else {}
        for r in sel:
            seq = read_contig(fna, r['contig'])
            cds = ''.join(seq[s - 1:e] for s, e in sorted(r['_cds']))
            r['gene_gc'] = round(gc(cds), 2)
            gs = bc.get(r['contig'], [])
            ov = [(min(G['end'], r['end']) - max(G['start'], r['start']) + 1, G) for G in gs if G['start'] <= r['end'] and G['end'] >= r['start']]
            if ov:
                o, G = max(ov, key=lambda x: x[0]); r.update(annot_gene=G['name'], annot_prod=G['prod'], annot_protein=G['prot'], annot_exons=G['ex'], annot_pseudo=int(G['pseudo']), annot_overlap_bp=o, annot_aa=G.get('aa', ''))
            if r['status'] == 'pseudogene':
                r['note'] = (f"in-frame stop(s) at seed aa {r['stop_pos_in_primary_seed'] or 'NA'}" if r['stops'] else '') + ('; ' if r['stops'] and r['frameshifts'] else '') + (f"{r['frameshifts']} frameshift(s) vs seed" if r['frameshifts'] else '')
                r['note'] += f"; NCBI model {r['annot_exons']} exon(s), {r['annot_aa']} aa" if r['annot_gene'] else '; no NCBI model overlaps'
            ids = {G['id'] for _, G in ov}
            left = [G for G in gs if G['end'] < r['start'] and G['id'] not in ids][-10:][::-1]
            right = [G for G in gs if G['start'] > r['end'] and G['id'] not in ids][:10]
            for sgn, lst in ((-1, left), (1, right)):
                for i, G in enumerate(lst, 1):
                    dist = r['start'] - G['end'] if sgn < 0 else G['start'] - r['end']
                    fn.write('\t'.join(map(str, [gca, r['locus_id'], sgn * i, G['name'], G['prot'], G['prod'], G['ex'], dist, int(bool(MOB.search(G['prod'] + ' ' + G['btype'])))])) + '\n')
    for r in rows: fo.write('\t'.join(str(r[h]) for h in H) + '\n')
    best_st = [r['status'] for r in sel if r['cls'] == 'chs_like']
    st = 'absent' if not best_st else ('intact' if 'intact' in best_st else 'truncated' if 'truncated' in best_st else 'pseudogene')
    gen.append(dict(gca=gca, organism=org, species=sp, complex=cx, level=gl[gca]['level'], n_contigs=len(ctg), assembly_len=tot, assembly_gc=round(agc, 2), n_chs_like=len(best_st), n_borderline=sum(r['cls'] == 'borderline' for r in sel),
                    n_ags_like=sum(r['cls'] == 'ags_like_or_other' for r in rows), genome_status=st if not (st == 'absent' and any(r['cls'] == 'borderline' for r in sel)) else 'absent(borderline_only)'))
fo.close(); fn.close()
with open(f'{S}/genomes.tsv', 'w') as f:
    ks = list(gen[0]); f.write('\t'.join(ks) + '\n')
    for g in gen: f.write('\t'.join(str(g[k]) for k in ks) + '\n')
cnt = collections.defaultdict(collections.Counter)
for g in gen: cnt[g['species']]['total'] += 1; cnt[g['species']][g['genome_status']] += 1
with open(f'{S}/species_summary.tsv', 'w') as f:
    f.write('species\tcomplex\tgenomes\tintact\ttruncated\tpseudogene\tabsent\tabsent(borderline_only)\n')
    cxm = {g['species']: g['complex'] for g in gen}
    for s in sorted(cnt, key=lambda s: -cnt[s]['total']):
        c = cnt[s]; f.write(f"{s}\t{cxm[s]}\t{c['total']}\t{c['intact']}\t{c['truncated']}\t{c['pseudogene']}\t{c['absent']}\t{c['absent(borderline_only)']}\n")
print('genomes parsed', len(gen), 'missing mp', nmiss)
