import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
#!/usr/bin/env python3
# Build <set>/<set>.cds.fna + .pep.faa for the extra dN/dS sets (intact CDS only).
# nos-ODH Boletales is built separately by nos_extract.py (miniprot on locus regions).
# chs / odh-island: single-exon miniprot models -> CDS cut from the assembly by catalogue coords.
# cus-like non-Boletales / pKerr14a acs-like: CDS from loci_pipeline round2 *_models.fna (as prep.py), fallback composition/cds_loci.fna.
import csv, os, glob
from Bio import SeqIO
from Bio.Seq import Seq
G = f'{GMO_ROOT}'; A = f'{G}/analysis'; X = f'{A}/dnds/extra'
T = list(csv.DictReader(open(f'{G}/paper/tables/TableS1_all_lineages.tsv'), delimiter='\t'))
def genome(gca):
    for f in [f'{A}/chs_fusarium/stage2/genomes_files.tsv', f'{A}/fusarium_odh_island/scan/genomes_files.tsv']:
        for l in open(f):
            t = l.rstrip('\n').split('\t')
            if t[0] == gca: return t[1]
    h = glob.glob(f'{G}/data/fungal_genomes/*/ncbi_dataset/data/{gca}/*_genomic.fna')
    return h[0]
def contig(fna, name):
    s = []; on = False
    for l in open(fna):
        if l[0] == '>':
            if on: break
            on = l[1:].split()[0] == name
        elif on: s.append(l.strip())
    return ''.join(s).upper()
def check(i, s):
    if len(s) % 3: return None, f'len%3={len(s)%3}'
    p = str(Seq(s).translate())
    if p.endswith('*'): p = p[:-1]; s = s[:-3]
    if '*' in p: return None, f'internal stop x{p.count("*")}'
    return (s, p), ('ok' if p[0] == 'M' else 'ok_noATG')
def by_coord(rows):
    out = []
    for r in rows:
        c = contig(genome(r['gca']), r['contig']); s = c[int(r['start'])-1:int(r['end'])]
        if r['strand'] == '-': s = str(Seq(s).reverse_complement())
        # if the model lacks the stop codon, add it when it follows directly
        ext = c[int(r['end']):int(r['end'])+3] if r['strand'] == '+' else str(Seq(c[max(0,int(r['start'])-4):int(r['start'])-1]).reverse_complement())
        res, msg = check(r['locus_id'], s)
        if res is None and len(s) % 3 == 0 and False: pass
        out.append((r['locus_id'], res, msg, r['species'], r['status']))
    return out
def from_models(rows):
    fam = {}
    for f in ['ACS', 'CUS']:
        fam.update({x.id.split('|')[0]: str(x.seq).upper() for x in SeqIO.parse(f'{A}/loci_pipeline/round2/proteins/{f}_models.fna', 'fasta')})
    prot = {}
    for f in ['ACS', 'CUS']:
        prot.update({x.id.split('|')[0]: str(x.seq).upper() for x in SeqIO.parse(f'{A}/loci_pipeline/round2/proteins/{f}_models.faa', 'fasta')})
    alt = {x.id.split('|')[0]: str(x.seq).upper() for x in SeqIO.parse(f'{A}/composition/cds_loci.fna', 'fasta')}
    out = []
    for r in rows:
        s = fam.get(r['locus_id']) or alt.get(r['locus_id'])
        res, msg = check(r['locus_id'], s) if s else (None, 'no CDS')
        if s and res is None and r['locus_id'] in prot:  # model CDS with 1-2 extra nt (contig-edge partial): pick the frame that reproduces the model protein
            q = prot[r['locus_id']].rstrip('*')
            for o in range(3):
                for e in range(3):
                    t = s[o:len(s)-e]
                    if len(t) % 3 == 0 and str(Seq(t).translate()).rstrip('*') == q:
                        res, m2 = check(r['locus_id'], t); msg = f'{m2} (trimmed {o}nt 5prime,{e}nt 3prime)'
        out.append((r['locus_id'], res, msg + ('' if r['locus_id'] in fam else ' [cds_loci.fna]'), r['species'], r['status']))
    return out
L = lambda pred: [r for r in T if pred(r) and r['status'] == 'intact']
sets = {
 'chs': by_coord(L(lambda r: r['lineage'].startswith('chs-like'))),
 'odh': by_coord(L(lambda r: r['lineage'] == 'odh island: opine dehydrogenase')),
 'pepP': by_coord(L(lambda r: r['lineage'] == 'odh island: Xaa-Pro aminopeptidase')),
 'g3': by_coord(L(lambda r: r['lineage'] == 'odh island: third gene')),
 'cus_nonBol': from_models(L(lambda r: r['lineage'] == 'cus-like (T-DNA type)' and r['order'] != 'Boletales')),
 'acs_pKerr14a': from_models(L(lambda r: r['lineage'] == 'pKerr14a-type acs-like')),
}
# odh island: Fp133 (GCA_046056315.1) + three F. langsethiae assemblies; the other F. poae assembly (GCA_019609905.1) is kept out
for k in ['odh', 'pepP', 'g3']: sets[k] = [x for x in sets[k] if not x[0].startswith('GCA_019609905.1')]
log = open(f'{X}/prep_log.tsv', 'w'); log.write('set\tlocus_id\tspecies\tcatalogue_status\tcheck\tcds_bp\n')
for k, v in sets.items():
    os.makedirs(f'{X}/{k}', exist_ok=True)
    with open(f'{X}/{k}/{k}.cds.fna', 'w') as fc, open(f'{X}/{k}/{k}.pep.faa', 'w') as fp:
        for i, res, msg, sp, st in v:
            log.write(f'{k}\t{i}\t{sp}\t{st}\t{msg}\t{len(res[0]) if res else ""}\n')
            if res: fc.write(f'>{i}\n{res[0]}\n'); fp.write(f'>{i}\n{res[1]}\n')
    print(k, len(v), 'intact;', sum(1 for x in v if x[1]), 'usable')
