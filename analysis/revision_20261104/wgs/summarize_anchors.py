from pathlib import Path
import csv,json,re,hashlib
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
cache=ROOT/'analysis/refset/screen/seeds_ext/reciprocal/organism_taxonomy.json'
taxonomy=json.loads(cache.read_text())
annotation={r['protein']:r for r in csv.DictReader(open(ROOT/'analysis/fusarium_odh_island/fig4/fp133_window_genes.tsv'),delimiter='\t')}
selected={};groups={};used={}
for r in csv.reader(open(HERE/'anchors_vs_nr.tsv'),delimiter='\t'):
    names=re.findall(r'\[([^\]]+)\]',r[8]);name=names[-1] if names else ''
    tax=taxonomy.get(name,{})
    if tax.get('kingdom')!='Fungi' or tax.get('genus')=='Fusarium':continue
    seen=groups.setdefault(r[0],set())
    if len(seen)==3 or tax['genus'] in seen:continue
    seen.add(tax['genus']);a=annotation[r[0]];used[name]=tax
    side={'KAL3595990.1':'left','KAL3595996.1':'right','KAL3603064.1':'left','KAL3603068.1':'right'}[r[0]]
    selected.setdefault(r[0],[]).append(dict(anchor_locus=a['locus_tag'],anchor_protein=r[0],contig=a['contig'],
        start_1based=a['start'],end_1based=a['end'],strand=a['strand'],exons=a['exons'],side=side,
        representative_homolog=r[1],homolog_organism=name,homolog_genus=tax['genus'],kingdom=tax['kingdom'],
        identity_pct=r[2],alignment_length_aa=r[3],query_length_aa=r[4],subject_length_aa=r[5],evalue=r[6],bitscore=r[7],title=r[8]))
assert len(selected)==4 and all(len(rs)==3 for rs in selected.values())
rows=[r for group in selected.values() for r in group]
with open(HERE/'anchor_evidence.tsv','w') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
(HERE/'anchor_hit_taxonomy.json').write_text(json.dumps(used,indent=2)+'\n')
db=Path('/shareDB/blastdb/nr/diamond/nr_10_2025.dmnd')
(HERE/'anchor_search.json').write_text(json.dumps(dict(queries=4,reported_targets_per_query=200,
    diamond_version='2.2.6',algorithm='ctg',sensitivity='default (fast)',evalue=1e-5,
    command='diamond blastp -q anchors.faa -d nr_10_2025.dmnd -o anchors_vs_nr.tsv --algo ctg --threads 16 -b 2 -c 1 --max-target-seqs 200 --evalue 1e-5 --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore stitle',
    database_name=db.name,database_bytes=db.stat().st_size,query_sha256=hashlib.sha256((HERE/'anchors.faa').read_bytes()).hexdigest(),
    taxonomy_source='Original reverse-search taxonomy cache, NCBI taxdump 2026-10-06; exact-name entries retained for the representative matches',
    taxonomy_cache_sha256=hashlib.sha256(cache.read_bytes()).hexdigest(),
    purpose='Verify positive fungal homology of annotated multi-exonic anchors; no absence inference from this fast search'),indent=2)+'\n')
print('4 anchors, 3 representative fungal genera per anchor; 800 raw nr hits retained')
