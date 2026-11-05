from pathlib import Path
import csv,json,hashlib,collections,gzip
from Bio import SeqIO
HERE=Path(__file__).resolve().parent
GMO=HERE.parents[2]
P=GMO/'analysis/refset/plasmids'
rows=list(csv.DictReader(open(P/'cds_table.tsv'),delimiter='\t'))
selected=[r for r in rows if r['in_TDNA']=='yes']
legacy=list(SeqIO.parse(P/'tdna_proteins_all.faa','fasta'))
assert len(legacy)==len(selected)==302
source={}
for acc in {r['plasmid'] for r in selected}:
    for rec in SeqIO.parse(gzip.open(HERE/'source_records'/(acc+'.gb.gz'),'rt'),'genbank'):
        for f in rec.features:
            if f.type=='CDS' and 'translation' in f.qualifiers:
                source[(acc,f.qualifiers.get('protein_id',[''])[0],str(int(f.location.start)+1),str(int(f.location.end)))]=f.qualifiers['translation'][0]
mapping=[];corrected=[]
for n,(r,s) in enumerate(zip(selected,legacy),1):
    assert r['protein_id']==s.id.split('|')[0]
    aa=source[(r['plasmid'],r['protein_id'],r['start'],r['end'])]
    assert aa==str(s.seq),(n,r)
    mapping.append(dict(input_row=n,source_accession=r['plasmid'],protein_accession=r['protein_id'],
        start_1based=r['start'],end_1based=r['end'],strand=r['strand'],gene=r['gene'],product=r['product'],
        legacy_header=s.description,legacy_accession=s.id.split('|')[1],sequence_sha256=hashlib.sha256(aa.encode()).hexdigest(),
        matches_source_translation=True,header_accession_correct=r['plasmid']==s.id.split('|')[1]))
    corrected.append(f'>{r["protein_id"]}|{r["plasmid"]}|{r["start"]}-{r["end"]}|row{n}\n{aa}\n')
def write(path,records):
    with open(path,'w') as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter='\t');w.writeheader();w.writerows(records)
write(HERE/'Table_S19b_protein_provenance.tsv',mapping)
(HERE/'tdna_proteins_provenance_corrected.faa').write_text(''.join(corrected))
meta=list(csv.DictReader(open(HERE/'Table_S19a_plasmid_sources.tsv'),delimiter='\t'))
counts=collections.Counter(r['plasmid'] for r in selected)
catalog={r['plasmid']:r for r in csv.DictReader(open(P/'plasmids.tsv'),delimiter='\t')}
for r in meta:
    x=catalog[r['accession']]
    r.update(n_border_hits=x['n_borders'],n_candidate_intervals=x['n_TDNA_intervals'],
             candidate_intervals_0based=x['TDNA_intervals'],n_selected_proteins=counts[r['accession']])
write(HERE/'Table_S19a_plasmid_sources.tsv',meta)
result=dict(records=len(meta),all_cds=len(rows),records_with_intervals=sum(int(r['n_TDNA_intervals'])>0 for r in catalog.values()),
            records_contributing_proteins=len(counts),selected_proteins=len(selected),
            legacy_distinct_accessions=len({r['legacy_accession'] for r in mapping}),
            mislabelled_headers=sum(not r['header_accession_correct'] for r in mapping),
            all_amino_acid_sequences_verified=True)
(HERE/'plasmid_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(result)
