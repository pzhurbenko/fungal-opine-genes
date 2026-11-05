"""Revision-table overlay used by the repository supplementary builder."""
from pathlib import Path
import csv,io,tarfile,re
ROOT=Path(__file__).resolve().parents[2]
REV=ROOT/'analysis/revision_20261104'

def revised_tables(original):
    """Use audited revision exports while preserving original column definitions."""
    result=list(original);indexes={r[0]:i for i,r in enumerate(result)}
    for p in sorted((REV/'supplement_tables').glob('S*.tsv')):
        data=list(csv.reader(open(p),delimiter='\t'));tid=p.stem
        title=data[0][0].split('. ',1)[1];legend=data[1][0]
        fields=data[2];rows=[]
        for row in data[3:]:
            vals=[]
            for v in row:
                if re.fullmatch(r'-?\d+',v):v=int(v)
                elif re.fullmatch(r'-?\d+\.\d+(?:[Ee][+-]?\d+)?',v):v=float(v)
                vals.append(v)
            rows.append(dict(zip(fields,vals)))
        oldcols=dict(result[indexes[tid]][3]) if tid in indexes else {}
        cols=[(f,oldcols.get(f,'See table legend; '+f.replace('_',' '))) for f in fields]
        record=(tid,title,legend,cols,rows)
        if tid in indexes:result[indexes[tid]]=record
        else:result.append(record)
    return result

def extend_data4(path):
    """Append the two completed M0 sensitivity runs to the original archive."""
    path=Path(path);buf=path.read_bytes()
    with tarfile.open(fileobj=io.BytesIO(buf),mode='r:gz') as old,tarfile.open(path,'w:gz') as new:
        manifest=old.extractfile('MANIFEST.tsv').read().decode()
        for member in old.getmembers():
            if member.name=='MANIFEST.tsv' or member.name.startswith('revision_20261104/'):continue
            new.addfile(member,old.extractfile(member) if member.isfile() else None)
        manifest='\n'.join(line for line in manifest.splitlines() if '\trevision_20261104/' not in line)+'\n'
        for p in sorted((REV/'dnds').rglob('*')):
            if not p.is_file() or '__pycache__' in p.parts:continue
            name='revision_20261104/'+str(p.relative_to(REV/'dnds'))
            new.add(p,arcname=name);manifest+='revision\t'+name+'\tintegration-filtered M0 input or output\n'
        data=manifest.encode();info=tarfile.TarInfo('MANIFEST.tsv');info.size=len(data);new.addfile(info,io.BytesIO(data))
