import sys
def get(path, ids):
    out={};h=None
    for l in open(path):
        if l.startswith('>'):
            a=l[1:].split()[0]; h=a if a in ids else None
            if h: out[h]=[l.rstrip()]
        elif h: out[h].append(l.strip())
    return {k:''.join(v[1:]) for k,v in out.items()}, 
