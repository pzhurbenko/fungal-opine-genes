import re,sys,collections
def parse(s):
    s=s.strip().rstrip(';'); i=0
    def node():
        nonlocal i
        if s[i]=='(':
            i+=1; ch=[]
            while True:
                ch.append(node())
                if s[i]==',': i+=1; continue
                if s[i]==')': i+=1; break
            m=re.compile(r'[^,():;]*').match(s,i); i=m.end()
            m=re.compile(r':[^,():;]*').match(s,i)
            if m: i=m.end()
            return ch
        m=re.compile(r'[^,():;]+').match(s,i); name=m.group(0); i=m.end()
        m=re.compile(r':[^,():;]*').match(s,i)
        if m: i=m.end()
        return name
    return node()
def leaves(n): return [n] if isinstance(n,str) else [l for c in n for l in leaves(c)]
def cls(l):
    l=re.sub(r':.*','',l)
    return re.split(r'__|_',l)[0] if '__' in l or True else l
def groups(path):
    t=parse(open(path).read())
    # constrained clades = non-trivial internal children (nested)
    out=[]
    def walk(n,top):
        if isinstance(n,str): return
        if not top and len(leaves(n))>1: out.append(leaves(n))
        for c in n: walk(c,False)
    walk(t,True)
    return leaves(t),out
if __name__=="__main__":
    for p in sys.argv[1:]:
        al,g=groups(p)
        print(p.split('/')[-1],len(al),[len(x) for x in g])
