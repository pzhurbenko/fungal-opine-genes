"""Minimal Markdown -> DOCX for README_supplementary (headings, paragraphs, pipe tables, **bold**, *italic*, `code`).
usage: md_to_docx.py IN.md OUT.docx"""
import re,sys
from docx import Document
from docx.shared import Cm,Pt
RX=re.compile(r'(\*\*[^*]+\*\*|\*[^*\s][^*]*?\*|`[^`]+`)')
def runs(p,text,size=None):
    for part in RX.split(text):
        if not part:continue
        if part.startswith('**'):r=p.add_run(part[2:-2]);r.bold=True
        elif part.startswith('`'):r=p.add_run(part[1:-1]);r.font.name='Courier New'
        elif part.startswith('*') and len(part)>2:r=p.add_run(part[1:-1]);r.italic=True
        else:r=p.add_run(part)
        if size:r.font.size=Pt(size)
def cells(l):return [c.strip() for c in l.strip().strip('|').split('|')]
def main(src,dst):
    d=Document();s=d.sections[0];s.page_width=Cm(21);s.page_height=Cm(29.7)
    for m in('left_margin','right_margin','top_margin','bottom_margin'):setattr(s,m,Cm(2.5))
    st=d.styles['Normal'];st.font.name='Palatino Linotype';st.font.size=Pt(10)
    L=open(src).read().split('\n');i=0
    while i<len(L):
        l=L[i]
        if not l.strip():i+=1;continue
        m=re.match(r'(#+)\s+(.*)',l)
        if m:
            h=d.add_heading(level=min(len(m[1]),4));runs(h,m[2]);i+=1;continue
        if l.startswith('|'):
            rows=[]
            while i<len(L) and L[i].startswith('|'):
                if not re.fullmatch(r'\|[\s:|-]+\|?',L[i].strip()):rows.append(cells(L[i]))
                i+=1
            n=max(len(r) for r in rows);t=d.add_table(rows=len(rows),cols=n);t.style='Table Grid'
            for a,r in enumerate(rows):
                for b,c in enumerate(r):
                    p=t.cell(a,b).paragraphs[0];runs(p,c,9)
                    if a==0:
                        for x in p.runs:x.bold=True
            d.add_paragraph();continue
        m=re.match(r'[-*]\s+(.*)',l)
        if m:runs(d.add_paragraph(style='List Bullet'),m[1]);i+=1;continue
        buf=[l];i+=1
        while i<len(L) and L[i].strip() and not L[i].startswith(('|','#')):buf.append(L[i]);i+=1
        runs(d.add_paragraph(),' '.join(buf))
    d.save(dst)
if __name__=='__main__':main(*sys.argv[1:3])
