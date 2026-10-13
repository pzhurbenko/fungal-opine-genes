#!/usr/bin/env python
"""
build_paper_html.py — version 2 (2026-10-12y).

Builds a single self-contained HTML page of the manuscript for reading and commenting:
  source text   : paper/manuscript_v2/manuscript_v2.md  (assembled by manuscript_v2/build_manuscript_v2.py;
                  open marks are visible as [CHECK: …] / [REF] / [PENDING: …] and are highlighted)
  main figures  : paper/figures/Fig{1..5}_*.png, embedded as base64 (downscaled to MAXW px) before each caption
  supplement    : paper/supplementary/README_supplementary.md (descriptions of Tables S1–S17, captions of
                  Figures S1–S11 and Data S1–S3) with the Figure S PNGs embedded and each TSV shown as a table
                  (first MAX_TSV_ROWS rows for the long ones)
  output        : paper/manuscript_v2/manuscript_v2.html

Run in env PlantCad2 (needs the `markdown` package; Pillow optional, for downscaling):
    python build_paper_html.py [-i manuscript_v2.md] [-o out.html] [--no-supp-tables]

Version 1 of this script (draft_v1 sources, service notes, docx parsing) is kept as build_paper_html.py.bak_20261012y.
"""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import argparse, base64, csv, datetime, html, io, re
from html.parser import HTMLParser
from pathlib import Path

import markdown

try:
    from PIL import Image
except ImportError:  # pragma: no cover
    Image = None

PAPER = Path(f'{GMO_ROOT}/paper')
MS = PAPER / 'manuscript_v2'
SUPP = PAPER / 'supplementary'
SUPP_FIG = SUPP / 'figures'
MAXW = 1600          # main figures
MAXW_SUPP = 1100     # supplementary figures (page size; the vector files are in supplementary/figures)
MAX_TSV_ROWS = 25    # rows shown for long supplementary tables

FIGURES = {
    1: PAPER / 'figures/Fig1_workflow.png',
    2: PAPER / 'figures/Fig2_screen.png',
    3: PAPER / 'figures/Fig3_trees.png',   # AU panel moved to Supplementary Figure S8 (2026-10-12q)
    4: PAPER / 'figures/Fig4_islands.png',
    5: PAPER / 'figures/Fig5_boletales_genera.png',   # former Figure 2 (2026-10-12aa); dN/dS is now Figure S11 (2026-10-13m)
}

_warn, _stats = [], {'imgs': 0, 'resized': 0, 'tables': 0, 'supp_tables': 0}


def esc(s):
    return html.escape(s, quote=False)


# ---------------- images ----------------
def img_data_uri(path, maxw, quantize=False):
    """PNG as data URI; downscaled to maxw px and, for the supplement, reduced to a 256-colour palette
    (about 3x smaller; the line-art trees lose nothing visible — the vector files are in supplementary/figures)."""
    if not path.exists():
        _warn.append(f'missing figure file: {path}')
        return None
    raw = path.read_bytes()
    if Image is not None:
        im = Image.open(path)
        changed = False
        if im.width > maxw:
            h = round(im.height * maxw / im.width)
            im = im.convert('RGB').resize((maxw, h), Image.LANCZOS)
            changed = True
        if quantize:
            im = im.convert('RGB').quantize(256, method=Image.MEDIANCUT)
            changed = True
        if changed:
            buf = io.BytesIO()
            im.save(buf, 'PNG', optimize=True)
            raw = buf.getvalue()
            _stats['resized'] += 1
    _stats['imgs'] += 1
    return 'data:image/png;base64,' + base64.b64encode(raw).decode()


def img_tag(path, alt, maxw=MAXW, quantize=False):
    uri = img_data_uri(path, maxw, quantize)
    if uri is None:
        return f'<div class="figmissing">Missing file: {esc(str(path))}</div>'
    return f'<img src="{uri}" alt="{esc(alt)}" class="zoomable" loading="lazy">'


# ---------------- marks ----------------
MARK_RE = re.compile(r'\[(CHECK|REF|PENDING|TO FILL|VERIFY WITH AUTHOR)((?::[^\]]*)?)\]')


def mark_sub(m):
    kind = m.group(1)
    cls = {'CHECK': 'f-chk', 'REF': 'f-ref', 'PENDING': 'f-todo', 'TO FILL': 'f-todo', 'VERIFY WITH AUTHOR': 'f-chk'}[kind]
    return f'<mark class="flag {cls}">[{esc(kind + m.group(2))}]</mark>'


def highlight_marks(h):
    return MARK_RE.sub(mark_sub, h)


# ---------------- markdown -> html ----------------
def md_to_html(text):
    return markdown.markdown(text, extensions=['tables', 'sane_lists'])


def slug(s, used):
    base = re.sub(r'[^a-z0-9]+', '-', s.lower()).strip('-')[:60] or 'sec'
    sid, i = base, 2
    while sid in used:
        sid, i = f'{base}-{i}', i + 1
    used.add(sid)
    return sid


HEAD_RE = re.compile(r'<h([23])>(.*?)</h\1>', re.S)
FIGCAP_RE = re.compile(r'<p><strong>Figure (\d+)\.</strong>')


def add_heading_ids(h, toc, used, prefix=''):
    def rep(m):
        lvl, inner = int(m.group(1)), m.group(2)
        plain = re.sub(r'<[^>]+>', '', inner)
        sid = slug(prefix + plain, used)
        toc.append((lvl, sid, plain))
        return f'<h{lvl} id="{sid}">{inner}</h{lvl}>'
    return HEAD_RE.sub(rep, h)


def insert_figures(h):
    def rep(m):
        n = int(m.group(1))
        p = FIGURES.get(n)
        img = img_tag(p, f'Figure {n}') if p else f'<div class="figmissing">No file mapped for Figure {n}</div>'
        return f'<figure id="fig{n}"><div class="fpanels"><div class="fpanel">{img}</div></div></figure>' + m.group(0)
    return FIGCAP_RE.sub(rep, h)


def wrap_tables(h):
    n = len(re.findall(r'<table>', h))
    _stats['tables'] += n
    return h.replace('<table>', '<div class="tablewrap"><table>').replace('</table>', '</table></div>')


# ---------------- supplementary ----------------
TSV_FILE_RE = re.compile(r'<p>File: <code>([^<]+\.tsv)</code>')
SUPP_FIG_HEAD_RE = re.compile(r'<h3 id="([^"]+)">(Figure S(\d+)\.[^<]*)</h3>')


def tsv_table_html(path):
    if not path.exists():
        _warn.append(f'missing table file: {path}')
        return f'<div class="figmissing">Missing file: {esc(path.name)}</div>'
    with path.open(encoding='utf-8', newline='') as fh:
        rows = list(csv.reader(fh, delimiter='\t'))
    if not rows:
        return ''
    head, body = rows[0], rows[1:]
    shown = body[:MAX_TSV_ROWS]
    out = ['<details class="tsv"><summary>Show table (', esc(path.name), f'; {len(body)} rows × {len(head)} columns',
           f'; first {len(shown)} rows shown' if len(body) > len(shown) else '', ')</summary><div class="tablewrap"><table><thead><tr>']
    out += [f'<th>{esc(c)}</th>' for c in head]
    out.append('</tr></thead><tbody>')
    for r in shown:
        out.append('<tr>' + ''.join(f'<td>{esc(c)}</td>' for c in r) + '</tr>')
    out.append('</tbody></table></div></details>')
    _stats['supp_tables'] += 1
    return ''.join(out)


def supplement_html(toc, used, with_tables=True):
    readme = SUPP / 'README_supplementary.md'
    if not readme.exists():
        _warn.append('README_supplementary.md not found')
        return '<p>Supplementary materials not found.</p>'
    text = readme.read_text(encoding='utf-8')
    text = re.sub(r'^# .*\n', '', text, count=1)  # drop the README title (we add our own h2)
    # demote: README '## Table S..' -> h3, '### Figure S..' -> h3 too (both are items of the supplement)
    text = re.sub(r'^## ', '### ', text, flags=re.M)
    h = md_to_html(text)
    h = add_heading_ids(h, toc, used, prefix='supp-')
    if with_tables:
        h = TSV_FILE_RE.sub(lambda m: m.group(0) + '</p>' + tsv_table_html(SUPP / m.group(1)) + '<p>', h)
        h = h.replace('<p></p>', '')

    def figrep(m):
        n = int(m.group(3))
        cands = sorted(SUPP_FIG.glob(f'FigureS{n}_*.png')) + sorted(SUPP.glob(f'FigureS{n}_*.png'))
        img = img_tag(cands[0], m.group(2), MAXW_SUPP, quantize=True) if cands else f'<div class="figmissing">No PNG for Figure S{n}</div>'
        return m.group(0) + f'<figure id="figS{n}"><div class="fpanels"><div class="fpanel">{img}</div></div></figure>'
    h = SUPP_FIG_HEAD_RE.sub(figrep, h)
    h = wrap_tables(h)
    files = [('Supplementary_Tables.xlsx', SUPP / 'Supplementary_Tables.xlsx'),
             ('Supplementary_Materials.docx', SUPP / 'Supplementary_Materials.docx'),
             ('README_supplementary.md', readme),
             ('Data_S1_proteins_of_loci_of_bacterial_origin.faa', SUPP / 'Data_S1_proteins_of_loci_of_bacterial_origin.faa'),
             ('Data_S3_screened_assemblies.tsv', SUPP / 'Data_S3_screened_assemblies.tsv'),
             ('Supplementary_Data_trees.tar.gz', SUPP / 'data/Supplementary_Data_trees.tar.gz')]
    lst = ''.join(f'<li><code>{esc(n)}</code>{"" if p.exists() else " (missing)"}</li>' for n, p in files)
    return (f'<div class="svc">Supplementary files: <ul class="svc">{lst}</ul></div>' + h)  # 2026-10-13: no local path; div (ul inside p is invalid)


# ---------------- validation ----------------
class _Checker(HTMLParser):
    def __init__(self):
        super().__init__()
        self.counts = {}
        self.imgs = 0
        self.bad_img = 0
        self.stack = []
        self.errors = []

    def handle_starttag(self, tag, attrs):
        self.counts[tag] = self.counts.get(tag, 0) + 1
        if tag == 'img':
            self.imgs += 1
            src = dict(attrs).get('src', '')
            if not src.startswith('data:image/png;base64,'):
                self.bad_img += 1
        if tag not in ('img', 'meta', 'br', 'hr', 'input', 'link'):
            self.stack.append(tag)

    def handle_endtag(self, tag):
        if tag in ('img', 'meta', 'br', 'hr', 'input', 'link'):
            return
        if self.stack and self.stack[-1] == tag:
            self.stack.pop()
        elif tag in self.stack:
            self.errors.append(f'misnested </{tag}>')
            while self.stack and self.stack[-1] != tag:
                self.stack.pop()
            self.stack.pop()
        else:
            self.errors.append(f'stray </{tag}>')


def validate(page):
    c = _Checker()
    c.feed(page)
    c.close()
    return c


def build_toc(items):
    out, open_sub = ['<ul>'], False
    for lvl, sid, text in items:
        t = esc(text if len(text) < 70 else text[:67] + '...')
        if lvl == 2:
            if open_sub:
                out.append('</ul></li>')
                open_sub = False
            elif out[-1] != '<ul>':
                out.append('</li>')
            out.append(f'<li class="l2"><a href="#{sid}">{t}</a>')
        else:
            if not open_sub:
                out.append('<ul>')
                open_sub = True
            out.append(f'<li class="l3"><a href="#{sid}">{t}</a></li>')
    if open_sub:
        out.append('</ul>')
    out.append('</li></ul>')
    return ''.join(out)


# ---------------- main ----------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-i', '--inp', default=str(MS / 'manuscript_v2.md'))
    ap.add_argument('-o', '--out', default=str(MS / 'manuscript_v2.html'))
    ap.add_argument('--no-supp-tables', action='store_true', help='do not embed the TSV tables of the supplement')
    a = ap.parse_args()
    src, out = Path(a.inp), Path(a.out)

    text = src.read_text(encoding='utf-8')
    m = re.match(r'# (.+)\n', text)
    title = m.group(1).strip() if m else 'Manuscript'
    if m:
        text = text[m.end():]
    n_marks = {k: len(re.findall(r'\[' + k + r'(?::[^\]]*)?\]', text)) for k in ('CHECK', 'REF', 'PENDING', 'TO FILL', 'VERIFY WITH AUTHOR')}

    toc, used = [], set()
    body = md_to_html(text)
    body = add_heading_ids(body, toc, used)
    body = insert_figures(body)
    body = wrap_tables(body)
    body = highlight_marks(body)

    supp_id = slug('Supplementary tables, figures and data', used)
    toc.append((2, supp_id, 'Supplementary tables, figures and data'))
    supp = supplement_html(toc, used, with_tables=not a.no_supp_tables)
    supp = highlight_marks(supp)

    today = datetime.date.today().isoformat()
    sub = (f'Manuscript v2 — assembled {today} from manuscript_v2.md — target journal: Biology (MDPI). '
           f'Open marks: {n_marks["CHECK"]} [CHECK], {n_marks["REF"]} [REF], {n_marks["PENDING"]} [PENDING], {n_marks["TO FILL"]} [TO FILL], {n_marks["VERIFY WITH AUTHOR"]} [VERIFY WITH AUTHOR] '
           f'(hidden by “clean reading”).')
    page = PAGE.format(title=esc(title), sub=esc(sub), toc=build_toc(toc), body=body,
                       supp_id=supp_id, supp=supp, css=CSS, js=JS)
    out.write_text(page, encoding='utf-8')

    chk = validate(page)
    missing = [w for w in _warn if w.startswith('missing')]
    print(f'{out}  {out.stat().st_size/1e6:.2f} MB; h2={chk.counts.get("h2",0)} h3={chk.counts.get("h3",0)}; '
          f'tables={chk.counts.get("table",0)} (main {_stats["tables"]-_stats["supp_tables"]}, supp tsv {_stats["supp_tables"]}); '
          f'images={chk.imgs} (embedded {_stats["imgs"]}, resized {_stats["resized"]}, non-embedded {chk.bad_img}); '
          f'parser errors={len(chk.errors)}; unclosed={chk.stack[:5]}; PIL={"yes" if Image else "NO"}')
    for w in _warn:
        print('WARNING:', w)
    for e in chk.errors[:10]:
        print('PARSER:', e)
    return 1 if (missing or chk.errors) else 0


CSS = r"""
*{box-sizing:border-box}
body{margin:0;font-family:Georgia,'Times New Roman',serif;line-height:1.55;color:#222;background:#fafaf8}
#toc{position:fixed;left:0;top:0;bottom:0;width:17em;overflow:auto;padding:1em .8em;background:#f0efea;border-right:1px solid #d8d6cc;font:13px/1.35 system-ui,sans-serif}
#toc ul{list-style:none;margin:0;padding:0}#toc ul ul{padding-left:1em}
#toc a{color:#234;text-decoration:none;display:block;padding:2px 0}#toc a:hover{background:#e2e0d6}
#toc .l2>a{font-weight:600;margin-top:.4em}
main{margin-left:18em;padding:1em 2em 4em;max-width:calc(50em + 4em)}
@media(max-width:1000px){#toc{position:static;width:auto;border:0;border-bottom:1px solid #ccc}main{margin-left:0;padding:1em}}
h1{font-size:1.7em;line-height:1.25;margin:.4em 0 .1em}.sub{font:14px system-ui,sans-serif;color:#555;margin-bottom:1em}
h2{font-size:1.4em;margin:1.8em 0 .5em;border-bottom:1px solid #ccc;padding-bottom:.15em}
h3{font-size:1.12em;margin:1.4em 0 .4em}
p{margin:.6em 0}
.flag{padding:0 .2em;border-radius:3px;font-style:normal}
.f-chk{background:#fff2a8;outline:1px solid #e0c200}.f-ref{background:#ffd9c9;outline:1px solid #d9774f}.f-todo{background:#ffc9c9;outline:1px solid #c33}
.tablewrap{overflow-x:auto;margin:.5em 0}
table{border-collapse:collapse;font:13px/1.35 system-ui,sans-serif;min-width:100%}
th,td{border:1px solid #ccc;padding:3px 6px;vertical-align:top}th{background:#ece9df;position:sticky;top:0}
figure{margin:1.5em 0 .3em;padding:.5em;background:#fff;border:1px solid #ddd}
.fpanels{display:flex;flex-direction:column;gap:.6em}.fpanel{position:relative}
figure img{width:100%;height:auto;cursor:zoom-in;display:block}
.figmissing{color:#c00;border:1px dashed #c00;padding:1em}
details.tsv{margin:.4em 0 1em;background:#f4f1e6;border:1px solid #d9d3b8;padding:.3em .8em;font:13px/1.4 system-ui,sans-serif}
details.tsv summary{cursor:pointer;font-weight:600}
.svc{font:13px system-ui,sans-serif;color:#555}
#top{font:13px system-ui,sans-serif;background:#fff;border:1px solid #d8d6cc;padding:.5em .8em;margin:.8em 0;position:sticky;top:0;z-index:5}
#top label{float:right;cursor:pointer;font-weight:600}
body.clean .flag{background:none;outline:0;font-size:0}
#lb{display:none;position:fixed;inset:0;background:#000d;z-index:99;overflow:auto;cursor:zoom-out;text-align:center}
#lb img{max-width:none;width:auto;max-height:none;margin:1em auto;background:#fff}
@media print{#toc,#top,#lb{display:none!important}main{margin:0;max-width:none;padding:0}.flag{background:none;outline:0}
figure,.tablewrap{break-inside:avoid}}
"""

JS = r"""
(function(){
var cb=document.getElementById('clean');
cb.addEventListener('change',function(){document.body.classList.toggle('clean',cb.checked);});
var lb=document.getElementById('lb'),li=lb.querySelector('img');
document.querySelectorAll('img.zoomable').forEach(function(im){im.addEventListener('click',function(){li.src=im.src;lb.style.display='block';});});
lb.addEventListener('click',function(){lb.style.display='none';});
document.addEventListener('keydown',function(e){if(e.key==='Escape')lb.style.display='none';});
})();
"""

PAGE = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{title}</title><style>{css}</style></head><body>
<nav id="toc"><b>Contents</b>{toc}</nav>
<main>
<h1>{title}</h1><div class="sub">{sub}</div>
<div id="top"><label><input type="checkbox" id="clean"> clean reading (hide marks)</label>
<span class="legend"><b>Marks for the co-author:</b> <mark class="flag f-chk">[CHECK]</mark> verify
<mark class="flag f-ref">[REF]</mark> reference needed <mark class="flag f-todo">[PENDING]</mark> not done <mark class="flag f-todo">[TO FILL]</mark> authors to fill in <mark class="flag f-chk">[VERIFY WITH AUTHOR]</mark> reference to verify</span></div>
{body}
<h2 id="{supp_id}">Supplementary tables, figures and data</h2>
{supp}
</main>
<div id="lb"><img alt=""></div>
<script>{js}</script></body></html>
"""

if __name__ == '__main__':
    raise SystemExit(main())
