#!/usr/bin/env python3
"""Assemble the clean manuscript (manuscript_v2.md) from the draft_v2_* files.

Usage (login node, any env with Python 3.8+):
    python build_manuscript_v2.py            # keeps [CHECK: …]/[PENDING: …]/[REF] marks visible
    python build_manuscript_v2.py --clean    # removes those marks too (submission copy)
    python build_manuscript_v2.py -o other.md

What it does
  * takes each draft file, drops its service header (everything up to the first '---' line)
    and its trailing '## Перенос: что куда ушло' section;
  * inserts top-level headings where the draft files lack them (2. Materials and Methods,
    3. Results) and demotes their '## 2.x' / '## 3.x' subsections to '###';
  * inserts Table 1 and Table 2 as markdown tables from paper/tables/*.tsv right after
    their caption paragraphs (Tables 3 and 4 are already markdown in the results draft);
  * removes service comments written as italic brackets, e.g. *[colleague, edited; R2:31]*,
    *[new: …]*, *[v1 new]*, *[edited: …]*;
  * keeps the open marks **[check…]**, **[pending…]**, **[ref]**, **[REF]**, [user…] — they are
    rewritten as [CHECK: …], [PENDING: …], [REF], [USER: …] so that a reader sees them — and
    writes all of them to open_marks.tsv (section, source file, line, type, text, context);
    with --clean the marks (and a Russian explanation that directly follows a mark) are removed;
  * appends References (references_v2.md, blocks I–III as they are; renumbering is a separate step)
    and a Supplementary Materials pointer;
  * checks that no '*[' remains, reports Cyrillic lines left in the output, and prints
    word counts per section.
"""
import argparse
import csv
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parent
TABLES = PAPER / "tables"

TITLE = "Opine metabolism genes of rhizobial origin in fungal genomes"  # decided 2026-10-12 (docs/decisions.md)
# The colleague's docx has no author block; keep a visible mark until it is provided.
AUTHORS = ("[TO FILL: First Author 1,*, Second Author 2 and Third Author 1]  \n"
           "[TO FILL: 1 Department, Institute, City, Country; 2 Department, Institute, City, Country]  \n"
           "[TO FILL: * Correspondence: email@institution.org]")

# (file, top heading to insert or None, demote '## N.x' subsections to '###')
SECTIONS = [
    ("draft_v2_summary_abstract.md", None, False),
    ("draft_v2_introduction.md", None, False),
    ("draft_v2_methods.md", "## 2. Materials and Methods", True),
    ("draft_v2_results.md", "## 3. Results", True),
    ("draft_v2_discussion_conclusions.md", None, False),
]
REFERENCES = HERE / "references_v2.md"

# Tables inserted after their caption paragraph: caption prefix -> (tsv, columns or None=all)
TABLE_TSV = {
    "**Table 1.**": TABLES / "Table1_funnel.tsv",
    "**Table 2.**": TABLES / "Table2_main.tsv",  # 2026-10-12x: seven lineages, 8 columns; full summary = Table S1c
}

KEEP_MARK = re.compile(
    r"(\*\*|\*)?\[(?P<kind>check|pending|ref|REF|user)(?P<rest>[^\]]*)\](?:\*\*|\*)?"
)
ITALIC_COMMENT = re.compile(r"\s?\*\[[^\]\n]*\]\*")
CYR = re.compile(r"[А-Яа-яЁё]")


def body_of(path: Path) -> list[str]:
    """Lines between the first '---' line and the '## Перенос' section."""
    lines = path.read_text(encoding="utf-8").splitlines()
    start = 0
    for i, l in enumerate(lines):
        if l.strip() == "---":
            start = i + 1
            break
    end = len(lines)
    for i, l in enumerate(lines):
        if l.startswith("## Перенос"):
            end = i
            break
    body = lines[start:end]
    while body and not body[0].strip():
        body.pop(0)
    while body and not body[-1].strip():
        body.pop()
    return body


def tsv_to_md(path: Path) -> list[str]:
    with path.open(encoding="utf-8", newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    head = [h.replace("_", " ") for h in rows[0]]
    out = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    for r in rows[1:]:
        r = r + [""] * (len(head) - len(r))
        out.append("| " + " | ".join(c.replace("|", "/") for c in r) + " |")
    return out



SMALL = {"a", "an", "and", "the", "of", "in", "for", "on", "to", "by", "as", "at", "or", "vs.", "with", "from", "among", "into", "against", "between", "under", "over"}
def title_case_heading(line: str) -> str:
    """MDPI headings in Title Case (decision 2026-10-12): capitalise lower-case words except small words
    (never the first word); words with any capital letter (dN/dS, RNA-seq, T-DNA, AU) and italic tokens are left as they are."""
    hashes, text = line.split(" ", 1)
    m = re.match(r"^((?:\d+\.)+\s*)?(.*)$", text)
    num, rest = (m.group(1) or ""), m.group(2)
    words = rest.split(" ")
    out = []
    for i, w in enumerate(words):
        if not w or w.startswith("*") or any(c.isdigit() for c in w):
            out.append(w); continue
        core = w.strip("(),:;")
        if i > 0 and core in SMALL:
            out.append(w); continue
        def cap(part):
            if not part or any(c.isupper() for c in part) or part in SMALL: return part
            k = 0
            while k < len(part) and not part[k].isalpha(): k += 1   # keep leading "(" etc.
            return part[:k] + part[k:k + 1].upper() + part[k + 1:]
        out.append("-".join(cap(part) for part in w.split("-")))
    return f"{hashes} {num}{' '.join(out)}"

def process_section(name: str, lines: list[str], clean: bool, marks: list[dict],
                    src: Path, offset: int) -> list[str]:
    out = []
    section = name
    for i, line in enumerate(lines):
        if line.startswith("## ") or line.startswith("### "):
            line = title_case_heading(line)
            section = line.lstrip("# ").strip()

        # 1. collect + rewrite the open marks
        def repl(m):
            kind = m.group("kind")
            rest = m.group("rest").strip()
            rest = rest[1:].strip() if rest.startswith(":") else rest
            tail = line[m.end():m.end() + 160]
            marks.append({"section": section, "file": src.name, "line": offset + i + 1,
                          "type": kind.lower(), "text": rest, "context": tail.strip()})
            if clean:
                return ""
            label = {"check": "CHECK", "pending": "PENDING", "ref": "REF", "REF": "REF",
                     "user": "USER"}[kind]
            return f"[{label}: {rest}]" if rest else f"[{label}]"

        new = KEEP_MARK.sub(repl, line)
        # 2. drop italic service comments
        new = ITALIC_COMMENT.sub("", new)
        if clean:
            # drop a Russian explanation that directly follows a removed mark: ' — …' up to
            # the sentence end, only when it contains Cyrillic
            new = re.sub(r"\s—\s(?=[^.\n]*[А-Яа-яЁё])[^.\n]*(\.|$)", "", new)
            new = re.sub(r"\s*\((?=[^)]*[А-Яа-яЁё])[^)]*\)", "", new)
        new = re.sub(r"[ \t]+$", "", new)
        new = re.sub(r"  +", " ", new) if not new.startswith("|") else new
        out.append(new)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clean", action="store_true", help="remove CHECK/PENDING/REF marks too")
    ap.add_argument("-o", "--out", default=str(HERE / "manuscript_v2.md"))
    ap.add_argument("--marks", default=str(HERE / "open_marks.tsv"))
    args = ap.parse_args()

    marks: list[dict] = []
    doc: list[str] = []
    stats: list[tuple[str, int]] = []

    doc += [f"# {TITLE}", "", AUTHORS, ""]   # placeholders kept in the clean copy too (marked [TO FILL])

    for fname, top, demote in SECTIONS:
        src = PAPER / fname
        raw = src.read_text(encoding="utf-8").splitlines()
        body = body_of(src)
        offset = next(i for i, l in enumerate(raw) if l.strip() == "---") + 1
        # skip leading blank lines removed in body_of
        while raw[offset].strip() == "":
            offset += 1
        if demote:
            body = [("#" + l) if re.match(r"^## \d+\.\d+", l) else l for l in body]
        if top:
            body = [top, ""] + body
        # insert Tables 1–2 after their caption paragraphs
        with_tables = []
        for l in body:
            with_tables.append(l)
            for cap, tsv in TABLE_TSV.items():
                if l.startswith(cap):
                    with_tables += [""] + tsv_to_md(tsv) + [""]
        body = process_section(fname, with_tables, args.clean, marks, src, offset)
        words = sum(len(re.findall(r"\S+", l)) for l in body if not l.startswith("|"))
        stats.append((fname, words))
        doc += body + ["", ""]

    # Supplementary Materials before References (MDPI order; 2026-10-12x)
    _par = PAPER / "submission" / "supplementary_materials_paragraph.md"   # generated by supplementary/build_supplementary.py
    _supp = _par.read_text(encoding="utf-8").strip() if _par.exists() else "**Supplementary Materials:** [PENDING: paragraph not generated]"
    doc += ["## Supplementary Materials", "", _supp.replace("**Supplementary Materials:** ", ""), "", ""]
    # Back matter (MDPI Biology): placeholders for the authors, marked [TO FILL] (author decision 2026-10-13)
    doc += ["## Author Contributions", "",
            "Conceptualization, [TO FILL: initials]; methodology, [TO FILL: initials]; software, [TO FILL: initials]; validation, [TO FILL: initials]; formal analysis, [TO FILL: initials]; investigation, [TO FILL: initials]; data curation, [TO FILL: initials]; writing—original draft preparation, [TO FILL: initials]; writing—review and editing, [TO FILL: initials]; visualization, [TO FILL: initials]; supervision, [TO FILL: initials]; funding acquisition, [TO FILL: initials]. All authors have read and agreed to the published version of the manuscript.", "",
            "## Funding", "",
            "This research was funded by [TO FILL: funder and grant number], or: This research received no external funding.", "",
            "## Institutional Review Board Statement", "", "Not applicable.", "",
            "## Informed Consent Statement", "", "Not applicable.", "",
            "## Data Availability Statement", "",
            "All genome assemblies and proteomes analysed are public NCBI records; their accessions are listed in Tables S6 and S7 and Data S3. The protein sequences of the loci of bacterial origin (Data S1), the alignments, trees and IQ-TREE reports (Data S2) and all Supplementary Tables are provided with this article. The code of the screening pipeline and of all analyses, together with the small input files needed to regenerate the tables and figures, is available at [TO FILL: GitHub URL https://github.com/pzhurbenko/fungal-opine-genes] under the MIT licence (code) and CC BY 4.0 (data).", "",
            "## Acknowledgments", "",
            "[TO FILL: e.g. We thank the sequencing consortia whose public genome assemblies made this study possible, and the Weizmann Institute HPC team for computing resources.]", "",
            "## Conflicts of Interest", "", "The authors declare no conflicts of interest.", ""]
    # References: single list renumbered by first citation (build_references.py, 2026-10-13)
    from build_references import renumber
    doc, refs_entries, ref_num, ref_unused = renumber(doc)
    refs = ["## References", ""] + refs_entries + [""]
    doc += refs
    print(f"references: {len(refs_entries)} entries; unused keys dropped: {ref_unused}")
    stats.append(("references (renumbered)", sum(len(l.split()) for l in refs)))

    text = "\n".join(doc).rstrip() + "\n"
    text = re.sub(r"\n{3,}", "\n\n", text)
    Path(args.out).write_text(text, encoding="utf-8")

    with open(args.marks, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["section", "file", "line", "type", "text", "context"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(marks)

    # checks
    leftovers = [(i + 1, l[:120]) for i, l in enumerate(text.splitlines()) if "*[" in l]
    cyr = [(i + 1, l[:120]) for i, l in enumerate(text.splitlines())
           if CYR.search(l) and not l.startswith("## References") and i < len(text.splitlines())]
    from collections import Counter
    kinds = Counter(m["type"] for m in marks)

    print(f"written: {args.out}  ({len(text.split())} words, clean={args.clean})")
    print("words per section:")
    for n, w in stats:
        print(f"  {n:45s} {w:6d}")
    print("open marks:", dict(kinds), "->", args.marks)
    print("leftover '*[' lines:", len(leftovers))
    for ln, l in leftovers[:10]:
        print("   ", ln, l)
    # Cyrillic inside the References block is expected (notes of references_v2.md)
    ref_start = next(i for i, l in enumerate(text.splitlines()) if l == "## References")
    cyr_body = [c for c in cyr if c[0] <= ref_start]
    print("Cyrillic lines in body (outside References):", len(cyr_body))
    for ln, l in cyr_body[:15]:
        print("   ", ln, l)
    return 0 if not leftovers else 1


if __name__ == "__main__":
    sys.exit(main())
