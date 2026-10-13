#!/usr/bin/env python3
"""Renumber citations by order of first appearance and build the single MDPI reference list.

Used by build_manuscript_v2.py (2026-10-13). Sources:
  references_v2.md   block I  : lines "N. entry"            (keys "N")
                     block II : lines "- [Mn] entry (cited …)" (keys "Mn")
  submission/refs_verify_1-34.tsv, refs_verify_35-63_M.tsv : columns key, status, corrected_entry_MDPI
Rules: entries with status OK/FIX use corrected_entry_MDPI (if non-empty); UNVERIFIABLE keep the
original plus " [VERIFY WITH AUTHOR]"; keys without a TSV row (64–70) keep the original text with
italic service notes *[…]* removed. Citations: "[3–5]" for runs of consecutive numbers, "[3,7]" otherwise
(MDPI: no space after the comma). Uncited entries are dropped. Writes reference_map.tsv.
"""
import csv, re
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parent
CIT = re.compile(r"\[((?:M?\d+)(?:\s*(?:,|–|-)\s*M?\d+)*)\]")
NOTE = re.compile(r"\s*\*\[[^\]]*\]\*")

def parse_sources():
    entries = {}
    blk = None
    for line in (HERE / "references_v2.md").read_text(encoding="utf-8").splitlines():
        if line.startswith("## (I)"): blk = "I"; continue
        if line.startswith("## (II)"): blk = "II"; continue
        if line.startswith("## (III)"): blk = "III"; continue
        if blk == "I":
            m = re.match(r"^(\d+)\.\s+(.*)$", line)
            if m: entries[m.group(1)] = NOTE.sub("", m.group(2)).strip()
        elif blk == "II":
            m = re.match(r"^- \[(M\d+)\]\s+(.*)$", line)
            if m:
                t = re.sub(r"\s*\(cited in v2:[^)]*\)", "", m.group(2))
                entries[m.group(1)] = NOTE.sub("", t).strip()
    verified = {}
    for f in ("refs_verify_1-34.tsv", "refs_verify_35-63_M.tsv"):
        with open(PAPER / "submission" / f, encoding="utf-8") as fh:
            rows = [l for l in fh if not l.startswith("#")]
        for r in csv.DictReader(rows, delimiter="\t"):
            k = (r.get("key") or r.get("n") or "").strip()   # first file uses column "n"
            if k: verified[k] = (r.get("status", "").strip(), (r.get("corrected_entry_MDPI") or "").strip())
    return entries, verified

def expand(spec):
    keys = []
    for part in re.split(r"\s*,\s*", spec.strip()):
        m = re.match(r"^(M?)(\d+)\s*(?:–|-)\s*(M?)(\d+)$", part)
        if m:
            a, b = int(m.group(2)), int(m.group(4)); pre = m.group(1)
            keys += [f"{pre}{i}" for i in range(a, b + 1)]
        else:
            keys.append(part.strip())
    return keys

def fmt(nums):
    nums = sorted(set(nums)); out = []; i = 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1: j += 1
        out.append(f"{nums[i]}–{nums[j]}" if j - i >= 2 else (f"{nums[i]},{nums[j]}" if j == i + 1 else str(nums[i])))
        i = j + 1
    return "[" + ",".join(out) + "]"

def mdpi_style(t):
    """*Journal* 2020, 8, 229 -> *Journal* **2020**, *8*, 229 (issue in parentheses dropped); leave entries that already have bold."""
    if "**" in t: return t
    return re.sub(r"(\*[^*]+\*\.?)\s+(\d{4}),\s*(\d+)(?:\(\d+\))?,\s*", r"\1 **\2**, *\3*, ", t, count=1)

def renumber(body_lines, in_refs_section=None):
    """body_lines: list of str (manuscript body, before References). Returns (new_lines, ref_lines, mapping)."""
    entries, verified = parse_sources()
    order = []
    for line in body_lines:
        for m in CIT.finditer(line):
            for k in expand(m.group(1)):
                if k in entries and k not in order: order.append(k)
                elif k not in entries: raise SystemExit(f"citation key without entry: [{k}] in: {line[:100]}")
    num = {k: i + 1 for i, k in enumerate(order)}
    def repl(m):
        ks = expand(m.group(1))
        return fmt([num[k] for k in ks])
    new_lines = [CIT.sub(repl, l) for l in body_lines]
    refs = []
    for k in order:
        st, corr = verified.get(k, ("", ""))
        text = corr if (st.startswith("OK") or st.startswith("FIX")) and corr else entries[k]
        text = NOTE.sub("", text).strip()
        text = re.sub(r"\s*\[[^\]]*[\u0400-\u04FF][^\]]*\]", "", text).strip()   # drop Cyrillic bracket notes
        text = mdpi_style(text)
        if st.startswith("UNVERIFIABLE"): text += " [VERIFY WITH AUTHOR: entry could not be verified in Crossref/PubMed]"
        refs.append(f"{num[k]}. {text}")
    with open(HERE / "reference_map.tsv", "w", encoding="utf-8") as fh:
        fh.write("key\tnumber\tstatus\n")
        for k in order: fh.write(f"{k}\t{num[k]}\t{verified.get(k, ('n/a',''))[0]}\n")
    unused = sorted(set(entries) - set(order), key=lambda x: (x[0] == 'M', int(x.lstrip('M'))))
    return new_lines, refs, num, unused
