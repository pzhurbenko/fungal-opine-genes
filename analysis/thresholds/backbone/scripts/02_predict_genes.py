#!/usr/bin/env python3
"""Predict protein-coding genes on every test contig with pyrodigal (Prodigal algorithm, metagenomic mode,
same procedure for both samples). Fungal contigs contain introns, so prodigal yields gene fragments there;
genes/kb is therefore a prokaryote-style coding-density measure, comparable between samples.
Outputs (work/): proteins.faa (ids '<sample>__<contig>__g<n>'), genes.tsv (contig_key, n_genes, coding_bp)
Requires PYTHONPATH=work/pylib (pyrodigal 3.7.1 installed there with pip --target).
"""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import sys
import pyrodigal

W = f'{GMO_ROOT}/analysis/thresholds/backbone/work'


def read_fasta(path):
    name, seq = None, []
    for l in open(path):
        if l.startswith('>'):
            if name: yield name, ''.join(seq)
            name, seq = l[1:].split()[0], []
        else: seq.append(l.strip())
    if name: yield name, ''.join(seq)


gf = pyrodigal.GeneFinder(meta=True)
with open(f'{W}/proteins.faa', 'w') as fp, open(f'{W}/genes.tsv', 'w') as ft:
    ft.write('contig_key\tn_genes\tcoding_bp\n')
    for name, s in read_fasta(f'{W}/contigs.fna'):
        genes = gf.find_genes(s.encode())
        cbp = 0
        for i, g in enumerate(genes, 1):
            fp.write(f'>{name}__g{i}\n{g.translate(include_stop=False)}\n'); cbp += g.end - g.begin + 1
        ft.write(f'{name}\t{len(genes)}\t{cbp}\n')
print('done', file=sys.stderr)
