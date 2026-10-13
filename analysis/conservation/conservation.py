import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
#!/usr/bin/env python3
# B2a: positions invariant among T-DNA references vs conservation in fungal intact genes.
from Bio import AlignIO
A = f'{GMO_ROOT}/analysis'
SIX = ['JBBPKL010003317', 'JBBPKL010018021', 'JBBPNB010000203', 'JBBPNB010000771', 'JBBPOG010012756', 'JBBPRY010001778']
out = open('conservation_summary.tsv', 'w')
out.write('family\tgroup\tn_refs\tn_fungal\tref_invariant_positions\tconserved_all_fungal\tconserved_ge90_fungal\tmean_identity_invariant\tcontrol_identity_variable\n')


def run(fam, ref_sel, groups):
    aln = AlignIO.read(f'{A}/trees_final/{fam}/{fam}.mafft.aln', 'fasta')
    recs = {r.id: str(r.seq).upper() for r in aln}
    L = len(next(iter(recs.values())))
    refs = [s for n, s in recs.items() if ref_sel(n)]
    core = [i for i in range(L) if sum(1 for s in refs if s[i] != '-') >= 0.9 * len(refs)]
    inv = [i for i in core if len({s[i] for s in refs}) == 1]
    var = [i for i in core if i not in inv]
    for gname, gsel in groups:
        fun = [s for n, s in recs.items() if n.startswith('FUNGI_') and gsel(n)]
        if not fun:
            continue
        cons_all = 0; cons90 = 0; ids = []; ids_var = []
        for i in inv:
            r = refs[0][i]
            m = sum(1 for s in fun if s[i] == r)
            ids.append(m / len(fun))
            cons_all += (m == len(fun)); cons90 += (m >= 0.9 * len(fun))
        for i in var:
            col = [s[i] for s in refs]
            top = max(set(col), key=col.count)
            ids_var.append(sum(1 for s in fun if s[i] == top) / len( fun))
        mi = sum(ids) / max(1, len(ids)); mv = sum(ids_var) / max(1, len(ids_var))
        out.write('\t'.join(map(str, [fam, gname, len(refs), len(fun), len(inv), cons_all, cons90, f'{mi:.3f}', f'{mv:.3f}'])) + '\n')
        print(fam, gname, 'refs', len(refs), 'fungal', len(fun), 'invariant', len(inv), 'all', cons_all, 'ge90', cons90, f'id {mi:.3f} vs control {mv:.3f}')

intact = lambda n: '_intact__' in n
run('SUS', lambda n: n.startswith('REF_') and '_sus__' in n, [('SUS T-DNA clade intact', lambda n: n.startswith('FUNGI_sus') and intact(n))])
run('CUS', lambda n: n.startswith('REF_') and '_cus__' in n, [('CUS cus intact', lambda n: n.startswith('FUNGI_cus_') and 'actino' not in n and intact(n) and 'JBBPKX010002805' not in n)])
run('ACS', lambda n: n.startswith('REF_') and '_acs__' in n and 'PON40589' not in n, [('ACS six lineage intact', lambda n: any(s in n for s in SIX)), ('ACS gdpd_meso all', lambda n: n.startswith('FUNGI_gdpd_meso'))])
run('ACS', lambda n: n.startswith('REF_') and '_gdpd_meso__' in n, [('ACS gdpd_meso vs acc refs', lambda n: n.startswith('FUNGI_gdpd_meso'))])
out.close()
