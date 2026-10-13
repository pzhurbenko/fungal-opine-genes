#!/usr/bin/env python3
"""Sensitivity of the protein-screen rules (contamination / A / B / C) of analysis/funnel/build_funnel.py.
The rule logic is COPIED from build_funnel.py (importing it would rewrite paper/tables/); the original is not modified.
Baseline: contamination o_bact_pid>=90; A: o_top50_bact>=35 & o_h>=50; B: pident_seed>=35 & cov_seed>=0.4 & pident_seed>=0.9*o_bact_pid
(propagated within gene x genus); C: A & not B & Rhizobiales best bacterial hit with o_bact_pid>=50 (propagated within gene x genus).
Outputs (analysis/thresholds/): funnel_grid.tsv (full factorial grid, contamination fixed at 90), funnel_oat.tsv (one-at-a-time sweeps incl. contamination 80-95),
funnel_flips.tsv (proteins whose fate changes vs baseline), margins_116.tsv (distance of each of the 116 proteins from every threshold),
contamination_gap.tsv (item 2: sorted best bacterial identity, Rhizobiales flag).
Env: PlantCad2 (pandas). Run: python sensitivity_funnel.py
"""
import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
import itertools, re, numpy as np, pandas as pd
R = f"{GMO_ROOT}/analysis/refset/screen/seeds_ext/reciprocal/reciprocal.tsv"
OUT = f"{GMO_ROOT}/analysis/thresholds/"
RHIZ = re.compile(r"^(Agrobacterium|Rhizobium|Allorhizobium|Neorhizobium|Pararhizobium|Sinorhizobium|Ensifer|Mesorhizobium|Bradyrhizobium|Devosia|Labrys|Bosea|Roseibium|"
                  r"Phyllobacterium|Ochrobactrum|Brucella|Methylobacterium|Aminobacter|Shinella|Gillisella|Rhizobiaceae|Hyphomicrobiales|Rhizobiales|Xanthobacter|Azorhizobium|"
                  r"Nitrobacter|Rhodopseudomonas|Pseudorhizobium|Chelativorans|Pelagibacterium|Rhizobium/Agrobacterium|Afipia|Microvirga|Hoeflea|Aquamicrobium|Nitratireductor)", re.I)
BASE = dict(contam=90, A_hits=35, A_margin=50, B_ident=35, B_cov=0.4, B_ratio=0.9, C_ident=50)
GRID = dict(A_hits=[25, 30, 35, 40, 45], A_margin=[20, 30, 50, 70, 80], B_ident=[30, 35, 40, 45], B_cov=[0.3, 0.4, 0.5, 0.6], B_ratio=[0.8, 0.9, 1.0], C_ident=[40, 45, 50, 55, 60])
OAT = dict(GRID, A_margin=[10, 15, 20, 30, 50, 70, 80], contam=[75, 80, 85, 90, 95])  # 10/15 added to reach the Beauveria margins
FATES = ["fungal", "contamination", "marginal", "B", "C", "single"]  # single = A, not B/C (unrelated family, one species)

r = pd.read_csv(R, sep="\t")
r["genus"] = r.organism.str.split().str[0]
grp, _ = pd.factorize(pd.Series(list(zip(r.gene, r.genus))))
ng = grp.max() + 1
nf = (r.o_h > 0).values
bp, h, t50, ps, cs = (r[c].values.astype(float) for c in ["o_bact_pid", "o_h", "o_top50_bact", "pident_seed", "cov_seed"])
rhiz = r.o_bact_org.fillna("").apply(lambda o: bool(RHIZ.match(o))).values


def classify(contam, A_hits, A_margin, B_ident, B_cov, B_ratio, C_ident):
    con = nf & (bp >= contam)
    A = nf & ~con & (t50 >= A_hits) & (h >= A_margin)
    b0 = A & (ps >= B_ident) & (cs >= B_cov) & (ps >= B_ratio * bp)
    c0 = A & ~b0 & (bp >= C_ident) & rhiz
    B = A & (np.bincount(grp, weights=b0, minlength=ng) > 0)[grp]
    C = A & ~B & (np.bincount(grp, weights=c0, minlength=ng) > 0)[grp]
    f = np.full(len(r), "single", dtype=object)
    f[~nf] = "fungal"; f[con] = "contamination"; f[nf & ~con & ~A] = "marginal"; f[B] = "B"; f[C] = "C"
    return f


def summary(f):
    d = {k: int((f == k).sum()) for k in FATES}
    d["ruleA_total"] = d["B"] + d["C"] + d["single"]
    for k in ["B", "C"]:
        d[f"pairs_{k}"] = len(set(zip(r.gene[f == k], r.genus[f == k])))
    return d


base = classify(**BASE)
assert summary(base)["B"] == 44 and summary(base)["C"] == 6 and summary(base)["contamination"] == 11, summary(base)  # = funnel_numbers.txt
lab = (r.genus + " " + r.gene + " " + r.protein).values

# ---- full factorial grid (contamination fixed at 90)
rows, diffcount, altfates = [], np.zeros(len(r), int), [set() for _ in range(len(r))]
keys = list(GRID)
for vals in itertools.product(*GRID.values()):
    p = dict(BASE, **dict(zip(keys, vals)))
    f = classify(**p); d = (f != base)
    diffcount += d
    for i in np.where(d)[0]: altfates[i].add(f[i])
    ch = "; ".join(f"{lab[i]}:{base[i]}->{f[i]}" for i in np.where(d)[0])
    rows.append(dict(p, **summary(f), n_changed=int(d.sum()), changed=ch))
grid = pd.DataFrame(rows); grid.to_csv(OUT + "funnel_grid.tsv", sep="\t", index=False)
NG = len(grid)

# ---- one-at-a-time sweeps
orows, flips = [], {}
for k, vs in OAT.items():
    for v in vs:
        p = dict(BASE, **{k: v}); f = classify(**p); d = np.where(f != base)[0]
        orows.append(dict(param=k, value=v, is_base=v == BASE[k], **summary(f), n_changed=len(d),
                          changed="; ".join(f"{lab[i]}:{base[i]}->{f[i]}" for i in d)))
        for i in d: flips.setdefault(i, []).append(f"{k}={v}:{f[i]}")
pd.DataFrame(orows).to_csv(OUT + "funnel_oat.tsv", sep="\t", index=False)

fl = []
for i in sorted(set(flips) | set(np.where(diffcount > 0)[0])):
    fl.append(dict(protein=r.protein[i], genus=r.genus[i], organism=r.organism[i], gene=r.gene[i], base_fate=base[i],
                   oat_flips="; ".join(flips.get(i, [])) or "none (only in multi-parameter combinations)",
                   grid_frac_changed=round(diffcount[i] / NG, 3), grid_alt_fates=",".join(sorted(altfates[i]))))
pd.DataFrame(fl).sort_values(["grid_frac_changed"], ascending=False).to_csv(OUT + "funnel_flips.tsv", sep="\t", index=False)

# ---- margins of the 116 proteins with non-fungal nearest relatives
m = r[nf].copy(); m["base_fate"] = base[nf]
m["rhizobial_best_bact"] = rhiz[nf]
m["m_contam_pid"] = m.o_bact_pid - BASE["contam"]
m["m_A_hits"] = m.o_top50_bact - BASE["A_hits"]
m["m_A_margin_bits"] = m.o_h - BASE["A_margin"]
m["m_B_ident"] = m.pident_seed - BASE["B_ident"]
m["m_B_cov"] = (m.cov_seed - BASE["B_cov"]).round(3)
m["m_B_ratio"] = (m.pident_seed / m.o_bact_pid - BASE["B_ratio"]).round(3)
m["m_C_ident"] = np.where(m.rhizobial_best_bact, m.o_bact_pid - BASE["C_ident"], np.nan)
own_B = m.o_top50_bact.ge(0) & (m.pident_seed >= 35) & (m.cov_seed >= 0.4) & (m.pident_seed >= 0.9 * m.o_bact_pid)
m["B_own_pass"] = own_B  # whether the protein itself passes B criteria (vs. via gene x genus propagation)
m["grid_frac_changed"] = diffcount[nf] / NG
cols = ["protein", "organism", "genus", "gene", "base_fate", "o_bact_org", "rhizobial_best_bact", "o_bact_pid", "o_top50_bact", "o_h", "pident_seed", "cov_seed",
        "m_contam_pid", "m_A_hits", "m_A_margin_bits", "m_B_ident", "m_B_cov", "m_B_ratio", "B_own_pass", "m_C_ident", "grid_frac_changed"]
m.sort_values(["base_fate", "gene", "organism"])[cols].to_csv(OUT + "margins_116.tsv", sep="\t", index=False)

# ---- item 2: contamination gap
c = m[["protein", "organism", "gene", "o_bact_org", "o_bact_pid", "rhizobial_best_bact", "base_fate"]].sort_values("o_bact_pid", ascending=False)
c.to_csv(OUT + "contamination_gap.tsv", sep="\t", index=False)
v = np.sort(m.o_bact_pid.values)[::-1]
gaps = -np.diff(v); j = gaps.argmax()
print(f"largest gap in best bacterial identity: {v[j+1]} -> {v[j]} ({gaps[j]:.1f} points)")
for t in [80, 85, 90, 95]: print(f"contamination at {t}%: {(m.o_bact_pid >= t).sum()}")
nr = m[~m.rhizobial_best_bact].o_bact_pid.sort_values(ascending=False).values
print("non-Rhizobiales best bacterial hit, top values:", nr[:8])
print("grid points:", NG, "| points with no change vs baseline:", int((grid.n_changed == 0).sum()))
# ---- exact stability interval per threshold (fine one-at-a-time scan; other thresholds at baseline)
FINE = dict(contam=np.arange(60, 100.01, 0.1), A_hits=np.arange(0, 51, 1), A_margin=np.arange(0, 201, 1), B_ident=np.arange(20, 65.01, 0.1),
            B_cov=np.arange(0.0, 1.001, 0.01), B_ratio=np.arange(0.5, 1.3001, 0.01), C_ident=np.arange(25, 80.01, 0.1))
srows = []
for k, vs in FINE.items():
    ok = np.array([(classify(**dict(BASE, **{k: v})) == base).all() for v in vs])
    i0 = int(np.argmin(np.abs(vs - BASE[k]))); lo = hi = i0
    while lo > 0 and ok[lo - 1]: lo -= 1
    while hi < len(vs) - 1 and ok[hi + 1]: hi += 1
    below = "" if lo == 0 else "; ".join(f"{lab[i]}:{base[i]}->{f[i]}" for f in [classify(**dict(BASE, **{k: vs[lo - 1]}))] for i in np.where(f != base)[0])
    above = "" if hi == len(vs) - 1 else "; ".join(f"{lab[i]}:{base[i]}->{f[i]}" for f in [classify(**dict(BASE, **{k: vs[hi + 1]}))] for i in np.where(f != base)[0])
    srows.append(dict(param=k, base=BASE[k], stable_from=round(vs[lo], 3), stable_to=round(vs[hi], 3), scan_min=round(vs[0], 3), scan_max=round(vs[-1], 3),
                      first_change_below=f"at {vs[lo-1]:.3g}: {below}" if lo else "none in scan", first_change_above=f"at {vs[hi+1]:.3g}: {above}" if hi < len(vs) - 1 else "none in scan"))
st = pd.DataFrame(srows); st.to_csv(OUT + "funnel_stability.tsv", sep="\t", index=False)
print(st.to_string(index=False, max_colwidth=110))
print(pd.DataFrame(orows)[["param", "value", "contamination", "marginal", "B", "C", "single", "pairs_B", "pairs_C", "n_changed"]].to_string(index=False))
