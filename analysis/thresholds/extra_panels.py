# executed by make_thresholds_figure.py (shares its namespace: panel, PANELS, T, colours)
import os
ci = pd.read_csv(T + "chs/chs_identity.tsv", sep="\t")
dp = ci[ci.source == "diamond_protein"]
lin = dp.cls.str.startswith("chs_lineage") & (dp.pid_to_chs_seed >= 0.4)  # 2 Dactylonectria paralogues at 24-25% are native homologues of a chs-lineage genus
PANELS.append(lambda ax: panel(ax, [(100 * dp.pid_to_chs_seed[~lin] if dp.pid_to_chs_seed.max() <= 1 else dp.pid_to_chs_seed[~lin], "native homologues", GREY),
                                    (100 * dp.pid_to_chs_seed[lin] if dp.pid_to_chs_seed.max() <= 1 else dp.pid_to_chs_seed[lin], "chs lineage", VERM)],
                               45, "identity to chs seed, %", "h  chs-like proteins:\n≥45% identity to seed", np.arange(10, 70, 2), (37.4, 44.0)))
vi = pd.read_csv(T + "chs/vanettenii_identity.tsv", sep="\t")
dep = vi.depends_on_70pct_rule == 1
PANELS.append(lambda ax: panel(ax, [(100 * vi.pid_to_vanettenii_EEU41297[~dep], "pass via chs seed", GREY), (100 * vi.pid_to_vanettenii_EEU41297[dep], "pass only via 70% criterion", VERM)],
                               70, "identity to F. vanettenii, %", "i  chs-like loci:\n≥70% to F. vanettenii", np.arange(50, 101, 2), (50, 75.2)))
fi = pd.read_csv(T + "pkerr/fragment_identity.tsv", sep="\t")
fr = fi[fi.is_six == 0]; asg = fr.basis == "similarity"
PANELS.append(lambda ax: panel(ax, [(fr.pid1[~asg].dropna(), "not assigned", GREY), (fr.pid1[asg], "assigned fragments", VERM)],
                               60, "identity to nearest intact, %", "j  pKerr14a-type fragments:\n≥60% identity", np.arange(20, 101, 2.5), (49.5, 62.7)))
bb = T + "backbone/backbone_contigs.tsv"
if os.path.exists(bb):
    b = pd.read_csv(bb, sep="\t")
    b["best_backbone_pid"] = b.best_backbone_pid.fillna(0)
    PANELS.append(lambda ax: panel(ax, [(b.best_backbone_pid[b["sample"] == "fungal"], "fungal locus contigs", GREY), (b.best_backbone_pid[b["sample"] == "bacterial"], "contaminant contigs", VERM)],
                                   60, "best id. to vir/tra/rep protein, %", "k  vir/tra/rep, test contigs:\n≥60% identity", np.arange(0, 101, 4), None))
bo = pd.read_csv(T + "backbone/screen_boletales_genome_maxbb.tsv", sep="\t", header=None, names=["gca", "organism", "pid", "contig", "len", "label"])
hb = bo.organism.str.contains("harrisonii")
PANELS.append(lambda ax: panel(ax, [(bo.pid[~hb], "other Boletales", GREY), (bo.pid[hb], "B. harrisonii", VERM)],
                               60, "max vir/tra/rep id., %", "l  vir/tra/rep, Boletales:\n≥60% identity", np.arange(20, 101, 2.5), (55.2, 89.6)))
