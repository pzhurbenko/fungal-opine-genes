import sys as _sys, os as _os; _sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..')))
from gmo_paths import GMO_ROOT  # noqa: E402
#!/usr/bin/env python3
import re, os
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
D = f'{GMO_ROOT}/analysis/dnds'
fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
for ax, fam in zip(axes, ['SUS', 'CUS']):
    f = f'{D}/{fam}.yn00.out'
    if not os.path.exists(f): ax.set_title(f'{fam}: нет данных'); continue
    vals = []
    for l in open(f):
        m = re.match(r'^\s*\d+\s+\d+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+([\d.]+)\s+[-\d.]+\s+\+-\s+[\d.]+\s+([-\d.]+)', l)
        if m and float(m.group(2)) > 0.01 and float(m.group(1)) < 5: vals.append(float(m.group(1)))
    ax.hist(vals, bins=40, color='#4a90d9', edgecolor='white'); ax.axvline(1, color='k', ls='--'); ax.set_xlabel('dN/dS (yn00, попарно)'); ax.set_ylabel('число пар'); ax.set_title(f'{fam}: {len(vals)} пар, медиана {sorted(vals)[len(vals)//2]:.2f}')
plt.tight_layout(); plt.savefig('fig/FB2_dnds.png', dpi=150); plt.savefig('fig/FB2_dnds.svg'); print('ok')
