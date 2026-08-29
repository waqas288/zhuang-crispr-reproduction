import xenaPython as xena
import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

HUB = 'https://tcga.xenahubs.net'

def num(s):
    return pd.to_numeric(s, errors='coerce')

def fetch(expr_ds, surv_ds, genes):
    samples = xena.dataset_samples(HUB, expr_ds, None)
    vals = xena.dataset_probe_values(HUB, expr_ds, samples, genes)[1]
    expr = pd.DataFrame(dict(zip(genes, vals)), index=samples)
    for g in genes:
        expr[g] = num(expr[g])
    fields = xena.dataset_field(HUB, surv_ds)
    ssamp = xena.dataset_samples(HUB, surv_ds, None)
    sv = xena.dataset_fetch(HUB, surv_ds, ssamp, fields)
    surv = pd.DataFrame(dict(zip(fields, sv)), index=ssamp)
    surv['time'] = num(surv['OS.time'])
    surv['event'] = num(surv['OS'])
    return expr, surv[['time','event']]

def km(sub, ax, title):
    kmf = KaplanMeierFitter()
    lo = sub[sub.grp == 'IFNGR2 low']
    hi = sub[sub.grp == 'IFNGR2 high']
    for d_, lab in [(lo,'IFNGR2 low'), (hi,'IFNGR2 high')]:
        kmf.fit(d_['time'], d_['event'], label=f'{lab} (n={len(d_)})')
        kmf.plot_survival_function(ax=ax, ci_show=False)
    r = logrank_test(lo['time'], hi['time'], lo['event'], hi['event'])
    ax.set_title(f'{title}\nP = {r.p_value:.4f}', fontsize=10)
    ax.set_xlabel('days'); ax.set_ylabel('overall survival')
    return r.p_value, len(lo), len(hi)

results = {}
fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))

# ---------- LAML ----------
e, s = fetch('TCGA.LAML.sampleMap/HiSeqV2', 'survival/LAML_survival.txt', ['IFNGR2','NKG7'])
d = e.join(s, how='inner').dropna(subset=['IFNGR2','NKG7','time','event'])
print('LAML samples with expression + survival:', len(d))

q33, q66 = d.NKG7.quantile([1/3, 2/3])
d = d.assign(nkg7=np.where(d.NKG7 <= q33, 'low',
                    np.where(d.NKG7 >= q66, 'high', 'mid')))

for strat, ax in zip(['low','high'], axes[:2]):
    sub = d[d.nkg7 == strat].copy()
    q40, q60 = sub.IFNGR2.quantile([0.40, 0.60])
    sub = sub[(sub.IFNGR2 <= q40) | (sub.IFNGR2 >= q60)].copy()
    sub['grp'] = np.where(sub.IFNGR2 <= q40, 'IFNGR2 low', 'IFNGR2 high')
    results[f'LAML_NKG7-{strat}'] = km(sub, ax, f'TCGA-LAML, NKG7-{strat}')

# ---------- KIRC ----------
e2, s2 = fetch('TCGA.KIRC.sampleMap/HiSeqV2', 'survival/KIRC_survival.txt', ['IFNGR2'])
d2 = e2.join(s2, how='inner').dropna(subset=['IFNGR2','time','event'])
print('KIRC samples with expression + survival:', len(d2))

a33, a66 = d2.IFNGR2.quantile([1/3, 2/3])
k = d2[(d2.IFNGR2 <= a33) | (d2.IFNGR2 >= a66)].copy()
k['grp'] = np.where(k.IFNGR2 <= a33, 'IFNGR2 low', 'IFNGR2 high')
results['KIRC'] = km(k, axes[2], 'TCGA-KIRC')

fig.tight_layout()
fig.savefig('figures/fig6_tcga_survival.png', dpi=200)

paper = {'LAML_NKG7-low': (0.2680, '22 vs 22'),
         'LAML_NKG7-high': (0.2117, '22 vs 23'),
         'KIRC': (0.0001, '200 vs 201')}

print()
print('=== REPRODUCTION SCORECARD ===')
print(f"{'comparison':18s} {'our P':>9s} {'paper P':>9s}   {'our n':>12s}   {'paper n':>12s}")
for key, (p, nlo, nhi) in results.items():
    pp, pn = paper[key]
    print(f'{key:18s} {p:9.4f} {pp:9.4f}   {nlo:>4d} vs {nhi:<4d}   {pn:>12s}')
