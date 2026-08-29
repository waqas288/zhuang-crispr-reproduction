"""
Reproduces the visual format of Zhuang et al. 2019 Figures 1D, 2A and 2B.

Fig 1D/2A: gene index (x) vs normalised score (y), top hits labelled.
           Paper plots enriched and depleted separately, one panel per
           selection pressure. Genes found in BOTH pressures are red.
Fig 2B:    paired sgRNA counts, Ctrl -> Surv, one line per guide.
"""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from adjustText import adjust_text

BOTH = ['NCR3LG1','IFNGR2','B2M','HLA-E','JAK1','JAK2','STAT1','ABL1','PSMB5']

low  = pd.read_csv('low_rra.gene_summary.txt',  sep='\t')
high = pd.read_csv('high_rra.gene_summary.txt', sep='\t')


def scatter(df, direction, ax, title, n_label=14):
    """direction: 'pos' (enriched) or 'neg' (depleted)"""
    lfc = df['pos|lfc'] if direction == 'pos' else df['neg|lfc']
    rank = df[f'{direction}|rank']
    d = df.assign(score=lfc, rank=rank).sort_values('id').reset_index(drop=True)
    d['is_top'] = d['rank'] <= n_label
    d['idx'] = np.arange(len(d))

    ax.scatter(d.idx, d.score, s=6, c='#999999', rasterized=True)

    top = d[d.is_top]
    texts = []
    for _, r in top.iterrows():
        c = '#d62728' if r['id'] in BOTH else '#555555'
        ax.scatter(r.idx, r.score, s=28, c=c, zorder=5)
        texts.append(ax.text(r.idx, r.score, r['id'], fontsize=7, color=c))
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', lw=.4, color='grey'))

    ax.set_xlabel('Genes'); ax.set_ylabel('Score')
    ax.set_title(title, fontsize=10)


# ---------- Figure 1D style: ENRICHED ----------
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
scatter(low,  'pos', axes[0], 'Enriched sgRNAs\n(Low selection pressure)')
scatter(high, 'pos', axes[1], 'Enriched sgRNAs\n(High selection pressure)')
fig.tight_layout(); fig.savefig('figures/fig8_paperstyle_enriched.png', dpi=200)

# ---------- Figure 2A style: DEPLETED ----------
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
scatter(low,  'neg', axes[0], 'Depleted sgRNAs\n(Low selection pressure)')
scatter(high, 'neg', axes[1], 'Depleted sgRNAs\n(High selection pressure)')
fig.tight_layout(); fig.savefig('figures/fig9_paperstyle_depleted.png', dpi=200)

# ---------- Figure 2B style: paired guide counts ----------
cnt = pd.read_csv('zhuang.count.txt', sep='\t')
cols = ['Low_Ctrl_1','Low_Ctrl_2','Low_Surv_1','Low_Surv_2']
norm = cnt[cols] / cnt[cols].sum() * cnt[cols].sum().mean()
d2 = cnt[['sgRNA','Gene']].join(norm)
d2['Ctrl'] = d2[['Low_Ctrl_1','Low_Ctrl_2']].mean(axis=1)
d2['Surv'] = d2[['Low_Surv_1','Low_Surv_2']].mean(axis=1)

genes = ['IFNGR2','B2M','HLA-E']
fig, axes = plt.subplots(1, 3, figsize=(11, 4))
for g, ax in zip(genes, axes):
    sub = d2[d2.Gene == g]
    for _, r in sub.iterrows():
        ax.plot([0, 1], [r.Ctrl, r.Surv], marker='o', ms=6, lw=1.2)
    ax.set_xticks([0, 1]); ax.set_xticklabels(['Ctrl', 'Surv.'])
    ax.set_xlim(-0.3, 1.3); ax.set_ylabel('sgRNA count')
    ax.set_title(g, fontsize=11)
fig.tight_layout(); fig.savefig('figures/fig10_paperstyle_guidecounts.png', dpi=200)

print('wrote fig8, fig9, fig10')
