"""
QC and analysis figures (not reproductions of the paper's panels --
see 08_paper_figures.py for those).

fig1  volcano, low-pressure MLE, with the 1,000 non-targeting controls overlaid
fig2  the empirical noise floor -- KEY FIGURE
fig3  per-guide CPM, controls vs survivors
fig4  RRA vs MLE rank-rank scatter
fig5  library representation QC

FIGURE ERRORS MADE AND CORRECTED (recorded deliberately)
--------------------------------------------------------
v1 volcano clipped FDR at 1e-6, so six genes at FDR=0 stacked on one
   horizontal line with overlapping labels. Misleading -- implied identical
   significance.
v2 plotted z-score against a p-value DERIVED FROM z. Every point fell on a
   smooth curve because that just draws the mapping function. Redundant axes.
v3 (this version) uses beta vs Wald p -- two genuinely different quantities.

fig3 v1 used symlog, producing a meaningless negative CPM axis, and
   interleaved Ctrl/Surv so the direction was unreadable. Now grouped, log.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from adjustText import adjust_text

TARGETS = [("IFNGR2", "#2a9d8f"), ("JAK1", "#2a9d8f"), ("B2M", "#2a9d8f"),
           ("STAT1", "#2a9d8f"), ("HLA-E", "#2a9d8f"), ("JAK2", "#2a9d8f"),
           ("IRF1", "#2a9d8f"), ("NCR3LG1", "#d1495b")]

ml  = pd.read_csv("low_mle.gene_summary.txt",  sep="\t")
mh  = pd.read_csv("high_mle.gene_summary.txt", sep="\t")
rra = pd.read_csv("low_rra.gene_summary.txt",  sep="\t")
cnt = pd.read_csv("zhuang.count.txt",          sep="\t")

SAMPLES = ["Low_Ctrl_1", "Low_Surv_1", "Low_Ctrl_2", "Low_Surv_2",
           "High_Ctrl", "High_1", "High_2"]


# ---------------------------------------------------------------- fig 1
# Volcano: effect size (beta) vs significance (Wald p).
# The blue cloud is the 1,000 non-targeting controls -- they form a band
# below -log10(p) ~ 1.5, and every labelled gene sits above it.
fig, ax = plt.subplots(figsize=(7.5, 6))

x = ml["NKselection|beta"]
y = -np.log10(ml["NKselection|wald-p-value"].clip(lower=1e-12))
ax.scatter(x, y, s=4, c="#d9d9d9", rasterized=True)

nt = ml[ml.Gene.str.startswith("NonTargetingControl")]
ax.scatter(nt["NKselection|beta"],
           -np.log10(nt["NKselection|wald-p-value"].clip(lower=1e-12)),
           s=8, c="#4a90d9", alpha=.6, label="non-targeting (n=1000)")

labels = []
for gene, colour in TARGETS:
    row = ml[ml.Gene == gene]
    if not len(row):
        continue
    x0 = row["NKselection|beta"].iloc[0]
    y0 = -np.log10(max(row["NKselection|wald-p-value"].iloc[0], 1e-12))
    ax.scatter(x0, y0, s=55, c=colour, edgecolors="k", linewidths=.6, zorder=5)
    labels.append(ax.text(x0, y0, gene, fontsize=9))
adjust_text(labels, ax=ax, arrowprops=dict(arrowstyle="-", lw=.5))

ax.axvline(0, lw=.5, c="k")
ax.set_xlabel("MLE beta (NK selection)")
ax.set_ylabel("-log10(Wald p)")
ax.set_title("Low selection pressure — MAGeCK MLE")
ax.legend(fontsize=8, loc="upper center")
fig.tight_layout()
fig.savefig("figures/fig1_volcano_low_mle.png", dpi=200)


# ---------------------------------------------------------------- fig 2
# THE KEY FIGURE.
# The 1,000 non-targeting guides cut nothing, so their beta distribution IS
# the noise. Low arm SD 0.210, high arm SD 1.189 -- 5.7x wider.
#
# This resolves an apparent paradox: high-pressure effect sizes are 2-3x
# larger, yet only 1 gene reaches FDR<0.05 there vs 7 in the low arm.
# Normalised to noise, IFNGR2 falls from 4.4 SD-units to 1.8.
# The high-pressure arm has LARGER EFFECTS AND LESS POWER.
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)

for ax, (d, name) in zip(axes, [(ml, "Low pressure"), (mh, "High pressure")]):
    beta = d[d.Gene.str.startswith("NonTargetingControl")]["NKselection|beta"]
    ax.hist(beta, bins=60, color="#4a90d9", alpha=.75,
            label=f"controls: SD={beta.std():.3f}")
    for gene, colour in [("IFNGR2", "#2a9d8f"), ("B2M", "#e76f51"),
                         ("NCR3LG1", "#d1495b")]:
        row = d[d.Gene == gene]
        if len(row):
            ax.axvline(row["NKselection|beta"].iloc[0], c=colour, lw=2, label=gene)
    ax.set_title(name)
    ax.set_xlabel("MLE beta")
    ax.legend(fontsize=8)

axes[0].set_ylabel("number of non-targeting guides")
fig.suptitle("Empirical noise floor: guides that cut nothing", y=1.02)
fig.tight_layout()
fig.savefig("figures/fig2_noise_floor.png", dpi=200, bbox_inches="tight")


# ---------------------------------------------------------------- fig 3
# Per-guide CPM. Controls grouped left of the dotted line, survivors right,
# so direction reads left-to-right. Log scale (NOT symlog -- negative CPM
# is meaningless).
LOW  = ["Low_Ctrl_1", "Low_Ctrl_2", "Low_Surv_1", "Low_Surv_2"]
HIGH = ["High_Ctrl", "High_1", "High_2"]

cpm = cnt[LOW + HIGH] / cnt[LOW + HIGH].sum() * 1e6
d2 = cnt[["sgRNA", "Gene"]].join(cpm)

fig, axes = plt.subplots(2, 2, figsize=(11, 8))
for row_i, (gene, label) in enumerate([("NCR3LG1", "NCR3LG1 — enriched"),
                                       ("IFNGR2", "IFNGR2 — depleted")]):
    sub = d2[d2.Gene == gene]
    for col_i, (cols, arm) in enumerate([(LOW, "Low pressure"),
                                         (HIGH, "High pressure")]):
        ax = axes[row_i][col_i]
        for _, guide in sub.iterrows():
            ax.plot(range(len(cols)),
                    np.maximum(guide[cols].values, 0.05),
                    marker="o", ms=4, lw=1, alpha=.75)
        ax.set_yscale("log")
        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels(cols, rotation=30, ha="right", fontsize=8)
        ax.set_ylabel("CPM (log)")
        ax.set_title(f"{label} — {arm}", fontsize=10)
        ax.axvline(1.5 if col_i == 0 else 0.5, ls=":", c="grey", lw=.8)

fig.suptitle("Guide-level counts, controls left of dotted line", y=1.0)
fig.tight_layout()
fig.savefig("figures/fig3_guide_level.png", dpi=200, bbox_inches="tight")


# ---------------------------------------------------------------- fig 4
# Where the two methods disagree. RRA rewards CONSISTENCY OF DIRECTION across
# guides; MLE rewards EFFECT SIZE with variance accounted for. Neither is
# wrong -- reporting both and explaining the difference is stronger than
# picking one and hiding the other.
merged = rra[["id", "neg|rank"]].merge(
    ml[["Gene", "NKselection|beta"]], left_on="id", right_on="Gene")
merged["mle_rank"] = merged["NKselection|beta"].rank()

fig, ax = plt.subplots(figsize=(6, 6))
ax.scatter(np.log10(merged["neg|rank"]), np.log10(merged["mle_rank"]),
           s=3, c="#cccccc", rasterized=True)

for gene, _ in TARGETS:
    row = merged[merged.id == gene]
    if not len(row):
        continue
    x0 = np.log10(row["neg|rank"].iloc[0])
    y0 = np.log10(row["mle_rank"].iloc[0])
    ax.scatter(x0, y0, s=45, c="#2a9d8f", edgecolors="k", linewidths=.5, zorder=5)
    ax.annotate(gene, (x0, y0), fontsize=8, xytext=(4, 4),
                textcoords="offset points")

ax.plot([0, 4.4], [0, 4.4], ls="--", lw=.8, c="k")
ax.set_xlabel("log10 RRA rank (depleted)")
ax.set_ylabel("log10 MLE rank (by beta)")
ax.set_title("RRA vs MLE — where the methods disagree")
fig.tight_layout()
fig.savefig("figures/fig4_rra_vs_mle.png", dpi=200)


# ---------------------------------------------------------------- fig 5
# Library representation. The high-pressure samples spike at zero --
# 26,238 guides entirely absent from High_1. That is SELECTION, not drift:
# NCR3LG1 sits at 1683 CPM there against 7.7 in its control (~200x).
# Three rounds of 90% killing left a population dominated by resistant clones.
GINI = [0.090, 0.088, 0.092, 0.098, 0.226, 0.386, 0.366]
ZERO = [487, 385, 536, 1001, 11620, 26238, 23427]

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

for s in SAMPLES:
    axes[0].hist(np.log10(cnt[s] + 1), bins=80, histtype="step", lw=1.2, label=s)
axes[0].set_xlabel("log10(count+1)")
axes[0].set_ylabel("guides")
axes[0].legend(fontsize=7)
axes[0].set_title("Count distribution per sample")

axes[1].bar(range(7), GINI, color="#4a90d9")
axes[1].set_xticks(range(7))
axes[1].set_xticklabels(SAMPLES, rotation=45, ha="right", fontsize=8)
axes[1].set_ylabel("Gini index")
axes[1].set_title("Library evenness (higher = more skewed)")
for i, z in enumerate(ZERO):
    axes[1].text(i, GINI[i] + .008, f"{z:,}", ha="center", fontsize=7)

fig.tight_layout()
fig.savefig("figures/fig5_library_qc.png", dpi=200)

print("wrote fig1 - fig5")
