"""
Reproduces Zhuang et al. 2019 Figure 2C (pathway analysis).

The paper used Ingenuity Pathway Analysis (IPA), which is commercial.
Documented substitution: STRING enrichment API (free).

Two analyses:
  A. Unbiased  - top 100 depleted genes from our own MLE result.
                 This is the real reproduction: does the pathway signal
                 emerge without being told what to look for?
  B. Confirmatory - the paper's own hand-picked gene sets (Fig 2D/2E).
                 Circular by construction, but shows the sets are coherent.
"""
import pandas as pd
import requests

API = "https://string-db.org/api/tsv/enrichment"
SPECIES = 9606
CALLER = "zhuang_repro"
CATS = ["Process", "KEGG", "RCTM"]


def enrich(genes, label, n_show=20):
    r = requests.post(API, data={
        "identifiers": "\r".join(genes),
        "species": SPECIES,
        "caller_identity": CALLER,
    })
    r.raise_for_status()
    rows = [l.split("\t") for l in r.text.strip().split("\n")]
    df = pd.DataFrame(rows[1:], columns=rows[0])
    df["fdr"] = df["fdr"].astype(float)
    df["number_of_genes"] = df["number_of_genes"].astype(int)
    df = df[df.category.isin(CATS)].nsmallest(n_show, "fdr")

    print("=" * 70)
    print(f"{label}  (n = {len(genes)} genes)")
    print("=" * 70)
    print(df[["category", "description", "fdr", "number_of_genes"]].to_string(index=False))
    print()
    df.to_csv(f"results_enrichment_{label.replace(' ', '_').lower()}.csv", index=False)
    return df


# ---------- A. unbiased: our own top hits ----------
ml = pd.read_csv("low_mle.gene_summary.txt", sep="\t")
real = ml[~ml.Gene.str.startswith(("NonTargetingControl", "hsa-mir"))]
top100 = real.nsmallest(100, "NKselection|beta")["Gene"].tolist()
a = enrich(top100, "unbiased top100 depleted")

# did the paper's two themes appear on their own?
hits = a[a.description.str.contains("nterferon|ntigen|MHC", regex=True, case=False)]
print(">>> paper's themes recovered unsupervised:", len(hits), "terms")
if len(hits):
    print(hits[["description", "fdr"]].to_string(index=False))
print()

# ---------- B. confirmatory: paper's own sets ----------
enrich(["B2M", "HLA-A", "HLA-C", "HLA-E", "TAP1", "TAPBP", "PSMB5"],
       "paper set antigen presentation", n_show=8)
enrich(["JAK1", "JAK2", "IFNGR1", "IFNGR2", "STAT1", "IRF1", "IFITM1", "TAP1"],
       "paper set interferon signaling", n_show=8)
