import requests, pandas as pd, time

API = "https://string-db.org/api"
SPECIES = 9606          # Homo sapiens
CALLER  = "zhuang_repro"

sets = {
 "antigen_presentation": ["B2M","HLA-A","HLA-C","HLA-E","TAP1","TAPBP","PSMB5"],
 "interferon_signaling": ["JAK1","JAK2","IFNGR1","IFNGR2","STAT1","IRF1","IFITM1","TAP1"],
}

for name, genes in sets.items():
    print("="*60)
    print(name, "|", len(genes), "genes")

    # --- network image (PNG), same content as paper Fig 2D/2E ---
    r = requests.post(f"{API}/highres_image/network", data={
        "identifiers": "\r".join(genes),
        "species": SPECIES,
        "caller_identity": CALLER,
        "network_flavor": "confidence",
    })
    r.raise_for_status()
    out = f"figures/fig7_string_{name}.png"
    open(out, "wb").write(r.content)
    print("  wrote", out)
    time.sleep(1)

    # --- interaction table ---
    r = requests.post(f"{API}/tsv/network", data={
        "identifiers": "\r".join(genes),
        "species": SPECIES,
        "caller_identity": CALLER,
    })
    r.raise_for_status()
    lines = [l.split("\t") for l in r.text.strip().split("\n")]
    df = pd.DataFrame(lines[1:], columns=lines[0])
    df["score"] = df["score"].astype(float)
    df = df[["preferredName_A","preferredName_B","score"]].sort_values("score", ascending=False)
    df.to_csv(f"results_string_{name}.csv", index=False)
    print(f"  {len(df)} interactions, median score {df.score.median():.3f}")
    print(df.head(12).to_string(index=False))
    time.sleep(1)

    # --- enrichment: does STRING call these pathways itself? ---
    r = requests.post(f"{API}/tsv/enrichment", data={
        "identifiers": "\r".join(genes),
        "species": SPECIES,
        "caller_identity": CALLER,
    })
    r.raise_for_status()
    lines = [l.split("\t") for l in r.text.strip().split("\n")]
    e = pd.DataFrame(lines[1:], columns=lines[0])
    e["fdr"] = e["fdr"].astype(float)
    top = e[e.category.isin(["Process","KEGG","RCTM"])].nsmallest(8, "fdr")
    print("\n  top enriched terms:")
    print(top[["category","description","fdr"]].to_string(index=False))
    print()
    time.sleep(1)
