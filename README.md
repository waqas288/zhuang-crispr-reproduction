# Reproducing a genome-wide CRISPR screen: Zhuang et al. 2019

An independent reproduction of the computational analyses in **Zhuang X, Veltri DP, Long EO (2019), *Frontiers in Immunology* 10:2879**, starting from raw sequencing reads.

**Paper:** [doi:10.3389/fimmu.2019.02879](https://doi.org/10.3389/fimmu.2019.02879) (open access, CC BY)
**Data:** GEO [GSE139313](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139313) · SRA SRP226783

---

## What this is

A pooled CRISPR knockout screen in the leukaemia cell line K562, selected by co-incubation with primary human NK cells. The paper asks which gene knockouts make cancer cells resistant to immune killing, and which make them more vulnerable.

I started from the seven raw FASTQ files and reran the analysis independently — guide counting, QC, hit calling, pathway enrichment, protein network analysis, and the TCGA survival analysis.

**This is a reproduction, not new research.** The purpose is to demonstrate that I can take a published screen from raw reads to biological conclusions and check my answers against a known result.

---

## Result

**Four of the paper's five main figures reproduced in full.** Figures 3 and 4 are wet-lab experiments and are not reproducible computationally.

| Paper figure | Analysis | Status |
|---|---|---|
| Fig 1D | Enriched sgRNAs, both selection pressures | ✅ |
| Fig 2A | Depleted sgRNAs, both selection pressures | ✅ |
| Fig 2B | Paired sgRNA counts | ✅ |
| Fig 2C | Pathway analysis | ✅ (tool substituted — see below) |
| Fig 2D, 2E | STRING protein networks | ✅ |
| Fig 5 | TCGA survival analysis | ✅ |
| Fig 1B, 1C, 3, 4, S1 | Cytotoxicity, flow cytometry, arrayed validation | ❌ wet lab |

### The paper's headline findings, recovered independently

**Resistance — NCR3LG1**

| | Rank | FDR | Effect |
|---|---|---|---|
| Low pressure | **1** of 21,697 | 0.005 | LFC +1.20 |
| High pressure | **1** of 21,697 | 0.002 | LFC +7.81 |

**Vulnerability — IFN-γ signalling and antigen presentation** (MAGeCK MLE, low pressure)

| Gene | Rank | beta | FDR |
|---|---|---|---|
| JAK1 | 1 | −0.992 | 0.000 |
| IFNGR2 | 2 | −0.933 | 0.000 |
| B2M | 8 | −0.662 | 0.000 |
| STAT1 | 12 | −0.604 | 0.000 |
| HLA-E | 14 | −0.598 | 0.000 |
| JAK2 | 21 | −0.527 | 0.200 |

**Pathway enrichment — the strongest result.** Taking my own top 100 depleted genes, with nothing indicating what to look for, **interferon gamma signalling was the top enriched term of everything tested** (FDR 7.3e-05), followed by antigen presentation (4.0e-04). Seven of the top fourteen terms are the paper's two themes.

The paper's own IPA analysis was run on a hand-selected gene list. This one was not, and the same pathways emerged.

**TCGA survival**

| Comparison | Mine | Paper | My n | Paper n |
|---|---|---|---|---|
| LAML, NKG7-low | P = 0.2866 | 0.2680 | 22 vs 22 | 22 vs 22 |
| LAML, NKG7-high | P = 0.1788 | 0.2117 | 23 vs 22 | 22 vs 23 |
| KIRC | **P = 0.0001** | 0.0001 | 204 vs 203 | 200 vs 201 |

Arriving at 22 vs 22 independently — from percentile cutoffs applied to a stratified subset — indicates the filtering matched step for step.

---

## What did not reproduce

Reported rather than omitted.

| Gene | Paper | Mine (low pressure) |
|---|---|---|
| PSMB5 | Highlighted in Fig 2A | rank 17,580 (RRA) / 16,055 (MLE) |
| ABL1 | Highlighted in Fig 2A | ~4,500 both methods |
| IFNGR1 | In the interferon network | 1,510 (RRA) / 3,885 (MLE) |
| HLA-A | Highlighted | 14,477 low — but **rank 68 under high pressure** |

---

## Findings not in the paper

**The empirical noise floor.** The GeCKO v2 library contains 1,000 non-targeting control guides, which cut nothing. Their spread measures the noise directly:

| Selection arm | Control beta SD |
|---|---|
| Low pressure | **0.210** |
| High pressure | **1.189** |

Noise is **5.7× wider** under high selection. This explains an apparent paradox: high-pressure effect sizes are 2–3× larger, yet only one gene reaches FDR < 0.05 there versus seven in the low arm. Normalised to the noise, IFNGR2 falls from 4.4 SD-units to 1.8. **The high-pressure arm has larger effects and less power.**

**A false positive caught by inspecting raw counts.** GPSM1 reached FDR 0.000 but carried a z-score of only 1.70. Its guide-level data showed one guide jumping 6× in one replicate while the other five moved ~0.3. Excluded.

**A candidate the paper does not mention.** ZNF474 shows **6/6 guides depleted in both replicates** — more internally consistent than IFNGR2 (5/6) — and ranks 3rd by Wald p-value. It carries a permutation FDR of 0.41 and comes from a single screen with two replicates, so it is an observation, not a finding.

**A discrepancy between MAGeCK's two p-value columns.** ZNF474 has a Wald p 14× smaller than HLA-E but a far worse FDR, because the `fdr` column derives from the permutation test while `wald-p-value` is parametric. The two rank genes differently and the choice has to be deliberate.

---

## Method

```
7 FASTQ (SRA)
   → mageck count          119,461 guides, 21,697 genes, ~53% mapped
   → QC                    Gini, zero counts, replicate correlation
   → label verification    NCR3LG1 CPM check before any statistics
   → mageck test (RRA)     both arms, --paired for the low arm
   → mageck mle            both arms — matches the paper's method
   → STRING                networks + enrichment
   → xenaPython/lifelines  TCGA survival
```

### Deviations from the paper, and why

| Aspect | Paper | Here | Reason |
|---|---|---|---|
| Tool version | MAGeCK-VISPR 0.5.4 | MAGeCK 0.5.9.5 | Current bioconda build |
| Statistic | MLE | RRA **and** MLE | RRA first, MLE added to match |
| Trim | not reported | `--trim-5 auto` | Staggered library measured at offsets 28–37 |
| Pairing | not reported | `--paired`, low arm | Controls correlate better with their own survivors (0.843) than with each other (0.765) |
| Pathway tool | Ingenuity IPA | STRING | IPA is commercial |
| Fig 2A x-axis | ~1,600 genes | all 21,697 | They applied an unstated threshold |

Recovering the paper's top hits at rank 1 using a *different* statistic is arguably stronger evidence than matching exactly — it shows the result is not an artifact of one method.

---

## Data notes

**Library file.** GeCKO v2 (Addgene #1000000049) ships as two CSVs with three problems, none of which raise errors: classic Mac `\r` line endings, columns in the order `gene, sgRNA_id, sequence` where MAGeCK requires `sgRNA_id, sequence, gene`, and a trailing empty column. Fed in unmodified, MAGeCK produces a wrong answer rather than a failure. Merged and reformatted to 123,411 guides — matching the published figure exactly.

**Guide count.** 123,411 in, 119,461 out. MAGeCK collapses duplicate sequences; GeCKO v2 contains genuinely identical guides where genes are near-identical (PCDHA/PCDHG clusters, POTEB, PPIAL4, RGPD).

**Mapping rate.** The MAGeCK log reports ~69% and ~53%. These are different measurements: 69% is the trim test on the first 100,001 reads, 53% is the rate across all reads. 53% is below the usual 65% guideline and probably means `--trim-5 auto` did not capture all ten stagger offsets. 13–27M mapped reads per sample remains ample for a 120k library.

---

## Reproducing this

```bash
conda env create -f environment.yml
conda activate crispr

bash   scripts/01_download.sh          # ~10 GB from SRA
python scripts/02_prepare_library.py
bash   scripts/03_count.sh
bash   scripts/04_test.sh
bash   scripts/05_mle.sh
python scripts/06_figures.py
python scripts/07_tcga_survival.py
python scripts/08_paper_figures.py
python scripts/09_string_networks.py
python scripts/10_pathway_enrichment.py
```

FASTQ and count files are **not** committed. `data/README.md` has the accessions and download instructions.

---

## Repository

```
├── scripts/             numbered, runnable
├── figures/             10 figures
├── results/             gene summaries, enrichment tables
└── data/README.md       accessions and download instructions
```

---

## Software

MAGeCK 0.5.9.5 · Python 3.11 · pandas · scipy · statsmodels · matplotlib · lifelines · xenaPython · STRING REST API

## References

Zhuang X, Veltri DP, Long EO (2019) *Front Immunol* 10:2879. doi:10.3389/fimmu.2019.02879

Li W et al. (2014) MAGeCK. *Genome Biol* 15:554. doi:10.1186/s13059-014-0554-8

Li W et al. (2015) MAGeCK-VISPR. *Genome Biol* 16:281. doi:10.1186/s13059-015-0843-6

Szklarczyk D et al. STRING v11. *Nucleic Acids Res*

Goldman MJ et al. (2020) UCSC Xena. *Nat Biotechnol* 38:675–678
