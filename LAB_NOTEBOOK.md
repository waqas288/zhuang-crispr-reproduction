# Lab notebook — Zhuang et al. 2019 CRISPR screen reproduction

**Purpose:** portfolio artifact demonstrating end-to-end pooled CRISPR screen analysis, from raw FASTQ to hit lists, benchmarked against a published paper's own results.

**Rule for this document:** every step gets *what was done* and *why*. If a step has no reason, it should not be in the repo.

---

## 0. Project definition

### 0.1 The paper

**Zhuang X, Veltri DP, Long EO.** Genome-Wide CRISPR Screen Reveals Cancer Cell Resistance to NK Cells Induced by NK-Derived IFN-γ. *Frontiers in Immunology* 2019;10:2879.
**DOI:** 10.3389/fimmu.2019.02879 — open access (CC BY)
**Data:** GEO **GSE139313** · BioProject PRJNA579227 · SRA **SRP226783**

### 0.2 Why this paper and not another

Candidates were checked against four criteria. Most failed on the first.

| Criterion | Why it matters | Zhuang |
|---|---|---|
| Raw FASTQ actually deposited | Without reads you can't demonstrate `mageck count`, which is the step a real lab needs | ✅ SRA SRP226783 |
| Open access full text | Methods must be readable to reproduce them | ✅ CC BY |
| Same toolchain | Reproducing with the same software isolates *your* execution as the variable | ✅ MAGeCK-VISPR v0.5.4 |
| Specific checkable findings | Vague conclusions can't be scored as reproduced or not | ✅ see §0.3 |
| Laptop-tractable | No cluster available | ✅ ~10 GB compressed |

**Rejected candidates and why — this list is itself part of the record:**

- **Steinhart et al. 2016** (Wnt/FZD5 pancreatic, Moffat & Angers). Strong findings, ideal author list. **No raw reads locatable** — CRISPRAnalyzeR repackaged only read counts, and no GEO/SRA accession surfaced. Would have meant starting mid-pipeline.
- **Yin et al. 2024** (BCL2L1 + radiotherapy, *Life Science Alliance*). Recent, open access, MAGeCK+DrugZ. **Data availability deposits only RNA-seq (GSE236331); CRISPR screen raw reads not deposited.** Only processed DrugZ output in Table S1.
- **CLD demo data** (Heigwer et al. 2016). Has FASTQ, but it's a *software* paper — the screen validates a guide-design tool, so reproduction targets are thin.
- **2025–26 papers** (Bradu, Zhu, Burrell). Mostly Perturb-seq and CRISPRa — different pipeline, single-cell scale, not laptop-tractable.

> **Lesson to carry forward:** data availability is the binding constraint, not scientific interest. Check the Data Availability statement *before* reading the paper properly.

### 0.3 Reproduction targets

From the paper's own results. These are scored pass/fail at the end.

| Direction | Genes | Meaning |
|---|---|---|
| **Enriched** (resistance) | **NCR3LG1** | Encodes B7H6, ligand for NK activating receptor NKp30. Losing it hides the cell from NK killing |
| **Depleted** (vulnerability) | **B2M, HLA-E, HLA-A, HLA-C, TAP1, TAPBP, PSMB5** | Antigen presentation pathway |
| **Depleted** (vulnerability) | **IFNGR1, IFNGR2, JAK1, JAK2, STAT1, IRF1, IFITM1** | IFNGR-JAK-STAT pathway |
| **Depleted** | **ABL1** | K562 is BCR-ABL positive (Philadelphia chromosome) |

**Internal control built into the design:** genes appearing in *both* low and high selection pressure are highlighted in red in the paper's Figure 2A. Two independent selection conditions agreeing is stronger evidence than either alone — and it gives us a consistency check that doesn't depend on matching their exact numbers.

### 0.4 What this is NOT

- Not novel research. ARID1A-style novelty checks do not apply — reproduction is the stated claim.
- Not a paper. Do not attempt to publish.
- Its job is to be a link in supervisor emails for **Canada, Australia, China** (Italy and Erasmus Mundus don't need it).

---

## 1. The experiment being reproduced (wet lab, from the paper)

Not done by us — but must be understood to defend any analysis of it.

| Parameter | Value | Why it matters |
|---|---|---|
| Cell line | K562 (human CML), Cas9-expressing, blasticidin-selected | Standard NK target line; low MHC-I makes it NK-sensitive |
| Library | GeCKO v2 human, Addgene **#1000000049** (2-plasmid) | 6 sgRNA/gene; 2-plasmid format needs pre-existing Cas9 |
| MOI | **0.3** | Low on purpose → ~1 guide per cell, so effects are attributable |
| Selection | Puromycin, 7 days | Removes uninfected cells |
| Cells | 50 × 10⁶ transduced | Library representation (~400×) |
| Effector:Target | 0.3:1 IL-2-activated primary NK cells | The selection pressure |
| Low pressure | 1 round, until only 10% of K562 survive | Milder — catches subtler effects |
| High pressure | 3 rounds of co-incubation | Stronger — enriches the strongest resistors |
| Control | K562 cultured in parallel, no NK exposure | **Critical**: controls for growth in culture, so hits reflect NK selection not fitness |
| Replicates | Low: 2 biological · High: 2 technical | Biological > technical for generalisability |
| Sequencing | 75-bp single-end, NextSeq 500 | Our reads measure 80 bp — extra bases are the stagger |
| Analysis | MAGeCK-VISPR v0.5.4 | Same tool we're using |

**Screen type:** this is a **survival/selection screen**, not a dropout screen. That means **both directions are informative** — `pos|` columns for resistance, `neg|` for vulnerability. DepMap-style dropout screens only read one direction.

---

## 2. Environment setup

### 2.1 Why WSL

MAGeCK is distributed via bioconda, which has no Windows build. `pip install mageck` fails — it isn't on PyPI. Nextflow and Docker also need Linux. One WSL install unblocks the entire field.

### 2.2 Why a separate conda environment

Base was pinned to Python 3.13; MAGeCK's newest build caps at 3.12. Rather than downgrade base and risk breaking other work, isolate.

```bash
conda create -n crispr -c bioconda -c conda-forge python=3.11 mageck -y
conda activate crispr
conda install -c conda-forge matplotlib-base pandas numpy scipy statsmodels -y
conda install -c conda-forge openjdk=21 -y      # Nextflow needs Java 11+
conda install -c bioconda -c conda-forge sra-tools -y
```

**Gotchas hit, recorded so they aren't hit again:**

- **`matplotlib` not `matplotlib-base`** pulls the full Qt GUI stack (~70 MB, `qt6-main` alone is 58 MB). It timed out on download and is useless in a headless WSL session. `matplotlib-base` does everything except open interactive windows.
- **Java is per-environment.** The Java visible in `(base)` was JetBrains' bundled runtime and vanished on `conda activate crispr`. Install openjdk into the env.
- **conda environments don't persist across terminal sessions.** Every new WSL window starts in `(base)`. `conda activate crispr` is the first command, every time.
- **bioconda is a channel, not a program.** Without `-c bioconda`, `conda install mageck` returns "not found" — same failure mode as pip.

### 2.3 Docker

Installed via **Docker Desktop on Windows**, with Settings → Resources → **WSL Integration** enabled for the Ubuntu distro. Verified with `docker run hello-world`.

Docker is not a second Linux to work inside. Nextflow asks it to pull a small prebuilt image per pipeline step, run it, discard it. Avoids installing 30 tools with conflicting dependencies, and makes the analysis reproducible on any machine.

---

## 3. Data acquisition

### 3.1 The seven samples

| SRA run | Label | Condition |
|---|---|---|
| SRR10338906 | Low_Ctrl_1 | Low pressure, control, bio rep 1 |
| SRR10338907 | Low_Surv_1 | Low pressure, survivors, bio rep 1 |
| SRR10338908 | Low_Ctrl_2 | Low pressure, control, bio rep 2 |
| SRR10338909 | Low_Surv_2 | Low pressure, survivors, bio rep 2 |
| SRR10338910 | High_Ctrl | High pressure, control |
| SRR10338911 | High_1 | High pressure, survivors, tech rep 1 |
| SRR10338912 | High_2 | High pressure, survivors, tech rep 2 |

⚠️ **Order assumed from the GEO GSM listing. Must be confirmed against the SRA Run Table before trusting any result.** Mislabelling controls and survivors would invert every conclusion and produce a plausible-looking wrong answer.

### 3.2 Download

```bash
mkdir -p ~/zhuang && cd ~/zhuang
for SRR in SRR10338906 SRR10338907 SRR10338908 SRR10338909 \
           SRR10338910 SRR10338911 SRR10338912; do
  echo ">>> $SRR $(date)"
  fastq-dump --gzip $SRR
done 2>&1 | tee download.log
```

**Why `fastq-dump --gzip` and not `fasterq-dump`:** `fasterq-dump` is faster but writes *uncompressed* FASTQ — ~6 GB per sample, ~45 GB total. Disk was the binding constraint. `--gzip` lands at ~1.5 GB per sample, ~10 GB total. MAGeCK reads `.gz` directly, so nothing is lost but time.

**Why `~/zhuang` and not `/mnt/c/...`:** file I/O across the WSL↔Windows boundary is markedly slower. Work in WSL's native filesystem; copy small result files back at the end.

**Result:** 41–53 million reads per sample. Sizes 775 MB – 1.8 GB.

> **Storage note:** WSL reports a 1 TB virtual disk regardless of the physical drive. The real constraint is free space on the Windows side. `df -h ~` will lie to you. Check `Get-PSDrive C` in PowerShell. Also: the WSL virtual disk does **not** shrink when you delete files inside it.

### 3.3 Read structure — the stagger

First read:

```
@SRR10338906.1 1 length=80
ATCGANAGTTGCTTGCTTTATAAATCTTNTGGAAAGGACGAAACACCTNAAGCAGTTCCAACTGTTACGGTTTTAGAGCN
```

`AAAGGACGAAAC` is the end of the vector sequence immediately preceding the guide. In two reads it appeared at slightly different positions — but two reads is not evidence. Tested properly across 100,000 reads:

```bash
zcat SRR10338906.fastq.gz | awk 'NR%4==2' | head -100000 | \
  awk '{n=index($0,"AAAGGACGAAAC"); if(n>0) print n}' | sort -n | uniq -c
```

**Result:**

```
  12749  28      6328  34
  11476  29      6296  35
  11423  30      6841  36
  11648  31      6187  37
  10828  32
  12413  33
```

Ten distinct offsets, positions 28–37, ~98% of reads. A fixed-position library would show one spike.

**Why staggers exist:** Illumina estimates base quality by comparing across all clusters at each cycle. If every read carries the same base at position 1, 2, 3…, base diversity is too low to calibrate and quality degrades. A variable-length filler before the constant region offsets reads from each other and restores diversity.

**Consequence:** a fixed `--trim-5 N` cannot work. Use `--trim-5 auto`.

**Also note:** the `N`s scattered through the reads are low-quality base calls. Normal, not a problem.

---

## 4. The library file

### 4.1 Source

Addgene GeCKO v2 human, **#1000000049** (two-plasmid). Sequences distributed as two CSVs — half-libraries A and B, 3 guides per gene each, 6 combined. The paper used the full library, so both are needed.

*(The plasmid itself is a paid physical product. The sequence CSVs are free — a distinction worth being clear about.)*

### 4.2 Three problems, all silent

**Problem 1 — line endings.** Files use `\r` only (classic Mac). Linux tools read the whole file as one record; `head` dumps everything onto one line.

**Problem 2 — column order.** Files are `gene_id, UID, seq`. MAGeCK requires **`sgRNA_id, sequence, gene`**. Feeding it unswapped produces garbage, not an error.

**Problem 3 — trailing empty column** on every row.

> **This is the most dangerous class of bioinformatics bug: input that is wrong but parseable.** Nothing raises an exception. You get an answer, and it's wrong. Same category as the DepMap `IsDefaultEntryForModel == True` filter that silently returned zero rows.

### 4.3 Fix

Normalise line endings → strip empties → reorder to `UID, seq, gene_id` → concatenate A and B → keep only 20 bp ACGT guides.

**Result: 123,411 guides.** Matches the published GeCKO v2 figure exactly — that agreement is the validation that the merge was correct.

| Property | Value |
|---|---|
| Total guides | 123,411 |
| Gene labels | 21,915 |
| Genes with 6 guides | 18,940 |
| miRNA entries (`hsa-mir-*`) | 1,853 |
| **Non-targeting controls** | **1,000** (`NonTargetingControlGuideForHuman_0001`–`_1000`) |

**Why the controls matter:** these guides cut nothing. They define what "no effect" looks like *in this dataset*, rather than in theory. MAGeCK can normalise against them via `--control-sgrna`, and they serve as the internal negative control. Their spread is also a direct read on technical noise.

Output: `GeCKOv2_library_mageck.csv`

---

## 5. Counting — DONE

```bash
cd ~/zhuang
conda activate crispr
mageck count -l GeCKOv2_library_mageck.csv -n zhuang \
  --sample-label "Low_Ctrl_1,Low_Surv_1,Low_Ctrl_2,Low_Surv_2,High_Ctrl,High_1,High_2" \
  --fastq SRR10338906.fastq.gz SRR10338907.fastq.gz SRR10338908.fastq.gz \
          SRR10338909.fastq.gz SRR10338910.fastq.gz SRR10338911.fastq.gz \
          SRR10338912.fastq.gz \
  --trim-5 auto
```

**What this does:** for each read, find the constant vector sequence, take the following 20 bases, look them up in the library, tally. Output is a counts table — guides as rows, samples as columns.

**Error encountered:** `FileNotFoundError: GeCKOv2_library_mageck.csv` — the file was downloaded to Windows, not copied into WSL. Fix:

```bash
cp /mnt/c/Users/WAQAS/Downloads/GeCKOv2_library_mageck.csv ~/zhuang/
wc -l ~/zhuang/GeCKOv2_library_mageck.csv    # expect 123411
```

### 5.1 Results

```
119,461 guides · 21,697 genes
```

Library entered as 123,411 guides, came out as **119,461** — a loss of 3,950. **Not an error:** MAGeCK collapses identical sequences. GeCKO v2 contains genuinely duplicated guides where genes are near-identical (the PCDHA/PCDHG protocadherin clusters, POTEB/POTEB2, PPIAL4 family, RGPD family). Verifiable in the library file itself.

### 5.2 Mapping rates

| Sample | Total reads | Mapped | % |
|---|---|---|---|
| Low_Ctrl_1 | 49,933,559 | 26,522,030 | 53.1 |
| Low_Surv_1 | 52,864,160 | 27,641,345 | 52.3 |
| Low_Ctrl_2 | 41,279,851 | 23,580,548 | 57.1 |
| Low_Surv_2 | 43,850,660 | 25,002,636 | 57.0 |
| High_Ctrl | 30,633,514 | 15,696,597 | 51.2 |
| High_1 | 26,389,723 | 13,637,139 | 51.7 |
| High_2 | 24,897,799 | 13,013,857 | 52.3 |

⚠️ **Two different numbers appear in the log and they are not the same measurement.** The log first reports ~69% — that is MAGeCK testing trim lengths on the first 100,001 reads only. The countsummary reports ~53%, which is the real rate across all reads.

**53% is below the ~65% rule of thumb.** Probable cause: `--trim-5 auto` did not capture all ten stagger offsets measured in §3.3, so reads at the uncaptured offsets failed to match.

**Not fatal** — 13–27M mapped reads per sample is ample for a 120k library. But it is a real loss and must be stated in the README rather than hidden.

**Open improvement:** supply the offsets explicitly instead of `auto` and check whether mapping recovers. Would be a concrete, documentable enhancement.

### 5.3 Library representation QC

| Sample | Zero-count guides | Gini index |
|---|---|---|
| Low_Ctrl_1 | 487 | 0.090 |
| Low_Surv_1 | 385 | 0.088 |
| Low_Ctrl_2 | 536 | 0.092 |
| Low_Surv_2 | 1,001 | 0.098 |
| **High_Ctrl** | **11,620** | **0.226** |
| **High_1** | **26,238** | **0.386** |
| **High_2** | **23,427** | **0.366** |

Low-pressure arm: excellent. Gini ~0.09, under 1% of guides lost.

High-pressure arm: High_1 lost **22% of the library entirely**.

**Initial interpretation (WRONG, recorded deliberately):** attributed mainly to random bottlenecking from three rounds of 90% killing.

**Corrected interpretation:** the §5.5 CPM check shows NCR3LG1 at **1,683 CPM in High_1 against 7.7 in High_Ctrl** — roughly 200-fold. That is not drift, it is ferocious selection. Most guides disappeared *because those cells were killed*, which is the experiment working exactly as designed. A Gini of 0.39 in a three-round selection arm is expected, not a defect.

**Consequence for interpretation:** the high-pressure arm has lower statistical power but larger effect sizes. Weight the low-pressure arm for significance; use the high arm for confirmation.

### 5.4 Correlation structure — changed the analysis

log10(count+1) Pearson:

|  | Low_Ctrl_1 | Low_Surv_1 | Low_Ctrl_2 | Low_Surv_2 |
|---|---|---|---|---|
| **Low_Ctrl_1** | 1.000 | **0.843** | 0.765 | 0.738 |
| **Low_Surv_1** | 0.843 | 1.000 | 0.777 | 0.746 |
| **Low_Ctrl_2** | 0.765 | 0.777 | 1.000 | 0.776 |
| **Low_Surv_2** | 0.738 | 0.746 | 0.776 | 1.000 |

High_1 vs High_2 = **0.760** (low for *technical* replicates — same bottleneck cause).

> **A control correlates better with its own treated sample (0.843) than with the other control (0.765).** Replicate identity dominates the treatment effect.

**Action taken:** ran `mageck test` with **`--paired`**, so Surv_1 is compared against Ctrl_1 and Surv_2 against Ctrl_2, rather than pooling all controls against all survivors. Pooling would have let batch differences leak into the treatment contrast.

**This is a decision driven by observed data, not a default.** Exactly the kind of choice to be able to defend.

### 5.5 Label verification — the check that had to pass first

§3.1 flagged the sample-to-condition mapping as **assumed** from the GEO listing. If controls and survivors were swapped, every downstream conclusion would invert and still look plausible. Tested directly using the paper's known biology, before any statistics:

**Mean CPM per gene per sample**

| Gene | Low_Ctrl_1 | Low_Surv_1 | Low_Ctrl_2 | Low_Surv_2 | High_Ctrl | High_1 | High_2 |
|---|---|---|---|---|---|---|---|
| **NCR3LG1** ↑ | 9.3 | **25.5** | 8.7 | **18.3** | 7.7 | **1682.9** | **1368.5** |
| IFNGR2 ↓ | 9.1 | 3.2 | 6.3 | 2.7 | 7.5 | 0.1 | 0.5 |
| B2M ↓ | 10.3 | 5.6 | 6.9 | 4.4 | 8.4 | 0.4 | 0.3 |
| JAK1 ↓ | 8.9 | 3.2 | 8.5 | 4.2 | 5.3 | 1.1 | 0.0 |
| STAT1 ↓ | 6.1 | 2.5 | 5.7 | 3.9 | 4.6 | 1.0 | 0.5 |

**PASS.** All five move in the direction the paper reports, in every comparison, in both arms. Labels are correct.

**Note this was visible in raw normalised counts before any statistical test.** When an effect is real and strong, it shows up before the statistics — the statistics quantify confidence, they don't manufacture the signal.

---

## 6. Hit calling — RRA

```bash
mageck test -k zhuang.count.txt \
  -t Low_Surv_1,Low_Surv_2 -c Low_Ctrl_1,Low_Ctrl_2 -n low_rra --paired

mageck test -k zhuang.count.txt \
  -t High_1,High_2 -c High_Ctrl -n high_rra
```

High arm cannot be paired — one control, two survivors.

### 6.1 Low pressure

**Top depleted (vulnerability):**

| Rank | Gene | FDR | goodsgrna | LFC |
|---|---|---|---|---|
| 1 | **JAK1** | 0.0025 | 11 | −1.60 |
| 2 | **IFNGR2** | 0.0025 | 11 | −1.46 |
| 3 | SEC62 | 0.43 | 9 | −0.81 |
| 4 | **IRF1** | 0.43 | 10 | −0.43 |
| 5 | **STAT1** | 0.43 | 10 | −0.68 |
| 8 | **HLA-E** | 0.43 | 10 | −1.13 |
| 14 | **B2M** | 0.43 | 9 | −1.55 |
| 15 | **ERAP1** | 0.43 | 8 | −0.62 |

**Top enriched (resistance):**

| Rank | Gene | FDR | goodsgrna | LFC |
|---|---|---|---|---|
| 1 | **NCR3LG1** | 0.0050 | 10 | **+1.20** |
| 2 | MRPL41 | 0.012 | 5 | +0.18 |

### 6.2 High pressure

| Rank | Gene | FDR | LFC |
|---|---|---|---|
| **pos 1** | **NCR3LG1** | 0.0025 | **+7.81** |
| neg 10 | **B2M** | 0.45 | −3.07 |
| neg 20 | **JAK2** | 0.50 | −4.39 |
| neg 38 | **STAT1** | 0.69 | −4.06 |
| neg 50 | **TAP1** | 0.72 | −3.92 |
| neg 60 | **IFNGR2** | 0.72 | −4.75 |
| neg 68 | **HLA-A** | 0.72 | −3.78 |
| neg 97 | **HLA-E** | 0.72 | −3.35 |

⚠️ Top depleted gene in the high arm is `hsa-mir-6859-2`, LFC −8.77, **`goodsgrna` = 1**. One guide. That is an artifact, not biology — and it is precisely what the `goodsgrna` column exists to expose. Must be flagged in the write-up, not silently reported.

---

## 7. Hit calling — MLE (matches the paper's method)

The paper used MAGeCK-VISPR and reported normalized log fold change per gene — that is MLE output. RRA was our substitution; MLE closes the gap.

```bash
printf "Samples\tbaseline\tNKselection\nLow_Ctrl_1\t1\t0\nLow_Surv_1\t1\t1\nLow_Ctrl_2\t1\t0\nLow_Surv_2\t1\t1\n" > design_low.txt
mageck mle -k zhuang.count.txt -d design_low.txt -n low_mle
```

**Design matrix:** `baseline` = 1 everywhere (intercept). `NKselection` = 1 only for survivors. The beta for `NKselection` is the effect of interest — negative = depleted, positive = enriched.

**Why MLE has more power than RRA:** RRA treats the samples as lists to rank. MLE fits one model across all four, so replicate variance is estimated properly instead of collapsed.

Runtime ~14 minutes (RRA was ~1 minute).

### 7.1 Results — six genes at FDR 0.000

| Gene | beta | FDR | rank (of 21,697) |
|---|---|---|---|
| **JAK1** | −0.992 | **0.000** | 1 |
| **IFNGR2** | −0.933 | **0.000** | 2 |
| **B2M** | −0.662 | **0.000** | 8 |
| **STAT1** | −0.604 | **0.000** | 12 |
| **HLA-E** | −0.598 | **0.000** | 14 |
| JAK2 | −0.527 | 0.200 | 21 |
| IRF1 | −0.418 | 0.455 | 71 |
| TAPBP | −0.400 | 0.569 | 91 |
| HLA-C | −0.338 | 0.772 | 189 |
| ERAP1 | −0.330 | 0.792 | 216 |
| TAP1 | −0.186 | 0.950 | 1596 |
| IFNGR1 | −0.112 | 0.992 | 3885 |
| ABL1 | −0.108 | 0.992 | 4054 |
| HLA-A | −0.088 | 0.992 | 5024 |
| PSMB5 | +0.083 | 0.992 | 16055 |
| **NCR3LG1** | **+0.819** | **0.000** | **21,696** (last = most positive) |

Under RRA only 3 genes cleared FDR 0.05, most sitting at 0.43. Under MLE, **six hit 0.000**. That is the power gain, demonstrated rather than asserted.

### 7.2 The non-targeting controls — the most instructive result in the project

Top 15 negative betas contained **7 `NonTargetingControlGuideForHuman_*` entries**. Top 15 positive betas contained **6**.

These guides cut nothing. Their true effect is zero. Yet several carry larger |beta| than B2M or STAT1:

| Entry | beta | FDR |
|---|---|---|
| NonTargeting_0849 | **+0.845** | 0.863 |
| NCR3LG1 | +0.819 | **0.000** |
| NonTargeting_0576 | −0.811 | 0.871 |
| B2M | −0.662 | **0.000** |

> **This is the empirical noise floor, measured in this dataset rather than assumed.** Any gene whose beta falls inside the range spanned by the controls cannot be distinguished from noise on effect size alone.

**And it demonstrates why FDR is reported rather than fold change.** MAGeCK separates them correctly — controls carry large betas but FDR 0.86–0.94, while B2M at −0.66 gets FDR 0.000. The model accounts for guide-level variance; the raw magnitude does not.

miRNA entries cluster in the same region — consistent, since most have few guides and low counts.

**Open improvement:** re-run with `--control-sgrna` pointing at the 1,000 non-targeting guides, so normalisation uses the measured null rather than the global median. Compare against the current run.

### 7.3 RRA vs MLE disagreement

| Gene | RRA rank | MLE rank | Note |
|---|---|---|---|
| JAK1 | 1 | 1 | agree |
| IFNGR2 | 2 | 2 | agree |
| SEC62 | 3 | not top | RRA only |
| IRF1 | 4 | 71 | RRA favours |
| B2M | 14 | 8 | MLE favours |
| JAK2 | 106 | 21 | MLE favours |

**Why they differ:** RRA rewards *consistency of direction across guides*. MLE rewards *effect size with variance accounted for*. Neither is wrong. Reporting both, and explaining the difference, is stronger than picking one and hiding the other.

---

## 8. Reproduction scorecard

### Reproduced

| Paper claim | Evidence |
|---|---|
| **NCR3LG1 loss → NK resistance** | rank 1 in *both* arms, RRA *and* MLE. FDR 0.000. LFC +7.8 high arm |
| **IFN-γ pathway loss → vulnerability** | JAK1 rank 1, IFNGR2 rank 2, STAT1 rank 12, JAK2 rank 21 — all FDR 0.000 except JAK2 (0.20) |
| **Antigen presentation loss → vulnerability** | B2M rank 8 (FDR 0.000), HLA-E rank 14 (FDR 0.000), HLA-C 189, TAPBP 91 |
| **Overlap between selection pressures** | B2M, IFNGR2, JAK1, JAK2, STAT1, HLA-E, HLA-C, TAPBP all appear in both arms — the paper's own internal control |

### Not reproduced — report, do not hide

| Gene | Result | Note |
|---|---|---|
| **PSMB5** | RRA 17,580 low / 6,348 high · MLE 16,055 | Highlighted in the paper's Fig 2A. Clear divergence |
| **ABL1** | ~4,000–4,500 both methods | Paper reports it sensitises K562 (BCR-ABL positive) |
| **IFNGR1** | 1,510 RRA / 3,885 MLE | Its partner IFNGR2 is rank 2 — asymmetry worth commenting on |
| **HLA-A** | 14,477 low / **68 high** | Only surfaces under strong selection |

### Found independently, not on the target list

**ERAP1** — RRA rank 15 low, MLE rank 216. Trims peptides for MHC class I loading, so it belongs with TAP1 and TAPBP in the antigen presentation pathway. The pipeline recovered a pathway member that was not searched for.

---

## 9. Methodological deviations from the paper — state these explicitly

| Aspect | Paper | This reproduction | Why |
|---|---|---|---|
| Tool version | MAGeCK-VISPR v0.5.4 | MAGeCK v0.5.9.5 | Current bioconda build |
| Statistic | MLE (normalized LFC) | RRA **and** MLE | RRA run first; MLE added to match |
| Trim | not reported | `--trim-5 auto` | Stagger measured at offsets 28–37; not stated in the paper, so this is our parameter |
| Pairing | not reported | `--paired`, low arm | Driven by observed correlation structure (§5.4) |
| Mapping rate | not reported | ~53% | Our number; disclose it |

> **Recovering the paper's top hits at rank 1 using a *different* statistic (RRA) is arguably stronger evidence than matching exactly — it shows the result is not an artifact of one method.**

---

## 10. Control-sgRNA normalisation

```bash
grep "NonTargetingControlGuideForHuman" zhuang.count.txt | cut -f1 > nt_controls.txt   # 1000 guides
mageck test -k zhuang.count.txt -t Low_Surv_1,Low_Surv_2 -c Low_Ctrl_1,Low_Ctrl_2 \
  -n low_rra_ctrlnorm --paired --control-sgrna nt_controls.txt --norm-method control
```

**Rationale:** default normalisation scales samples to a common median across all guides. That assumes most guides are unaffected — reasonable, but it is an assumption. Normalising to the 1,000 non-targeting guides uses the *measured* null instead.

**Result: the paper's genes barely move.**

| Gene | RRA rank | RRA-ctrl rank | RRA FDR | RRA-ctrl FDR |
|---|---|---|---|---|
| IFNGR2 | 2 | 1 | 0.0025 | 0.0025 |
| JAK1 | 1 | 2 | 0.0025 | 0.0025 |
| STAT1 | 5 | 7 | 0.432 | 0.279 |
| HLA-E | 8 | 10 | 0.432 | 0.279 |
| B2M | 14 | 17 | 0.432 | 0.279 |
| JAK2 | 106 | 83 | 0.791 | 0.583 |
| TAP1 | 1057 | 909 | 0.974 | 0.723 |

Genes at FDR<0.05: **2 neg / 2 pos** (default) → **2 neg / 4 pos** (control-normalised). FDRs improve modestly; ranks shift by 1–3 places.

**Interpretation:** the result is robust to normalisation choice. That is a cheap, genuine robustness check.

---

## 11. High-pressure MLE

```bash
printf "Samples\tbaseline\tNKselection\nHigh_Ctrl\t1\t0\nHigh_1\t1\t1\nHigh_2\t1\t1\n" > design_high.txt
mageck mle -k zhuang.count.txt -d design_high.txt -n high_mle
```

| Gene | MLE-low beta | MLE-low FDR | MLE-high beta | MLE-high FDR |
|---|---|---|---|---|
| **NCR3LG1** | +0.819 | **0.000** | **+5.774** | **0.000** |
| IFNGR2 | −0.933 | **0.000** | −2.133 | 0.378 |
| B2M | −0.662 | **0.000** | −2.163 | 0.375 |
| HLA-E | −0.598 | **0.000** | −2.306 | 0.375 |
| JAK1 | −0.992 | **0.000** | −1.357 | 0.681 |
| STAT1 | −0.604 | **0.000** | −1.124 | 0.787 |
| JAK2 | −0.527 | 0.200 | −1.553 | 0.595 |

**Genes at FDR<0.05: MLE-low = 7 · MLE-high = 1.**

Counter-intuitive at first glance: high-pressure betas are **2–3× larger**, yet almost nothing is significant. §12 explains why.

---

## 12. The empirical noise floor — the central result of this project

The 1,000 non-targeting control guides cut nothing. Their true effect is zero. So their beta distribution *is* the noise.

| Arm | n | mean | **SD** | min | max |
|---|---|---|---|---|---|
| Low pressure | 1000 | 0.016 | **0.210** | −0.811 | +0.845 |
| High pressure | 1000 | 0.252 | **1.189** | −3.718 | +4.773 |

> **Noise is 5.7× wider in the high-pressure arm.**

**This resolves the §11 paradox.** Effect sizes grew, but noise grew faster:

| Gene | Low beta | Low SD-units | High beta | High SD-units |
|---|---|---|---|---|
| IFNGR2 | −0.933 | **4.4** | −2.133 | **1.8** |
| B2M | −0.662 | **3.2** | −2.163 | **1.8** |
| NCR3LG1 | +0.819 | **3.9** | +5.774 | **4.9** |

Only NCR3LG1's signal grows faster than the noise — which is why it is the sole survivor at FDR<0.05 in the high arm.

**Three conclusions:**
1. **The high-pressure arm is less powerful, not more.** Three rounds of 90% killing bottlenecked the library (26,238 guides lost in High_1) and inflated variance beyond the gain in effect size.
2. **Report FDR, not fold change.** Several control guides carry larger |beta| than B2M (e.g. NonTargeting_0849 at +0.845 vs NCR3LG1 at +0.819) — but FDR 0.86 vs 0.000. The model separates them; raw magnitude does not.
3. The high-arm control mean drifts to +0.25 (from 0.016), a mild positive bias from the bottleneck. Worth noting.

**Figure:** `figures/fig2_noise_floor.png`

---

## 13. Guide-level validation — where the summary statistic misleads

MLE flagged **GPSM1** at FDR 0.000, which is not in the paper. But its z-score was only 1.70 (vs JAK1 −4.42). Checked the raw counts.

**GPSM1 — ARTIFACT**

| guide | Ctrl_1 | Ctrl_2 | Surv_1 | Surv_2 | LFC r1 | LFC r2 |
|---|---|---|---|---|---|---|
| HGLibA_20155 | 7.99 | 11.92 | 8.54 | **73.55** | 0.08 | **2.53** |
| others (×5) | ~6–9 | ~4–10 | ~7–11 | ~9–11 | 0.23–0.35 | −0.08–1.01 |

One guide, one replicate, a 6× jump. Five other guides move ~0.3. **Excluded.** The low z-score was the warning.

**ZNF474 — the cleanest signal in the dataset**

| | concordant guides | Wald p | FDR |
|---|---|---|---|
| **ZNF474** | **6/6 down, both replicates** | 0.00007 | **0.410** |
| IFNGR2 (paper hit) | 5/6 down | 0.00001 | 0.000 |
| NCR3LG1 (paper hit) | 5/6 up | — | 0.000 |
| GPSM1 (artifact) | 4/6, one outlier | — | 0.000 |

**ZNF474 is more internally consistent than a published finding**, ranks 3rd by Wald p, and is not mentioned in the paper. Uncharacterised zinc finger protein.

**Stated honestly:** one screen, two replicates. Could be off-target, could be a general fitness effect. Not a discovery — an observation with the evidence attached.

### The FDR / Wald discrepancy

| Gene | Wald p | FDR |
|---|---|---|
| ZNF474 | **0.00007** | 0.410 |
| HLA-E | 0.00099 | **0.000** |

A gene with a **14× smaller p-value** gets a **far worse FDR**. Not a bug: MAGeCK's `fdr` column derives from the **permutation** test, which accounts for guide count and consistency; `wald-p-value` is the parametric test on the beta estimate alone. The two rank genes differently and the choice must be deliberate.

---

## 14. Figures

| File | Content |
|---|---|
| `fig1_volcano_low_mle.png` | Volcano: beta vs −log10(Wald p). 1,000 controls in blue form the noise band; all 8 labelled genes sit above it |
| `fig2_noise_floor.png` | **Key figure.** Control beta distributions, low vs high, with SD annotated and IFNGR2/B2M/NCR3LG1 marked |
| `fig3_guide_level.png` | Per-guide CPM, controls vs survivors, grouped by arm. Shows guide-to-guide variation directly |
| `fig4_rra_vs_mle.png` | Rank-rank scatter. JAK1/IFNGR2 on the diagonal; IRF1/ERAP1 favoured by RRA; JAK2/TAPBP by MLE |
| `fig5_library_qc.png` | Count distributions + Gini. High-pressure samples spike at zero |
| `fig6_tcga_survival.png` | Kaplan-Meier, reproduction of paper Figure 5 |

**Figure-making errors made and corrected (recorded deliberately):**
- v1 volcano clipped FDR at 1e-6 → six genes stacked on one horizontal line, labels overlapping. Misleading.
- v2 plotted z-score against a p-value *derived from* z → produced a smooth V, not a scatter. Redundant axes.
- v3 (final) uses beta vs Wald p — two genuinely different quantities.
- fig3 v1 used `symlog`, producing a meaningless negative CPM axis, and interleaved Ctrl/Surv so the direction was unreadable.

---

## 15. TCGA survival reproduction (paper Figure 5)

**This is a separate analysis from the screen** — clinical outcome data, not CRISPR.

**Source:** UCSC Xena, hub `https://tcga.xenahubs.net`
- `TCGA.LAML.sampleMap/HiSeqV2` + `survival/LAML_survival.txt`
- `TCGA.KIRC.sampleMap/HiSeqV2` + `survival/KIRC_survival.txt`

Fetched via `xenaPython` API (queries specific genes rather than downloading full matrices — necessary given limited disk).

**Method, per the paper:** LAML split by NKG7 at 33rd/66th percentiles, then each stratum split by IFNGR2 at 40th/60th. KIRC split by IFNGR2 at 33rd/66th. Overall survival, log-rank (Mantel-Cox).

### Scorecard

| Comparison | Our P | Paper P | Our n | Paper n |
|---|---|---|---|---|
| LAML, NKG7-low | 0.2866 | 0.2680 | **22 vs 22** | **22 vs 22** |
| LAML, NKG7-high | 0.1788 | 0.2117 | **23 vs 22** | 22 vs 23 |
| **KIRC** | **0.0001** | **0.0001** | 204 vs 203 | 200 vs 201 |

**Reproduced.** The sample sizes are the strongest evidence — 22 and 22 arrived at independently, from percentile cutoffs on a stratified subset, means the filtering matched step for step.

**Why LAML p-values differ slightly:** with 44 patients, one sample crossing a percentile boundary moves the log-rank statistic materially. The 23/22 flip in the high stratum is one patient at the boundary, assigned differently by ties handling. Conclusion unchanged — non-significant either way.

**Version discrepancy, documented not chased:** 606 KIRC samples had both expression and survival here; the paper's groups total 401, implying ~601 available to them. Xena dataset versions have changed since 2019.

**Honest assessment of the finding itself:** this is the weakest analysis in the paper. Two non-significant results, and the KIRC result is a correlation between one gene's expression and survival across 400 patients — suggestive, not causal, and IFNGR2 is one of ~20,000 genes that could have been tested. Faithfully reproduced; that is separate from what it means.

---

## 16. STRING protein networks (paper Figures 2D, 2E)

`string_networks.py` — STRING v11+ REST API, free, no key required. The paper also used STRING, so this is a direct method match rather than a substitution.

### Antigen presentation (Fig 2D)
Genes: B2M, HLA-A, HLA-C, HLA-E, TAP1, TAPBP, PSMB5

**17 interactions · median confidence 0.989**

| A | B | score |
|---|---|---|
| HLA-E | B2M | 0.999 |
| TAPBP | B2M | 0.999 |
| HLA-A | B2M | 0.999 |
| HLA-C | B2M | 0.999 |
| HLA-A | TAPBP | 0.998 |
| TAP1 | TAPBP | 0.995 |

**B2M is the hub** — 0.999 edges to all four HLA/TAPBP partners. Matches the paper's Fig 2D topology, where B2M sits central.

### Interferon signalling (Fig 2E)
Genes: JAK1, JAK2, IFNGR1, IFNGR2, STAT1, IRF1, IFITM1, TAP1

**20 interactions · median confidence 0.955**

Ten edges at 0.999 form a dense IFNGR1–IFNGR2–JAK1–JAK2–STAT1 core. IRF1–STAT1 at 0.999. **IFITM1 attaches peripherally at 0.845** — matching the thinner line drawn in the paper's figure.

Output: `figures/fig7_string_antigen_presentation.png`, `figures/fig7_string_interferon_signaling.png`, plus `results_string_*.csv`

---

## 17. Pathway enrichment (paper Figure 2C) — the strongest result

`pathway_enrichment.py`

**Method substitution, documented:** the paper used **Ingenuity Pathway Analysis (IPA)**, which is commercial and unavailable. Substituted **STRING enrichment API** (GO Process, KEGG, Reactome).

Two analyses were run deliberately, because they answer different questions.

### A. Unbiased — the real reproduction

Input: **our own top 100 depleted genes** by MLE beta, excluding non-targeting controls and miRNAs. Nothing tells the analysis what to look for.

| Rank | Category | Term | FDR | genes |
|---|---|---|---|---|
| **1** | RCTM | **Interferon gamma signaling** | **7.3e-05** | 8 |
| 2 | RCTM | NoRC negatively regulates rRNA expression | 1.9e-04 | 7 |
| **3** | Process | **Interferon-gamma-mediated signaling pathway** | **2.6e-04** | 5 |
| **4** | RCTM | **Antigen Presentation: Folding, assembly and peptide loading of class I MHC** | **4.0e-04** | 5 |
| **5** | Process | **Antigen processing and presentation of endogenous peptide antigen** | **6.9e-04** | 5 |
| **6** | RCTM | **Regulation of IFNG signaling** | **7.8e-04** | 4 |
| **7** | Process | **Antigen processing and presentation via MHC class I** | **2.4e-03** | 5 |

> **Interferon gamma signalling is the single top term out of everything STRING tested. Seven of the top fourteen terms are the paper's two themes.**

**This is a genuine reproduction of the paper's biological conclusion**, obtained from raw FASTQ through independent counting, independent statistics, and an unbiased gene list.

**Arguably stronger than the paper's own version:** their IPA analysis was run on a hand-selected list. This one was not, and the same pathways emerged.

### B. Confirmatory — the paper's own gene sets

Running enrichment on the paper's hand-picked Fig 2D/2E lists returns the expected terms at FDR 1e-13 to 6e-13. **This is close to circular by construction** — confirming that a curated interferon gene list is enriched for interferon signalling. Included because it is what the paper's figure shows; flagged as circular in the script docstring.

### Honest observation the paper does not report

The unbiased top-100 list also returns unrelated terms at high rank:

| Term | FDR |
|---|---|
| NoRC negatively regulates rRNA expression | 1.9e-04 |
| Chromatin modifying enzymes | 5.0e-03 |
| DNA methylation | 5.1e-03 |
| Deposition of CENPA nucleosomes at centromere | 3.0e-03 |

These are chromatin and transcription complexes — almost certainly **general fitness genes**. Cells that grow slowly for any reason deplete under selection, regardless of NK biology.

**Consequence:** an unbiased top-100 list mixes fitness effects with the NK-specific signal. That is a real limitation of the unbiased approach, and it makes the paper's targeted selection defensible even though it is less rigorous. Worth stating in both directions.

---

## 18. Reproduction status — complete

All four of the paper's computational analyses have been reproduced.

| Paper analysis | Figure | Status | Evidence |
|---|---|---|---|
| Screen hit calling | 1D, 2A | ✅ | NCR3LG1 rank 1 both arms; JAK1 1, IFNGR2 2, B2M 8, STAT1 12, HLA-E 14 |
| Pathway analysis | 2C | ✅ (tool substituted) | Interferon gamma signalling = top term, unbiased, FDR 7.3e-05 |
| STRING networks | 2D, 2E | ✅ | 17 and 20 interactions, median confidence 0.99 and 0.96 |
| TCGA survival | 5 | ✅ | KIRC P=0.0001 exact; LAML n = 22 vs 22 exact |

**Not reproducible** (wet-lab): Figures 1B, 1C, 3, 4 and S1 — cytotoxicity assays, degranulation, flow cytometry, arrayed gRNA validation.

**Paper claims not recovered:** PSMB5, ABL1, IFNGR1, HLA-A (low-pressure arm). Reported, not hidden.

**Found independently, not in the paper:** ZNF474 (6/6 concordant guides, both replicates, Wald p 7e-05).

---

## 19. Paper-style figures (Figures 1D, 2A, 2B)

`paper_style_figures.py`

### Design decision made and corrected

**v1 sorted genes by rank on the x-axis.** All labelled hits collapsed to x ≈ 0 in a vertical stack with leader lines fanning out — unreadable, and not what the paper shows.

**v2 sorts alphabetically by gene ID**, marking the top 14 by rank. This spreads labels across the plot and reproduces the paper's scatter format.

| Output | Reproduces | Content |
|---|---|---|
| `fig8_paperstyle_enriched.png` | Fig 1D | Enriched sgRNAs, low and high pressure. **NCR3LG1 isolated above the cloud in both panels** (+1.2 low, +7.8 high) |
| `fig9_paperstyle_depleted.png` | Fig 2A | Depleted sgRNAs. **JAK1, IFNGR2, B2M, HLA-E in red** at the bottom of the low panel |
| `fig10_paperstyle_guidecounts.png` | Fig 2B | Paired guide counts Ctrl→Surv for IFNGR2, B2M, HLA-E |

Genes found in **both** selection pressures are coloured red, following the paper's convention.

### Deviations from the paper's figures — state these

**x-axis range.** The paper's Fig 2A runs to ~1,600 genes (low) and ~7,000 (high). Ours runs to 21,697. **They plotted only genes passing an unstated threshold; we plot all of them.** Showing the full distribution is more honest about where hits sit relative to the bulk, but it is a difference.

**Guides per gene in Fig 2B.** The paper shows 3–5 lines per gene. We show all 6. Theirs was a subset.

**Visible artifacts, deliberately not hidden:**
- `hsa-mir-6859-2` sits alone at −8.8 in the depleted high panel. Single-guide artifact (§13). Left labelled so the failure mode is visible.
- The high-pressure panels show a curved arc of points — guides with zero counts hitting a floor value, a consequence of the 26,238 guides lost in High_1.
- In fig10, one B2M guide and one HLA-E guide **rise** rather than fall. Consistent with the guide-level data in §13 and with the general point that guide efficacy varies.

---

## 20. Figure inventory — complete

### Reproducing the paper

| File | Paper figure |
|---|---|
| `fig6_tcga_survival.png` | Fig 5A/B/C |
| `fig7_string_antigen_presentation.png` | Fig 2D |
| `fig7_string_interferon_signaling.png` | Fig 2E |
| `fig8_paperstyle_enriched.png` | Fig 1D |
| `fig9_paperstyle_depleted.png` | Fig 2A |
| `fig10_paperstyle_guidecounts.png` | Fig 2B |
| (`results_enrichment_*.csv`) | Fig 2C |

### Additional — not in the paper

| File | Content |
|---|---|
| `fig1_volcano_low_mle.png` | Volcano with the 1,000 non-targeting controls overlaid |
| `fig2_noise_floor.png` | **Key figure.** Control beta distributions, low vs high, SD annotated |
| `fig3_guide_level.png` | Per-guide CPM, controls vs survivors, grouped by arm |
| `fig4_rra_vs_mle.png` | Rank-rank scatter showing where the two methods disagree |
| `fig5_library_qc.png` | Count distributions and Gini index per sample |

### Not reproducible — wet lab

Figures **1A** (schematic), **1B** (cytotoxicity assay), **1C** (degranulation), **3** (flow cytometry, MHC-I), **4** (arrayed gRNA validation) and **S1**. These require primary NK cells, a flow cytometer and a laboratory.

> **Claim to make: four of the five main figures reproduced in full. Figures 3 and 4 are experimental.** That is more precise, and more credible, than "reproduced the paper."

---

## 21. Remaining work — packaging only

- [ ] Scripts folder: `01_download.sh` … `08_paper_figures.py`
- [ ] README
- [ ] Push to GitHub

**Optional, not blocking:** explicit trim offsets instead of `--trim-5 auto`, to test whether mapping exceeds 53%.

**nf-core status:** attempted, deferred. Nextflow 26.04.6 fails on the pipeline config (`Variable declarations cannot be mixed with config statements`); `NXF_VER=24.10.5` loads it, then the schema requires a `.tsv`/`.txt` library and a samplesheet that was never built. The analysis was already complete, so the pipeline would add packaging rather than results.

---

## 22. Defensibility checklist

1. Why `--trim-5 auto`? *(stagger measured at offsets 28–37 across 100k reads; the paper reports no trim setting)*
2. Why 123,411 → 119,461 guides? *(MAGeCK collapses duplicate sequences: PCDHA/PCDHG, POTEB, PPIAL4, RGPD)*
3. Why is mapping 53% when the log says 69%? *(69% = trim test on the first 100k reads; 53% = full-file rate)*
4. Why `--paired`? *(each control correlates better with its own survivor, 0.843, than with the other control, 0.765)*
5. Why is the high arm's Gini 0.39? *(three rounds of 90% killing; NCR3LG1 at 1683 vs 7.7 CPM proves selection, not drift)*
6. Why did MLE beat RRA on FDR? *(joint model estimates replicate variance instead of collapsing to ranks)*
7. **Why does the high arm have bigger effects but fewer hits?** *(noise SD 1.189 vs 0.210; IFNGR2 falls from 4.4 to 1.8 SD-units)*
8. What do the 1,000 non-targeting controls establish? *(the empirical noise floor)*
9. Why is `hsa-mir-6859-2` not a real hit? *(`goodsgrna` = 1)*
10. Why was GPSM1 excluded despite FDR 0.000? *(one guide, one replicate, 6× jump; z = 1.70)*
11. Why does ZNF474 have a smaller Wald p than HLA-E but a worse FDR? *(permutation vs parametric test)*
12. **Why run pathway enrichment twice?** *(the paper's own gene sets are circular; the top-100 list is not)*
13. Why do chromatin terms appear in the unbiased enrichment? *(general fitness genes deplete under any selection)*
14. Why do the LAML survival p-values differ slightly? *(n=44; one boundary sample moves the log-rank statistic)*
15. Why substitute STRING for IPA? *(IPA is commercial; STRING was already used by the paper for Fig 2D/E)*
16. Why does your Fig 2A x-axis run to 21,697 when theirs stops at 1,600? *(they applied an unstated threshold; we plot everything)*
17. Why was Steinhart rejected as the reproduction target despite better authors? *(no locatable raw reads)*

> **Corollary: do not let an AI write this code.** An artifact you cannot defend is worse than no artifact.

---

## Running log

| Step | Status |
|---|---|
| Paper selected, data availability verified | ✅ |
| WSL + conda + Docker | ✅ |
| 7 FASTQ downloaded (~10 GB) | ✅ |
| Stagger characterised (offsets 28–37) | ✅ |
| Library reformatted → 123,411 guides | ✅ |
| `mageck count` → 119,461 guides, 53% mapped | ✅ |
| QC: Gini, zero counts, correlation | ✅ |
| **Label verification via NCR3LG1 CPM** | ✅ **PASS** |
| RRA, both arms | ✅ |
| Control-sgRNA normalised RRA | ✅ |
| MLE, both arms | ✅ |
| Noise floor quantified | ✅ |
| Guide-level validation | ✅ |
| **TCGA survival (Fig 5)** | ✅ |
| **STRING networks (Fig 2D/E)** | ✅ |
| **Pathway enrichment (Fig 2C)** | ✅ |
| **Paper-style figures (Fig 1D, 2A, 2B)** | ✅ |
| **ANALYSIS COMPLETE** | ✅ |
| Scripts folder | ⏳ |
| README | ⏳ |
