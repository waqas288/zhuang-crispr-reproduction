"""
Prepare the GeCKO v2 human sgRNA library for MAGeCK.

Source: Addgene pooled library #1000000049 (GeCKO v2 human, two-plasmid).
Ships as two CSVs -- half-libraries A and B, 3 guides per gene each, 6 combined.
The paper used the full library, so both are needed.

Sequence CSVs are a free download. The plasmid itself is a paid physical product.

THREE PROBLEMS, none of which raise errors
------------------------------------------
1. LINE ENDINGS. The files use classic Mac '\r' only. Linux tools read the
   whole file as one record -- `head` dumps everything onto a single line.

2. COLUMN ORDER. Files are `gene_id, UID, seq`.
   MAGeCK requires `sgRNA_id, sequence, gene`.
   Fed in unmodified, MAGeCK produces a WRONG ANSWER rather than a failure.

3. A trailing empty column on every row.

This is the most dangerous class of bioinformatics bug: input that is wrong
but parseable. Nothing throws. You get an answer, and it is incorrect.
(Same category as the DepMap `IsDefaultEntryForModel == True` filter that
silently returns zero rows because the column holds 'Yes'/'No' strings.)

VALIDATION
----------
Output should be 123,411 guides. That matches the published GeCKO v2 figure
exactly -- the agreement is the check that the merge was correct.

Expected composition:
  21,915 gene labels
  18,940 genes with the full 6 guides
   1,853 miRNA entries (hsa-mir-*)
   1,000 non-targeting controls (NonTargetingControlGuideForHuman_0001..1000)

The non-targeting controls matter later: they cut nothing, so their spread
measures the noise floor empirically. Used in 04_test.sh for control
normalisation and in 06_figures.py to establish what "no effect" looks like.
"""

import sys

INPUTS = [
    "human_geckov2_library_a_09mar2015.csv",
    "human_geckov2_library_b_09mar2015.csv",
]
OUTPUT = "GeCKOv2_library_mageck.csv"

rows = []

for path in INPUTS:
    try:
        raw = open(path, "rb").read().decode("utf-8", "replace")
    except FileNotFoundError:
        sys.exit(f"ERROR: {path} not found. Download both half-libraries from "
                 f"Addgene #1000000049 and place them in this directory.")

    # normalise line endings -- files use bare \r
    raw = raw.replace("\r\n", "\n").replace("\r", "\n")
    lines = [l for l in raw.split("\n") if l.strip()]
    print(f"{path}: {len(lines)} lines")

    for line in lines:
        parts = [p.strip() for p in line.split(",") if p.strip() != ""]
        if len(parts) < 3:
            continue
        gene, sgid, seq = parts[0], parts[1], parts[2]
        # reorder to MAGeCK's expected sgRNA_id, sequence, gene
        rows.append((sgid, seq, gene))

# drop header and anything that is not a clean 20bp ACGT guide
clean = [r for r in rows
         if len(r[1]) == 20 and set(r[1]) <= set("ACGT")]

print(f"\nparsed:  {len(rows)}")
print(f"kept:    {len(clean)}  (20bp ACGT only)")
print(f"dropped: {len(rows) - len(clean)}  (header rows, malformed entries)")

with open(OUTPUT, "w") as fh:
    for sgid, seq, gene in clean:
        fh.write(f"{sgid},{seq},{gene}\n")

# --- validation ---
genes = {}
for _, _, g in clean:
    genes[g] = genes.get(g, 0) + 1

from collections import Counter
per_gene = Counter(genes.values())

print(f"\nwrote {OUTPUT}")
print(f"  guides:            {len(clean)}   (expect 123,411)")
print(f"  gene labels:       {len(genes)}")
print(f"  genes w/ 6 guides: {per_gene.get(6, 0)}")
print(f"  miRNA entries:     {sum(1 for g in genes if g.startswith('hsa-mir'))}")
print(f"  non-targeting:     {sum(1 for g in genes if g.startswith('NonTargetingControl'))}")

if len(clean) != 123411:
    print("\n  WARNING: guide count does not match the published figure. "
          "Check that both half-libraries were present and complete.")
