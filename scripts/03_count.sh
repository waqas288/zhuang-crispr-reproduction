#!/usr/bin/env bash
set -euo pipefail
#
# FASTQ -> counts table.
#
# --trim-5 auto is REQUIRED. GeCKO v2 is a STAGGERED library: a variable-length
# filler sits between the sequencing primer and the constant vector region, so
# the guide does not start at a fixed position. Measured directly:
#
#   zcat SRR10338906.fastq.gz | awk 'NR%4==2' | head -100000 | \
#     awk '{n=index($0,"AAAGGACGAAAC"); if(n>0) print n}' | sort -n | uniq -c
#
# -> ten distinct offsets, positions 28-37, ~98% of reads. A fixed --trim-5
#    would discard most of them.
#
# Staggers exist because Illumina calibrates base quality across clusters at
# each cycle; without offsetting, every read would be identical for the first
# ~30 bases and base diversity would be too low to call accurately.
#
# Sample order follows the GEO GSM listing. VERIFY against the SRA Run Table
# before trusting results -- swapping control and survivor labels would invert
# every conclusion while still producing a plausible-looking answer.

mageck count \
  -l GeCKOv2_library_mageck.csv \
  -n zhuang \
  --sample-label "Low_Ctrl_1,Low_Surv_1,Low_Ctrl_2,Low_Surv_2,High_Ctrl,High_1,High_2" \
  --fastq data/SRR10338906.fastq.gz data/SRR10338907.fastq.gz \
          data/SRR10338908.fastq.gz data/SRR10338909.fastq.gz \
          data/SRR10338910.fastq.gz data/SRR10338911.fastq.gz \
          data/SRR10338912.fastq.gz \
  --trim-5 auto

echo
echo ">>> QC summary"
cat zhuang.countsummary.txt
echo
echo "NOTE: the .log file reports TWO mapping rates. ~69% is MAGeCK's trim test"
echo "      on the first 100,001 reads only. ~53% is the real full-file rate."
