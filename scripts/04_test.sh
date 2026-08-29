#!/usr/bin/env bash
set -euo pipefail
#
# Hit calling with MAGeCK RRA (robust rank aggregation).
#
# RRA rather than a mean fold change because guide-level variation exceeds the
# biological signal. Demonstrated on practice data: BUB3's 10 guides spanned
# -4.216 to +0.987, SD 1.62 on a mean of -1.37. Two guides called it the most
# essential gene in the dataset; two called it enriched.
#
# RRA ranks all guides, then asks whether a gene's guides cluster at one end.
# It tests CONSISTENCY OF DIRECTION, not magnitude of average. The goodsgrna
# column reports how many guides individually support the call.
#
# Read BOTH directions: this is a selection screen, not a dropout screen.
#   neg| = vulnerability (depleted in survivors)
#   pos| = resistance    (enriched in survivors)

# --- low selection pressure, PAIRED ---
#
# --paired because the correlation structure demands it. log10 count correlation:
#   Low_Ctrl_1 vs Low_Surv_1 (its own survivor) = 0.843
#   Low_Ctrl_1 vs Low_Ctrl_2 (the other control) = 0.765
# Replicate identity dominates the treatment effect. Pooling all controls
# against all survivors would let batch differences leak into the contrast.
mageck test -k zhuang.count.txt \
  -t Low_Surv_1,Low_Surv_2 -c Low_Ctrl_1,Low_Ctrl_2 \
  -n low_rra --paired

# --- high selection pressure ---
# Cannot be paired: one control, two technical replicates.
mageck test -k zhuang.count.txt \
  -t High_1,High_2 -c High_Ctrl \
  -n high_rra

# --- control-normalised variant ---
#
# Default normalisation scales samples to a common median across all guides,
# assuming most guides are unaffected. Reasonable, but an assumption.
# Normalising to the 1,000 non-targeting guides uses the MEASURED null instead.
# Result: the paper's genes shift by only 1-3 rank positions -> robust.
grep "NonTargetingControlGuideForHuman" zhuang.count.txt | cut -f1 > nt_controls.txt
wc -l nt_controls.txt          # expect 1000

mageck test -k zhuang.count.txt \
  -t Low_Surv_1,Low_Surv_2 -c Low_Ctrl_1,Low_Ctrl_2 \
  -n low_rra_ctrlnorm --paired \
  --control-sgrna nt_controls.txt --norm-method control
