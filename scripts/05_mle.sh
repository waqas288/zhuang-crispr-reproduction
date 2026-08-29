#!/usr/bin/env bash
set -euo pipefail
#
# Hit calling with MAGeCK MLE. This matches the paper's method -- MAGeCK-VISPR
# reports normalised log fold change per gene, which is MLE output.
#
# MLE has more power than RRA here: RRA treats the samples as lists to rank,
# while MLE fits one model across all of them, estimating replicate variance
# properly instead of collapsing it. Result: 7 genes at FDR<0.05 vs 3 under RRA.
#
# Design matrix:
#   baseline    = 1 for every sample (the intercept)
#   NKselection = 1 only for survivors
# The NKselection beta is the effect of interest.
#   negative -> depleted under NK pressure (vulnerability)
#   positive -> enriched (resistance)
#
# Runtime ~14-15 min per arm. RRA takes ~1 min.

printf "Samples\tbaseline\tNKselection\nLow_Ctrl_1\t1\t0\nLow_Surv_1\t1\t1\nLow_Ctrl_2\t1\t0\nLow_Surv_2\t1\t1\n" > config/design_low.txt
printf "Samples\tbaseline\tNKselection\nHigh_Ctrl\t1\t0\nHigh_1\t1\t1\nHigh_2\t1\t1\n" > config/design_high.txt

echo ">>> design matrices"
cat config/design_low.txt; echo; cat config/design_high.txt; echo

mageck mle -k zhuang.count.txt -d config/design_low.txt  -n low_mle
mageck mle -k zhuang.count.txt -d config/design_high.txt -n high_mle

echo
echo "NOTE: the gene_summary 'fdr' column comes from the PERMUTATION test,"
echo "      'wald-p-value' from the parametric test on the beta estimate."
echo "      They rank genes differently. ZNF474 has a Wald p 14x smaller than"
echo "      HLA-E but a much worse FDR. Choose deliberately."
