# Known correctness issue in checkpoint results

The VL-BuRR rows for C. elegans and rice are provisional and must not be used
in the paper yet. The `autocsf-eval` frequency model accumulates class counts
in IEEE-754 single-precision `float` values. Incrementing such a counter stops
working reliably above `2^24 = 16,777,216`.

The majority classes in C. elegans (about 57.4 million keys) and rice (about
144 million keys) exceed that limit. VL-BuRR therefore receives incorrect
class probabilities and emits an empty filter (`Max length filter: 0`) for both
datasets. E. coli and SRR10211353 do not exceed the limit and are unaffected.

Before finalizing results, change the frequency model to accumulate counts in
an integer type, normalize after counting, add a regression test, rebuild, and
rerun VL-BuRR on C. elegans and rice. Remove this file only after the corrected
results have replaced the provisional rows.
