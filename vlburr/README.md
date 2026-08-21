# VL-BuRR evaluation adapter

The build starts from the commit in `UPSTREAM_COMMIT` and applies the patches
in this directory. `integer-frequency-counts.patch` fixes the `autocsf-eval`
frequency model's loss of unit increments above `2^24` observations per class.
`release-filter-input.patch` releases a completed construction buffer before
allocating the next O(N) buffer; it changes peak memory, not the algorithm.

This directory vendors the small CSF evaluation layer from
`vihan-lakshman/autocsf-bench` at commit
`76eed05857c8a1a911abffcb2da4174c275ec068`. It is not a dependency.

`patches/nofilter.patch` adds the upstream no-filter comparison to
LearnedStaticFunction's `ribbon_learned_bench`. The implementation remains the
GPLv3 LearnedStaticFunction submodule. The scripts here convert the public TSV
tables to its `.lrbin` input and parse its reported storage/build/query metrics.
