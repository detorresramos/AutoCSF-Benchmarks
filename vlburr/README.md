# VL-BuRR evaluation adapter

This directory vendors the small CSF evaluation layer from
`vihan-lakshman/autocsf-bench` at commit
`76eed05857c8a1a911abffcb2da4174c275ec068`. It is not a dependency.

`patches/nofilter.patch` adds the upstream no-filter comparison to
LearnedStaticFunction's `ribbon_learned_bench`. The implementation remains the
GPLv3 LearnedStaticFunction submodule. The scripts here convert the public TSV
tables to its `.lrbin` input and parse its reported storage/build/query metrics.
