# VL-BuRR evaluation adapter

The fourth method in the comparison. Unlike the rules in
`methods/decision_rules.py`, VL-BuRR is not a policy over the shared Caramel
implementation — it is a separate binary with its own CSF, its own filter, and
its own input format. That is why it gets a subtree rather than a module.

- `UPSTREAM_COMMIT` / `build.sh`: clone GPLv3
  [LearnedStaticFunction](https://github.com/gvinciguerra/LearnedStaticFunction)
  at the pinned commit into `data/cache/`, apply the patches below, and build
  `ribbon_learned_bench`. Upstream source is never vendored into this
  repository; only the patches are.
- `gen_datasets.py`: writes the synthetic comparison's values as `.lrbin`, using
  the same `common.data_generation.gen_alpha_values` and the same seed as the
  local methods, so every method sees an identical value sequence.
- `table_to_lrbin.py` / `run_genomics.sh` / `parse_genomics.py`: convert a
  processed k-mer table to `.lrbin`, run the benchmark, and parse its reported
  storage/build/query metrics.

## Patches

- `nofilter.patch` adds the no-filter comparison to `ribbon_learned_bench`, so
  VL-BuRR can be measured against its own unfiltered baseline.
- `integer-frequency-counts.patch` fixes the `autocsf-eval` frequency model's
  loss of unit increments above `2^24` observations per class. `make verify`
  carries a regression check for this.
- `release-filter-input.patch` releases a completed construction buffer before
  allocating the next O(N) buffer. It changes peak memory, not the algorithm.
- `macos-cstdint.patch` is applied only on Darwin.
