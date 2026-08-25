# AutoCSF Benchmarks

## Introduction

This repository implements the experiments in [**AutoCSF: A Provably Safe
Indexing Framework for Filter-Augmented Compressed Static
Functions**][paper] (CIKM '26). The command `make reproduce` regenerates the
experimental results from the paper from source and checks it against the published
numbers in `baselines/`.

**The problem.** A compressed static function (CSF) maps a fixed key set to
values in space proportional to the entropy of the value sequence, without
storing the keys. It has a floor: no CSF uses less than one bit per key. When a
single value dominates — 72.7% of rice's 198 million 15-mers share one count —
that floor, not the entropy, is what binds, and the standard remedy is to put an
approximate membership filter in front of the CSF. Only minority keys go into
the filter, so a negative answer means "the dominating value" and never reaches
the CSF, letting the index break the one-bit barrier.

The filter is not free. It occupies space of its own, and the keys it reports as
false positives still have to be stored in the CSF. Whether the trade pays off
depends on how skewed the values are, and the paper's starting observation is
that every prior method decides using an idealized cost model — so each has a
**dead zone**, a band of dominating-value fractions $\alpha$ where it recommends
a filter that makes the index larger.

**What the experiments test.** AutoCSF (Algorithm 1) bounds the *difference* in
space between the filtered and unfiltered designs instead of modelling either
one, and filters only when that lower bound (Theorem 4.2) is positive. The
guarantee is safety: with high probability it does not recommend a filter that
increases space. The three experiments test whether the bounds are tight enough
to be useful, whether the decision rule avoids the dead zone that prior rules
fall into, and whether both hold on real workloads.

| Experiment | Question | In the paper | Main output |
|---|---|---|---|
| Bound validation | Do the theoretical bounds match measured savings as $\alpha$ and $\varepsilon$ vary? | §6.1–6.2, Figures 3 and 4 | Nine alpha/epsilon sweep plots |
| Synthetic comparison | Do HKP, BCSF, or GFT recommend filters that increase space? | §6.3, Figure 2 | The four-method comparison plot |
| Genomics comparison | Does the same behavior occur on real 15-mer count tables? | §6.4, Table 1 | The four-dataset savings table |

[paper]: https://doi.org/10.1145/3799682.3841067

Everything else in this repository supports those experiments: `methods/`
contains the four decision frameworks, `common/` contains shared math and data
generation, and `datasets/` prepares and validates genomics inputs.

## Running it

```bash
git clone --recursive https://github.com/detorresramos/AutoCSF-Benchmarks.git
cd AutoCSF-Benchmarks
sudo ./setup.sh          # builds CaramelDB, VL-BuRR, and the Python venv
make reproduce           # runs all three experiments and checks the results
```

That is the whole thing. `make reproduce` writes every number and figure the
paper uses into `results/`, then verifies them and writes a reproduction
receipt. Individual experiments:

```bash
make bound-validation      # experiment 1: the nine alpha/epsilon sweep plots
make datasets              # download and validate the genomics tables
make synthetic-comparison  # experiment 2: the four-method comparison plot
make genomics-comparison   # experiment 3: the savings table
```

`make datasets` is only needed for the genomics experiment. The four processed
tables occupy about 500 MB; allow roughly 10 GB of working space for downloaded
assemblies and build products. The synthetic comparison uses one million keys
and the genomics comparison includes rice with 198 million distinct 15-mers, so
a full `make reproduce` is not a quick smoke test.

### Where it runs

Everything above works natively on Ubuntu 24.04 x86-64 and on macOS arm64, with
one exception: **VL-BuRR is x86-64 only**, because upstream uses SSE/AVX
intrinsics that do not compile on ARM. `setup.sh` detects this and skips it
rather than failing, so on a Mac you get bound validation, genomics, and three
of the four methods in the synthetic comparison.

To get the VL-BuRR arm on a Mac, prefix any command with `docker_run.sh`:

```bash
./scripts/docker_run.sh make synthetic-comparison
```

That builds a `linux/amd64` image and runs the target inside it. The container
is also the reference environment — it is where the committed reproduction
receipts were produced, and the only place all four methods have run together —
so `./scripts/docker_run.sh make reproduce` is the most faithful reproduction
regardless of host.

On macOS, `setup.sh` additionally needs `libomp` (`brew install libomp`), and
`cmake` and `zstd` on `PATH`. Use `./setup.sh --no-system` to skip the
`apt-get` step if the system packages are already installed.

## 1. Bound validation

Code: `experiments/bound_validation/`

This experiment sweeps the dominating value fraction $\alpha$ and filter false
positive rate $\varepsilon$ across three minority distributions and five filter
configurations. It compares empirical bits/key savings with AutoCSF's lower and
upper bounds.

```bash
make bound-validation
```

Outputs:

- Numerical measurements: `results/bound_validation/data/`
- Nine paper plots: `results/bound_validation/figures/paper/`
- Individual diagnostic plots: `results/bound_validation/figures/alpha_sweep/`
  and `results/bound_validation/figures/epsilon_sweep/`

## 2. Synthetic method comparison

Code: `experiments/synthetic_comparison/`

This experiment compares the filter decisions made by HKP, BCSF, GFT,
and AutoCSF on Uniform-100 and Zipfian values while sweeping $\alpha$. Its
metric is bits/key saved relative to each method's own unfiltered CSF; negative
values mean the selected filter made the index larger.

```bash
make synthetic-comparison
```

Outputs:

- Measurements: `results/synthetic_comparison/data.csv`
- Paper plot: `results/synthetic_comparison/method-comparison.png`

## 3. Real-world genomics comparison

Code: `experiments/genomics_comparison/`

This runs the same four decision frameworks on E. coli Sakai, SRR10211353,
C. elegans, and rice. The camera-ready result is the filter-versus-no-filter
savings table.

```bash
make datasets
make genomics-comparison
```

Outputs:

- Raw records: `results/genomics_comparison/genomics.json`
- Measurements: `results/genomics_comparison/genomics.csv`
- Readable table, with the plain-CSF baseline and the filter each method
  selected: `results/genomics_comparison/genomics.md`

The primary comparison is bits/key saved relative to each method's own plain
CSF. HKP, BCSF, and AutoCSF share the Caramel/GOV CSF and Bloom-filter
implementation; only their decision rule differs. GFT uses VL-BuRR's native
Huffman CSF and filter, so absolute end-to-end sizes do not share a baseline.

## How each method is run

Four decision rules are compared. Three are policies over one shared
implementation; the fourth is a separate system built from its authors' code.

### HKP, BCSF, and AutoCSF

These three build a Caramel/GOV CSF over the same keys and values with the same
Bloom filter implementation, and differ only in `select_filter()`
(`methods/decision_rules.py`). It takes four summary statistics — $\alpha$,
$n/N$, the zeroth-order entropy $H_0$, and how many keys the filter would hold —
and returns a Bloom configuration, or `None` when the rule declines to filter:

- **HKP** (Hreinsson, Krøyer and Pagh, 2009) filters when $\alpha > 0.63$ and
  targets $\varepsilon = 1 - \alpha$. That is a false positive rate, not a
  filter, so the harness realizes it as the Bloom configuration nearest to it in
  log space — searching only configurations that no other beats on both size and
  $\varepsilon$, so the rule cannot be charged for a wasteful filter.
- **BCSF** (Shibuya, Belazzougui and Kucherov, 2022) uses their
  $\varepsilon^* = C_{BF}(1-\alpha) / (C_{CSF}\,\alpha \ln 2)$, with
  $C_{BF} = 1.44$ and their piecewise data-driven fit for $C_{CSF}(H_0)$, both
  reproduced verbatim in `common/bcsf.py`.
- **AutoCSF** evaluates the Theorem 4.2 lower bound at every discrete filter
  configuration available and takes the largest, declining to filter when even
  that is not positive.

Because the CSF and the filter are identical across all three, any difference in
measured savings is attributable to the decision rule alone.

### GFT (VL-BuRR)

The fourth method is the generalised filter trick of Hermann et al. (2025). It
is not a policy over Caramel — it is a separate system with its own CSF, its own
filter, and its own input format — so it is built and run as one.

`methods/vlburr/build.sh` clones the authors' own repository,
[LearnedStaticFunction](https://github.com/gvinciguerra/LearnedStaticFunction),
at the commit pinned in `methods/vlburr/UPSTREAM_COMMIT`, applies the patches
below, and builds `ribbon_learned_bench`. Upstream source is never vendored
here; only the patches are.

It runs through the authors' own entry point, `ribbon_learned_bench -c CSF` —
the configuration their paper reports as "Ours (CSF)": their learned static
function with the model replaced by one that emits global value frequencies, and
a single global Huffman code shared by all keys. The filtered arm uses their
`FilterLengthStrategyOpt` unmodified, choosing an integer filter weight per node
of the code tree. The unfiltered arm it is measured against uses their own
`FilterLengthStrategyNoFilter`.

Four patches are applied, all listed in `methods/vlburr/README.md`. One is
material to the reported numbers: upstream's frequency model accumulates class
counts in a `float`, which stops incrementing above $2^{24}$. On C. elegans and
rice this understates the dominating fraction badly enough to suppress filtering
altogether, so the counts are made exact before normalizing; `make verify`
carries a regression check for it. The other three register the unfiltered arm
as a benchmark, release a build buffer to lower peak memory, and fix a macOS
include.

The rest of `methods/vlburr/` is adapters. `gen_datasets.py` and
`table_to_lrbin.py` convert the synthetic and genomics inputs to VL-BuRR's
`.lrbin` format — the synthetic path uses the same generator and seed as the
other three methods, so every method sees an identical value sequence — and
`parse_genomics.py` reads back the metrics the binary reports.

## Genomics data format

Processed datasets are UTF-8, headerless, k-mer-sorted TSV compressed with
Zstandard:

```text
AAAAAAAAAAAAAAA\t3
AAAAAAAAAAAAAAC\t1
```

Each row is a distinct forward-strand (non-canonical) 15-mer and its exact
count. Lowercase bases are normalized to uppercase and windows containing other
characters are omitted. See `datasets/README.md` for sources, checksums,
generation, validation, and Hugging Face staging.

## Implementations and reproducibility

All four methods live under `methods/`, as described in
[How each method is run](#how-each-method-is-run). The rest of the layout:

- `deps/CaramelDB/`: the pinned CSF/filter implementation, used unmodified.
  `setup.sh` is the single build definition — it builds the static library the
  genomics harness links, installs the `carameldb` Python module every other
  experiment imports, and builds VL-BuRR. The Dockerfile just calls it.
- `baselines/`: the accepted numbers, one CSV per experiment. Every run is
  checked against these.
- `reference/`: the figures and genomics table as they appear in the accepted
  paper, for side-by-side comparison with a fresh run.
- `results/`: everything a run produces. Not tracked in git -- `make reproduce`
  regenerates it, and `results/reproduction/receipt.json` records what the last
  `make verify` found.

### Checking a run against the accepted results

`make verify` runs at the end of `make reproduce`, and can be run on its own. It
checks that each experiment produced a complete set of outputs -- dataset
checksums, the 45 bound-validation sweeps and their figures, the four-method
synthetic table, the 4x4 genomics matrix, the ten paper plots -- then compares
every measurement in `results/` against `baselines/` and writes
`results/reproduction/receipt.json`.

Comparison is numeric with a 0.02 bits/key tolerance rather than exact. Most
measurements are bit-stable, but binary-fuse CSF construction is randomized:
repeated builds of the same input differ by a few thousandths of a bit per key,
which is why the nine `binary_fuse` datasets never reproduce exactly.

If a run changes the accepted numbers on purpose, `make baselines` promotes what
is in `results/` to `baselines/`.
