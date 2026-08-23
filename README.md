# AutoCSF Benchmarks

Reproduction artifact for **AutoCSF: A Provably Safe Indexing Framework for
Filter-Augmented Compressed Static Functions**.

The paper contains three experiments:

| Experiment | Question | Main output |
|---|---|---|
| Bound validation | Do the theoretical bounds match measured savings as $\alpha$ and $\varepsilon$ vary? | Nine alpha/epsilon sweep plots |
| Synthetic comparison | Do HKP, BCSF, or VL-BuRR recommend filters that increase space? | The four-method comparison plot |
| Genomics comparison | Does the same behavior occur on real 15-mer count tables? | The four-dataset savings table |

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

This experiment compares the filter decisions made by HKP, BCSF, VL-BuRR,
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
- Camera-ready table: `results/genomics_comparison/genomics-paper.tex`
- Human-readable table: `results/genomics_comparison/genomics-paper.md`
- Audit details: `results/genomics_comparison/genomics-audit.md`

The primary comparison is bits/key saved relative to each method's own plain
CSF. HKP, BCSF, and AutoCSF share the Caramel/GOV CSF and Bloom-filter
implementation; only their decision rule differs. VL-BuRR uses its native
Huffman CSF and filter, so absolute end-to-end sizes do not share a baseline.

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

All four methods live under `methods/`:

- `methods/decision_rules.py`: AutoCSF, HKP, and BCSF selection rules over the
  same Caramel implementation. They differ only in how they pick a filter.
- `methods/vlburr/`: VL-BuRR is an external tool rather than a decision rule, so
  it gets a subtree: a pinned upstream commit, four patches (see
  `methods/vlburr/README.md`), a build script, and the adapters that convert
  data into its input format and parse its reported metrics.
- `deps/CaramelDB/`: the pinned CSF/filter implementation, used unmodified.
  `setup.sh` is the single build definition — it builds the static library the
  genomics harness links, installs the `carameldb` Python module every other
  experiment imports, and builds VL-BuRR. The Dockerfile just calls it.
- `baselines/`: the accepted numbers, one CSV per experiment. Every run is
  checked against these.
- `reference/`: the figures and genomics table as they appear in the accepted
  paper, for side-by-side comparison with a fresh run.
- `results/reproduction/receipt.json`: what the last `make verify` found.

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
