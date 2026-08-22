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

## Quick start

On Ubuntu 24.04 x86-64:

```bash
git clone --recursive https://github.com/detorresramos/AutoCSF-Benchmarks.git
cd AutoCSF-Benchmarks
./setup.sh
```

Run an experiment:

```bash
make bound-validation
make synthetic-comparison
```

The genomics experiment first needs the processed datasets, which `make
datasets` downloads and validates:

```bash
make datasets
make genomics-comparison
```

The four processed tables occupy about 500 MB. Allow roughly 10 GB of working
space for downloaded assemblies, build products, and processed data.

`make reproduce` runs all three experiments and verifies their outputs. The
synthetic comparison uses one million keys and the genomics comparison includes
rice with 198 million distinct 15-mers, so the complete run is not a quick
smoke test.

On Apple Silicon, the checked-in Dockerfile provides the required
`linux/amd64` environment for VL-BuRR. `scripts/docker_run.sh` is the only
container entry point; it builds the image and runs any target inside it:

```bash
./scripts/docker_run.sh make bound-validation
```

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
  it gets a subtree: a pinned upstream commit, small correctness/memory patches,
  a build script, and the adapters that convert data into its input format and
  parse its reported metrics.
- `deps/CaramelDB/`: the pinned CSF/filter implementation. `setup.sh` is the
  single build definition — it patches CaramelDB, builds the static library that
  the genomics harness links, installs the `carameldb` Python module that every
  other experiment imports, and builds VL-BuRR. The Dockerfile just calls it.
- `reference/`: the figures and genomics table as they appear in the accepted
  paper, for side-by-side comparison with a fresh run. Nothing verifies their
  contents, so treat them as a record rather than as ground truth.
- `results/reproduction/`: reproduction receipts.

### Checking a run against the accepted results

`make verify` checks that the committed results are complete and internally
consistent: dataset checksums and profiles, the VL-BuRR integer-frequency
regression, the expected 45 bound-validation datasets and their figures, the
four-method synthetic table, the 4x4 genomics matrix, the ten paper plots, the
presence of the `reference/` copies, and the receipts.

`make repro-clean` goes further. It builds a fresh Linux container that never
sees `results/`, recomputes all 45 bound-validation datasets and all 88
synthetic-comparison rows, renders the ten plots, and compares the fresh numbers
against the committed ones before promoting anything. `SCOPE=small make
repro-clean` adds E. coli genomics on all four methods.

Comparison is numeric with a 0.02 bits/key tolerance rather than exact. Most
measurements are bit-stable, but binary-fuse CSF construction is randomized:
repeated builds of the same input differ by a few thousandths of a bit per key,
which is why the nine `binary_fuse` datasets never reproduce exactly.
