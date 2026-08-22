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

Run an experiment with one plainly named command:

```bash
make bound-validation
make synthetic-comparison
make genomics-comparison
```

The genomics experiment first needs the processed datasets:

```bash
make datasets
make validate
make genomics-comparison
```

The four processed tables occupy about 500 MB. Allow roughly 10 GB of working
space for downloaded assemblies, build products, and processed data.

`make reproduce` runs all three experiments, verifies their outputs, and marks
the repository's completion hook satisfied. The
synthetic comparison uses one million keys and the genomics comparison includes
rice with 198 million distinct 15-mers, so the complete run is not a quick
smoke test. Use `make repro-small` for E. coli on all four methods plus the nine
bound-validation plots and the synthetic-comparison plot.

On Apple Silicon, the checked-in Dockerfile provides the required
`linux/amd64` environment for VL-BuRR:

```bash
./scripts/docker_run.sh make repro-small
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

- `methods/decision_rules.py`: AutoCSF, HKP, and BCSF selection rules over the
  same Caramel implementation.
- `vlburr/`: pinned upstream commit, small correctness/memory patches, build,
  and runners for VL-BuRR.
- `deps/CaramelDB/`: the pinned CSF/filter implementation.
- `reference/accepted-paper/`: immutable historical PNGs used as visual
  references, protected by checksums.
- `results/reproduction/`: bounded and full synthetic reproduction receipts.

`make synthetic-clean` builds a fresh Linux container, removes bundled
synthetic results inside it, recomputes all 45 bound-validation datasets and
all 88 comparison rows, renders the ten plots, and validates the fresh numbers
before promoting them.

`make verify` validates dataset checksums and profiles, checks the VL-BuRR
integer-frequency regression, verifies the expected tables/figures/receipts,
and marks the Codex stop hook complete.
