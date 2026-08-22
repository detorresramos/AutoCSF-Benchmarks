# AutoCSF Benchmarks

Reproduction artifact for **AutoCSF: A Provably Safe Indexing Framework for
Filter-Augmented Compressed Static Functions**.

This repository contains the AutoCSF implementation used in the experiments,
the HKP and BCSF decision rules, the VL-BuRR baseline, dataset preparation
tools, and scripts that generate the paper's figures and tables.

## Reproduce the paper

An AI coding agent should read [`AGENTS.md`](AGENTS.md) first. On Ubuntu
24.04 x86-64, the intended human workflow is:

```bash
git clone --recursive https://github.com/detorresramos/AutoCSF-Benchmarks.git
cd AutoCSF-Benchmarks
./setup.sh

# Download the published, processed genomics tables.
python datasets/manage.py download
python datasets/manage.py validate

# Reproduce synthetic experiments and existing figures.
make experiments plots

# Or force a clean linux/amd64 Docker rebuild and validate fresh numbers.
make synthetic-clean

# Re-run the camera-ready genomics comparison.
./run_genomics.sh
```

`run_genomics.sh` writes raw JSON results and generated Markdown/LaTeX tables
under `results/`. One run is the default; pass `--repetitions 3` only when the
machine has enough time and memory. It never edits paper source files.

On Apple Silicon, use the checked-in `linux/amd64` Dockerfile for VL-BuRR,
whose upstream implementation uses x86 intrinsics. The bounded verification
command is `make repro-small`: E. coli on all four methods plus all ten plots.
Rice needs about 14 GB for VL-BuRR; a memory failure should be reported rather
than replaced with a downsampled result.

`make synthetic-clean` deletes every bundled synthetic result inside a fresh
container before measurement, reruns all 45 theory experiments and all 88
four-method comparison rows, renders the ten paper plots from those fresh
numbers, and writes `results/reproduction/synthetic-receipt.json`. The command
only promotes the fresh artifacts into the checkout after the numerical
comparison passes.

## Genomics data

The public dataset artifact contains both original source files and processed
count tables. Processed files use one transparent format:

```text
AAAAAAAAAAAAAAA\t3
AAAAAAAAAAAAAAC\t1
```

They are UTF-8 TSV files without a header, sorted by k-mer and compressed as
`.tsv.zst`. Each row is a distinct forward-strand (non-canonical) 15-mer and
its exact occurrence count. ASCII `acgt` is normalized to uppercase; windows
containing other characters are omitted. See [`datasets/README.md`](datasets/README.md).

By default, experiments download the processed files. Regenerating them from
the original FASTA/FASTQ files is optional:

```bash
python datasets/manage.py generate --dataset celegans
python datasets/manage.py validate --dataset celegans
```

Set `AUTOCSF_DATASET_REPO` to the Hugging Face dataset repository identifier.
Until the public repository is created, the default is the proposed
`detorresramos/autocsf-genomics`.

## Methods

The camera-ready table compares four decision frameworks:

VL-BuRR is built with a documented integer frequency-counter correction for
classes larger than `2^24`. C. elegans and rice were rerun with that fix. Rice
requires roughly 16 GB or more for the container; the recorded run used a
21.2 GB Docker memory allocation.

- **HKP** — idealized Bloom-filter cost model from Hreinsson et al.
- **BCSF** — Shibuya et al.'s entropy-based heuristic.
- **VL-BuRR** — generalized filter trick using its native implementation.
- **AutoCSF** — the provably safe lower-bound criterion.

HKP, BCSF, and AutoCSF use the same CaramelDB CSF and Bloom-filter
implementation; only their parameter-selection rule differs. VL-BuRR is built
from a pristine, commit-pinned LearnedStaticFunction/BuRR checkout by
`vlburr/build.sh`. The separate `autocsf-bench` evaluation repository is not a
dependency; only its small evaluation patch was vendored here.

## Repository layout

```text
datasets/              download, generation, validation, and provenance
shared/                AutoCSF, HKP, and BCSF decision rules
theory_validation/     bound-validation experiments
shibuya_comparison/    BCSF comparison experiments
baselines/             end-to-end benchmark runner
deps/CaramelDB/        pinned C++ implementation
deps/LearnedStaticFunction/ upstream VL-BuRR/LSF implementation
results/               generated raw results and tables (ignored)
```

## Development setup

For a quick local build without installing system packages:

```bash
./setup.sh --no-system
source .venv/bin/activate
python -m unittest discover
```

The canonical performance environment is Ubuntu x86-64. The committed camera-
ready table uses native Apple Silicon for the three Caramel methods and the
linux/amd64 container for x86-only VL-BuRR. The primary comparison is bits per
key saved relative to each method's own plain CSF: HKP, BCSF, and AutoCSF use
the Caramel/GOV plain CSF, while VL-BuRR uses `Filtered-HuffmanCSF_No`. Absolute
sizes from unlike CSF implementations do not share a baseline. Structure size
is platform-independent; those mixed-host timings are diagnostic,
not a fair cross-method speed comparison. The JSON records both environments.

The generated genomics outputs separate these two questions:

- `results/genomics-paper.md` reports bits/key saved relative to each method's
  own plain CSF and isolates the filter-selection decision.
- `results/genomics-total.md` reports final end-to-end bits/key, including the
  different underlying CSF implementations.
- `results/genomics-audit.md` shows plain size, filtered size, savings, and the
  selected filter parameters for every dataset/method pair.
