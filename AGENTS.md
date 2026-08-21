# Reproduction instructions

The goal of this repository is a fresh-clone reproduction of every AutoCSF
paper experiment. Prefer simple scripts and explicit files over framework
machinery.

## Safe default workflow

1. Read `README.md` and `datasets/README.md`.
2. Initialize submodules with `git submodule update --init --recursive`.
3. Run `./setup.sh --no-system` when system dependencies already exist.
4. Download processed genomics data with `python datasets/manage.py download`.
5. Validate checksums/statistics before running experiments.
6. Write generated results only beneath `results/` or existing `figures/data/`.

Do not download hundreds of gigabytes or launch the full benchmark unless the
user explicitly asks. `python datasets/manage.py validate --manifest-only` is
the lightweight dataset check.

## Reproducibility rules

- Never hand-edit generated tables or figures.
- Do not mix timing results from different hosts.
- Record git commits, compiler versions, command lines, and dataset SHA-256s.
- The authoritative genomics format is headerless `kmer<TAB>count` in `.tsv.zst`.
- Genomics k-mers are `k=15`, forward/non-canonical; ASCII `acgt` is normalized
  to uppercase and other bases invalidate the overlapping window.
- Treat files beneath `data/cache/` as disposable method-specific derivatives.
- Preserve unrelated working-tree changes, especially inside submodules.

## Paper methods

The camera-ready genomics comparison is HKP, BCSF, VL-BuRR, and AutoCSF.
HKP/BCSF/AutoCSF must share the same CSF and Bloom-filter implementation so
that only the decision rule differs. VL-BuRR uses its native implementation.

The VL-BuRR evaluation path is based on `vihan-lakshman/autocsf-bench`, but the
small required patch and runner live directly in this repository. Do not add a
runtime or submodule dependency on that separate evaluation repository.
