# Genomics datasets

The Hugging Face **dataset repository** (informally called the data space in
project discussion) is the public home for raw inputs, processed tables, and
their manifests. It is separate from this code repository so a fresh clone
does not contain several gigabytes of data.

Proposed repository: `detorresramos/autocsf-genomics`.

## Datasets

| name | source | expected role |
|---|---|---|
| `ecoli_sakai` | E. coli Sakai source used by Shibuya et al. | extreme skew |
| `srr10211353` | SRA run SRR10211353 | low skew |
| `celegans` | RefSeq GCF_000002985.6 (WBcel235) | moderate skew |
| `rice` | RefSeq GCF_001433935.1 (IRGSP-1.0) | dead-zone, alpha about 0.727 |

## Published layout

```text
README.md
raw/<source files>
processed/<name>_k15.tsv.zst
manifests/<name>_k15.json
```

Each manifest records source accession and checksum, generation command,
generator commit, processed checksum, number of distinct k-mers, modal count,
majority fraction, and count vocabulary size.

## Counting semantics

- `k = 15`
- forward strand only (`canonical = false`)
- overlapping k-mers
- records do not cross FASTA/FASTQ sequence boundaries
- uppercase `A`, `C`, `G`, and `T` only
- any window containing another character is omitted
- minimum count is one

These semantics reproduce the `locom.py count` path from Shibuya et al.'s
artifact. `locom` is provenance, not a runtime dependency.

## Commands

```bash
python datasets/manage.py list
python datasets/manage.py download
python datasets/manage.py validate
python datasets/manage.py generate --dataset rice
```

Generation expects the corresponding raw files under `data/raw/`, a C++17
compiler, and `zstd`. `datasets/count_kmers.cpp` uses a dense 2-bit k=15 table,
counts forward-strand windows, and emits keys in lexical order. The generator
uses about 4 GiB of RAM regardless of genome size.
