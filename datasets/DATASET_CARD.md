---
pretty_name: AutoCSF Genomics k-mer Count Tables
license: other
task_categories:
  - tabular-classification
---

# AutoCSF genomics datasets

This artifact contains the original public FASTA/FASTQ inputs and deterministic
processed tables used by the AutoCSF genomics experiments. Each processed file
is a Zstandard-compressed, headerless UTF-8 TSV containing
`15-mer<TAB>occurrence-count`, sorted bytewise by k-mer.

The bundled dense 2-bit counter counts overlapping forward-strand 15-mers,
normalizes ASCII `acgt` to uppercase, keeps counts of one, and omits windows
containing other bases.
`manifests/` records accessions, raw and processed SHA-256 checksums,
row counts, modal count, alpha, and the number of distinct count values.

Sources are NCBI RefSeq for E. coli Sakai, C. elegans, and rice, and ENA for
SRR10211353. The older BCSF paper called the Sakai dataset `B000007`; its
resolved RefSeq assembly is `GCF_000008865.2`.
