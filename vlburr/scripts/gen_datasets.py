#!/usr/bin/env python3
"""Generate AutoCSF-style synthetic datasets in the .lrbin format consumed by
ribbon_learned_bench.

Value distributions follow the `gen_alpha_values` logic from the companion
filter-optimization code: a single majority value with frequency `alpha`, and
the remaining keys drawn from a configurable `minority_dist`:

  unique       - each minority key gets a unique value (worst case)
  zipfian      - Zipf(s=1.5)
  geometric    - Geometric(p=0.3)
  uniform_10   - uniform over 10 distinct values
  uniform_100  - uniform over 100 distinct values
  two_value    - single other value (best case)

ribbon_learned_bench reads uint16 labels (see dataset_reader.hpp), so the
uint32 values produced by gen_alpha_values are remapped to contiguous uint16
labels. This is lossless for this benchmark: ModelFreq and the CSF depend only
on the multiset of class frequencies, not on value identities. Distributions
that produce more than 65535 distinct values (e.g. `unique` at N=1M) cannot be
represented and are rejected.

File format (matches include/lsf/dataset_reader.hpp):
  <dataset>_X.lrbin: uint64 N, uint64 num_features, then N*num_features float32
  <dataset>_y.lrbin: uint16 n_classes, then N uint16 labels
"""
import argparse, struct
from pathlib import Path
import numpy as np

UINT16_MAX_CLASSES = 65535

MINORITY_DISTRIBUTIONS = [
    "unique",
    "zipfian",
    "geometric",
    "uniform_10",
    "uniform_100",
    "two_value",
]


def _generate_minority_values(num_minority: int, minority_dist: str,
                              rng: np.random.Generator) -> np.ndarray:
    # Values start at 1 to avoid collision with majority_value=2^32-1.
    if minority_dist == "unique":
        return np.arange(1, num_minority + 1, dtype=np.uint32)
    elif minority_dist == "zipfian":
        return rng.zipf(1.5, size=num_minority).astype(np.uint32)
    elif minority_dist == "geometric":
        return rng.geometric(0.3, size=num_minority).astype(np.uint32)
    elif minority_dist == "uniform_10":
        return rng.integers(1, 11, size=num_minority, dtype=np.uint32)
    elif minority_dist == "uniform_100":
        return rng.integers(1, 101, size=num_minority, dtype=np.uint32)
    elif minority_dist == "two_value":
        return np.ones(num_minority, dtype=np.uint32)
    else:
        raise ValueError(f"Unknown minority distribution: {minority_dist}")


def gen_alpha_values(n: int, alpha: float, seed: int,
                     minority_dist: str) -> np.ndarray:
    """Faithful port of the companion code's gen_alpha_values: the most common
    value (2^32-1) has frequency exactly floor(n*alpha); the rest follow
    minority_dist. Values are shuffled."""
    rng = np.random.default_rng(seed)
    num_most_common = int(n * alpha)
    num_minority = n - num_most_common
    majority_value = np.uint32(2**32 - 1)
    values = np.full(n, majority_value, dtype=np.uint32)
    if num_minority > 0:
        values[num_most_common:] = _generate_minority_values(num_minority, minority_dist, rng)
    rng.shuffle(values)
    return values


def remap_to_uint16(values: np.ndarray) -> tuple[np.ndarray, int]:
    """Factorize arbitrary uint32 values into contiguous uint16 labels,
    preserving the frequency multiset. Returns (labels, n_classes)."""
    _, labels = np.unique(values, return_inverse=True)
    n_classes = int(labels.max()) + 1
    if n_classes > UINT16_MAX_CLASSES:
        raise ValueError(
            f"distribution produces {n_classes} distinct values, exceeding the "
            f"uint16 limit of {UINT16_MAX_CLASSES} that ribbon_learned_bench "
            f"supports (label_type=uint16_t). Use a smaller --n or a different "
            f"--minority-dist (e.g. 'unique' is infeasible at this N).")
    return labels.astype(np.uint16), n_classes


def write_dataset(out_dir: Path, name: str, labels: np.ndarray, n_classes: int,
                  features: np.ndarray) -> None:
    N, num_features = features.shape
    (out_dir / f"{name}_X.lrbin").write_bytes(
        struct.pack('<QQ', N, num_features) + features.tobytes()
    )
    (out_dir / f"{name}_y.lrbin").write_bytes(
        struct.pack('<H', n_classes) + labels.tobytes()
    )


def dist_token(minority_dist: str) -> str:
    """Underscore-free token for dataset names (so they stay easy to parse)."""
    return minority_dist.replace("_", "")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, required=True,
                    help="Output directory (will be created if needed).")
    ap.add_argument("--n", type=int, default=1_000_000,
                    help="Number of examples per dataset.")
    ap.add_argument("--minority-dist", nargs="+", default=["uniform_100"],
                    choices=MINORITY_DISTRIBUTIONS,
                    help="Minority value distribution(s) to generate.")
    ap.add_argument("--p-max", "--alpha", type=float, nargs="+", dest="p_max",
                    default=[0.5, 0.6, 0.7, 0.8, 0.9],
                    help="Majority-value frequencies (alpha) to sweep.")
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3, 4],
                    help="Random seeds.")
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    count = 0
    for dist in args.minority_dist:
        for alpha in args.p_max:
            for seed in args.seeds:
                # Match the original per-dataset seeding so runs are reproducible.
                ds_seed = seed * 100 + int(alpha * 100)
                values = gen_alpha_values(args.n, alpha, ds_seed, dist)
                labels, n_classes = remap_to_uint16(values)
                # Single dummy feature: ModelFreq ignores features but the reader expects them.
                features = np.random.default_rng(ds_seed).standard_normal((args.n, 1)).astype(np.float32)
                name = f"acsf_{dist_token(dist)}_p{int(alpha*100):02d}_s{seed}"
                write_dataset(args.out, name, labels, n_classes, features)
                count += 1
    print(f"Wrote {count} datasets to {args.out} "
          f"(N={args.n}, dists={args.minority_dist}, alpha={args.p_max}, seeds={args.seeds})")


if __name__ == "__main__":
    main()
