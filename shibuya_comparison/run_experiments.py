"""Data generation for prior methods comparison. Builds CSFs and saves results to JSON.

For each (alpha, distribution), three methods decide:
  1. Whether to use a prefilter at all
  2. If yes, what Bloom parameters (bits_per_element, num_hashes) to use

We build each recommendation and measure actual bits/key.

Methods compared:
  - AutoCSF (ours): maximizes lower bound over discrete Bloom configs
  - BCSF (Shibuya et al.): heuristic CSF cost model with idealized Bloom constant
  - LSF (Hermann et al.): per-node integer-r optimization for VL-BuRR fingerprints,
    mapped to the smallest Bloom config achieving the target FPR. Note this gives
    LSF a slight disadvantage relative to its native VL-BuRR instantiation, since
    Bloom filters have ~44% higher per-key overhead than VL-BuRR at the same FPR.
    This isolates LSF's decision rule from its native filter backend.

Usage:
    python shibuya_comparison/run_experiments.py
"""
import argparse
import json
import math
import os
import sys

import carameldb
import numpy as np
from carameldb import BloomFilterConfig
from tqdm import tqdm

_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_dir, ".."))

from shared.data_gen import compute_actual_alpha, gen_alpha_values, gen_keys
from shared.measure import measure_csf
from shared.shibuya import C_BF, empirical_entropy, shibuya_bloom_params
from shared.theory import best_discrete_bloom_all_k

SEED = 42
ALPHAS = list(np.arange(0.50, 1.00, 0.01))
DEFAULT_DISTS = ["unique", "zipfian", "uniform_100"]

FIGURES_DIR = os.path.join(_dir, "figures")
DATA_DIR = os.path.join(FIGURES_DIR, "data")


def build_and_measure_bpk(keys, values, prefilter=None):
    structure = carameldb.Caramel(keys, values, prefilter=prefilter, verbose=False)
    return structure.get_stats().in_memory_bytes * 8 / len(keys)


def our_recommendation(alpha, n_over_N):
    """AutoCSF: pick the (bits_per_element, num_hashes) maximizing the lower bound.

    Returns None if no config has positive lower bound.
    """
    bpe, k, lb = best_discrete_bloom_all_k(alpha, n_over_N)
    if lb <= 0:
        return None
    return bpe, k


def hkp_recommendation(alpha):
    """HKP: activate filter iff alpha > 0.63, with target FPR = 1 - alpha.

    Uses the canonical BCSF Bloom construction:
        bits_per_element = ceil(C_BF * log2(1/eps))
        num_hashes = round(bits_per_element * ln(2))
    """
    if alpha <= 0.63:
        return None
    eps = 1.0 - alpha
    bits_per_element = max(1, math.ceil(C_BF * math.log2(1 / eps)))
    num_hashes = max(1, round(bits_per_element * math.log(2)))
    return bits_per_element, num_hashes


def lsf_recommendation(alpha):
    """LSF (Hermann et al. 2025): per-node integer-r optimization.

    LSF activates the filter when the rare-branch probability p = 1 - alpha
    satisfies p < 1/3, i.e., when alpha > 2/3. Given activation, it picks
    integer r >= 1 minimizing space(p, r) = p*r + p + (1-p)*2^(-r), where r is
    the number of filter bits per stored key and epsilon = 2^(-r). This rule
    is derived for VL-BuRR fingerprint filters with cost b(epsilon) =
    log_2(1/epsilon) bits per stored key.

    To compare on the same filter substrate as BCSF and AutoCSF, we map the
    target FPR = 2^(-r) to the smallest Bloom configuration (bits_per_element,
    num_hashes) achieving FPR <= 2^(-r). This isolates LSF's decision rule
    from its native VL-BuRR backend.

    Returns None if LSF would not activate the filter, else (bpe, num_hashes).
    """
    p = 1.0 - alpha
    # Activation threshold from LSF: alpha > 2/3 equivalently p < 1/3.
    if p >= 1.0 / 3.0:
        return None

    # Find optimal integer r >= 1 minimizing space(p, r).
    # space(p, r) is unimodal in r, so we can stop once cost starts increasing.
    best_r = None
    best_cost = float("inf")
    for r in range(1, 20):
        cost = p * r + p + (1 - p) * (2 ** -r)
        if cost < best_cost:
            best_cost = cost
            best_r = r
        else:
            break

    target_epsilon = 2.0 ** -best_r

    # Map target FPR to smallest Bloom config achieving FPR <= target_epsilon.
    # Bloom FPR with k hashes and bpe bits per element: (1 - exp(-k/bpe))^k.
    # We minimize total filter cost (1-alpha) * bpe over feasible configs.
    best_bloom = None
    best_bloom_size = float("inf")
    for bpe in range(1, 16):
        for k in range(1, 8):
            fpr = (1.0 - math.exp(-k / bpe)) ** k
            if fpr <= target_epsilon:
                size_bpk = (1.0 - alpha) * bpe
                if size_bpk < best_bloom_size:
                    best_bloom_size = size_bpk
                    best_bloom = (bpe, k)

    return best_bloom


def run_experiments(n, dists, skip_lsf):
    results = {}
    for dist in dists:
        rows = []
        for alpha in tqdm(ALPHAS, desc=f"prior methods {dist}"):
            keys = gen_keys(n)
            values = gen_alpha_values(n, alpha, seed=SEED, minority_dist=dist)
            actual_alpha = compute_actual_alpha(values)
            H0 = empirical_entropy(values)

            baseline = measure_csf(keys, values, "none", minority_dist=dist)
            n_over_N = baseline.huffman_num_symbols / n
            baseline_bpk = baseline.bits_per_key

            our_rec = our_recommendation(actual_alpha, n_over_N)
            shib_rec = shibuya_bloom_params(actual_alpha, H0)
            hkp_rec = hkp_recommendation(actual_alpha)
            lsf_rec = None if skip_lsf else lsf_recommendation(actual_alpha)

            if our_rec is None:
                our_bpk = baseline_bpk
            else:
                config = BloomFilterConfig(bits_per_element=our_rec[0], num_hashes=our_rec[1])
                our_bpk = build_and_measure_bpk(keys, values, prefilter=config)

            if shib_rec is None:
                shib_bpk = baseline_bpk
            else:
                config = BloomFilterConfig(bits_per_element=shib_rec[0], num_hashes=shib_rec[1])
                shib_bpk = build_and_measure_bpk(keys, values, prefilter=config)

            if hkp_rec is None:
                hkp_bpk = baseline_bpk
            else:
                config = BloomFilterConfig(bits_per_element=hkp_rec[0], num_hashes=hkp_rec[1])
                hkp_bpk = build_and_measure_bpk(keys, values, prefilter=config)

            if lsf_rec is None:
                lsf_bpk = baseline_bpk
            else:
                config = BloomFilterConfig(bits_per_element=lsf_rec[0], num_hashes=lsf_rec[1])
                lsf_bpk = build_and_measure_bpk(keys, values, prefilter=config)

            rows.append(
                {
                    "alpha": round(actual_alpha, 4),
                    "requested_alpha": round(alpha, 4),
                    "n_over_N": n_over_N,
                    "H0": H0,
                    "baseline_bpk": baseline_bpk,
                    "our_bpk": our_bpk,
                    "our_params": list(our_rec) if our_rec else None,
                    "shib_bpk": shib_bpk,
                    "shib_params": list(shib_rec) if shib_rec else None,
                    "hkp_bpk": hkp_bpk,
                    "hkp_params": list(hkp_rec) if hkp_rec else None,
                    "lsf_bpk": lsf_bpk,
                    "lsf_params": list(lsf_rec) if lsf_rec else None,
                }
            )

            print(
                f"  {dist} a={actual_alpha:.4f}: baseline={baseline_bpk:.2f}  "
                f"ours={our_bpk:.2f} {our_rec}  "
                f"shibuya={shib_bpk:.2f} {shib_rec}  "
                f"hkp={hkp_bpk:.2f} {hkp_rec}  "
                f"lsf={lsf_bpk:.2f} {lsf_rec}"
            )

        results[dist] = rows
    return results


def save_json(data, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"Saved: {path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, default=100_000, help="Number of keys")
    parser.add_argument(
        "--dists",
        default=",".join(DEFAULT_DISTS),
        help="Comma-separated distributions to run",
    )
    parser.add_argument(
        "--output",
        default="prior_methods_comparison.json",
        help="Output JSON filename (under figures/data/)",
    )
    parser.add_argument("--skip-lsf", action="store_true", help="Skip LSF builds")
    args = parser.parse_args()

    dists = args.dists.split(",")
    os.makedirs(DATA_DIR, exist_ok=True)
    results = run_experiments(args.n, dists, args.skip_lsf)
    data = {
        "N": args.n,
        "seed": SEED,
        "distributions": results,
    }
    save_json(data, os.path.join(DATA_DIR, args.output))


if __name__ == "__main__":
    main()
