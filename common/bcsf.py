"""Shibuya's cost model for CSF+Bloom filter epsilon derivation.

Implements the epsilon selection from Shibuya et al., "Space-efficient
representation of genomic k-mer count tables" (shibuya2022space):
    epsilon* = C_BF * (1 - alpha) / (C_CSF * alpha * ln2)

where C_BF = 1.44 (idealized Bloom constant) and C_CSF is their data-driven
estimate of CSF cost as a function of the empirical entropy H0. Both the
constant and the piecewise C_CSF fit below are taken from that paper; they are
reproduced here so BCSF's decisions can be compared against AutoCSF's on
identical inputs, not because this repository endorses the model.

The canonical Bloom filter construction for a target epsilon is:
    bits_per_element = ceil(log2(1/eps) * log2(e))
    num_hashes = round(bits_per_element * ln(2))
"""

import math

import numpy as np

C_BF = 1.44  # idealized Bloom filter constant (log2(e))


def empirical_entropy(values):
    """Compute H0 (zeroth-order empirical entropy) from value distribution."""
    _, counts = np.unique(values, return_counts=True)
    probs = counts / counts.sum()
    return -np.sum(probs * np.log2(probs))


def shibuya_csf_cost(H0):
    """Piecewise C_CSF cost function from Shibuya's model.

    Reproduced verbatim from shibuya2022space:

        C_CSF = 0.22 H0^2 + 0.18 H0 + 1.16   if H0 < 2
                1.1 H0 + 0.2                 otherwise

    These coefficients are an empirical fit reported in that paper, not a
    derivation; that they are a heuristic rather than a principled model is
    precisely what AutoCSF sets out to improve on.
    """
    if H0 < 2:
        return 0.22 * H0**2 + 0.18 * H0 + 1.16
    return 1.1 * H0 + 0.2


def shibuya_optimal_epsilon(alpha, H0):
    """epsilon* = C_BF * (1 - alpha) / (C_CSF * alpha * ln2)"""
    C_CSF = shibuya_csf_cost(H0)
    return C_BF * (1 - alpha) / (C_CSF * alpha * math.log(2))


def shibuya_bloom_params(alpha, H0):
    """Compute Bloom filter parameters from Shibuya's optimal epsilon.

    Uses the canonical Bloom filter construction:
        bits_per_element = ceil(1.44 * log2(1/eps))
        num_hashes = round(bits_per_element * ln(2))

    Returns (bits_per_element, num_hashes) or None if eps* >= 1 (no filter).
    """
    eps_star = shibuya_optimal_epsilon(alpha, H0)
    if eps_star >= 1:
        return None
    bits_per_element = max(1, math.ceil(C_BF * math.log2(1 / eps_star)))
    num_hashes = max(1, round(bits_per_element * math.log(2)))
    return bits_per_element, num_hashes
