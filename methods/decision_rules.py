from collections import namedtuple
import os
import sys
import tempfile

import numpy as np

Profile = namedtuple("Profile", "alpha n_over_N entropy n_filter")

_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_dir, ".."))


def _import_carameldb():
    import carameldb
    return carameldb


def _import_shared():
    from common import theory
    from common.data_generation import compute_actual_alpha
    from common.measurement import create_filter_config
    from common.bcsf import empirical_entropy, shibuya_bloom_params
    return theory, compute_actual_alpha, create_filter_config, empirical_entropy, shibuya_bloom_params


def _csf_stats_to_dict(stats):
    d = {
        "in_memory_bytes": stats.in_memory_bytes,
        "solution_bytes": stats.solution_bytes,
        "filter_bytes": stats.filter_bytes,
        "metadata_bytes": stats.metadata_bytes,
    }
    bs = stats.bucket_stats
    d["bucket_stats"] = {
        "num_buckets": bs.num_buckets,
        "total_solution_bits": bs.total_solution_bits,
        "avg_solution_bits": bs.avg_solution_bits,
        "min_solution_bits": bs.min_solution_bits,
        "max_solution_bits": bs.max_solution_bits,
    }
    hs = stats.huffman_stats
    d["huffman_stats"] = {
        "num_unique_symbols": hs.num_unique_symbols,
        "max_code_length": hs.max_code_length,
        "avg_bits_per_symbol": hs.avg_bits_per_symbol,
        "code_length_distribution": list(hs.code_length_distribution),
    }
    fs = stats.filter_stats
    if fs is not None:
        d["filter_stats"] = {
            "type": fs.type,
            "size_bytes": fs.size_bytes,
            "num_elements": fs.num_elements,
            "num_hashes": fs.num_hashes,
            "size_bits": fs.size_bits,
            "fingerprint_bits": fs.fingerprint_bits,
        }
    else:
        d["filter_stats"] = None
    return d


# The experiments name the rules after their papers; CSFFilter kept the older
# strategy names, so translate rather than rename the public argument.
METHOD_FOR_STRATEGY = {"optimal": "autocsf", "shibuya": "bcsf", "hkp": "hkp"}
EPSILON_STRATEGIES = tuple(METHOD_FOR_STRATEGY)

# Hreinsson, Kroyer and Pagh (2009), Section 5.1. Not alpha + 0.086 = 1.82(1 - alpha):
# Gallager's bound reduces to r <= alpha above 0.5, which is the only regime in which
# the decision arises.
HKP_CROSSOVER_ALPHA = 0.63


def profile(values):
    """The summary statistics every decision rule needs.

    The rules depend on the value distribution only through these four numbers,
    which is what lets the synthetic experiments (which hold the values in
    memory) and the genomics harness (which reads them from a dataset manifest)
    share one implementation.
    """
    _, compute_actual_alpha, _, empirical_entropy, _ = _import_shared()
    alpha = compute_actual_alpha(values)
    n = len(values)
    return Profile(
        alpha=alpha,
        n_over_N=len(np.unique(values)) / n,
        entropy=empirical_entropy(values),
        n_filter=int(n * (1 - alpha)),
    )


def select_filter(method, stats, filter_type="bloom"):
    """Each method's filter choice, or None when it declines to filter.

    This is the only place the three decision rules are implemented; the
    genomics harness and the synthetic comparison both come through here.
    """
    theory, _, _, _, shibuya_bloom_params = _import_shared()

    if method == "hkp":
        # HKP fixes a false positive rate, not a filter; turning that rate into
        # (bits, hashes) is ours to do. Searching the full grid picks
        # configurations that are dominated -- larger *and* less accurate than
        # an available alternative -- so restrict to the frontier.
        if stats.alpha <= HKP_CROSSOVER_ALPHA:
            return None
        bpe, k = theory.nearest_bloom_config(1.0 - stats.alpha)
        return {"bloom_bits_per_element": bpe, "bloom_num_hashes": k}

    if method == "bcsf":
        chosen = shibuya_bloom_params(stats.alpha, stats.entropy)
        if chosen is None:
            return None
        return {"bloom_bits_per_element": chosen[0], "bloom_num_hashes": chosen[1]}

    if method != "autocsf":
        raise ValueError(f"Unknown method: {method}")

    if filter_type == "xor":
        bits, bound = theory.best_discrete_xor(stats.alpha, stats.n_over_N)
        return {"fingerprint_bits": bits} if bound > 0 else None
    if filter_type == "binary_fuse":
        bits, bound = theory.best_discrete_binary_fuse(
            stats.alpha, stats.n_over_N, stats.n_filter
        )
        return {"fingerprint_bits": bits} if bound > 0 else None
    if filter_type == "bloom":
        bpe, k, bound = theory.best_discrete_bloom_all_k(stats.alpha, stats.n_over_N)
        return {"bloom_bits_per_element": bpe, "bloom_num_hashes": k} if bound > 0 else None
    raise ValueError(f"Unknown filter type: {filter_type}")


class CSFFilter:
    def __init__(self, filter_type="binary_fuse", epsilon_strategy="optimal"):
        if epsilon_strategy not in EPSILON_STRATEGIES:
            raise ValueError(
                f"Unknown epsilon_strategy: {epsilon_strategy}, "
                f"must be one of {EPSILON_STRATEGIES}"
            )
        self.epsilon_strategy = epsilon_strategy
        if epsilon_strategy == "shibuya":
            self.filter_type = "bloom"
        else:
            self.filter_type = filter_type
        self.name = f"csf_filter_{epsilon_strategy}_{self.filter_type}"
        self._params = None

    def construct(self, keys, values):
        carameldb = _import_carameldb()
        _, _, create_filter_config, _, _ = _import_shared()
        self._params = select_filter(
            METHOD_FOR_STRATEGY[self.epsilon_strategy], profile(values), self.filter_type
        )
        config = None if self._params is None else create_filter_config(self.filter_type, **self._params)
        return carameldb.Caramel(keys, values, prefilter=config, verbose=False)

    @staticmethod
    def query(structure, key):
        return structure.query(key)

    @staticmethod
    def measure_memory(keys, values):
        return None

    def measure_memory_from_structure(self, structure):
        stats = structure.get_stats()
        with tempfile.NamedTemporaryFile(suffix=".csf", delete=False) as f:
            tmp_path = f.name
        try:
            structure.save(tmp_path)
            serialized_bytes = os.path.getsize(tmp_path)
        finally:
            os.unlink(tmp_path)
        return {
            "serialized": serialized_bytes,
            "csf_stats": _csf_stats_to_dict(stats),
        }

    def get_params(self):
        return self._params
