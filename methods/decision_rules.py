import os
import sys
import tempfile

import numpy as np

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


EPSILON_STRATEGIES = ("optimal", "shibuya", "hkp")

# Hreinsson, Kroyer and Pagh, "Storing a compressed function with constant time
# access" (2009), bound the redundancy of the filtered construction by
#     r < min(alpha + 0.086, 1.82 * (1 - alpha))
# where the first term is Gallager's Huffman bound (the no-filter regime) and the
# second is the filtered regime. HKP's decision criterion is to filter beyond the
# crossover of the two. Solving alpha + 0.086 == 1.82 * (1 - alpha) gives 0.6149;
# the paper rounds this to "roughly 0.63".
HKP_CROSSOVER_ALPHA = (1.82 - 0.086) / (1.0 + 1.82)


def _find_hkp_params(keys, values):
    """Realize HKP's epsilon=1-alpha rule with the nearest Bloom config."""
    _, compute_actual_alpha, _, _, _ = _import_shared()
    alpha = compute_actual_alpha(values)
    if alpha <= HKP_CROSSOVER_ALPHA:
        return None
    target = max(1e-12, 1.0 - alpha)
    candidates = []
    for bits_per_element in range(1, 17):
        for num_hashes in range(1, 9):
            epsilon = (1.0 - np.exp(-num_hashes / bits_per_element)) ** num_hashes
            candidates.append((abs(np.log(epsilon) - np.log(target)), bits_per_element, num_hashes))
    _, bits_per_element, num_hashes = min(candidates)
    return {"bloom_bits_per_element": bits_per_element, "bloom_num_hashes": num_hashes}


def _find_optimal_params(filter_type, keys, values):
    theory, compute_actual_alpha, _, _, _ = _import_shared()
    n = len(keys)
    alpha = compute_actual_alpha(values)
    n_over_N = len(np.unique(values)) / n
    n_filter = int(n * (1 - alpha))

    if filter_type == "xor":
        bits, bound = theory.best_discrete_xor(alpha, n_over_N)
        return {"fingerprint_bits": bits} if bound > 0 else None
    elif filter_type == "binary_fuse":
        bits, bound = theory.best_discrete_binary_fuse(alpha, n_over_N, n_filter)
        return {"fingerprint_bits": bits} if bound > 0 else None
    elif filter_type == "bloom":
        bpe, k, bound = theory.best_discrete_bloom_all_k(alpha, n_over_N)
        return {"bloom_bits_per_element": bpe, "bloom_num_hashes": k} if bound > 0 else None
    else:
        raise ValueError(f"Unknown filter type: {filter_type}")


def _find_shibuya_params(keys, values):
    _, compute_actual_alpha, _, empirical_entropy, shibuya_bloom_params = _import_shared()
    alpha = compute_actual_alpha(values)
    H0 = empirical_entropy(values)
    result = shibuya_bloom_params(alpha, H0)
    if result is None:
        return None
    bits_per_element, num_hashes = result
    return {"bloom_bits_per_element": bits_per_element, "bloom_num_hashes": num_hashes}


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
        if self.epsilon_strategy == "optimal":
            self._params = _find_optimal_params(self.filter_type, keys, values)
        elif self.epsilon_strategy == "shibuya":
            self._params = _find_shibuya_params(keys, values)
        else:
            self._params = _find_hkp_params(keys, values)
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
