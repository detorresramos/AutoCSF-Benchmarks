#!/usr/bin/env python3
"""Benchmark the four camera-ready methods on processed genomics tables."""

import argparse
import csv
import json
import os
from pathlib import Path
import platform
import math
import statistics
import subprocess
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from methods.decision_rules import CSFFilter
from common.bcsf import shibuya_bloom_params
from common.theory import best_discrete_bloom_all_k, nearest_bloom_config

METHODS = ("hkp", "bcsf", "vlburr", "autocsf")
DATASETS = ("ecoli_sakai", "srr10211353", "celegans", "rice")


def rows(path):
    process = subprocess.Popen(["zstd", "-dc", str(path)], stdout=subprocess.PIPE, text=True)
    assert process.stdout is not None
    for line in process.stdout:
        key, value = line.rstrip("\n").split("\t")
        yield key, int(value)
    if process.wait():
        raise RuntimeError(f"zstd failed for {path}")


def load_table(name):
    path = ROOT / "data" / "processed" / f"{name}_k15.tsv.zst"
    if not path.exists():
        raise FileNotFoundError(path)
    keys, values = zip(*rows(path))
    return list(keys), np.asarray(values, dtype=np.uint32)


def modal_fraction(values):
    _, counts = np.unique(values, return_counts=True)
    return float(counts.max() / len(values))


def local_method(name):
    return {
        "hkp": lambda: CSFFilter("bloom", "hkp"),
        "bcsf": lambda: CSFFilter("bloom", "shibuya"),
        "autocsf": lambda: CSFFilter("bloom", "optimal"),
    }[name]


def query_batch(method, structure, keys, seed=42, n=10000):
    rng = np.random.default_rng(seed)
    sample = [keys[i] for i in rng.integers(0, len(keys), size=n)]
    for key in sample[:100]:
        method.query(structure, key)
    started = time.perf_counter_ns()
    for key in sample:
        method.query(structure, key)
    return (time.perf_counter_ns() - started) / n


def run_local(name, keys, values, repeats):
    builds, queries, sizes, params = [], [], [], None
    for repetition in range(repeats):
        method = local_method(name)()
        started = time.perf_counter()
        structure = method.construct(keys, values)
        builds.append(time.perf_counter() - started)
        memory = method.measure_memory_from_structure(structure)
        sizes.append(memory["serialized"])
        queries.append(query_batch(method, structure, keys, seed=42 + repetition))
        params = method.get_params()
        for key, expected in zip(keys[:1000], values[:1000]):
            if structure.query(key) != int(expected):
                raise AssertionError(f"{name} returned a wrong value for {key}")
    return {
        "serialized_bytes": int(statistics.median(sizes)),
        "build_seconds": statistics.median(builds),
        "query_ns": statistics.median(queries),
        "parameters": params,
        "repetitions": repeats,
    }


def run_vlburr(dataset, repeats):
    command = ROOT / "methods" / "vlburr" / "run_genomics.sh"
    result = subprocess.run(
        [str(command), dataset, str(repeats)], capture_output=True, text=True, check=True
    )
    parsed = [json.loads(line[5:]) for line in result.stdout.splitlines() if line.startswith("JSON ")]
    if len(parsed) != 1:
        raise RuntimeError(f"expected one VL-BuRR JSON record:\n{result.stdout}\n{result.stderr}")
    return parsed[0]


def dataset_profile(dataset):
    manifest = json.loads((ROOT / "data/manifests" / f"{dataset}_k15.json").read_text())
    histogram = {int(value): int(frequency) for value, frequency in manifest["histogram"].items()}
    n = manifest["records"]
    entropy = -sum((frequency / n) * math.log2(frequency / n) for frequency in histogram.values())
    return manifest, entropy


def bloom_parameters(method, manifest, entropy):
    alpha = manifest["alpha"]
    if method == "hkp":
        if alpha <= 0.63:
            return None
        return nearest_bloom_config(1 - alpha)
    if method == "bcsf":
        return shibuya_bloom_params(alpha, entropy)
    bpe, hashes, bound = best_discrete_bloom_all_k(
        alpha, manifest["distinct_values"] / manifest["records"]
    )
    return (bpe, hashes) if bound > 0 else None


def plain_table(dataset):
    compressed = ROOT / "data/processed" / f"{dataset}_k15.tsv.zst"
    output = ROOT / "data/cache/tables" / f"{dataset}_k15.tsv"
    output.parent.mkdir(parents=True, exist_ok=True)
    if not output.exists():
        subprocess.run(["zstd", "-d", "-f", str(compressed), "-o", str(output)], check=True)
    return output


def run_native_caramel(dataset, method, repeats, manifest, entropy):
    binary = ROOT / "data/cache/bin/caramel_bench"
    if not binary.exists():
        subprocess.run([str(ROOT / "experiments/genomics_comparison/build_native.sh")], check=True)
    params = None if method == "plain" else bloom_parameters(method, manifest, entropy)
    kind, bpe, hashes = ("none", 0, 0) if params is None else ("bloom", *params)
    result = subprocess.run(
        [str(binary), str(plain_table(dataset)), kind, str(bpe), str(hashes), str(repeats)],
        capture_output=True, text=True, check=True,
    )
    measured = json.loads(result.stdout.strip().splitlines()[-1])
    measured["parameters"] = None if params is None else {"bloom_bits_per_element": bpe, "bloom_num_hashes": hashes}
    return measured


def environment():
    def output(command):
        try:
            return subprocess.check_output(command, text=True, stderr=subprocess.DEVNULL).strip()
        except (OSError, subprocess.CalledProcessError):
            return "unavailable"
    return {
        "platform": platform.platform(),
        "python": platform.python_version(),
        "git_commit": output(["git", "rev-parse", "HEAD"]),
        "carameldb_commit": output(["git", "-C", "deps/CaramelDB", "rev-parse", "HEAD"]),
        "vlburr_commit": (ROOT / "methods/vlburr/UPSTREAM_COMMIT").read_text().strip(),
        "docker": os.environ.get("AUTOCSF_DOCKER_DESCRIPTION", "not recorded"),
    }


def write_tables(records, output):
    """Write genomics.csv, which verify.py and export_baseline.py read, and genomics.md."""
    output.mkdir(parents=True, exist_ok=True)
    fields = ["dataset", "n", "alpha", "method", "bits_per_key", "bits_saved_vs_plain", "build_seconds", "query_ns"]
    with (output / "genomics.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows({key: row.get(key) for key in fields} for row in records)

    def number(row, key, digits):
        value = row.get(key)
        return f"{value:.{digits}f}" if value is not None else row.get("status", "N/A")

    def selected_filter(row):
        params = row.get("parameters")
        if not params:
            return "no filter"
        if row["method"] == "vlburr":
            return "no filter (Opt chose empty)" if row["bits_saved_vs_plain"] == 0 else "VL-BuRR optimized filter"
        return f"Bloom b={params['bloom_bits_per_element']}, k={params['bloom_num_hashes']}"

    md = [
        "Bits per key saved by each method's selected filter, relative to that",
        "method's own unfiltered CSF. VL-BuRR uses its Huffman CSF; the other",
        "three share Caramel/GOV, so only their decision rule differs.",
        "",
        "| Dataset | N | alpha | Method | Plain bpk | Filtered bpk | Saved bpk | Selected filter | build (s) | query (ns) |",
        "|---|---:|---:|---|---:|---:|---:|---|---:|---:|",
    ]
    for row in records:
        plain = row.get("baseline_bits_per_key")
        if plain is None and row.get("bits_per_key") is not None:
            plain = row["bits_per_key"] + row["bits_saved_vs_plain"]
        md.append(
            f"| {row['dataset']} | {row['n']:,} | {row['alpha']:.4f} | {row['method']} | "
            f"{plain:.4f} | {number(row, 'bits_per_key', 4)} | {number(row, 'bits_saved_vs_plain', 4)} | "
            f"{selected_filter(row)} | {number(row, 'build_seconds', 4)} | {number(row, 'query_ns', 2)} |"
        )
    (output / "genomics.md").write_text("\n".join(md) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", nargs="+", choices=DATASETS, default=list(DATASETS))
    parser.add_argument("--method", nargs="+", choices=METHODS, default=list(METHODS))
    parser.add_argument("--repetitions", type=int, default=1)
    parser.add_argument("--output", type=Path, default=ROOT / "results/genomics_comparison")
    parser.add_argument("--plain-results", type=Path,
                        help="reuse plain-table byte counts from an earlier JSON run")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    reused_plain = {}
    if args.plain_results:
        prior = json.loads(args.plain_results.read_text())
        reused_plain = {
            row["dataset"]: row["plain_serialized_bytes"]
            for row in prior["records"]
        }
    records = []
    for dataset in args.dataset:
        manifest, entropy = dataset_profile(dataset)
        if dataset in reused_plain:
            plain_bytes = reused_plain[dataset]
        else:
            plain = run_native_caramel(dataset, "plain", args.repetitions, manifest, entropy)
            plain_bytes = plain["serialized_bytes"]
        for method in args.method:
            measured = run_vlburr(dataset, args.repetitions) if method == "vlburr" else run_native_caramel(dataset, method, args.repetitions, manifest, entropy)
            bpk = measured.get("bits_per_key", measured["serialized_bytes"] * 8 / manifest["records"])
            baseline_bytes = measured.get("plain_serialized_bytes", plain_bytes) if method == "vlburr" else plain_bytes
            savings = measured.get("bits_saved_vs_plain", (baseline_bytes - measured["serialized_bytes"]) * 8 / manifest["records"])
            records.append({
                "dataset": dataset, "n": manifest["records"], "alpha": manifest["alpha"], "method": method,
                "bits_per_key": bpk, "bits_saved_vs_plain": savings,
                "build_seconds": measured["build_seconds"], "query_ns": measured["query_ns"],
                "serialized_bytes": measured["serialized_bytes"], "plain_serialized_bytes": baseline_bytes,
                "baseline_bits_per_key": measured.get("baseline_bits_per_key", baseline_bytes * 8 / manifest["records"]),
                "parameters": measured.get("parameters"),
                "repetitions": measured.get("repetitions", args.repetitions),
            })
    payload = {"schema_version": 1, "environment": environment(), "records": records}
    (args.output / "genomics.json").write_text(json.dumps(payload, indent=2) + "\n")
    write_tables(records, args.output)


if __name__ == "__main__":
    main()
