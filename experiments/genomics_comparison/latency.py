#!/usr/bin/env python3
"""Construction and query cost of AutoCSF's filter against a plain CSF."""

import argparse
import json
from pathlib import Path
import statistics
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from experiments.genomics_comparison.run import (
    dataset_profile, environment, plain_table,
)
from methods.decision_rules import select_filter

BINARY = ROOT / "data/cache/bin/caramel_bench"


def run(args, dataset, bpe, hashes):
    if not BINARY.exists():
        subprocess.run([str(ROOT / "experiments/genomics_comparison/build_native.sh")], check=True)
    result = subprocess.run(
        [str(BINARY), str(plain_table(dataset)), "compare", str(bpe), str(hashes),
         str(args.rounds), str(args.batch)],
        capture_output=True, text=True, check=True,
    )
    return json.loads(result.stdout.strip().splitlines()[-1])


def build_seconds(dataset, kind, bpe, hashes):
    result = subprocess.run(
        [str(BINARY), str(plain_table(dataset)), kind, str(bpe), str(hashes), "1", "1"],
        capture_output=True, text=True, check=True,
    )
    return json.loads(result.stdout.strip().splitlines()[-1])["build_seconds"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", nargs="+", default=["celegans", "rice"])
    parser.add_argument("--rounds", type=int, default=31)
    parser.add_argument("--batch", type=int, default=100_000)
    parser.add_argument("--output", type=Path, default=ROOT / "results/genomics_comparison")
    args = parser.parse_args()

    records = []
    for dataset in args.dataset:
        manifest, stats = dataset_profile(dataset)
        selected = select_filter("autocsf", stats)
        if selected is None:
            raise SystemExit(f"AutoCSF declines to filter on {dataset}; nothing to compare")
        bpe, hashes = selected["bloom_bits_per_element"], selected["bloom_num_hashes"]
        n = manifest["records"]
        paired = run(args, dataset, bpe, hashes)
        records.append({
            "dataset": dataset,
            "n": n,
            "alpha": manifest["alpha"],
            "filter": {"bloom_bits_per_element": bpe, "bloom_num_hashes": hashes},
            "plain_bits_per_key": paired["plain_serialized_bytes"] * 8 / n,
            "filter_bits_per_key": paired["filter_serialized_bytes"] * 8 / n,
            "plain_build_seconds": build_seconds(dataset, "none", 0, 0),
            "filter_build_seconds": build_seconds(dataset, "bloom", bpe, hashes),
            "query": paired,
        })

    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "latency.json").write_text(
        json.dumps({"schema_version": 2, "environment": environment(), "records": records},
                   indent=2) + "\n")

    md = [
        "AutoCSF's selected filter against the plain CSF it is built from.",
        "Query latency is a paired measurement: both indexes are built in one",
        "process and answer the same key sample each round.",
        "",
        "| Dataset | N | Construction (s) | | Query (ns) | | Filter faster |",
        "|---|---:|---:|---:|---:|---:|---:|",
        "| | | plain | filter | plain | filter | (of rounds) |",
    ]
    for row in records:
        q = row["query"]
        md.append(
            f"| {row['dataset']} | {row['n']:,} | {row['plain_build_seconds']:.1f} | "
            f"{row['filter_build_seconds']:.1f} | {q['plain_query_ns']:.1f} | "
            f"{q['filter_query_ns']:.1f} | {q['filter_faster_rounds']}/{q['rounds']} |"
        )
    (args.output / "latency.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))
    for row in records:
        q = row["query"]
        print(f"  {row['dataset']:10s} paired delta {q['paired_delta_ns']:+.2f} ns "
              f"({100 * q['paired_delta_ns'] / q['plain_query_ns']:+.2f}%), "
              f"{row['plain_bits_per_key']:.4f} -> {row['filter_bits_per_key']:.4f} bits/key")


if __name__ == "__main__":
    main()
