#!/usr/bin/env python3
"""Merge native Caramel and x86-container VL-BuRR measurements."""

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from baselines.run_genomics import write_tables
from vlburr.parse_genomics import parse as parse_vlburr


def correct_vlburr(row):
    """Recompute VL-BuRR savings against its own no-filter construction."""
    if row["method"] != "vlburr":
        return
    log = ROOT / "data" / "cache" / "vlburr" / row["dataset"] / "benchmark.log"
    with log.open() as handle:
        measured = parse_vlburr(handle, row["n"])
    for key in ("serialized_bytes", "plain_serialized_bytes", "bits_per_key",
                "baseline_bits_per_key", "bits_saved_vs_plain", "build_seconds",
                "query_ns", "parameters", "repetitions"):
        row[key] = measured[key]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+", help="JSON path, optionally suffixed :REPETITIONS")
    parser.add_argument("--output", type=Path, default=ROOT / "results")
    args = parser.parse_args()
    records = []
    environments = []
    for specification in args.inputs:
        raw_path, separator, raw_repetitions = specification.rpartition(":")
        path = Path(raw_path if separator and raw_repetitions.isdigit() else specification)
        repetitions = int(raw_repetitions) if separator and raw_repetitions.isdigit() else None
        payload = json.loads(path.read_text())
        environment_index = len(environments)
        environment = payload["environment"]
        environment["source_result"] = str(path)
        if repetitions is not None:
            environment["repetitions"] = repetitions
        environments.append(environment)
        for row in payload["records"]:
            correct_vlburr(row)
            row["environment_index"] = environment_index
            if repetitions is not None:
                row["repetitions"] = repetitions
            records.append(row)
    by_pair = {(row["dataset"], row["method"]): row for row in records}
    dataset_order = ("ecoli_sakai", "srr10211353", "celegans", "rice")
    method_order = ("hkp", "bcsf", "vlburr", "autocsf")
    records = [by_pair[(dataset, method)] for dataset in dataset_order for method in method_order]
    expected_datasets = {"ecoli_sakai", "srr10211353", "celegans", "rice"}
    expected_methods = {"hkp", "bcsf", "vlburr", "autocsf"}
    expected = {(dataset, method) for dataset in expected_datasets for method in expected_methods}
    if set(by_pair) != expected:
        missing = sorted(expected - set(by_pair))
        raise SystemExit(f"incomplete 4x4 matrix; missing {missing}")
    args.output.mkdir(parents=True, exist_ok=True)
    payload = {"schema_version": 1, "environments": environments, "records": records}
    (args.output / "genomics.json").write_text(json.dumps(payload, indent=2) + "\n")
    write_tables(records, args.output)
    print(f"Wrote {len(records)} rows to {args.output}")


if __name__ == "__main__":
    main()
