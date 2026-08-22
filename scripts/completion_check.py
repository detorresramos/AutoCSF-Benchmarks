#!/usr/bin/env python3
"""Single source of truth for `make verify` and the Codex Stop hook."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
STATUS = ROOT / ".codex" / "reproduction-status.json"
DATASETS = {"ecoli_sakai", "srr10211353", "celegans", "rice"}
METHODS = {"hkp", "bcsf", "vlburr", "autocsf"}
PLOTS = {
    "alpha_sweep_uniform_100.png", "alpha_sweep_unique.png", "alpha_sweep_zipfian.png",
    "epsilon_sweep_uniform_100_alpha0.7.png", "epsilon_sweep_uniform_100_alpha0.9.png",
    "epsilon_sweep_unique_alpha0.7.png", "epsilon_sweep_unique_alpha0.9.png",
    "epsilon_sweep_zipfian_alpha0.7.png", "epsilon_sweep_zipfian_alpha0.9.png",
    "method-comparison.png",
}
THEORY_FILTERS = {"xor", "binary_fuse", "bloom_k1", "bloom_k2", "bloom_k3"}
THEORY_DISTRIBUTIONS = {"unique", "zipfian", "uniform_100"}


def digest(path):
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def inspect():
    errors = []
    table = ROOT / "results" / "genomics.csv"
    if not table.exists():
        errors.append("missing results/genomics.csv")
    else:
        with table.open() as handle:
            rows = list(csv.DictReader(handle))
        pairs = {(row.get("dataset"), row.get("method")) for row in rows}
        expected = {(dataset, method) for dataset in DATASETS for method in METHODS}
        if len(rows) != 16 or pairs != expected:
            errors.append(f"genomics table is not the expected 4x4 matrix ({len(rows)} rows)")
        numeric = ("n", "alpha", "bits_per_key", "bits_saved_vs_plain", "build_seconds", "query_ns")
        for index, row in enumerate(rows, 2):
            try:
                if any(not float(row[key]) == float(row[key]) for key in numeric):
                    raise ValueError
            except (KeyError, ValueError):
                errors.append(f"invalid numeric value in genomics.csv row {index}")
        payload_path = ROOT / "results" / "genomics.json"
        if payload_path.exists():
            payload = json.loads(payload_path.read_text())
            for row in payload.get("records", []):
                if row.get("method") == "vlburr":
                    baseline = row.get("baseline_bits_per_key")
                    if baseline is None or abs((baseline - row["bits_per_key"]) - row["bits_saved_vs_plain"]) > 1e-9:
                        errors.append(f"VL-BuRR row uses an invalid baseline: {row.get('dataset')}")
    for suffix in ("json", "md", "tex"):
        if not (ROOT / "results" / f"genomics.{suffix}").exists():
            errors.append(f"missing results/genomics.{suffix}")
    for name in ("genomics-paper.md", "genomics-paper.tex",
                 "genomics-audit.md", "genomics-audit.tex",
                 "genomics-total.md", "genomics-total.tex"):
        if not (ROOT / "results" / name).exists():
            errors.append(f"missing results/{name}")
    actual_plots = {path.name for path in (ROOT / "results" / "figures").glob("*.png")}
    for name in sorted(PLOTS - actual_plots):
        errors.append(f"missing plot: {name}")
    reference = ROOT / "reference" / "accepted-paper"
    sums = reference / "SHA256SUMS"
    if not sums.exists():
        errors.append("missing accepted-paper plot checksums")
    else:
        for line in sums.read_text().splitlines():
            expected_digest, name = line.split(maxsplit=1)
            path = reference / name.strip().lstrip("*")
            if not path.exists() or digest(path) != expected_digest:
                errors.append(f"accepted-paper reference changed: {path.name}")
    expected_alpha = {
        f"alpha_sweep_{filter_name}_{distribution}"
        for filter_name in THEORY_FILTERS for distribution in THEORY_DISTRIBUTIONS
    }
    expected_epsilon = {
        f"epsilon_sweep_{filter_name}_{distribution}_alpha{alpha}"
        for filter_name in THEORY_FILTERS for distribution in THEORY_DISTRIBUTIONS
        for alpha in ("0.7", "0.9")
    }
    data_dir = ROOT / "theory_validation" / "figures" / "data"
    alpha_dir = ROOT / "theory_validation" / "figures" / "alpha_sweep"
    epsilon_dir = ROOT / "theory_validation" / "figures" / "epsilon_sweep"
    for stem in sorted(expected_alpha):
        if not (data_dir / f"{stem}.json").exists() or not (alpha_dir / f"{stem}.png").exists():
            errors.append(f"missing original alpha-sweep artifact: {stem}")
    for stem in sorted(expected_epsilon):
        if not (data_dir / f"{stem}.json").exists() or not (epsilon_dir / f"{stem}.png").exists():
            errors.append(f"missing original epsilon-sweep artifact: {stem}")
    for dataset in DATASETS:
        manifest = ROOT / "data" / "manifests" / f"{dataset}_k15.json"
        processed = ROOT / "data" / "processed" / f"{dataset}_k15.tsv.zst"
        if not manifest.exists() or not processed.exists():
            errors.append(f"missing dataset or manifest: {dataset}")
    receipt = ROOT / "results" / "reproduction" / "receipt.json"
    if not receipt.exists():
        errors.append("missing independent reproduction receipt")
    else:
        try:
            data = json.loads(receipt.read_text())
            if data.get("status") != "passed" or data.get("dataset") != "ecoli_sakai":
                errors.append("independent reproduction receipt did not pass E. coli")
            receipt_methods = {
                item.get("name") if isinstance(item, dict) else item
                for item in data.get("methods", [])
            }
            receipt_plots = {Path(item).name for item in data.get("plots", [])}
            if receipt_methods != METHODS or receipt_plots != PLOTS:
                errors.append("independent reproduction receipt has incomplete methods/plots")
            vlburr = next((item for item in data.get("methods", [])
                           if isinstance(item, dict) and item.get("name") == "vlburr"), None)
            if (not vlburr or vlburr.get("baseline_bits_per_key") is None or
                    abs((vlburr["baseline_bits_per_key"] - vlburr["bits_per_key"])
                        - vlburr["bits_saved_vs_plain"]) > 1e-9):
                errors.append("independent receipt uses an invalid VL-BuRR baseline")
        except (OSError, json.JSONDecodeError):
            errors.append("independent reproduction receipt is invalid JSON")
    synthetic_receipt = ROOT / "results" / "reproduction" / "synthetic-receipt.json"
    if not synthetic_receipt.exists():
        errors.append("missing clean synthetic reproduction receipt")
    else:
        try:
            synthetic = json.loads(synthetic_receipt.read_text())
            expected = {"data_files": 45, "method_rows": 88, "plots": 10}
            observed = synthetic.get("fresh_artifacts", {})
            if synthetic.get("status") != "passed" or any(
                observed.get(key) != value for key, value in expected.items()
            ):
                errors.append("clean synthetic reproduction receipt did not pass")
            if synthetic.get("comparison", {}).get("status") != "passed":
                errors.append("fresh synthetic results were not validated against accepted results")
        except (OSError, json.JSONDecodeError):
            errors.append("clean synthetic reproduction receipt is invalid JSON")
    known_issue = ROOT / "results" / "KNOWN_ISSUES.md"
    if known_issue.exists():
        errors.append("unresolved benchmark correctness issue; see results/KNOWN_ISSUES.md")
    return errors


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--hook", action="store_true")
    parser.add_argument("--mark-complete", action="store_true")
    args = parser.parse_args()
    errors = inspect()
    if args.mark_complete and not errors:
        STATUS.parent.mkdir(parents=True, exist_ok=True)
        STATUS.write_text(json.dumps({"status": "complete"}, indent=2) + "\n")
    if args.hook:
        try:
            hook_input = json.load(sys.stdin)
        except json.JSONDecodeError:
            hook_input = {}
        state = json.loads(STATUS.read_text()) if STATUS.exists() else {"status": "inactive"}
        if state.get("status") == "active" and errors:
            print(json.dumps({"decision": "block", "reason": "AutoCSF reproduction is incomplete: " + "; ".join(errors[:5]) + ". Continue the implementation and run make verify."}))
        else:
            print("{}")
        return
    if errors:
        print("Completion check failed:", file=sys.stderr)
        for error in errors:
            print(f"- {error}", file=sys.stderr)
        raise SystemExit(1)
    print("All AutoCSF reproduction artifacts are complete.")


if __name__ == "__main__":
    main()
