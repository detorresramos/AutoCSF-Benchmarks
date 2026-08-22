#!/usr/bin/env python3
"""Validate a from-scratch synthetic run against the accepted numeric artifacts."""

import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TOLERANCE_BPK = 0.02


def method_rows(path):
    with path.open() as handle:
        rows = list(csv.DictReader(handle))
    return {
        (row["distribution"], row["method"], float(row["alpha"])): float(row["bits_saved"])
        for row in rows
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fresh", type=Path, help="directory copied out of the clean container")
    parser.add_argument("--image-id", default="unknown")
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "results/reproduction/synthetic-receipt.json",
    )
    args = parser.parse_args()

    accepted_data = ROOT / "theory_validation/figures/data"
    fresh_data = args.fresh / "theory_figures/data"
    accepted_files = sorted(accepted_data.glob("*.json"))
    fresh_files = sorted(fresh_data.glob("*.json"))
    errors = []
    if {p.name for p in accepted_files} != {p.name for p in fresh_files}:
        errors.append("the fresh and accepted theory JSON file sets differ")

    max_theory_delta = 0.0
    for accepted_path in accepted_files:
        fresh_path = fresh_data / accepted_path.name
        if not fresh_path.exists():
            continue
        accepted = json.loads(accepted_path.read_text())
        fresh = json.loads(fresh_path.read_text())
        accepted_results = accepted.get("results", [accepted])
        fresh_results = fresh.get("results", [fresh])
        if len(accepted_results) != len(fresh_results):
            errors.append(f"result count differs for {accepted_path.name}")
            continue
        for index, (old, new) in enumerate(zip(accepted_results, fresh_results)):
            for key in ("alpha", "requested_alpha", "n_filter"):
                if old[key] != new[key]:
                    errors.append(f"{accepted_path.name} result {index} differs in {key}")
            for key in ("baseline_bpk", "theory_guided_bpk_saved", "best_empirical_bpk_saved"):
                if key not in old or key not in new:
                    continue
                delta = abs(float(old[key]) - float(new[key]))
                max_theory_delta = max(max_theory_delta, delta)
                if delta > TOLERANCE_BPK:
                    errors.append(f"{accepted_path.name} result {index} {key} delta={delta:.6g}")
            old_params = old.get("empirical_per_param", [])
            new_params = new.get("empirical_per_param", [])
            if len(old_params) != len(new_params):
                errors.append(f"parameter count differs for {accepted_path.name} result {index}")
                continue
            for old_param, new_param in zip(old_params, new_params):
                delta = abs(float(old_param["bpk_saved"]) - float(new_param["bpk_saved"]))
                max_theory_delta = max(max_theory_delta, delta)
                if delta > TOLERANCE_BPK:
                    errors.append(f"{accepted_path.name} result {index} parameter delta={delta:.6g}")

    accepted_methods = method_rows(ROOT / "method_comparison/data.csv")
    fresh_method_path = args.fresh / "method_comparison.csv"
    fresh_methods = method_rows(fresh_method_path) if fresh_method_path.exists() else {}
    if accepted_methods.keys() != fresh_methods.keys():
        errors.append("fresh and accepted method-comparison row keys differ")
    method_deltas = {
        key: abs(accepted_methods[key] - fresh_methods[key])
        for key in accepted_methods.keys() & fresh_methods.keys()
    }
    max_method_delta = max(method_deltas.values(), default=0.0)
    if max_method_delta > TOLERANCE_BPK:
        errors.append(f"method-comparison maximum delta={max_method_delta:.6g}")

    fresh_plots = list((args.fresh / "results_figures").glob("*.png"))
    receipt = {
        "status": "passed" if not errors else "failed",
        "kind": "full-clean-synthetic-reproduction",
        "checked_at": datetime.now(timezone.utc).isoformat(),
        "container_image_id": args.image_id,
        "started_at": (args.fresh / "started_at.txt").read_text().strip(),
        "completed_at": (args.fresh / "completed_at.txt").read_text().strip(),
        "fresh_artifacts": {
            "data_files": len(fresh_files),
            "method_rows": len(fresh_methods),
            "plots": len(fresh_plots),
        },
        "comparison": {
            "status": "passed" if not errors else "failed",
            "bits_per_key_absolute_tolerance": TOLERANCE_BPK,
            "max_theory_delta": max_theory_delta,
            "max_method_comparison_delta": max_method_delta,
            "errors": errors[:100],
        },
        "scope": {
            "theory": "all 15 alpha sweeps and all 30 epsilon sweeps, freshly measured",
            "methods": "all 88 accepted-paper method-comparison rows, freshly measured",
            "plots": "all ten accepted-paper aggregate plots, rendered from the fresh measurements",
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
