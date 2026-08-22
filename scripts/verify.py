#!/usr/bin/env python3
"""Single verification entry point for the AutoCSF artifact.

    verify.py                 check that the committed results are complete and
                              internally consistent
    verify.py --fresh DIR     additionally compare a from-scratch run copied out
                              of a clean container against the committed results,
                              and write a reproduction receipt

Numeric comparison covers all three experiments. Two runs never agree bit for
bit: binary-fuse CSF construction is randomized, so measured sizes move by a few
thousandths of a bit per key between runs. TOLERANCE_BPK is what separates that
from a real regression.
"""

import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import struct

ROOT = Path(__file__).resolve().parents[1]
TOLERANCE_BPK = 0.02

DATASETS = {"ecoli_sakai", "srr10211353", "celegans", "rice"}
METHODS = {"hkp", "bcsf", "vlburr", "autocsf"}
THEORY_FILTERS = {"xor", "binary_fuse", "bloom_k1", "bloom_k2", "bloom_k3"}
THEORY_DISTRIBUTIONS = {"unique", "zipfian", "uniform_100"}
EPSILON_ALPHAS = ("0.7", "0.9")
PLOTS = {
    "alpha_sweep_uniform_100.png", "alpha_sweep_unique.png", "alpha_sweep_zipfian.png",
    "epsilon_sweep_uniform_100_alpha0.7.png", "epsilon_sweep_uniform_100_alpha0.9.png",
    "epsilon_sweep_unique_alpha0.7.png", "epsilon_sweep_unique_alpha0.9.png",
    "epsilon_sweep_zipfian_alpha0.7.png", "epsilon_sweep_zipfian_alpha0.9.png",
    "method-comparison.png",
}
RECEIPTS = {"synthetic": "synthetic-receipt.json", "small": "receipt.json"}


def theory_stems():
    alpha = {
        f"alpha_sweep_{f}_{d}"
        for f in THEORY_FILTERS for d in THEORY_DISTRIBUTIONS
    }
    epsilon = {
        f"epsilon_sweep_{f}_{d}_alpha{a}"
        for f in THEORY_FILTERS for d in THEORY_DISTRIBUTIONS for a in EPSILON_ALPHAS
    }
    return alpha, epsilon


# ---------------------------------------------------------------------------
# Regression checks
# ---------------------------------------------------------------------------

def frequency_counter_regression():
    """VL-BuRR's float32 frequency model lost unit increments above 2^24."""
    def float32(value):
        return struct.unpack("f", struct.pack("f", value))[0]

    errors = []
    limit = 1 << 24
    if float32(float32(limit) + 1.0) != limit:
        errors.append("float32 saturation assumption no longer holds")
    records, majority = 69_680_812, 57_445_000
    if not float32(limit / records) < 0.25:
        errors.append("frequency-counter regression: pre-fix probability changed")
    if not 0.82 < float32(majority / records) < 0.83:
        errors.append("frequency-counter regression: post-fix probability changed")
    return errors


# ---------------------------------------------------------------------------
# Completeness of the committed artifact
# ---------------------------------------------------------------------------

def check_bound_validation():
    """The per-run JSON is not tracked; baseline.csv is the committed record."""
    errors = []
    alpha_dir = ROOT / "results/bound_validation/figures/alpha_sweep"
    epsilon_dir = ROOT / "results/bound_validation/figures/epsilon_sweep"
    alpha_stems, epsilon_stems = theory_stems()
    for stem in sorted(alpha_stems):
        if not (alpha_dir / f"{stem}.png").exists():
            errors.append(f"missing alpha-sweep figure: {stem}")
    for stem in sorted(epsilon_stems):
        if not (epsilon_dir / f"{stem}.png").exists():
            errors.append(f"missing epsilon-sweep figure: {stem}")

    baseline_path = ROOT / "results/bound_validation/baseline.csv"
    if not baseline_path.exists():
        errors.append("missing results/bound_validation/baseline.csv")
        return errors
    covered = {key[0] for key in baseline_measurements(baseline_path)}
    for stem in sorted((alpha_stems | epsilon_stems) - covered):
        errors.append(f"baseline.csv has no measurements for {stem}")
    return errors


def check_synthetic():
    path = ROOT / "results/synthetic_comparison/data.csv"
    if not path.exists():
        return ["missing results/synthetic_comparison/data.csv"]
    rows = synthetic_rows(path)
    present = {method for _, method, _ in rows}
    if present != METHODS:
        return [f"synthetic comparison is missing methods: {sorted(METHODS - present)}"]
    return []


def check_genomics():
    errors = []
    genomics_dir = ROOT / "results/genomics_comparison"
    table = genomics_dir / "genomics.csv"
    if not table.exists():
        return ["missing results/genomics_comparison/genomics.csv"]
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
    payload_path = genomics_dir / "genomics.json"
    if payload_path.exists():
        payload = json.loads(payload_path.read_text())
        for row in payload.get("records", []):
            if row.get("method") != "vlburr":
                continue
            baseline = row.get("baseline_bits_per_key")
            if baseline is None or abs((baseline - row["bits_per_key"]) - row["bits_saved_vs_plain"]) > 1e-9:
                errors.append(f"VL-BuRR row uses an invalid baseline: {row.get('dataset')}")
    for name in ("genomics.json", "genomics.md", "genomics.tex",
                 "genomics-paper.md", "genomics-paper.tex",
                 "genomics-audit.md", "genomics-audit.tex",
                 "genomics-total.md", "genomics-total.tex"):
        if not (genomics_dir / name).exists():
            errors.append(f"missing results/genomics_comparison/{name}")
    return errors


def check_plots():
    errors = []
    actual = {p.name for p in (ROOT / "results/bound_validation/figures/paper").glob("*.png")}
    if (ROOT / "results/synthetic_comparison/method-comparison.png").exists():
        actual.add("method-comparison.png")
    for name in sorted(PLOTS - actual):
        errors.append(f"missing plot: {name}")
    return errors


def check_reference():
    """The accepted-paper figures and table, kept for side-by-side comparison.

    Presence only. These are not checksummed, so nothing here detects a
    reference that has been edited or replaced.
    """
    reference = ROOT / "reference"
    errors = []
    for name in sorted(PLOTS | {"genomics-table.md"}):
        if not (reference / name).exists():
            errors.append(f"missing accepted-paper reference: {name}")
    return errors


def check_datasets():
    errors = []
    for dataset in sorted(DATASETS):
        manifest = ROOT / "data/manifests" / f"{dataset}_k15.json"
        processed = ROOT / "data/processed" / f"{dataset}_k15.tsv.zst"
        if not manifest.exists() or not processed.exists():
            errors.append(f"missing dataset or manifest: {dataset}")
    return errors


def check_receipts():
    errors = []
    for scope, name in sorted(RECEIPTS.items()):
        path = ROOT / "results/reproduction" / name
        if not path.exists():
            errors.append(f"missing {scope} reproduction receipt ({name})")
            continue
        try:
            receipt = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError):
            errors.append(f"{name} is invalid JSON")
            continue
        if receipt.get("status") != "passed":
            errors.append(f"{name} did not pass")
        comparison = receipt.get("comparison")
        if comparison is not None and comparison.get("status") != "passed":
            errors.append(f"{name} records an unvalidated comparison")
    return errors


def check_artifact(scope_datasets=True):
    errors = frequency_counter_regression()
    errors += check_bound_validation()
    errors += check_synthetic()
    errors += check_genomics()
    errors += check_plots()
    errors += check_reference()
    errors += check_receipts()
    if scope_datasets:
        errors += check_datasets()
    if (ROOT / "results/KNOWN_ISSUES.md").exists():
        errors.append("unresolved benchmark correctness issue; see results/KNOWN_ISSUES.md")
    return errors


# ---------------------------------------------------------------------------
# Numeric comparison of a fresh run against the committed results
# ---------------------------------------------------------------------------

def synthetic_rows(path):
    with path.open() as handle:
        rows = list(csv.DictReader(handle))
    return {
        (row["distribution"], row["method"], float(row["alpha"])): float(row["bits_saved"])
        for row in rows
    }


def genomics_rows(path):
    with path.open() as handle:
        rows = list(csv.DictReader(handle))
    return {
        (row["dataset"], row["method"]): (
            float(row["bits_per_key"]), float(row["bits_saved_vs_plain"])
        )
        for row in rows
    }


def baseline_measurements(path):
    """Read the committed bound-validation baseline into {key: value}."""
    with path.open() as handle:
        return {
            (row["sweep"], int(row["index"]), row["metric"], row["param"]): float(row["value"])
            for row in csv.DictReader(handle)
        }


def fresh_measurements(data_dir):
    """Read a fresh run's JSON into the same shape the baseline uses."""
    out = {}
    for path in sorted(data_dir.glob("*.json")):
        data = json.loads(path.read_text())
        for index, record in enumerate(data.get("results", [data])):
            for metric in ("baseline_bpk", "theory_guided_bpk_saved",
                           "best_empirical_bpk_saved"):
                if metric in record:
                    out[(path.stem, index, metric, "")] = float(record[metric])
            for entry in record.get("empirical_per_param", []):
                param = next(v for k, v in entry.items() if k != "bpk_saved")
                out[(path.stem, index, "bpk_saved", str(param))] = float(entry["bpk_saved"])
    return out


def compare_bound_validation(fresh_dir, errors):
    fresh_data = fresh_dir / "bound_validation/data"
    baseline_path = ROOT / "results/bound_validation/baseline.csv"
    if not fresh_data.is_dir():
        return 0.0, 0
    if not baseline_path.exists():
        errors.append("missing results/bound_validation/baseline.csv")
        return 0.0, 0

    accepted = baseline_measurements(baseline_path)
    fresh = fresh_measurements(fresh_data)
    missing = sorted(accepted.keys() - fresh.keys())
    extra = sorted(fresh.keys() - accepted.keys())
    if missing:
        errors.append(f"bound validation is missing accepted measurements: {missing[:3]}")
    if extra:
        errors.append(f"bound validation produced unexpected measurements: {extra[:3]}")

    max_delta = 0.0
    for key in accepted.keys() & fresh.keys():
        delta = abs(accepted[key] - fresh[key])
        max_delta = max(max_delta, delta)
        if delta > TOLERANCE_BPK:
            sweep, index, metric, param = key
            label = f"{sweep} result {index} {metric}" + (f" param={param}" if param else "")
            errors.append(f"{label} delta={delta:.6g}")
    return max_delta, len(list(fresh_data.glob("*.json")))


def compare_table(accepted_path, fresh_path, reader, label, errors, allow_subset=False):
    """Compare one CSV of measurements, returning (max delta, fresh row count).

    allow_subset is for bounded runs that deliberately cover only part of the
    matrix: rows the fresh run did not produce are not an error, but a row it
    produced that the accepted results do not contain always is.
    """
    if not fresh_path.exists():
        return 0.0, 0
    accepted = reader(accepted_path) if accepted_path.exists() else {}
    fresh = reader(fresh_path)
    missing = sorted(accepted.keys() - fresh.keys())
    extra = sorted(fresh.keys() - accepted.keys())
    if extra:
        errors.append(f"{label} produced rows absent from the accepted results: {extra[:4]}")
    if missing and not allow_subset:
        errors.append(f"{label} is missing accepted rows: {missing[:4]}")
    max_delta = 0.0
    for key in accepted.keys() & fresh.keys():
        old, new = accepted[key], fresh[key]
        pairs = zip(old, new) if isinstance(old, tuple) else [(old, new)]
        for old_value, new_value in pairs:
            delta = abs(old_value - new_value)
            max_delta = max(max_delta, delta)
            if delta > TOLERANCE_BPK:
                errors.append(f"{label} {key} delta={delta:.6g}")
    return max_delta, len(fresh)


def compare(fresh_dir, scope, image_id):
    errors = []
    theory_delta, data_files = compare_bound_validation(fresh_dir, errors)
    method_delta, method_rows = compare_table(
        ROOT / "results/synthetic_comparison/data.csv",
        fresh_dir / "synthetic_comparison/data.csv",
        synthetic_rows, "synthetic comparison", errors,
    )
    genomics_delta, genomics_records = compare_table(
        ROOT / "results/genomics_comparison/genomics.csv",
        fresh_dir / "genomics_comparison/genomics.csv",
        genomics_rows, "genomics comparison", errors,
        allow_subset=scope == "small",
    )

    plots = list((fresh_dir / "bound_validation/figures/paper").glob("*.png"))
    plots += list((fresh_dir / "synthetic_comparison").glob("method-comparison.png"))

    def stamp(name):
        path = fresh_dir / name
        return path.read_text().strip() if path.exists() else "unknown"

    return {
        "status": "passed" if not errors else "failed",
        "kind": f"clean-{scope}-reproduction",
        "checked_at": datetime.now(timezone.utc).isoformat(),
        "container_image_id": image_id,
        "started_at": stamp("started_at.txt"),
        "completed_at": stamp("completed_at.txt"),
        "fresh_artifacts": {
            "data_files": data_files,
            "method_rows": method_rows,
            "genomics_records": genomics_records,
            "plots": len(plots),
        },
        "comparison": {
            "status": "passed" if not errors else "failed",
            "bits_per_key_absolute_tolerance": TOLERANCE_BPK,
            "max_bound_validation_delta": theory_delta,
            "max_synthetic_comparison_delta": method_delta,
            "max_genomics_comparison_delta": genomics_delta,
            "errors": errors[:100],
        },
    }


# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fresh", type=Path, help="directory copied out of a clean container")
    parser.add_argument("--scope", choices=sorted(RECEIPTS), default="synthetic")
    parser.add_argument("--image-id", default="unknown")
    parser.add_argument(
        "--receipt", type=Path,
        help="write the receipt here instead of results/reproduction/ (for testing)",
    )
    parser.add_argument(
        "--skip-datasets", action="store_true",
        help="do not require the processed genomics tables to be present locally",
    )
    args = parser.parse_args()

    if args.fresh:
        receipt = compare(args.fresh, args.scope, args.image_id)
        output = args.receipt or ROOT / "results/reproduction" / RECEIPTS[args.scope]
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(receipt, indent=2) + "\n")
        print(json.dumps(receipt, indent=2))
        if receipt["status"] != "passed":
            raise SystemExit(1)
        return

    errors = check_artifact(scope_datasets=not args.skip_datasets)
    if errors:
        print("Verification failed:")
        for error in errors:
            print(f"- {error}")
        raise SystemExit(1)
    print("All AutoCSF reproduction artifacts are complete and consistent.")


if __name__ == "__main__":
    main()
