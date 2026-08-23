#!/usr/bin/env python3
"""Check what is in results/ against the accepted numbers in baselines/.

Verifies that every experiment produced a complete set of outputs, compares each
measurement against its accepted value, and writes results/reproduction/receipt.json.

Two runs never agree bit for bit: binary-fuse CSF construction is randomized, so
measured sizes move by a few thousandths of a bit per key between runs.
TOLERANCE_BPK is what separates that from a real regression.

Run scripts/export_baseline.py to move results/ into baselines/ after a run that
legitimately changes the accepted numbers.
"""

import argparse
import csv
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import platform
import struct
import subprocess

ROOT = Path(__file__).resolve().parents[1]
BASELINES = ROOT / "baselines"
RESULTS = ROOT / "results"
RECEIPT = RESULTS / "reproduction/receipt.json"
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


def environment():
    def output(command):
        try:
            return subprocess.check_output(command, text=True, stderr=subprocess.DEVNULL).strip()
        except (OSError, subprocess.CalledProcessError):
            return "unavailable"
    return {
        "platform": platform.platform(),
        "python": platform.python_version(),
        "git_commit": output(["git", "-C", str(ROOT), "rev-parse", "HEAD"]),
        "carameldb_commit": output(["git", "-C", str(ROOT / "deps/CaramelDB"), "rev-parse", "HEAD"]),
        "docker": os.environ.get("AUTOCSF_DOCKER_DESCRIPTION", "not recorded"),
    }


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
    """The per-run JSON is not tracked; baselines/bound_validation.csv is."""
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

    baseline_path = BASELINES / "bound_validation.csv"
    if not baseline_path.exists():
        errors.append("missing baselines/bound_validation.csv")
        return errors
    covered = {key[0] for key in baseline_measurements(baseline_path)}
    for stem in sorted((alpha_stems | epsilon_stems) - covered):
        errors.append(f"baselines/bound_validation.csv has no measurements for {stem}")
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


def check_artifact(scope_datasets=True):
    errors = frequency_counter_regression()
    errors += check_bound_validation()
    errors += check_synthetic()
    errors += check_genomics()
    errors += check_plots()
    errors += check_reference()
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


def compare_bound_validation(errors):
    fresh_data = RESULTS / "bound_validation/data"
    baseline_path = BASELINES / "bound_validation.csv"
    if not fresh_data.is_dir():
        return 0.0, 0
    if not baseline_path.exists():
        errors.append("missing baselines/bound_validation.csv")
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


def compare_table(accepted_path, fresh_path, reader, label, errors):
    """Compare one CSV of measurements, returning (max delta, fresh row count)."""
    if not fresh_path.exists():
        return 0.0, 0
    accepted = reader(accepted_path) if accepted_path.exists() else {}
    fresh = reader(fresh_path)
    missing = sorted(accepted.keys() - fresh.keys())
    extra = sorted(fresh.keys() - accepted.keys())
    if extra:
        errors.append(f"{label} produced rows absent from the accepted results: {extra[:4]}")
    if missing:
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


def compare(errors):
    theory_delta, data_files = compare_bound_validation(errors)
    method_delta, method_rows = compare_table(
        BASELINES / "synthetic_comparison.csv",
        RESULTS / "synthetic_comparison/data.csv",
        synthetic_rows, "synthetic comparison", errors,
    )
    genomics_delta, genomics_records = compare_table(
        BASELINES / "genomics_comparison.csv",
        RESULTS / "genomics_comparison/genomics.csv",
        genomics_rows, "genomics comparison", errors,
    )
    return {
        "bits_per_key_absolute_tolerance": TOLERANCE_BPK,
        "bound_validation": {"files": data_files, "max_delta": theory_delta},
        "synthetic_comparison": {"rows": method_rows, "max_delta": method_delta},
        "genomics_comparison": {"rows": genomics_records, "max_delta": genomics_delta},
    }


# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--skip-datasets", action="store_true",
        help="do not require the processed genomics tables to be present locally",
    )
    args = parser.parse_args()

    errors = check_artifact(scope_datasets=not args.skip_datasets)
    comparison = compare(errors)
    receipt = {
        "status": "passed" if not errors else "failed",
        "checked_at": datetime.now(timezone.utc).isoformat(),
        "environment": environment(),
        "comparison": comparison,
        "errors": errors[:100],
    }
    RECEIPT.parent.mkdir(parents=True, exist_ok=True)
    RECEIPT.write_text(json.dumps(receipt, indent=2) + "\n")

    if errors:
        print("Verification failed:")
        for error in errors:
            print(f"- {error}")
        raise SystemExit(1)
    print(f"All results match baselines/ within {TOLERANCE_BPK} bits/key. "
          f"Receipt: {RECEIPT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
