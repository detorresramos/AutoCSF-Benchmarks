#!/usr/bin/env python3
"""Condense the bound-validation JSON into the committed comparison baseline.

The per-run JSON under results/bound_validation/data/ is 39k lines of
pretty-printed detail and is not tracked. This writes the subset that
verify.py actually compares -- every bits-per-key number, nothing else -- as one
long-format CSV that is small enough to live in git and to read in a diff.

Regenerate it after any run that legitimately changes the accepted numbers.
"""

import csv
import glob
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "results/bound_validation/data"
BASELINE = ROOT / "results/bound_validation/baseline.csv"
FIELDS = ["sweep", "index", "alpha", "metric", "param", "value"]


def rows_for(path):
    data = json.loads(path.read_text())
    sweep = path.stem
    # Alpha sweeps hold a list of per-alpha records; epsilon sweeps are a single
    # record at the top level. Normalise both to an indexed list.
    records = data.get("results", [data])
    for index, record in enumerate(records):
        alpha = record.get("alpha", data.get("alpha"))
        for metric in ("baseline_bpk", "theory_guided_bpk_saved",
                       "best_empirical_bpk_saved"):
            if metric in record:
                yield {"sweep": sweep, "index": index, "alpha": alpha,
                       "metric": metric, "param": "", "value": repr(record[metric])}
        for entry in record.get("empirical_per_param", []):
            param = next(v for k, v in entry.items() if k != "bpk_saved")
            yield {"sweep": sweep, "index": index, "alpha": alpha,
                   "metric": "bpk_saved", "param": param,
                   "value": repr(entry["bpk_saved"])}


def main():
    paths = sorted(DATA_DIR.glob("*.json"))
    if not paths:
        raise SystemExit(f"no JSON under {DATA_DIR}; run the experiment first")
    rows = [row for path in paths for row in rows_for(path)]
    BASELINE.parent.mkdir(parents=True, exist_ok=True)
    with BASELINE.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"{BASELINE.relative_to(ROOT)}: {len(rows)} rows from {len(paths)} sweeps")


if __name__ == "__main__":
    main()
