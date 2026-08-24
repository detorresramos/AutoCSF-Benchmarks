#!/usr/bin/env python3
"""Promote the numbers in results/ to the accepted baselines in baselines/.

Bound validation is condensed to its bits-per-key numbers; the other two
experiments are copied as they stand.
"""

import csv
import json
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
BASELINES = ROOT / "baselines"
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


def export_bound_validation():
    paths = sorted((RESULTS / "bound_validation/data").glob("*.json"))
    if not paths:
        raise SystemExit("no bound-validation JSON; run the experiment first")
    rows = [row for path in paths for row in rows_for(path)]
    out = BASELINES / "bound_validation.csv"
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    return out, f"{len(rows)} rows from {len(paths)} sweeps"


def export_copy(source, name):
    if not source.exists():
        raise SystemExit(f"missing {source.relative_to(ROOT)}; run the experiment first")
    out = BASELINES / name
    shutil.copyfile(source, out)
    with out.open() as handle:
        count = sum(1 for _ in handle) - 1
    return out, f"{count} rows"


def main():
    BASELINES.mkdir(parents=True, exist_ok=True)
    exports = [
        export_bound_validation(),
        export_copy(RESULTS / "synthetic_comparison/data.csv", "synthetic_comparison.csv"),
        export_copy(RESULTS / "genomics_comparison/genomics.csv", "genomics_comparison.csv"),
    ]
    for path, detail in exports:
        print(f"{path.relative_to(ROOT)}: {detail}")


if __name__ == "__main__":
    main()
