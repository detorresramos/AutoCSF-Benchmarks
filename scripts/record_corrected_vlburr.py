#!/usr/bin/env python3
"""Update source result records from corrected VL-BuRR runs."""

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from vlburr.parse_genomics import parse


def update(path, dataset, measured):
    payload = json.loads(path.read_text())
    row = next(row for row in payload["records"] if row["dataset"] == dataset and row["method"] == "vlburr")
    row.update(measured)
    path.write_text(json.dumps(payload, indent=2) + "\n")


with (ROOT / "data/cache/vlburr/celegans/benchmark.log").open() as handle:
    corrected = parse(handle, 69_680_812)
update(ROOT / "results/full-vlburr/genomics.json", "celegans", corrected)

oom = {
    "bits_per_key": None,
    "bits_saved_vs_plain": None,
    "build_seconds": None,
    "query_ns": None,
    "serialized_bytes": None,
    "baseline_bits_per_key": 2.266399,
    "status": "OOM (8.2 GB limit)",
    "parameters": {
        "implementation": "Filtered-HuffmanCSF_Opt",
        "correctness_patch": "integer-frequency-counts.patch",
        "memory_patch": "release-filter-input.patch",
    },
    "repetitions": 1,
}
update(ROOT / "results/rice-vlburr/genomics.json", "rice", oom)
