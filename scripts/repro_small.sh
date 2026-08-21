#!/usr/bin/env bash
# Bounded clean-room check: E. coli on four methods plus all accepted plots.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

repetitions="${REPETITIONS:-1}"
if [[ ! -x .venv/bin/python ]]; then
  ./setup.sh
fi
.venv/bin/python datasets/manage.py fetch --dataset ecoli_sakai
.venv/bin/python datasets/manage.py generate --dataset ecoli_sakai
.venv/bin/python datasets/manage.py profile --dataset ecoli_sakai
.venv/bin/python datasets/manage.py validate --dataset ecoli_sakai
.venv/bin/python baselines/run_genomics.py --dataset ecoli_sakai \
  --repetitions "$repetitions" --output results/reproduction/ecoli
.venv/bin/python scripts/reproduce_plots.py

.venv/bin/python - <<'PY'
import csv
from pathlib import Path
root = Path.cwd()
rows = list(csv.DictReader((root / "results/reproduction/ecoli/genomics.csv").open()))
assert {row["method"] for row in rows} == {"hkp", "bcsf", "vlburr", "autocsf"}
assert len(list((root / "results/figures").glob("*.png"))) == 10
print("Bounded reproduction passed: E. coli 4 methods and 10 plots")
PY
