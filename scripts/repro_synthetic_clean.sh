#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
IMAGE=${AUTOCSF_SYNTHETIC_IMAGE:-autocsf-repro:synthetic-clean}
RUN_DIR=${AUTOCSF_SYNTHETIC_OUTPUT:-$(mktemp -d "${TMPDIR:-/tmp}/autocsf-synthetic.XXXXXX")}

mkdir -p "$RUN_DIR"
docker build --platform linux/amd64 -t "$IMAGE" "$ROOT"
docker run --rm --platform linux/amd64 -v "$RUN_DIR:/out" "$IMAGE" -lc '
set -euo pipefail
rm -rf theory_validation/figures method_comparison/data.csv results/figures data/cache/method-comparison
mkdir -p theory_validation/figures/data
date -u +%FT%TZ > /out/started_at.txt
.venv/bin/python theory_validation/run_experiments.py 2>&1 | tee /out/theory.log
.venv/bin/python method_comparison/run_experiments.py 2>&1 | tee /out/method.log
.venv/bin/python scripts/reproduce_plots.py 2>&1 | tee /out/plots.log
cp -a theory_validation/figures /out/theory_figures
cp method_comparison/data.csv /out/method_comparison.csv
cp -a results/figures /out/results_figures
date -u +%FT%TZ > /out/completed_at.txt
'

IMAGE_ID=$(docker image inspect "$IMAGE" --format '{{.Id}}')
"$ROOT/.venv/bin/python" "$ROOT/scripts/compare_synthetic_reproduction.py" \
  "$RUN_DIR" --image-id "$IMAGE_ID"

rm -rf "$ROOT/theory_validation/figures"
cp -a "$RUN_DIR/theory_figures" "$ROOT/theory_validation/figures"
cp "$RUN_DIR/method_comparison.csv" "$ROOT/method_comparison/data.csv"
rm -rf "$ROOT/results/figures"
cp -a "$RUN_DIR/results_figures" "$ROOT/results/figures"

echo "Clean synthetic reproduction passed. Raw run: $RUN_DIR"
