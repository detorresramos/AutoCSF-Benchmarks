#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
IMAGE=${AUTOCSF_SYNTHETIC_IMAGE:-autocsf-repro:synthetic-clean}
RUN_DIR=${AUTOCSF_SYNTHETIC_OUTPUT:-$(mktemp -d "${TMPDIR:-/tmp}/autocsf-synthetic.XXXXXX")}

mkdir -p "$RUN_DIR"
docker build --platform linux/amd64 -t "$IMAGE" "$ROOT"
docker run --rm --platform linux/amd64 -v "$RUN_DIR:/out" "$IMAGE" -lc '
set -euo pipefail
rm -rf results/bound_validation results/synthetic_comparison data/cache/method-comparison
mkdir -p results/bound_validation/data results/synthetic_comparison
date -u +%FT%TZ > /out/started_at.txt
.venv/bin/python experiments/bound_validation/run.py 2>&1 | tee /out/theory.log
.venv/bin/python experiments/synthetic_comparison/run.py 2>&1 | tee /out/method.log
.venv/bin/python scripts/reproduce_plots.py 2>&1 | tee /out/plots.log
cp -a results/bound_validation /out/bound_validation
cp -a results/synthetic_comparison /out/synthetic_comparison
date -u +%FT%TZ > /out/completed_at.txt
'

IMAGE_ID=$(docker image inspect "$IMAGE" --format '{{.Id}}')
"$ROOT/.venv/bin/python" "$ROOT/scripts/compare_synthetic_reproduction.py" \
  "$RUN_DIR" --image-id "$IMAGE_ID"

mkdir -p "$ROOT/results/bound_validation" "$ROOT/results/synthetic_comparison"
cp -a "$RUN_DIR/bound_validation/." "$ROOT/results/bound_validation/"
cp -a "$RUN_DIR/synthetic_comparison/." "$ROOT/results/synthetic_comparison/"

echo "Clean synthetic reproduction passed. Raw run: $RUN_DIR"
