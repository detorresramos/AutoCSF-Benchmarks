#!/usr/bin/env bash
# From-scratch reproduction in a clean container, compared against the results
# committed in this repository.
#
#   SCOPE=synthetic  (default)  all 45 bound-validation datasets, all 88
#                               synthetic-comparison rows, all ten paper plots
#   SCOPE=small                 the above plus E. coli genomics on four methods
#
# The container never sees results/; it computes everything itself and copies
# the output to /out. Nothing is promoted into results/ until the comparison
# passes.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SCOPE="${SCOPE:-synthetic}"
RUN_DIR="${AUTOCSF_REPRO_OUTPUT:-$(mktemp -d "${TMPDIR:-/tmp}/autocsf-repro.XXXXXX")}"
IMAGE="${AUTOCSF_IMAGE:-autocsf-benchmarks}"

case "$SCOPE" in
  synthetic) EXTRA="" ;;
  small)
    EXTRA='
.venv/bin/python datasets/manage.py fetch --dataset ecoli_sakai
.venv/bin/python datasets/manage.py generate --dataset ecoli_sakai
.venv/bin/python datasets/manage.py profile --dataset ecoli_sakai
.venv/bin/python datasets/manage.py validate --dataset ecoli_sakai
.venv/bin/python experiments/genomics_comparison/run.py --dataset ecoli_sakai \
  --repetitions "${REPETITIONS:-1}" --output results/genomics_comparison
cp -a results/genomics_comparison /out/genomics_comparison
'
    ;;
  *) echo "unknown SCOPE: $SCOPE (expected synthetic or small)" >&2; exit 1 ;;
esac

AUTOCSF_OUT_MOUNT="$RUN_DIR" "$ROOT/scripts/docker_run.sh" "
set -euo pipefail
rm -rf results data/cache/method-comparison
mkdir -p results/bound_validation/data results/synthetic_comparison
date -u +%FT%TZ > /out/started_at.txt
.venv/bin/python experiments/bound_validation/run.py 2>&1 | tee /out/bound_validation.log
.venv/bin/python experiments/synthetic_comparison/run.py 2>&1 | tee /out/synthetic.log
.venv/bin/python scripts/reproduce_plots.py 2>&1 | tee /out/plots.log
${EXTRA}
cp -a results/bound_validation /out/bound_validation
cp -a results/synthetic_comparison /out/synthetic_comparison
date -u +%FT%TZ > /out/completed_at.txt
"

IMAGE_ID="$(docker image inspect "$IMAGE" --format '{{.Id}}')"
"$ROOT/.venv/bin/python" "$ROOT/scripts/verify.py" \
  --fresh "$RUN_DIR" --scope "$SCOPE" --image-id "$IMAGE_ID"

for name in bound_validation synthetic_comparison genomics_comparison; do
  [[ -d "$RUN_DIR/$name" ]] || continue
  mkdir -p "$ROOT/results/$name"
  cp -a "$RUN_DIR/$name/." "$ROOT/results/$name/"
done

echo "Clean $SCOPE reproduction passed. Raw run: $RUN_DIR"
