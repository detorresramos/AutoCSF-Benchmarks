#!/usr/bin/env bash
# Run the Opt-vs-NoFilter head-to-head on every acsf_*_y.lrbin dataset in
# LRDATA_DIR, inside the autocsf-bench Docker image. Results are streamed to
# RESULTS_DIR/sweep.log (raw RESULT lines included).
#
# Usage:
#   bash run_sweep.sh [LRDATA_DIR] [RESULTS_DIR]
#
# Defaults: ./lrdata, ./results. Overridable via positional args or env vars.
set -euo pipefail

LRDATA_DIR="${1:-${LRDATA_DIR:-$PWD/lrdata}}"
RESULTS_DIR="${2:-${RESULTS_DIR:-$PWD/results}}"
IMAGE="${IMAGE:-autocsf-bench}"
DOCKER="${DOCKER:-sudo docker}"

if [[ ! -d "$LRDATA_DIR" ]]; then
  echo "error: LRDATA_DIR does not exist: $LRDATA_DIR" >&2
  exit 1
fi
mkdir -p "$RESULTS_DIR"
out="$RESULTS_DIR/sweep.log"

echo "image:   $IMAGE"
echo "lrdata:  $LRDATA_DIR"
echo "output:  $out"

# One docker invocation that loops over every dataset to amortize startup cost.
# taskset -c 0 pins the bench to a single core for stable timings.
$DOCKER run --rm \
  -v "$LRDATA_DIR:/lrdata" \
  --entrypoint bash \
  "$IMAGE" -c '
    set -e
    shopt -s nullglob
    files=(/lrdata/acsf_*_y.lrbin)
    if (( ${#files[@]} == 0 )); then
      echo "no acsf_*_y.lrbin files in /lrdata" >&2
      exit 1
    fi
    for f in "${files[@]}"; do
      ds=$(basename "$f" _y.lrbin)
      taskset -c 0 /lsf/build/ribbon_learned_bench -r /lrdata/ -c CSF -d "$ds"
    done
  ' | tee "$out"

n=$(grep -c '^RESULT' "$out" || true)
echo
echo "Done. $n RESULT lines in $out"
