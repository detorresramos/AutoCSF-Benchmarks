#!/usr/bin/env bash
# The only container entry point. Builds the linux/amd64 image (required for
# VL-BuRR) and runs a command in it.
#
#   scripts/docker_run.sh                 # make reproduce
#   scripts/docker_run.sh make bound-validation   # any other target
#
# By default results/ is bind-mounted so the container writes straight into the
# repository. Set AUTOCSF_OUT_MOUNT to a directory instead and results/ stays
# untouched: the container gets a writable /out and computes everything fresh.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
image="${AUTOCSF_IMAGE:-autocsf-benchmarks}"

mkdir -p "$ROOT/data/raw" "$ROOT/data/processed" "$ROOT/data/manifests" "$ROOT/results"

mounts=(
  --mount "type=bind,source=$ROOT/data/raw,target=/workspace/data/raw"
  --mount "type=bind,source=$ROOT/data/processed,target=/workspace/data/processed"
  --mount "type=bind,source=$ROOT/data/manifests,target=/workspace/data/manifests"
)
if [[ -n "${AUTOCSF_OUT_MOUNT:-}" ]]; then
  mkdir -p "$AUTOCSF_OUT_MOUNT"
  mounts+=(--mount "type=bind,source=$AUTOCSF_OUT_MOUNT,target=/out")
else
  mounts+=(--mount "type=bind,source=$ROOT/results,target=/workspace/results")
fi

docker build --platform linux/amd64 -t "$image" "$ROOT"
docker run --rm --platform linux/amd64 "${mounts[@]}" "$image" -lc "${*:-make reproduce}"
