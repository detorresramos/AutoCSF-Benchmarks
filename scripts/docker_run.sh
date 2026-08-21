#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
image="${AUTOCSF_IMAGE:-autocsf-benchmarks}"
docker build --platform linux/amd64 -t "$image" "$ROOT"
docker run --rm --platform linux/amd64 \
  -v "$ROOT/data:/workspace/data" \
  -v "$ROOT/results:/workspace/results" \
  "$image" -lc "${*:-make all}"
