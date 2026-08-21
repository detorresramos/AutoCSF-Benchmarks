#!/usr/bin/env bash
# Camera-ready real-genomics entry point.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

if [[ ! -x .venv/bin/python ]]; then
  echo "error: run ./setup.sh first" >&2
  exit 1
fi

# Benchmark startup checks manifests and file presence. Full byte-for-byte
# validation remains available as `make validate` and is intentionally not
# repeated before every timed subset run.
.venv/bin/python datasets/manage.py validate --manifest-only
.venv/bin/python baselines/run_genomics.py "$@"
