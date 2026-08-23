#!/usr/bin/env bash
# Camera-ready real-genomics entry point.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"

if [[ ! -x .venv/bin/python ]]; then
  echo "error: run ./setup.sh first" >&2
  exit 1
fi

# Manifests and file presence only; `make verify` does the full byte-for-byte
# validation.
.venv/bin/python datasets/manage.py validate --manifest-only
.venv/bin/python experiments/genomics_comparison/run.py "$@"
