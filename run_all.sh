#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
make experiments
make plots
make genomics
make stage
make repro-small
make verify
