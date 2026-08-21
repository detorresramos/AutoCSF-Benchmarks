#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
dataset="${1:?dataset name required}"
repeats="${2:-1}"
table="$ROOT/data/processed/${dataset}_k15.tsv.zst"
cache="$ROOT/data/cache/vlburr/$dataset"
binary="${VLBURR_BINARY:-$ROOT/data/cache/bin/ribbon_learned_bench}"
[[ -f "$table" ]] || { echo "missing $table" >&2; exit 1; }
[[ -x "$binary" ]] || { echo "missing VL-BuRR binary: $binary" >&2; exit 1; }
mkdir -p "$cache"
prefix="$cache/$dataset"
if [[ ! -f "${prefix}_y.lrbin" ]]; then
  python3 "$ROOT/vlburr/table_to_lrbin.py" "$table" "$prefix" \
    --manifest "$ROOT/data/manifests/${dataset}_k15.json"
fi
n="$(python3 - "$ROOT/data/manifests/${dataset}_k15.json" <<'PY'
import json,sys
print(json.load(open(sys.argv[1]))["records"])
PY
)"
log="$cache/benchmark.log"
: > "$log"
for ((i=0; i<repeats; i++)); do
  binary_dir="$(cd "$(dirname "$binary")" && pwd)"
  LD_LIBRARY_PATH="$binary_dir${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" \
    DYLD_LIBRARY_PATH="$binary_dir${DYLD_LIBRARY_PATH:+:$DYLD_LIBRARY_PATH}" \
    "$binary" -r "$cache/" -c CSF -d "$dataset" >> "$log"
done
python3 "$ROOT/vlburr/parse_genomics.py" "$n" < "$log"
