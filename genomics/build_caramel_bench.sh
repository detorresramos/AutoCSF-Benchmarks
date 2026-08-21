#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
build="$ROOT/deps/CaramelDB/pybind/build"
library="$build/libcaramel_lib.a"
[[ -f "$library" ]] || { echo "CaramelDB native library is missing; run ./setup.sh" >&2; exit 1; }
mkdir -p "$ROOT/data/cache/bin"
common=(-O3 -std=c++17 -I"$ROOT/deps/CaramelDB" -I"$ROOT/deps/CaramelDB/deps"
  -I"$ROOT/deps/CaramelDB/deps/cereal/include")
if [[ "$(uname -s)" == Darwin ]]; then
  c++ "${common[@]}" -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include \
    "$ROOT/genomics/caramel_bench.cpp" "$library" \
    -L/opt/homebrew/opt/libomp/lib -lomp -o "$ROOT/data/cache/bin/caramel_bench"
else
  # GNU mode provides libstdc++'s hash specialization for __uint128_t, which
  # CaramelDB's bucket implementation uses internally.
  common[1]=-std=gnu++17
  c++ "${common[@]}" -fopenmp "$ROOT/genomics/caramel_bench.cpp" "$library" \
    -o "$ROOT/data/cache/bin/caramel_bench"
fi
