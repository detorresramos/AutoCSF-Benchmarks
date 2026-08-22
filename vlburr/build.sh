#!/usr/bin/env bash
# Build from pristine upstream sources; never trust or modify a dirty submodule checkout.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
commit="$(cat "$ROOT/vlburr/UPSTREAM_COMMIT")"
source_dir="$ROOT/data/cache/lsf-source-$commit"
if [[ ! -d "$source_dir/.git" ]]; then
  mkdir -p "$ROOT/data/cache"
  git clone --recursive https://github.com/gvinciguerra/LearnedStaticFunction.git "$source_dir"
  git -C "$source_dir" checkout "$commit"
  git -C "$source_dir" submodule update --init --recursive
  patch -d "$source_dir" -p0 < "$ROOT/vlburr/patches/nofilter.patch"
fi
if [[ ! -f "$source_dir/.autocsf-integer-frequency-counts" ]]; then
  patch -d "$source_dir" -p0 < "$ROOT/vlburr/patches/integer-frequency-counts.patch"
  touch "$source_dir/.autocsf-integer-frequency-counts"
fi
if [[ ! -f "$source_dir/.autocsf-release-filter-input" ]]; then
  patch -d "$source_dir" -p0 < "$ROOT/vlburr/patches/release-filter-input.patch"
  touch "$source_dir/.autocsf-release-filter-input"
fi
if [[ "$(uname -s)" == Darwin && ! -f "$source_dir/.autocsf-macos-patched" ]]; then
  patch -d "$source_dir" -p0 < "$ROOT/vlburr/patches/macos-cstdint.patch"
  touch "$source_dir/.autocsf-macos-patched"
fi
compiler_flags="${CMAKE_CXX_FLAGS:-}"
if [[ "$(uname -s)" == Darwin && -z "$compiler_flags" ]]; then
  compiler_flags="-I/opt/homebrew/include"
fi
cmake -S "$source_dir" -B "$source_dir/build" -DCMAKE_BUILD_TYPE=Release \
  -DTFLITE_ENABLE_XNNPACK=OFF -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
  "-DCMAKE_CXX_FLAGS=$compiler_flags"
cmake --build "$source_dir/build" --target ribbon_learned_bench -j "${JOBS:-4}"
mkdir -p "$ROOT/data/cache/bin"
cp "$source_dir/build/ribbon_learned_bench" "$ROOT/data/cache/bin/ribbon_learned_bench"
find "$source_dir/build" -maxdepth 1 -type f \
  \( -name 'libRibbonSorter.so*' -o -name 'libRibbonSorter.dylib' \) \
  -exec cp {} "$ROOT/data/cache/bin/" \;
echo "Built $ROOT/data/cache/bin/ribbon_learned_bench from $commit"
