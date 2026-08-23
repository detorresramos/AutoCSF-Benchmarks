#!/usr/bin/env bash
# Install all dependencies for AutoCSF-Benchmarks on Ubuntu 24.04 x86_64.
#
# Usage:
#   sudo ./setup.sh           # install system packages + build everything
#   ./setup.sh --no-system    # skip apt-get, only build deps (already installed)
#   ./setup.sh --skip-vlburr  # skip VL-BuRR (it requires x86-64)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# ---------------------------------------------------------------------------
# Parse args
# ---------------------------------------------------------------------------
INSTALL_SYSTEM=true
INSTALL_PYTHON=true
BUILD_VLBURR=true
for arg in "$@"; do
    case "$arg" in
        --no-system) INSTALL_SYSTEM=false ;;
        --skip-python) INSTALL_PYTHON=false ;;
        --skip-vlburr) BUILD_VLBURR=false ;;
        *) echo "Unknown arg: $arg"; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# System packages
# ---------------------------------------------------------------------------
if $INSTALL_SYSTEM; then
    echo "=== Installing system packages ==="
    apt-get update && apt-get install -y \
        build-essential cmake git python3 python3-dev python3-pip python3-venv \
        zstd zlib1g-dev \
        llvm-17 clang-17 libc++-17-dev libc++abi-17-dev libtbb-dev
fi

# ---------------------------------------------------------------------------
# Python venv
# ---------------------------------------------------------------------------
echo "=== Setting up Python venv ==="
if $INSTALL_PYTHON; then
    if [ ! -d .venv ]; then
        python3 -m venv .venv
    fi
    source .venv/bin/activate
    pip install --upgrade pip
    pip install -r requirements.txt
else
    source .venv/bin/activate
fi

# ---------------------------------------------------------------------------
# Build CaramelDB: both the static library (linked by the genomics harness)
# and the carameldb Python module (imported by every other experiment).
# ---------------------------------------------------------------------------
CARAMEL_DIR="$SCRIPT_DIR/deps/CaramelDB"
if [ ! -d "$CARAMEL_DIR" ]; then
    echo "ERROR: CaramelDB not found at $CARAMEL_DIR" >&2
    echo "  Initialize the submodule: git submodule update --init deps/CaramelDB" >&2
    exit 1
fi

# This also builds deps/CaramelDB/build/libcaramel_lib.a as a side effect, which
# is the static library the genomics harness links; no separate cmake step.
echo "=== Installing the carameldb Python module ==="
CARAMEL_BUILD_JOBS="${JOBS:-4}" \
CMAKE_ARGS="-DCMAKE_INTERPROCEDURAL_OPTIMIZATION=OFF -DCMAKE_POLICY_VERSION_MINIMUM=3.5" \
  .venv/bin/pip install "$CARAMEL_DIR/cython"

# ---------------------------------------------------------------------------
# VL-BuRR from its pinned pristine upstream checkout. Upstream includes x86
# SSE/AVX intrinsics (immintrin.h) and does not compile on ARM, so on Apple
# Silicon this is skipped rather than failing the whole setup: everything
# except the VL-BuRR arm of the synthetic comparison still works natively.
# ---------------------------------------------------------------------------
case "$(uname -m)" in
    arm64|aarch64) SUPPORTS_VLBURR=false ;;
    *) SUPPORTS_VLBURR=true ;;
esac

if ! $BUILD_VLBURR; then
    echo "=== Skipping VL-BuRR (--skip-vlburr) ==="
elif ! $SUPPORTS_VLBURR; then
    echo ""
    echo "=== Skipping VL-BuRR: $(uname -m) is not supported ==="
    echo "  VL-BuRR uses x86 SSE/AVX intrinsics and cannot be built on ARM."
    echo "  Everything else is set up and runs natively:"
    echo "    make bound-validation"
    echo "    .venv/bin/python experiments/synthetic_comparison/run.py --component local"
    echo "  For VL-BuRR, and for the reference environment, use the container:"
    echo "    ./scripts/docker_run.sh make synthetic-comparison"
else
    echo "=== Building LSF (ribbon_learned_bench) ==="
    ./methods/vlburr/build.sh
fi

echo ""
echo "=== Setup complete ==="
echo "Activate the venv with: source .venv/bin/activate"
