#!/usr/bin/env bash
# Install all dependencies for AutoCSF-Benchmarks on Ubuntu 24.04 x86_64.
#
# Usage:
#   sudo ./setup.sh          # install system packages + build everything
#   ./setup.sh --no-system   # skip apt-get, only build deps (if packages already installed)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# ---------------------------------------------------------------------------
# Parse args
# ---------------------------------------------------------------------------
INSTALL_SYSTEM=true
INSTALL_PYTHON=true
for arg in "$@"; do
    case "$arg" in
        --no-system) INSTALL_SYSTEM=false ;;
        --skip-python) INSTALL_PYTHON=false ;;
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
# Build CaramelDB (expected as sibling directory)
# ---------------------------------------------------------------------------
CARAMEL_DIR="$SCRIPT_DIR/deps/CaramelDB"
if [ -d "$CARAMEL_DIR" ]; then
    echo "=== Building CaramelDB ==="
    # The genomics harness links the C++ library directly. Building only that
    # target avoids pulling the unrelated legacy Python bindings into the
    # camera-ready artifact.
    cmake -S "$CARAMEL_DIR" -B "$CARAMEL_DIR/pybind/build" \
      -DCMAKE_BUILD_TYPE=Release -DCMAKE_POLICY_VERSION_MINIMUM=3.5
    cmake --build "$CARAMEL_DIR/pybind/build" --target caramel_lib -j "${JOBS:-4}"
else
    echo "WARNING: CaramelDB not found at $CARAMEL_DIR"
    echo "  Initialize the submodule: git submodule update --init deps/CaramelDB"
fi

# ---------------------------------------------------------------------------
# VL-BuRR from its pinned pristine upstream checkout.
# ---------------------------------------------------------------------------
echo "=== Building LSF (ribbon_learned_bench) ==="
./vlburr/build.sh

echo ""
echo "=== Setup complete ==="
echo "Activate the venv with: source .venv/bin/activate"
