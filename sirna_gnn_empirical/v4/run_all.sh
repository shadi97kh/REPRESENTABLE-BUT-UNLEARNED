#!/usr/bin/env bash
set -euo pipefail
SIRNA_ROBUSTNESS_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 CUBLAS_WORKSPACE_CONFIG=:4096:8
exec /opt/miniforge3/bin/python "$SIRNA_ROBUSTNESS_DIR/runner.py" "$@"
