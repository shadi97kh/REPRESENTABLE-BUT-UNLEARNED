#!/usr/bin/env bash
set -euo pipefail
SIRNA_REVISION_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONDONTWRITEBYTECODE=1
exec /opt/miniforge3/bin/python "$SIRNA_REVISION_DIR/runner.py" "$@"
