#!/usr/bin/env bash
# Reproduced: 2026-09-25T12:59:25.846880
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="$(dirname "$SCRIPT_DIR")"
: "${CLAWBIO_ROOT:=C:\Users\manju\.claude\plugins\cache\clawbio\clawbio\0.7.1}"
if [ ! -d "$CLAWBIO_ROOT" ]; then
  echo "Invalid CLAWBIO_ROOT: $CLAWBIO_ROOT" >&2
  exit 1
fi

"${PYTHON:-python3}" "$CLAWBIO_ROOT/skills\skill-builder\skill_builder.py" \
  --input \
  "tumor_hot_cold_spec.json" \
  --output \
  "$OUTPUT_DIR"
