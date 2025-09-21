#!/usr/bin/env bash
# Source this script to load the proteomics analysis environment.
# Usage: source env/activate_proteomics.sh

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")"/.. && pwd)"
export MAMBA_ROOT_PREFIX="${PROJECT_ROOT}/env/.mamba"
MICROMAMBA_BIN="${PROJECT_ROOT}/env/bin/micromamba"

if [[ ! -x "${MICROMAMBA_BIN}" ]]; then
  echo "micromamba executable not found at ${MICROMAMBA_BIN}" >&2
  return 1
fi

# Initialise shell hook and activate environment
if ! command -v micromamba >/dev/null 2>&1; then
  eval "$(${MICROMAMBA_BIN} shell hook -s bash)"
fi

micromamba activate proteomics
