#!/usr/bin/env bash
# Convert Thermo RAW files to (indexed) mzML using ThermoRawFileParser with GNU parallel.
# Usage: scripts/sh/run_trfp_parallel.sh <PXD_ID> [jobs]
# Example: scripts/sh/run_trfp_parallel.sh PXD011796 4

set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <PXD_ID> [jobs]" >&2
  exit 1
fi

PXD="$1"
JOBS="${2:-4}"
PROJECT_ROOT="$(cd "$(dirname "$0")"/../.. && pwd)"
RAW_DIR="${PROJECT_ROOT}/data/raw/proteomics/${PXD}"
OUT_DIR="${PROJECT_ROOT}/data/interim/proteomics/${PXD}/mzML"
LOG_DIR="${PROJECT_ROOT}/logs/proteomics"
RUN_TS="$(date +%Y%m%d_%H%M%S)"
JOB_LOG="${LOG_DIR}/trfp_${PXD}_${RUN_TS}_parallel.tsv"
SUMMARY_LOG="${LOG_DIR}/trfp_${PXD}_${RUN_TS}.log"

mkdir -p "${OUT_DIR}" "${LOG_DIR}"

# shellcheck source=/dev/null
source "${PROJECT_ROOT}/env/activate_proteomics.sh"

if ! command -v ThermoRawFileParser >/dev/null 2>&1; then
  echo "ThermoRawFileParser not found in active environment" >&2
  exit 1
fi

if ! command -v parallel >/dev/null 2>&1; then
  echo "GNU parallel is required but was not found" >&2
  exit 1
fi

mapfile -t RAW_FILES < <(find "${RAW_DIR}" -maxdepth 1 -type f -name '*.raw' -print | sort)
TOTAL=${#RAW_FILES[@]}
if [[ ${TOTAL} -eq 0 ]]; then
  echo "No RAW files found in ${RAW_DIR}" >&2
  exit 1
fi

echo "[${RUN_TS}] Converting ${TOTAL} RAW files from ${RAW_DIR} -> ${OUT_DIR}" | tee -a "${SUMMARY_LOG}"
START_TIME=$(date +%s)

export OUT_DIR LOG_DIR SUMMARY_LOG PXD

parallel --eta --jobs "${JOBS}" --joblog "${JOB_LOG}" \
  'start_ts=$(date +%H:%M:%S);
   echo "[$start_ts] {/.}: starting" | tee -a "$SUMMARY_LOG";
   ThermoRawFileParser -i {} -o "$OUT_DIR" -f mzML -g -l 1 \
     > "$LOG_DIR/trfp_${PXD}_{/.}.log" 2>&1;
   status=$?;
   end_ts=$(date +%H:%M:%S);
   if [[ $status -eq 0 ]]; then
     echo "[$end_ts] {/.}: done" | tee -a "$SUMMARY_LOG";
   else
     echo "[$end_ts] {/.}: FAILED (exit $status)" | tee -a "$SUMMARY_LOG";
   fi;
   exit $status' ::: "${RAW_FILES[@]}"

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

printf '[%s] Finished conversion in %02dh:%02dm:%02ds\n' "$(date +%Y%m%d_%H%M%S)" \
  $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60)) | tee -a "${SUMMARY_LOG}"
