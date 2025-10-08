#!/usr/bin/env bash
set -euo pipefail

# Upload project data to Google Drive using rclone.
# Requires: rclone configured with a Drive remote (default: gdrive)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

REMOTE="gdrive"
DEST_ROOT="nets-fuction"
MODE="copy"     # copy (safe) or sync (mirror; deletes extraneous remote files)
INCLUDE_PROCESSED=false
INCLUDE_RESULTS=false
DRY_RUN=false

usage() {
  cat <<USAGE
Usage: $(basename "$0") [options]

Options:
  --remote NAME        rclone remote name (default: gdrive)
  --dest PATH          Remote base folder under Drive (default: nets-fuction)
  --mode copy|sync     Upload mode: copy (safe) or sync (mirror) (default: copy)
  --include-processed  Also upload data/processed (in addition to data/raw)
  --include-results    Also upload results/ (figures, tables)
  --dry-run            Print planned rclone commands without executing
  -h, --help           Show this help

Examples:
  # Copy data/raw to gdrive:nets-fuction/data/raw
  $(basename "$0")

  # Mirror (sync) raw + processed to a custom folder
  $(basename "$0") --mode sync --include-processed --dest MyProjectBackup
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --remote) REMOTE="$2"; shift 2 ;;
    --dest) DEST_ROOT="$2"; shift 2 ;;
    --mode) MODE="$2"; shift 2 ;;
    --include-processed) INCLUDE_PROCESSED=true; shift ;;
    --include-results) INCLUDE_RESULTS=true; shift ;;
    --dry-run) DRY_RUN=true; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[error] Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

if ! command -v rclone >/dev/null 2>&1; then
  echo "[error] rclone not found. Please install and configure it first." >&2
  exit 1
fi

LOCAL_TARGETS=("data/raw")
${INCLUDE_PROCESSED} && LOCAL_TARGETS+=("data/processed")
${INCLUDE_RESULTS} && LOCAL_TARGETS+=("results")

TS="$(date +%Y%m%d_%H%M%S)"
LOG_DIR="${REPO_ROOT}/logs/uploads"
mkdir -p "${LOG_DIR}"

rclone_cmd() {
  local src="$1"; shift
  local dst="$1"; shift
  local mode="$1"; shift
  local log_file="$1"; shift
  local base_opts=("--transfers" "8" "--checkers" "16" "--drive-chunk-size" "64M" "--drive-stop-on-upload-limit")
  if [[ "$mode" == "sync" ]]; then
    cmd=(rclone sync "${src}" "${dst}" "${base_opts[@]}")
  else
    cmd=(rclone copy "${src}" "${dst}" "${base_opts[@]}")
  fi
  echo "${cmd[*]} --log-file \"${log_file}\""
}

echo "[info] Remote: ${REMOTE}  Dest root: ${DEST_ROOT}  Mode: ${MODE}"
echo "[info] Will upload: ${LOCAL_TARGETS[*]}"

for local in "${LOCAL_TARGETS[@]}"; do
  abs_local="${REPO_ROOT}/${local}"
  if [[ ! -d "${abs_local}" ]]; then
    echo "[warn] Skip missing directory: ${abs_local}"
    continue
  fi
  remote_path="${REMOTE}:${DEST_ROOT}/${local}"
  log_file="${LOG_DIR}/upload_$(echo "${local}" | tr '/ ' '__')_${TS}.log"
  echo "[plan] ${abs_local} -> ${remote_path} (log: ${log_file})"
  if ${DRY_RUN}; then
    rclone_cmd "${abs_local}" "${remote_path}" "${MODE}" "${log_file}"
  else
    # Run in background to survive session; tail command hint provided
    nohup bash -c "$(rclone_cmd "${abs_local}" "${remote_path}" "${MODE}" "${log_file}")" >/dev/null 2>&1 &
    echo "[started] Background upload. Tail logs with: tail -f '${log_file}'"
  fi
done

echo "[done] Planned/started uploads. Logs under ${LOG_DIR}"

