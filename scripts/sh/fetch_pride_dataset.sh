#!/usr/bin/env bash
set -euo pipefail

# Fetch files from PRIDE (ProteomeXchange) using the public REST API and aria2c.
# Example:
#   scripts/sh/fetch_pride_dataset.sh PXD011796 --categories RAW,RESULT --dest data/raw/PXD011796
# Requires: curl, jq, aria2c

usage() {
  cat <<'USAGE'
Usage: fetch_pride_dataset.sh <PXD_ACCESSION> [options]

Options:
  --categories LIST   Comma-separated fileCategory values to include (default: RAW,RESULT,PEAK,SEARCH)
                      Use "ALL" to download everything exposed via FTP.
  --dest PATH         Output directory (default: data/raw/<accession>)
  --match REGEX       Only keep files whose basename matches REGEX (applied after category filter)
  --dry-run           Print URLs without invoking aria2c
  -h, --help          Show this message

Examples:
  # Download RAW/RESULT files for PXD046330 into data/raw/PXD046330
  fetch_pride_dataset.sh PXD046330

  # Fetch only PEAK files whose name contains "quant" into custom folder
  fetch_pride_dataset.sh PXD011796 --categories PEAK --match quant --dest tmp/PXD011796_peak
USAGE
}

if [[ $# -eq 0 ]]; then
  usage
  exit 1
fi

accession=""
dest=""
categories="RAW,RESULT,PEAK,SEARCH"
match_regex=""
dry_run=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)
      usage
      exit 0
      ;;
    --categories)
      if [[ $# -lt 2 ]]; then
        echo "[error] --categories requires a value" >&2
        exit 1
      fi
      categories="$2"
      shift 2
      ;;
    --dest)
      if [[ $# -lt 2 ]]; then
        echo "[error] --dest requires a path" >&2
        exit 1
      fi
      dest="$2"
      shift 2
      ;;
    --match)
      if [[ $# -lt 2 ]]; then
        echo "[error] --match requires a regex" >&2
        exit 1
      fi
      match_regex="$2"
      shift 2
      ;;
    --dry-run)
      dry_run=true
      shift
      ;;
    PXD*)
      if [[ -n "$accession" ]]; then
        echo "[error] accession already set to $accession" >&2
        exit 1
      fi
      accession="$1"
      shift
      ;;
    *)
      echo "[error] Unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
done

if [[ -z "$accession" ]]; then
  echo "[error] Missing PXD accession" >&2
  usage
  exit 1
fi

if [[ -z "$dest" ]]; then
  dest="data/raw/${accession}"
fi

if ! command -v curl >/dev/null 2>&1; then
  echo "[error] curl not found. Install it first." >&2
  exit 1
fi
if ! command -v jq >/dev/null 2>&1; then
  echo "[error] jq not found. Install it (e.g. brew install jq)." >&2
  exit 1
fi
if ! $dry_run && ! command -v aria2c >/dev/null 2>&1; then
  echo "[error] aria2c not found. Install it (e.g. brew install aria2)." >&2
  exit 1
fi

mkdir -p "$dest"

upper_categories=$(echo "$categories" | tr '[:lower:]' '[:upper:]')
declare -A include
if [[ "$upper_categories" != "ALL" ]]; then
  IFS=',' read -r -a cat_array <<< "$upper_categories"
  for cat in "${cat_array[@]}"; do
    trimmed=$(echo "$cat" | tr -d ' \t\r\n')
    [[ -z "$trimmed" ]] && continue
    include["$trimmed"]=1
  done
fi

limit=200
offset=0
manifest=$(mktemp)
trap 'rm -f "$manifest"' EXIT

fetch_page() {
  local url="https://www.ebi.ac.uk/pride/ws/archive/v2/projects/${accession}/files?offset=${offset}&limit=${limit}"
  curl -fsSL "$url"
}

page="$(fetch_page || true)"
if [[ -z "$page" ]]; then
  echo "[error] Empty response from PRIDE API" >&2
  exit 1
fi

while true; do
  count=$(printf '%s' "$page" | jq 'length')
  if [[ "$count" -eq 0 ]]; then
    break
  fi
  printf '%s' "$page" | jq -r '.[] | [.fileCategory.value, .fileName, (.publicFileLocations[] | select(.name == "FTP Protocol") | .value)] | @tsv' \
    >> "$manifest"
  offset=$((offset + limit))
  page="$(fetch_page || true)"
  if [[ -z "$page" ]]; then
    break
  fi
done

if [[ ! -s "$manifest" ]]; then
  echo "[error] No FTP URLs retrieved for ${accession}." >&2
  exit 1
fi

urls=()
while IFS=$'\t' read -r file_cat file_name ftp_url; do
  [[ -z "$ftp_url" ]] && continue
  cat_up=$(echo "$file_cat" | tr '[:lower:]' '[:upper:]')
  if [[ "$upper_categories" != "ALL" ]]; then
    if [[ -z "${include[$cat_up]:-}" ]]; then
      continue
    fi
  fi
  if [[ -n "$match_regex" ]]; then
    if [[ ! "$file_name" =~ $match_regex ]]; then
      continue
    fi
  fi
  urls+=("$ftp_url")
done < "$manifest"

if [[ ${#urls[@]} -eq 0 ]]; then
  echo "[warn] No files passed the filters for ${accession}." >&2
  exit 0
fi

printf 'Found %d files for %s\n' "${#urls[@]}" "$accession"
if $dry_run; then
  printf '%s\n' "${urls[@]}"
  exit 0
fi

list_file=$(mktemp)
trap 'rm -f "$manifest" "$list_file"' EXIT
printf '%s\n' "${urls[@]}" > "$list_file"

aria2c -c -x 16 -s 16 --max-tries=0 --retry-wait=30 -d "$dest" -i "$list_file"

echo "[done] Files saved to $dest"
