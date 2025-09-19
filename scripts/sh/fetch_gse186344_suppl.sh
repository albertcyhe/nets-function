#!/usr/bin/env bash
set -euo pipefail

# Download GSE186344 supplementary files directly from GEO (bypass corrupted RAW tar)
# Uses aria2c for robust, parallel downloads.
# Requires: aria2c (brew install aria2) and curl

DEST="data/raw/GSE186344_suppl"
URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186344/suppl/"

mkdir -p "${DEST}"

if ! command -v aria2c >/dev/null 2>&1; then
  echo "[error] aria2c not found. Please install it, e.g.: brew install aria2" >&2
  exit 1
fi
if ! command -v curl >/dev/null 2>&1; then
  echo "[error] curl not found. Please install curl." >&2
  exit 1
fi

echo "Listing files from ${URL} ..."

TMP_DIR=$(mktemp -d)
TOKENS="${TMP_DIR}/tokens.txt"
URLS="${TMP_DIR}/urls.txt"

# Fetch filelist.txt and extract URLs or filenames
curl -fsSL "${URL}/filelist.txt" -o "${TMP_DIR}/filelist.txt" || true
awk '{ for (i=1;i<=NF;i++) print $i }' "${TMP_DIR}/filelist.txt" \
  | awk '/^(https?|ftp):\/\/.*(csv\.gz|tsv\.gz|h5|mtx\.gz)$|^GSM[0-9]+.*(csv\.gz|tsv\.gz|h5|mtx\.gz)$/' \
  > "${TOKENS}" || true

> "${URLS}"
while IFS= read -r tok; do
  [ -z "${tok}" ] && continue
  if [[ "${tok}" =~ ^https?:// ]] || [[ "${tok}" =~ ^ftp:// ]]; then
    echo "${tok}" >> "${URLS}"
  else
    gsm=$(echo "${tok}" | grep -oE '^GSM[0-9]+')
    if [ -n "${gsm}" ]; then
      grp=$(echo "${gsm}" | sed -E 's/^GSM([0-9]{4}).*/GSM\1nnn/')
      echo "https://ftp.ncbi.nlm.nih.gov/geo/samples/${grp}/${gsm}/suppl/${tok}" >> "${URLS}"
    fi
  fi
done < "${TOKENS}"

if [ ! -s "${URLS}" ]; then
  echo "[error] Could not derive supplementary file URLs from GEO filelist.txt. Please provide one GSM ID to build URLs." >&2
  rm -rf "${TMP_DIR}"
  exit 1
fi

echo "Will download $(wc -l <"${URLS}") files to ${DEST}"
aria2c -c -x 16 -s 16 --max-tries=0 --retry-wait=30 -d "${DEST}" -i "${URLS}"

rm -rf "${TMP_DIR}"
echo "Done. Files are under ${DEST}"
