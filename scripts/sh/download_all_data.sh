#!/usr/bin/env bash
set -euo pipefail

# One-click retriever for core datasets and references used in this repo.
# Requires: curl, aria2c, jq, python (proteomics env optional for MET500 step)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

PRIDE_IDS=(
  PXD011796
  PXD005719
  PXD046330
  PXD051579
  PXD032767
  PXD018301
)

ensure_bin() {
  local name="$1"
  if ! command -v "$name" >/dev/null 2>&1; then
    echo "[error] Missing required tool: $name" >&2
    exit 1
  fi
}

download_file() {
  # download_file URL DEST
  local url="$1"; shift
  local dest="$1"; shift
  local dest_dir
  dest_dir="$(dirname "$dest")"
  mkdir -p "$dest_dir"
  # Use curl with resume and retry
  curl -L --fail --retry 3 --retry-delay 5 -C - -o "$dest" "$url"
}

fetch_pride() {
  local pxd="$1"
  echo "[info] Fetching PRIDE ${pxd} ..."
  bash "${REPO_ROOT}/scripts/sh/fetch_pride_dataset.sh" "$pxd" --categories RAW,RESULT --dest "${REPO_ROOT}/data/raw/proteomics/${pxd}" || true
}

fetch_geo_series_matrix() {
  # fetch_geo_series_matrix GSExxxxxx
  local gse="$1"
  local num="${gse#GSE}"
  if [[ ${#num} -lt 3 ]]; then
    echo "[warn] Unexpected GSE id: $gse" >&2
    return 0
  fi
  local head="${num::-3}"
  local family="GSE${head}nnn"
  local base="https://ftp.ncbi.nlm.nih.gov/geo/series/${family}/${gse}/matrix"
  local out_dir="${REPO_ROOT}/data/raw/transcriptomics/${gse}"
  mkdir -p "$out_dir"
  echo "[info] Downloading ${gse} series_matrix.txt.gz ..."
  download_file "${base}/${gse}_series_matrix.txt.gz" "${out_dir}/${gse}_series_matrix.txt.gz" || true
}

fetch_gpl_annot() {
  # fetch_gpl_annot GPLnnn
  local gpl="$1"
  local num="${gpl#GPL}"
  # Folders chunk by 3 digits for platforms too
  local head="${num%???}"
  local platform_base="https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL${head}nnn/${gpl}/annot"
  local out_dir="${REPO_ROOT}/data/raw/transcriptomics/${gpl}"
  mkdir -p "$out_dir"
  echo "[info] Downloading ${gpl}.annot.gz ..."
  download_file "${platform_base}/${gpl}.annot.gz" "${out_dir}/${gpl}.annot.gz" || true
}

main() {
  ensure_bin curl
  ensure_bin aria2c
  ensure_bin jq

  # PRIDE datasets
  for pxd in "${PRIDE_IDS[@]}"; do
    fetch_pride "$pxd"
  done

  # Metabolomics workbench studies
  echo "[info] Downloading Metabolomics Workbench studies ..."
  mkdir -p "${REPO_ROOT}/data/raw/metabolomics"
  aria2c -c -x 16 -s 16 --max-tries=0 --retry-wait=30 \
    -i "${REPO_ROOT}/resources/manifests/metabolomics_urls.txt" \
    -d "${REPO_ROOT}/data/raw/metabolomics" || true

  # Transcriptomics GEO series matrices + GPLs
  fetch_geo_series_matrix GSE12237
  fetch_gpl_annot GPL96
  fetch_geo_series_matrix GSE125989
  fetch_gpl_annot GPL571

  # Single-cell supplementary for GSE186344
  echo "[info] Downloading GSE186344 supplementary files ..."
  bash "${REPO_ROOT}/scripts/sh/fetch_gse186344_suppl.sh" || true

  # UCSC Xena MET500 (Python script handles URLs)
  echo "[info] Downloading MET500 (Xena) ..."
  python3 "${REPO_ROOT}/scripts/py/met500_fetch.py" || true

  # References: UniProt FASTA and Gencode GTF
  echo "[info] Downloading UniProt human proteome FASTA (UP000005640) ..."
  download_file \
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz" \
    "${REPO_ROOT}/data/raw/reference/UP000005640_9606.fasta.gz" || true

  echo "[info] Downloading Gencode v43 basic GTF ..."
  download_file \
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.basic.annotation.gtf.gz" \
    "${REPO_ROOT}/resources/annotation/gencode/gencode.v43.basic.annotation.gtf.gz" || true

  # Study supplements
  echo "[info] Downloading PXD032767 Supplementary Table S1 (labels) ..."
  download_file \
    "https://pmc.ncbi.nlm.nih.gov/articles/instance/9639356/bin/vdac161_suppl_supplementary_table_s1.xlsx" \
    "${REPO_ROOT}/data/interim/proteomics/PXD032767/metadata/vdac161_suppl_supplementary_table_s1.xlsx" || true

  echo "[note] PXD018301 Human512Reports.xlsx is not auto-downloaded. Place it under:"
  echo "       ${REPO_ROOT}/data/interim/proteomics/PXD018301/metadata/Human512Reports.xlsx"

  # Optional heavy downloads (instructions only)
  cat <<'NOTE'
[note] Optional: GSE96860 bigWig tracks are large and not fetched automatically.
       If needed, identify the bigWig URLs from GEO and place them under:
       data/raw/transcriptomics/GSE96860/ then run:
         python scripts/py/process_gse_bigwig.py \
           --bigwig-dir data/raw/transcriptomics/GSE96860 \
           --gtf resources/annotation/gencode/gencode.v43.basic.annotation.gtf.gz \
           --dataset GSE96860 \
           --output-dir data/processed/GSE96860
NOTE

  echo "[done] Data retrieval completed (with notes for optional/manual items)."
}

main "$@"

