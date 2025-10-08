#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
ENV_DIR="${REPO_ROOT}/env"
BIN_DIR="${ENV_DIR}/bin"
MICROMAMBA_TARBALL="${ENV_DIR}/micromamba.tar.bz2"
MICROMAMBA_ROOT="${ENV_DIR}/micromamba"
MICROMAMBA_BIN="${BIN_DIR}/micromamba"
MICROMAMBA_URL="https://micromamba.snakepit.net/api/micromamba/linux-64/latest"

ENV_NAME="nets-scrna"
ENV_SPEC="${ENV_DIR}/environment-scrna.yml"

command_exists() {
  command -v "$1" >/dev/null 2>&1
}

download() {
  local url="$1"; local dest="$2"
  if command_exists curl; then
    curl -L --fail --retry 3 --retry-delay 5 -o "$dest" "$url"
  elif command_exists wget; then
    wget -O "$dest" "$url"
  else
    echo "[error] Neither curl nor wget is available" >&2
    exit 1
  fi
}

ensure_micromamba() {
  if [[ -x "${MICROMAMBA_BIN}" ]]; then
    return
  fi
  mkdir -p "${ENV_DIR}" "${MICROMAMBA_ROOT}" "${BIN_DIR}"
  echo "[info] Downloading micromamba..."
  download "${MICROMAMBA_URL}" "${MICROMAMBA_TARBALL}"
  tar -xjf "${MICROMAMBA_TARBALL}" -C "${MICROMAMBA_ROOT}" --strip-components=1
  ln -sf ../micromamba/micromamba "${MICROMAMBA_BIN}"
  chmod +x "${MICROMAMBA_BIN}" "${MICROMAMBA_ROOT}/micromamba"
}

ensure_environment() {
  if [[ ! -f "${ENV_SPEC}" ]]; then
    echo "[error] Expected environment spec at ${ENV_SPEC}" >&2
    exit 1
  fi
  if "${MICROMAMBA_BIN}" env list | awk '{print $1}' | grep -Fxq "${ENV_NAME}"; then
    echo "[info] Updating existing micromamba env: ${ENV_NAME}"
    "${MICROMAMBA_BIN}" env update -n "${ENV_NAME}" -f "${ENV_SPEC}"
  else
    echo "[info] Creating micromamba env: ${ENV_NAME}"
    "${MICROMAMBA_BIN}" create -y -n "${ENV_NAME}" -f "${ENV_SPEC}"
  fi
}

main() {
  ensure_micromamba
  ensure_environment
  echo "[info] scRNA toolchain ready"
  echo "[info] Activate with: source <(\"${MICROMAMBA_BIN}\" shell hook --shell bash) && micromamba activate ${ENV_NAME}"
}

main "$@"

