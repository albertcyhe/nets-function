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

ENV_NAME="proteomics"
ENV_SPEC="${ENV_DIR}/proteomics.yml"

PHILOSOPHER_VERSION="v5.1.0"
PHILOSOPHER_ASSET="philosopher_${PHILOSOPHER_VERSION}_linux_amd64.zip"
PHILOSOPHER_URL="https://github.com/Nesvilab/philosopher/releases/download/${PHILOSOPHER_VERSION}/${PHILOSOPHER_ASSET}"
PHILOSOPHER_ARCHIVE="${ENV_DIR}/${PHILOSOPHER_ASSET}"
PHILOSOPHER_DIR="${ENV_DIR}/philosopher"
PHILOSOPHER_BIN="${PHILOSOPHER_DIR}/philosopher"

command_exists() {
  command -v "$1" >/dev/null 2>&1
}

download() {
  local url="$1"
  local dest="$2"

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

ensure_philosopher() {
  mkdir -p "${ENV_DIR}" "${BIN_DIR}"

  if [[ ! -x "${PHILOSOPHER_BIN}" ]]; then
    if [[ ! -f "${PHILOSOPHER_ARCHIVE}" ]]; then
      echo "[info] Downloading Philosopher ${PHILOSOPHER_VERSION}"
      download "${PHILOSOPHER_URL}" "${PHILOSOPHER_ARCHIVE}"
    fi

    mkdir -p "${PHILOSOPHER_DIR}"
    "${MICROMAMBA_BIN}" run -n "${ENV_NAME}" python -m zipfile -e "${PHILOSOPHER_ARCHIVE}" "${PHILOSOPHER_DIR}"
    chmod +x "${PHILOSOPHER_BIN}"
  fi

  ln -sf ../philosopher/philosopher "${BIN_DIR}/philosopher"
}

main() {
  ensure_micromamba
  ensure_environment
  ensure_philosopher

  echo "[info] Proteomics toolchain ready"
  echo "[info] Activate with: source <(\"${MICROMAMBA_BIN}\" shell hook --shell bash) && micromamba activate ${ENV_NAME}"
}

main "$@"
