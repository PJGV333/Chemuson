#!/usr/bin/env bash
set -euo pipefail

# Pinned official AppImageKit release asset (x86_64). The asset is verified
# before extraction or execution; a changed upstream payload fails closed.
readonly APPIMAGETOOL_URL="https://github.com/AppImage/AppImageKit/releases/download/continuous/appimagetool-x86_64.AppImage"
readonly APPIMAGETOOL_SHA256="b90f4a8b18967545fda78a445b27680a1642f1ef9488ced28b65398f2be7add2"
readonly APPIMAGETOOL_COMMIT="5735cc5"
readonly APPIMAGETOOL_ASSET_ID="98605504"

CACHE_BASE="${1:?usage: fetch_appimagetool.sh <fresh temporary directory>}"
mkdir -p "${CACHE_BASE}"
TOOL_PATH="${CACHE_BASE}/appimagetool-x86_64.AppImage"
EXTRACT_DIR="${CACHE_BASE}/squashfs-root"

if [[ ! -f "${TOOL_PATH}" ]] || ! printf '%s  %s\n' "${APPIMAGETOOL_SHA256}" "${TOOL_PATH}" | sha256sum --check --status; then
  DOWNLOAD_PATH="${TOOL_PATH}.download.$$"
  trap 'rm -f "${DOWNLOAD_PATH}"' EXIT
  curl --fail --location --silent --show-error --retry 3 \
    --proto '=https' --proto-redir '=https' --tlsv1.2 \
    --output "${DOWNLOAD_PATH}" "${APPIMAGETOOL_URL}"
  printf '%s  %s\n' "${APPIMAGETOOL_SHA256}" "${DOWNLOAD_PATH}" | sha256sum --check --status || {
    echo "ERROR: appimagetool asset SHA-256 mismatch; refusing to execute it." >&2
    exit 1
  }
  chmod 755 "${DOWNLOAD_PATH}"
  mv -f "${DOWNLOAD_PATH}" "${TOOL_PATH}"
  trap - EXIT
fi

printf '%s  %s\n' "${APPIMAGETOOL_SHA256}" "${TOOL_PATH}" | sha256sum --check --status
HEADER_MAGIC="$(od -An -tx1 -j8 -N3 "${TOOL_PATH}" | tr -d ' \n')"
if [[ "${HEADER_MAGIC}" != "414902" ]]; then
  echo "ERROR: pinned appimagetool is not an ELF Type 2 AppImage." >&2
  exit 1
fi

if [[ ! -x "${EXTRACT_DIR}/AppRun" ]]; then
  mkdir -p "${CACHE_BASE}/extract-work"
  (
    cd "${CACHE_BASE}/extract-work"
    timeout 60 "${TOOL_PATH}" --appimage-extract >/dev/null
  )
  mv "${CACHE_BASE}/extract-work/squashfs-root" "${EXTRACT_DIR}"
fi
APPIMAGETOOL_BIN="${EXTRACT_DIR}/AppRun"
if [[ ! -x "${APPIMAGETOOL_BIN}" ]]; then
  echo "ERROR: appimagetool AppRun missing after FUSE-free tool extraction." >&2
  exit 1
fi

TOOL_VERSION="$(timeout 30 "${APPIMAGETOOL_BIN}" --version 2>&1)"
if [[ "${TOOL_VERSION}" != *"commit ${APPIMAGETOOL_COMMIT}"* ]]; then
  printf 'ERROR: unexpected appimagetool version: %s\n' "${TOOL_VERSION}" >&2
  exit 1
fi
printf 'Verified appimagetool asset %s (SHA-256 %s): %s\n' \
  "${APPIMAGETOOL_ASSET_ID}" "${APPIMAGETOOL_SHA256}" "${TOOL_VERSION}" >&2
printf '%s\n' "${APPIMAGETOOL_BIN}"
