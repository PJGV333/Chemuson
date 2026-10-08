#!/usr/bin/env bash
set -euo pipefail

# Builder del ejecutable portable Linux de Chemuson.
# Conserva el sufijo histórico .AppImage por compatibilidad del updater; no
# crea un contenedor AppImage Type 2. También genera metadata de update.

VERSION="${1:?missing VERSION}"
DIST_DIR="${2:-dist}"
OUT_DIR="${3:-dist-appimage}"
OWNER="${4:-PJGV333}"
REPO="${5:-Chemuson}"
CHANNEL="${6:-stable}"
TAG="${7:-v${VERSION}}"
SOURCE_SHA="${8:-}"
BUILD_TYPE="${9:-release}"
SOURCE_BRANCH="${10:-}"

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
if [[ "${BUILD_TYPE}" == "preview" ]]; then
  python "${ROOT_DIR}/packaging/release/preview_build.py" verify \
    --version "${VERSION}" --git-sha "${SOURCE_SHA}" \
    --source-branch "${SOURCE_BRANCH}" >/dev/null
  APPIMAGE_NAME="Chemuson-v${VERSION}-preview-${SOURCE_SHA:0:8}-linux-x86_64.AppImage"
elif [[ "${BUILD_TYPE}" == "release" ]]; then
  python "${ROOT_DIR}/packaging/release/release_policy.py" \
    --version "${VERSION}" --channel "${CHANNEL}" --tag "${TAG}" >/dev/null
  APPIMAGE_NAME="Chemuson-v${VERSION}-linux-x86_64.AppImage"
  DOWNLOAD_BASE_URL="${APPIMAGE_DOWNLOAD_BASE_URL:-https://github.com/${OWNER}/${REPO}/releases/download/${TAG}}"
else
  echo "Unsupported build type: ${BUILD_TYPE}" >&2
  exit 2
fi
APPIMAGE_PATH="${OUT_DIR}/${APPIMAGE_NAME}"

mkdir -p "$OUT_DIR"

if compgen -G "${DIST_DIR}/*.AppImage" > /dev/null; then
  src="$(ls -1 "${DIST_DIR}"/*.AppImage | head -n1)"
  cp "$src" "$APPIMAGE_PATH"
  chmod +x "$APPIMAGE_PATH"
elif [[ -x "${DIST_DIR}/Chemuson" ]]; then
  cp "${DIST_DIR}/Chemuson" "$APPIMAGE_PATH"
  chmod +x "$APPIMAGE_PATH"
else
  echo "No Linux artifact found in ${DIST_DIR}" >&2
  exit 1
fi

# Preview outputs are Actions-only and must never contain public updater metadata.
if [[ "${BUILD_TYPE}" == "preview" ]]; then
  echo "Built preview portable executable ${APPIMAGE_PATH} (no public updater metadata)."
  exit 0
fi

# Metadata AppImageUpdate (best-effort para releases GitHub).
UPDATE_TRACK="latest"
if [[ "$CHANNEL" == "beta" ]]; then
  UPDATE_TRACK="prerelease"
fi
APPIMAGE_UPDATE_INFO="gh-releases-zsync|${OWNER}|${REPO}|${UPDATE_TRACK}|${APPIMAGE_NAME}.zsync"
printf '%s\n' "$APPIMAGE_UPDATE_INFO" > "${APPIMAGE_PATH}.updateinfo"

cat > "${APPIMAGE_PATH}.update.json" <<EOF
{
  "asset": "${APPIMAGE_NAME}",
  "owner": "${OWNER}",
  "repository": "${REPO}",
  "version": "${VERSION}",
  "channel": "${CHANNEL}",
  "tag": "${TAG}",
  "source_sha": "${SOURCE_SHA}",
  "download_base_url": "${DOWNLOAD_BASE_URL}",
  "appimage_update_information": "${APPIMAGE_UPDATE_INFO}"
}
EOF

# Si está disponible zsyncmake, genera sidecar para delta-updates.
if command -v zsyncmake >/dev/null 2>&1; then
  zsyncmake -u "${DOWNLOAD_BASE_URL}/${APPIMAGE_NAME}" \
    -o "${APPIMAGE_PATH}.zsync" \
    "$APPIMAGE_PATH" || true
fi

echo "Built ${APPIMAGE_PATH}"
