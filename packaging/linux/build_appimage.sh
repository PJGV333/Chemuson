#!/usr/bin/env bash
set -euo pipefail

# Builder AppImage de Chemuson.
# Si existe un AppImage previo lo reutiliza; en caso contrario
# envuelve el binario Linux portable y además genera metadata de update.

VERSION="${1:?missing VERSION}"
DIST_DIR="${2:-dist}"
OUT_DIR="${3:-dist-appimage}"
OWNER="${4:-PJGV333}"
REPO="${5:-Chemuson}"
CHANNEL="${6:-stable}"
TAG="${7:-v${VERSION}}"

APPIMAGE_NAME="Chemuson-v${VERSION}-linux-x86_64.AppImage"
APPIMAGE_PATH="${OUT_DIR}/${APPIMAGE_NAME}"
DOWNLOAD_BASE_URL="${APPIMAGE_DOWNLOAD_BASE_URL:-https://github.com/${OWNER}/${REPO}/releases/download/${TAG}}"

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
  "channel": "${CHANNEL}",
  "tag": "${TAG}",
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
