#!/usr/bin/env bash
set -euo pipefail

# Builder local para bundle Flatpak de Chemuson.

VERSION="${1:?missing VERSION}"
BRANCH="${2:-stable}"
OUT_DIR="${3:-dist-flatpak}"
MANIFEST_PATH="${4:-packaging/flatpak/io.github.PJGV333.Chemuson.yml}"
APP_ID="${APP_ID:-io.github.PJGV333.Chemuson}"
ARCH="${ARCH:-x86_64}"

BUILD_DIR="${OUT_DIR}/build-dir"
REPO_DIR="${OUT_DIR}/repo"
BUNDLE_PATH="${OUT_DIR}/Chemuson-v${VERSION}-linux-${ARCH}.flatpak"
FLATPAKREF_PATH="${OUT_DIR}/Chemuson-${BRANCH}.flatpakref"

mkdir -p "$OUT_DIR"

flatpak --user remote-add --if-not-exists flathub \
  https://dl.flathub.org/repo/flathub.flatpakrepo

flatpak-builder \
  --user \
  --force-clean \
  --default-branch="${BRANCH}" \
  --install-deps-from=flathub \
  --repo="${REPO_DIR}" \
  "${BUILD_DIR}" \
  "${MANIFEST_PATH}"

flatpak build-bundle "${REPO_DIR}" "${BUNDLE_PATH}" "${APP_ID}" "${BRANCH}"

REPO_URL="${CHEMUSON_FLATPAK_REPO_URL:-}"
if [[ -n "${REPO_URL}" ]]; then
  cat > "${FLATPAKREF_PATH}" <<EOF
[Flatpak Ref]
Title=Chemuson (${BRANCH})
Name=${APP_ID}
Branch=${BRANCH}
IsRuntime=false
Url=${REPO_URL}
RuntimeRepo=https://flathub.org/repo/flathub.flatpakrepo
EOF
fi

echo "Built ${BUNDLE_PATH}"
