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

if [[ ! -f "${MANIFEST_PATH}" ]]; then
  echo "ERROR: Flatpak manifest not found: ${MANIFEST_PATH}" >&2
  echo "Hint: pass a valid 4th argument or set FLATPAK_MANIFEST in build_linux_all.sh." >&2
  exit 1
fi

if ! flatpak --user remote-add --if-not-exists flathub \
  https://dl.flathub.org/repo/flathub.flatpakrepo; then
  echo "ERROR: Unable to configure flathub remote." >&2
  echo "Hint: check network/DNS and Flatpak user installation permissions." >&2
  exit 1
fi

BUILD_LOG="$(mktemp "${OUT_DIR}/flatpak-builder.XXXXXX.log")"

set +e
flatpak-builder \
  --user \
  --force-clean \
  --default-branch="${BRANCH}" \
  --install-deps-from=flathub \
  --repo="${REPO_DIR}" \
  "${BUILD_DIR}" \
  "${MANIFEST_PATH}" 2>&1 | tee "${BUILD_LOG}"
FLATPAK_BUILDER_EXIT="${PIPESTATUS[0]}"
set -e

if [[ "${FLATPAK_BUILDER_EXIT}" -ne 0 ]]; then
  if grep -Eq "Could not resolve hostname|Temporary failure in name resolution|Name or service not known|Unable to load summary from remote" "${BUILD_LOG}"; then
    echo "ERROR: Flatpak build failed due to network/DNS or sandbox egress restrictions." >&2
    echo "Hint: verify internet access from flatpak-builder sandbox and flathub reachability." >&2
  elif grep -Eq "No matching distribution found|Could not find a version that satisfies the requirement|ModuleNotFoundError: No module named|ERROR: .*pip.*failed" "${BUILD_LOG}"; then
    echo "ERROR: Flatpak build failed while installing Python dependencies inside flatpak-builder." >&2
    echo "Hint: verify pinned wheels/sha256 in the manifest and Python ABI compatibility with the selected runtime." >&2
  else
    echo "ERROR: flatpak-builder failed (exit ${FLATPAK_BUILDER_EXIT})." >&2
  fi
  echo "Build log: ${BUILD_LOG}" >&2
  exit "${FLATPAK_BUILDER_EXIT}"
fi

flatpak build-bundle "${REPO_DIR}" "${BUNDLE_PATH}" "${APP_ID}" "${BRANCH}"
rm -f "${BUILD_LOG}"

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
