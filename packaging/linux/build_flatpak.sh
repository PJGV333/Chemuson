#!/usr/bin/env bash
set -euo pipefail

# Builder local para bundle Flatpak de Chemuson.

VERSION="${1:?missing VERSION}"
BRANCH="${2:-stable}"
OUT_DIR="${3:-dist-flatpak}"
MANIFEST_PATH="${4:-packaging/flatpak/io.github.PJGV333.Chemuson.yml}"
APP_ID="${APP_ID:-io.github.PJGV333.Chemuson}"
ARCH="${ARCH:-x86_64}"
RUNTIME_REPO="${CHEMUSON_FLATPAK_RUNTIME_REPO:-https://dl.flathub.org/repo/flathub.flatpakrepo}"
REPO_TITLE="${CHEMUSON_FLATPAK_REPO_TITLE:-Chemuson (${BRANCH})}"
REPO_COMMENT="${CHEMUSON_FLATPAK_REPO_COMMENT:-Canal oficial Flatpak de Chemuson (${BRANCH}).}"
REPO_DESCRIPTION="${CHEMUSON_FLATPAK_REPO_DESCRIPTION:-Repositorio oficial Flatpak de Chemuson para el canal ${BRANCH}.}"
REPO_HOMEPAGE="${CHEMUSON_FLATPAK_HOMEPAGE:-https://github.com/PJGV333/Chemuson}"
REPO_ICON_URL="${CHEMUSON_FLATPAK_ICON_URL:-}"
REPO_GPG_KEY_ID="${CHEMUSON_FLATPAK_GPG_KEY_ID:-}"
REPO_GPG_HOMEDIR="${CHEMUSON_FLATPAK_GPG_HOMEDIR:-}"
REPO_PUBLIC_KEY_FILE="${CHEMUSON_FLATPAK_PUBLIC_KEY_FILE:-}"
REPO_CONFIG_BASENAME="${CHEMUSON_FLATPAK_CONFIG_BASENAME:-Chemuson}"
REPO_REMOTE_NAME="${CHEMUSON_FLATPAK_REMOTE_NAME:-chemuson-${BRANCH}}"

BUILD_DIR="${OUT_DIR}/build-dir"
REPO_DIR="${OUT_DIR}/repo"
BUNDLE_PATH="${OUT_DIR}/Chemuson-v${VERSION}-linux-${ARCH}.flatpak"

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

BUILDER_ARGS=(
  --user
  --force-clean
  --default-branch="${BRANCH}"
  --install-deps-from=flathub
  --repo="${REPO_DIR}"
)
if [[ -n "${REPO_GPG_KEY_ID}" ]]; then
  BUILDER_ARGS+=("--gpg-sign=${REPO_GPG_KEY_ID}")
fi
if [[ -n "${REPO_GPG_HOMEDIR}" ]]; then
  BUILDER_ARGS+=("--gpg-homedir=${REPO_GPG_HOMEDIR}")
fi

set +e
flatpak-builder "${BUILDER_ARGS[@]}" "${BUILD_DIR}" "${MANIFEST_PATH}" \
  2>&1 | tee "${BUILD_LOG}"
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

REPO_URL="${CHEMUSON_FLATPAK_REPO_URL:-}"
UPDATE_REPO_ARGS=(
  "--title=${REPO_TITLE}"
  "--comment=${REPO_COMMENT}"
  "--description=${REPO_DESCRIPTION}"
  "--homepage=${REPO_HOMEPAGE}"
  "--default-branch=${BRANCH}"
)
if [[ -n "${REPO_ICON_URL}" ]]; then
  UPDATE_REPO_ARGS+=("--icon=${REPO_ICON_URL}")
fi
if [[ -n "${REPO_GPG_KEY_ID}" ]]; then
  UPDATE_REPO_ARGS+=("--gpg-sign=${REPO_GPG_KEY_ID}")
fi
if [[ -n "${REPO_GPG_HOMEDIR}" ]]; then
  UPDATE_REPO_ARGS+=("--gpg-homedir=${REPO_GPG_HOMEDIR}")
fi
flatpak build-update-repo "${UPDATE_REPO_ARGS[@]}" "${REPO_DIR}"

BUNDLE_ARGS=()
if [[ -n "${REPO_URL}" ]]; then
  BUNDLE_ARGS+=("--repo-url=${REPO_URL}")
  BUNDLE_ARGS+=("--runtime-repo=${RUNTIME_REPO}")
fi
if [[ -n "${REPO_PUBLIC_KEY_FILE}" ]]; then
  BUNDLE_ARGS+=("--gpg-keys=${REPO_PUBLIC_KEY_FILE}")
fi
flatpak build-bundle "${BUNDLE_ARGS[@]}" "${REPO_DIR}" "${BUNDLE_PATH}" "${APP_ID}" "${BRANCH}"
rm -f "${BUILD_LOG}"

if [[ -n "${REPO_URL}" ]]; then
  REMOTE_FILE_ARGS=(
    "--repo-url" "${REPO_URL}"
    "--app-id" "${APP_ID}"
    "--branch" "${BRANCH}"
    "--out-dir" "${OUT_DIR}"
    "--basename" "${REPO_CONFIG_BASENAME}"
    "--repo-title" "${REPO_TITLE}"
    "--ref-title" "${REPO_TITLE}"
    "--homepage" "${REPO_HOMEPAGE}"
    "--comment" "${REPO_COMMENT}"
    "--description" "${REPO_DESCRIPTION}"
    "--default-branch" "${BRANCH}"
    "--runtime-repo" "${RUNTIME_REPO}"
    "--suggest-remote-name" "${REPO_REMOTE_NAME}"
  )
  if [[ -n "${REPO_ICON_URL}" ]]; then
    REMOTE_FILE_ARGS+=("--icon-url" "${REPO_ICON_URL}")
  fi
  if [[ -n "${REPO_PUBLIC_KEY_FILE}" ]]; then
    REMOTE_FILE_ARGS+=("--gpg-key-file" "${REPO_PUBLIC_KEY_FILE}")
  fi
  python packaging/release/generate_flatpak_remote_files.py "${REMOTE_FILE_ARGS[@]}"
  python packaging/release/validate_flatpak_remote_artifacts.py \
    --mode build-output \
    --root "${OUT_DIR}" \
    --channel "${BRANCH}" \
    --basename "${REPO_CONFIG_BASENAME}"
fi

if [[ ! -s "${REPO_DIR}/config" ]] || [[ ! -s "${REPO_DIR}/summary" ]]; then
  echo "ERROR: Flatpak repo metadata was not generated correctly in ${REPO_DIR}." >&2
  exit 1
fi

echo "Built ${BUNDLE_PATH}"
