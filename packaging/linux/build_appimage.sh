#!/usr/bin/env bash
set -euo pipefail

# Construye un AppImage Type 2 auténtico desde el ejecutable PyInstaller.
# El nombre oficial y el contrato de AppImageUpdate permanecen estables.

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
CHECKED_OUT_SHA="$(git -C "${ROOT_DIR}" rev-parse HEAD | tr '[:upper:]' '[:lower:]')"
SOURCE_SHA="${SOURCE_SHA,,}"
if [[ "${CHECKED_OUT_SHA}" != "${SOURCE_SHA}" ]]; then
  echo "ERROR: AppImage source SHA does not match the checked-out commit." >&2
  exit 2
fi
APP_ID="io.github.PJGV333.Chemuson"
ICON_SOURCE="${ROOT_DIR}/packaging/flatpak/${APP_ID}.svg"
DESKTOP_TEMPLATE="${ROOT_DIR}/packaging/linux/appimage/${APP_ID}.desktop.in"
METADATA_SOURCE="${ROOT_DIR}/packaging/flatpak/${APP_ID}.metainfo.xml"

if [[ ! "${SOURCE_SHA}" =~ ^([0-9a-fA-F]{40}|[0-9a-fA-F]{64})$ ]]; then
  echo "ERROR: AppImage builds require the full source Git SHA." >&2
  exit 2
fi

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

SOURCE_BINARY="${DIST_DIR}/Chemuson"
APPIMAGE_PATH="${OUT_DIR}/${APPIMAGE_NAME}"
mkdir -p "${OUT_DIR}"
if [[ ! -s "${SOURCE_BINARY}" || ! -x "${SOURCE_BINARY}" ]]; then
  echo "ERROR: PyInstaller executable missing or empty: ${SOURCE_BINARY}" >&2
  exit 1
fi
for required in "${ICON_SOURCE}" "${DESKTOP_TEMPLATE}" "${METADATA_SOURCE}"; do
  if [[ ! -s "${required}" ]]; then
    echo "ERROR: AppImage resource is missing or empty: ${required}" >&2
    exit 1
  fi
done

APPDIR="$(mktemp -d "${OUT_DIR}/.chemuson-appdir.XXXXXX")"
APPIMAGETOOL_WORK="$(mktemp -d "${TMPDIR:-/tmp}/chemuson-appimagetool.XXXXXX")"
cleanup() {
  rm -rf -- "${APPDIR}" "${APPIMAGETOOL_WORK}"
}
trap cleanup EXIT

mkdir -p \
  "${APPDIR}/usr/bin" \
  "${APPDIR}/usr/share/applications" \
  "${APPDIR}/usr/share/icons/hicolor/scalable/apps" \
  "${APPDIR}/usr/share/metainfo"
install -m755 "${ROOT_DIR}/packaging/linux/appimage/AppRun" "${APPDIR}/AppRun"
install -m755 "${SOURCE_BINARY}" "${APPDIR}/usr/bin/Chemuson"
sed "s/@APP_VERSION@/${VERSION}/g" "${DESKTOP_TEMPLATE}" > "${APPDIR}/${APP_ID}.desktop"
install -m644 "${APPDIR}/${APP_ID}.desktop" "${APPDIR}/usr/share/applications/${APP_ID}.desktop"
install -m644 "${ICON_SOURCE}" "${APPDIR}/${APP_ID}.svg"
install -m644 "${ICON_SOURCE}" "${APPDIR}/.DirIcon"
install -m644 "${ICON_SOURCE}" "${APPDIR}/usr/share/icons/hicolor/scalable/apps/${APP_ID}.svg"
install -m644 "${METADATA_SOURCE}" "${APPDIR}/usr/share/metainfo/${APP_ID}.appdata.xml"

if [[ -n "${APPIMAGETOOL_BIN:-}" ]]; then
  # Test/local override only. Official workflows do not set this; they use the
  # pinned downloader below and the subsequent full package validator.
  APPIMAGETOOL="${APPIMAGETOOL_BIN}"
else
  APPIMAGETOOL="$(bash "${ROOT_DIR}/packaging/linux/fetch_appimagetool.sh" "${APPIMAGETOOL_WORK}")"
fi
if [[ ! -x "${APPIMAGETOOL}" ]]; then
  echo "ERROR: verified appimagetool executable is unavailable." >&2
  exit 1
fi

if [[ "${BUILD_TYPE}" == "release" ]]; then
  UPDATE_TRACK="latest"
  if [[ "${CHANNEL}" == "beta" ]]; then
    UPDATE_TRACK="prerelease"
  fi
  UPDATE_INFO="gh-releases-zsync|${OWNER}|${REPO}|${UPDATE_TRACK}|${APPIMAGE_NAME}.zsync"
  ARCH=x86_64 "${APPIMAGETOOL}" --updateinformation "${UPDATE_INFO}" "${APPDIR}" "${APPIMAGE_PATH}"
  printf '%s\n' "${UPDATE_INFO}" > "${APPIMAGE_PATH}.updateinfo"

  if [[ -n "${APPIMAGETOOL_ZSYNCMAKE:-}" ]]; then
    ZSYNCMAKE_BIN="${APPIMAGETOOL_ZSYNCMAKE}"
  elif [[ -x "$(dirname "${APPIMAGETOOL}")/usr/bin/zsyncmake" ]]; then
    ZSYNCMAKE_BIN="$(dirname "${APPIMAGETOOL}")/usr/bin/zsyncmake"
  else
    ZSYNCMAKE_BIN="$(command -v zsyncmake || true)"
  fi
  if [[ -z "${ZSYNCMAKE_BIN}" || ! -x "${ZSYNCMAKE_BIN}" ]]; then
    echo "ERROR: pinned AppImageKit zsyncmake is required for the official update contract." >&2
    exit 1
  fi
  rm -f -- "${APPIMAGE_PATH}.zsync"
  "${ZSYNCMAKE_BIN}" -u "${DOWNLOAD_BASE_URL}/${APPIMAGE_NAME}" \
    -o "${APPIMAGE_PATH}.zsync" "${APPIMAGE_PATH}"
  if [[ ! -s "${APPIMAGE_PATH}.zsync" ]]; then
    echo "ERROR: official AppImage .zsync updater payload was not generated." >&2
    exit 1
  fi

  cat > "${APPIMAGE_PATH}.update.json" <<EOF
{
  "asset": "${APPIMAGE_NAME}",
  "owner": "${OWNER}",
  "repository": "${REPO}",
  "version": "${VERSION}",
  "channel": "${CHANNEL}",
  "tag": "${TAG}",
  "source_sha": "${SOURCE_SHA,,}",
  "download_base_url": "${DOWNLOAD_BASE_URL}",
  "appimage_update_information": "${UPDATE_INFO}"
}
EOF
else
  # A preview is a real Type 2 package but has no public updater track/info.
  ARCH=x86_64 "${APPIMAGETOOL}" "${APPDIR}" "${APPIMAGE_PATH}"
fi

python "${ROOT_DIR}/packaging/release/validate_appimage.py" \
  --header-only --appimage "${APPIMAGE_PATH}"

echo "Built genuine AppImage Type 2: ${APPIMAGE_PATH}"
