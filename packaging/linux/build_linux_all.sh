#!/usr/bin/env bash
set -euo pipefail

# Build helper: genera AppImage + Flatpak + checksums en un solo flujo.
# Uso:
#   bash packaging/linux/build_linux_all.sh <VERSION> [CHANNEL]
# Ejemplo:
#   bash packaging/linux/build_linux_all.sh 0.2.3-beta.2 stable

VERSION="${1:?missing VERSION (e.g. 0.2.3-beta.2)}"
CHANNEL="${2:-stable}"

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
VENV_DIR="${VENV_DIR:-${ROOT_DIR}/.venv}"
DIST_DIR="${DIST_DIR:-${ROOT_DIR}/dist}"
APPIMAGE_OUT_DIR="${APPIMAGE_OUT_DIR:-${ROOT_DIR}/dist-appimage}"
FLATPAK_OUT_DIR="${FLATPAK_OUT_DIR:-${ROOT_DIR}/dist-flatpak}"
OWNER="${OWNER:-PJGV333}"
REPO="${REPO:-Chemuson}"
TAG="${TAG:-v${VERSION}}"
FLATPAK_MANIFEST="${FLATPAK_MANIFEST:-${ROOT_DIR}/packaging/flatpak/io.github.PJGV333.Chemuson.yml}"

APPIMAGE_SCRIPT="${ROOT_DIR}/packaging/linux/build_appimage.sh"
FLATPAK_SCRIPT="${ROOT_DIR}/packaging/linux/build_flatpak.sh"
CHECKSUMS_SCRIPT="${ROOT_DIR}/packaging/release/generate_checksums.py"
SIGN_SCRIPT="${ROOT_DIR}/packaging/release/sign_hmac.py"

if ! command -v python3 >/dev/null 2>&1; then
  echo "python3 is required" >&2
  exit 1
fi
if ! command -v flatpak >/dev/null 2>&1; then
  echo "flatpak is required (install: flatpak flatpak-builder)" >&2
  exit 1
fi
if ! command -v flatpak-builder >/dev/null 2>&1; then
  echo "flatpak-builder is required (install: flatpak-builder)" >&2
  exit 1
fi
if [[ ! -f "${FLATPAK_MANIFEST}" ]]; then
  echo "Flatpak manifest not found: ${FLATPAK_MANIFEST}" >&2
  echo "Set FLATPAK_MANIFEST or pass a valid path to build_flatpak.sh." >&2
  exit 1
fi

cd "${ROOT_DIR}"

if [[ ! -d "${VENV_DIR}" ]]; then
  python3 -m venv "${VENV_DIR}"
fi
# shellcheck disable=SC1090
source "${VENV_DIR}/bin/activate"

python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m pip install pyinstaller

echo "[1/5] Building portable Linux binary with PyInstaller"
pyinstaller --clean --noconfirm chemuson.spec

echo "[2/5] Building AppImage + update metadata"
bash "${APPIMAGE_SCRIPT}" \
  "${VERSION}" \
  "${DIST_DIR}" \
  "${APPIMAGE_OUT_DIR}" \
  "${OWNER}" \
  "${REPO}" \
  "${CHANNEL}" \
  "${TAG}"

echo "[3/5] Building Flatpak bundle"
# build_flatpak.sh already classifies common failures:
# - manifest path
# - network/DNS/sandbox
# - python dependencies in flatpak-builder
bash "${FLATPAK_SCRIPT}" \
  "${VERSION}" \
  "${CHANNEL}" \
  "${FLATPAK_OUT_DIR}" \
  "${FLATPAK_MANIFEST}"

echo "[4/5] Generating checksums"
python "${CHECKSUMS_SCRIPT}" "${APPIMAGE_OUT_DIR}" "${FLATPAK_OUT_DIR}"

echo "[5/5] Optional signing (HMAC)"
if [[ -n "${CHEMUSON_SIGN_KEY:-}" ]]; then
  python "${SIGN_SCRIPT}" "${APPIMAGE_OUT_DIR}" "${FLATPAK_OUT_DIR}"
  echo "HMAC signatures generated (*.sig)."
else
  echo "CHEMUSON_SIGN_KEY is not set; skipping HMAC signing."
fi

echo "Done. Artifacts:"
echo "- AppImage: ${APPIMAGE_OUT_DIR}"
echo "- Flatpak:  ${FLATPAK_OUT_DIR}"
