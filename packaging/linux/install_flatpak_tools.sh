#!/usr/bin/env bash
set -euo pipefail

# Preparación determinista y acotada del toolchain Flatpak del anfitrión en los
# runners alojados por GitHub (Ubuntu 24.04).
#
# Por qué existe este script y no un `apt-get update && apt-get install -y`
# directo: ese cierre de dependencias descarga ~99 MB en ~123 paquetes. Con un
# mirror sano tarda 0,3-0,7 min, pero cuando el nodo `azure.archive.ubuntu.com`
# asignado al runner está degradado apt permanece en ese origen y el mismo paso
# ha tardado 37-44 min, consumiendo casi todo el `timeout-minutes` del job y
# dejando el empaquetado real sin tiempo (véanse actions/runner-images #14823,
# #14594 y #12949, además del incidente Flatpak registrado en
# docs/release/PREVIEW_BUILDS.md).
#
# Reglas de seguridad e integridad conservadas:
#   * apt sigue siendo la única fuente de paquetes; sólo se admiten los mirrors
#     oficiales listados en el lockfile. No se añaden PPA ni repos externos.
#   * nunca se desactiva la verificación de firmas: el índice firmado sigue
#     aplicando y apt comprueba el hash de cada .deb antes de instalarlo.
#   * la caché de .deb se reutiliza a través de apt (`--no-download`), así que
#     un archivo obsoleto o alterado simplemente se vuelve a descargar.
#   * cada fase de red tiene un límite de tiempo explícito y el script falla de
#     forma cerrada si ninguna fuente oficial responde.
#
# Las rutas y la identidad de la distribución pueden reescribirse por variables
# de entorno para un arnés local de simulación; en Actions usan los valores
# reales del runner y no cambian el comportamiento de producción.

# Sin valor: se usa sudo (runners de GitHub). Vacio deliberado: sin envoltorio
# de privilegios, solo para el arnés de simulación local.
SUDO="${CHEMUSON_APT_SUDO-sudo}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
LOCK_FILE="${CHEMUSON_FLATPAK_TOOLCHAIN_LOCK:-${REPO_ROOT}/packaging/flatpak/toolchain.lock}"
SOURCES_FILE="${CHEMUSON_APT_SOURCES_FILE:-/etc/apt/sources.list.d/ubuntu.sources}"
MIRROR_LIST="${CHEMUSON_APT_MIRROR_LIST:-/etc/apt/apt-mirrors.txt}"
APT_ARCHIVES="${CHEMUSON_APT_ARCHIVES_DIR:-/var/cache/apt/archives}"

lock_value() {
  local key="$1"
  [[ -f "${LOCK_FILE}" ]] || return 0
  sed -n "s/^${key}=//p" "${LOCK_FILE}" | head -n 1
}

if [[ ! -f "${LOCK_FILE}" ]]; then
  echo "ERROR: Flatpak toolchain lock file not found: ${LOCK_FILE}" >&2
  exit 2
fi

IFS=' ' read -r -a PACKAGES <<< "${CHEMUSON_FLATPAK_APT_PACKAGES:-$(lock_value apt_packages)}"
if [[ "${#PACKAGES[@]}" -eq 0 ]]; then
  echo "ERROR: no apt packages declared for the Flatpak toolchain." >&2
  exit 2
fi
IFS=' ' read -r -a ALLOWED_HOSTS <<< "${CHEMUSON_FLATPAK_APT_ALLOWED_HOSTS:-$(lock_value allowed_apt_hosts)}"
if [[ "${#ALLOWED_HOSTS[@]}" -eq 0 ]]; then
  ALLOWED_HOSTS=(archive.ubuntu.com security.ubuntu.com azure.archive.ubuntu.com us.archive.ubuntu.com)
fi

CACHE_DIR="${CHEMUSON_FLATPAK_APT_CACHE_DIR:-${HOME}/$(lock_value apt_cache_dir)}"
if [[ -z "${CACHE_DIR}" || "${CACHE_DIR}" == "${HOME}/" ]]; then
  CACHE_DIR="${HOME}/.cache/chemuson-flatpak-apt"
fi
mkdir -p "${CACHE_DIR}"

PROBE_SECONDS="$(lock_value apt_mirror_probe_seconds)"
UPDATE_TIMEOUT="$(lock_value apt_update_timeout_seconds)"
DOWNLOAD_TIMEOUT="$(lock_value apt_download_timeout_seconds)"
INSTALL_TIMEOUT="$(lock_value apt_install_timeout_seconds)"
ACQUIRE_RETRIES="$(lock_value apt_acquire_retries)"
ACQUIRE_TIMEOUT="$(lock_value apt_acquire_timeout_seconds)"
for numeric in PROBE_SECONDS UPDATE_TIMEOUT DOWNLOAD_TIMEOUT INSTALL_TIMEOUT ACQUIRE_RETRIES ACQUIRE_TIMEOUT; do
  if [[ ! "${!numeric}" =~ ^[0-9]+$ ]]; then
    echo "ERROR: invalid numeric value for ${numeric}: '${!numeric}'" >&2
    exit 2
  fi
done

# Opciones de adquisición acotadas: un mirror que sólo estrangula el ancho de
# banda puede tardar decenas de minutos sin devolver un error. Con Timeout bajo,
# apt agota el intento en segundos y el orden del mirror list decide el fallback.
APT_OPTS=(
  "-o" "Acquire::Retries=${ACQUIRE_RETRIES}"
  "-o" "Acquire::http::Timeout=${ACQUIRE_TIMEOUT}"
  "-o" "Acquire::https::Timeout=${ACQUIRE_TIMEOUT}"
  "-o" "Acquire::ForceIPv4=true"
)

OS_ID="${CHEMUSON_APT_OS_ID:-$(. /etc/os-release && printf '%s' "${ID:-}")}"
SUITE="${CHEMUSON_APT_OS_CODENAME:-$(. /etc/os-release && printf '%s' "${VERSION_CODENAME:-}")}"
if [[ "${OS_ID}" != "ubuntu" || ! "${SUITE}" =~ ^[a-z]+$ ]]; then
  echo "ERROR: expected an Ubuntu runner with a codename, got OS_ID='${OS_ID}' SUITE='${SUITE}'." >&2
  exit 2
fi

now_ms() { date +%s%3N; }
now_s() { date +%s; }
elapsed() { printf '%s' "$(( $(now_s) - $1 ))"; }
log() { printf '%s\n' "$*"; }

probe_host() {
  local host="$1" scheme url
  for scheme in https http; do
    url="${scheme}://${host}/ubuntu/dists/${SUITE}/Release"
    if curl -fsS --max-time "${PROBE_SECONDS}" -o /dev/null "${url}"; then
      printf '%s' "${scheme}"
      return 0
    fi
  done
  return 1
}

# --- 1. Elegir un origen oficial sano en vez de aceptar el nodo asignado ------
declare -a PROBES=()
PROBE_T0="$(now_s)"
for host in "${ALLOWED_HOSTS[@]}"; do
  attempt_t0="$(now_ms)"
  if scheme="$(probe_host "${host}")"; then
    latency=$(( $(now_ms) - attempt_t0 ))
    PROBES+=("${latency} ${host} ${scheme}")
    log "apt origin probe: ${host} OK in ${latency} ms (${scheme})"
  else
    log "apt origin probe: ${host} unreachable within ${PROBE_SECONDS}s"
  fi
done
probe_seconds="$(elapsed "${PROBE_T0}")"

if [[ "${#PROBES[@]}" -eq 0 ]]; then
  echo "ERROR: no official Ubuntu archive answered within ${PROBE_SECONDS}s; refusing to install from an unverified origin." >&2
  exit 3
fi

mapfile -t ORDERED_PROBES < <(printf '%s\n' "${PROBES[@]}" | sort -n -k1,1)
probe_entry() { printf '%s' "${1#* }"; }

write_mirror_list_from() {
  # Emite los origenes sanos restantes, empezando en el indice indicado.
  local start="$1" index=0 entry host scheme line=""
  for (( index=start; index<${#ORDERED_PROBES[@]}; index++ )); do
    entry="$(probe_entry "${ORDERED_PROBES[index]}")"
    host="${entry%% *}"
    scheme="${entry##* }"
    line+="${scheme}://${host}/ubuntu/	priority:$(( index - start + 1 ))"$'\n'
  done
  # El fallo de escritura no aborta: apply_origin comprueba el resultado y
  # cambia al modo URIs: directas sobre un origen oficial sondeado.
  printf '%s' "${line}" | ${SUDO} tee "${MIRROR_LIST}" >/dev/null || true
}

rewrite_sources_uri() {
  # Solo se reescriben las lineas URIs: cuyo host ya es un mirror oficial
  # admitido; las firmas (Signed-By) y los demas repositorios no se tocan.
  local host="$1" scheme="$2" direct="${3:-0}" allowed
  allowed="$(printf '%s|' "${ALLOWED_HOSTS[@]}" | sed 's/|$//')"
  local updated
  if ! updated="$(${SUDO} awk -v host="${host}" -v scheme="${scheme}" -v direct="${direct}" -v allowed="${allowed}" '
    /^URIs:/ {
      # Las actualizaciones de seguridad quedan en el archivo dedicado de
      # Canonical salvo que ese sea justamente el origen seleccionado.
      if ($2 ~ /security\.ubuntu\.com/ && host != "security.ubuntu.com") { print; next }
      # Fallback: si el runner declara mirror+file: y ese archivo no se pudo
      # reescribir, se apunta el stanza principal a un origen oficial sondeado.
      if (direct && $2 ~ /^mirror\+file:/) { $2 = scheme "://" host "/ubuntu/"; print; next }
      if ($2 ~ "://(" allowed ")") { $2 = scheme "://" host "/ubuntu/"; print; next }
      print; next
    }
    { print }
  ' "${SOURCES_FILE}")"; then
    echo "ERROR: could not read ${SOURCES_FILE}; refusing to install from an unverified origin." >&2
    exit 3
  fi
  if ! printf '%s\n' "${updated}" | ${SUDO} tee "${SOURCES_FILE}" >/dev/null; then
    echo "ERROR: ${SOURCES_FILE} is not writable; refusing to install from an unverified origin." >&2
    exit 3
  fi
}

CURRENT_HOST=""
CURRENT_SCHEME=""
ORIGIN_MODE="unchanged"
if [[ -r "${SOURCES_FILE}" ]] && grep -q "mirror+file:" "${SOURCES_FILE}"; then
  # El runner usa mirror+file: reordenar prioridades es la via menos invasiva.
  ORIGIN_MODE="mirror-list"
  declared_mirror="$(grep -oE 'mirror\+file:[^[:space:]]+' "${SOURCES_FILE}" | head -n 1 | cut -d: -f2-)"
  if [[ -n "${declared_mirror}" && -z "${CHEMUSON_APT_MIRROR_LIST:-}" ]]; then
    MIRROR_LIST="${declared_mirror}"
  fi
elif [[ -r "${SOURCES_FILE}" ]]; then
  ORIGIN_MODE="sources-uri"
fi

apply_origin() {
  local index="$1" entry
  entry="$(probe_entry "${ORDERED_PROBES[index]}")"
  CURRENT_HOST="${entry%% *}"
  CURRENT_SCHEME="${entry##* }"
  case "${ORIGIN_MODE}" in
    mirror-list)
      write_mirror_list_from "${index}"
      # Confirmar que la reescritura tuvo efecto: sin esa comprobacion el job
      # seguiria dependiendo del nodo asignado al mirror list del runner.
      if ! head -n 1 "${MIRROR_LIST}" 2>/dev/null | grep -q "://${CURRENT_HOST}/"; then
        log "WARNING: ${MIRROR_LIST} is not writable; pinning the main stanza to the selected official origin instead."
        ORIGIN_MODE="sources-direct"
        rewrite_sources_uri "${CURRENT_HOST}" "${CURRENT_SCHEME}" 1
      fi
      ;;
    sources-uri) rewrite_sources_uri "${CURRENT_HOST}" "${CURRENT_SCHEME}" ;;
    sources-direct) rewrite_sources_uri "${CURRENT_HOST}" "${CURRENT_SCHEME}" 1 ;;
    *) ;;
  esac
  log "apt origin applied: ${CURRENT_SCHEME}://${CURRENT_HOST}/ubuntu/ through ${ORIGIN_MODE}"
}

apply_origin 0
log "apt origin selected: ${CURRENT_SCHEME}://${CURRENT_HOST}/ubuntu/ (fastest of ${#PROBES[@]} healthy official mirrors)"

# --- 2. Indexes: acotados y con un unico rebote de mirror --------------------
stage_cached_debs() {
  local staged=0 deb
  ${SUDO} install -d -m 0755 "${APT_ARCHIVES}" 2>/dev/null || true
  shopt -s nullglob
  for deb in "${CACHE_DIR}"/*.deb; do
    ${SUDO} cp -f "${deb}" "${APT_ARCHIVES}/" 2>/dev/null || true
    ${SUDO} chmod 0644 "${APT_ARCHIVES}/$(basename "${deb}")" 2>/dev/null || true
    staged=$(( staged + 1 ))
  done
  shopt -u nullglob
  printf '%s' "${staged}"
}

STAGED_DEBS="$(stage_cached_debs)"
log "cached .deb archives staged for apt verification: ${STAGED_DEBS}"

UPDATE_T0="$(now_s)"
apt_update_ok=false
if timeout "${UPDATE_TIMEOUT}" ${SUDO} apt-get "${APT_OPTS[@]}" update -qq; then
  apt_update_ok=true
else
  log "WARNING: apt-get update exceeded ${UPDATE_TIMEOUT}s or failed on ${CURRENT_HOST}; failing over to the next official mirror."
  if [[ "${#ORDERED_PROBES[@]}" -gt 1 && "${ORIGIN_MODE}" != "unchanged" ]]; then
    apply_origin 1
    if timeout "${UPDATE_TIMEOUT}" ${SUDO} apt-get "${APT_OPTS[@]}" update -qq; then
      apt_update_ok=true
    fi
  fi
fi
update_seconds="$(elapsed "${UPDATE_T0}")"
if [[ "${apt_update_ok}" != true ]]; then
  echo "ERROR: apt-get update could not complete from the official Ubuntu mirrors within ${UPDATE_TIMEOUT}s each (elapsed ${update_seconds}s)." >&2
  exit 3
fi

# --- 3. Descargar una sola vez; la instalacion posterior es offline ---------
attempt_download() {
  timeout "${DOWNLOAD_TIMEOUT}" ${SUDO} apt-get "${APT_OPTS[@]}" install --yes --download-only "${PACKAGES[@]}"
}

DOWNLOAD_T0="$(now_s)"
download_ok=false
if attempt_download; then
  download_ok=true
else
  log "WARNING: bounded package download failed after ${DOWNLOAD_TIMEOUT}s on ${CURRENT_HOST}; failing over to the next official mirror."
  if [[ "${#ORDERED_PROBES[@]}" -gt 1 && "${ORIGIN_MODE}" != "unchanged" ]]; then
    apply_origin 1
    if attempt_download; then
      download_ok=true
    fi
  fi
fi
download_seconds="$(elapsed "${DOWNLOAD_T0}")"

# Publicar los .deb realmente descargados para que `actions/cache` los reutilice
# en el siguiente run; apt los volvera a verificar contra el indice firmado.
PUBLISHED_DEBS=0
shopt -s nullglob
for deb in "${APT_ARCHIVES}"/*.deb "${APT_ARCHIVES}"/partial/*.deb; do
  if [[ ! -f "${deb}" ]]; then
    continue
  fi
  ${SUDO} cp -f "${deb}" "${CACHE_DIR}/" 2>/dev/null || true
  ${SUDO} chmod 0644 "${CACHE_DIR}/$(basename "${deb}")" 2>/dev/null || true
  PUBLISHED_DEBS=$(( PUBLISHED_DEBS + 1 ))
done
shopt -u nullglob
log "cached .deb archives published for reuse: ${PUBLISHED_DEBS}"

INSTALL_T0="$(now_s)"
install_ok=false
if timeout "${INSTALL_TIMEOUT}" ${SUDO} apt-get "${APT_OPTS[@]}" install --yes --no-download "${PACKAGES[@]}"; then
  install_ok=true
else
  log "offline install from verified archives was incomplete; retrying with a bounded online install."
  if timeout "${INSTALL_TIMEOUT}" ${SUDO} apt-get "${APT_OPTS[@]}" install --yes "${PACKAGES[@]}"; then
    install_ok=true
  fi
fi
install_seconds="$(elapsed "${INSTALL_T0}")"
if [[ "${install_ok}" != true ]]; then
  echo "ERROR: Flatpak toolchain packages were not installed: ${PACKAGES[*]} (download_ok=${download_ok})" >&2
  exit 4
fi

# --- 4. Verificar que cada herramienta existe y supera el minimo fijado -----
version_at_least() {
  local actual="$1" minimum="$2"
  [[ "$(printf '%s\n%s\n' "${minimum}" "${actual}" | sort -V | head -n 1)" == "${minimum}" ]]
}

tool_binary() {
  case "$1" in
    gnupg) printf "gpg" ;;
    *) printf '%s' "$1" ;;
  esac
}

declare -a TOOL_REPORT=()
VERIFY_T0="$(now_s)"
for pkg in "${PACKAGES[@]}"; do
  dpkg_version="$(dpkg-query -W -f='${Version}' "${pkg}" 2>/dev/null || true)"
  if [[ -z "${dpkg_version}" ]]; then
    echo "ERROR: ${pkg} is not registered as installed after preparation." >&2
    exit 4
  fi
  binary="$(tool_binary "${pkg}")"
  raw_version=""
  if command -v "${binary}" >/dev/null 2>&1; then
    raw_version="$("${binary}" --version 2>/dev/null | head -n 3 | tr '\n' ' ' || true)"
  fi
  numeric="$(printf '%s' "${raw_version}" | grep -oE '[0-9]+(\.[0-9]+)+' | head -n 1 || true)"
  minimum="$(lock_value "min_${pkg}")"
  if [[ -n "${minimum}" ]]; then
    if [[ -z "${numeric}" ]]; then
      echo "ERROR: could not read the version of ${pkg} from '${raw_version}'." >&2
      exit 4
    fi
    if ! version_at_least "${numeric}" "${minimum}"; then
      echo "ERROR: ${pkg} ${numeric} is older than the locked minimum ${minimum}." >&2
      exit 4
    fi
  fi
  TOOL_REPORT+=("${pkg}=${dpkg_version}${numeric:+ (cli ${numeric})}")
done
verify_seconds="$(elapsed "${VERIFY_T0}")"

TOTAL_SECONDS=$(( update_seconds + download_seconds + install_seconds + verify_seconds ))
log "Flatpak toolchain ready: ${TOOL_REPORT[*]}"
log "Flatpak toolchain preparation timings (s): probe=${probe_seconds} update=${update_seconds} download=${download_seconds} install=${install_seconds} verify=${verify_seconds} total=${TOTAL_SECONDS}"

if [[ -n "${GITHUB_STEP_SUMMARY:-}" ]]; then
  {
    printf '### Flatpak toolchain preparation\n\n'
    printf -- '- apt origin used: `%s://%s/ubuntu/` (probe %s s, mode %s)\n' "${CURRENT_SCHEME}" "${CURRENT_HOST}" "${probe_seconds}" "${ORIGIN_MODE}"
    printf -- '- cached .deb staged / published: `%s` / `%s`\n' "${STAGED_DEBS}" "${PUBLISHED_DEBS}"
    printf -- '- packages: %s\n' "${TOOL_REPORT[*]}"
    printf -- '| phase | seconds |\n| --- | --- |\n'
    printf '| apt update | %s |\n| package download | %s |\n| offline install | %s |\n| tool verification | %s |\n| total | %s |\n' \
      "${update_seconds}" "${download_seconds}" "${install_seconds}" "${verify_seconds}" "${TOTAL_SECONDS}"
  } >> "${GITHUB_STEP_SUMMARY}"
fi
