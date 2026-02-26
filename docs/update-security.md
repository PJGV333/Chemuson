# Seguridad del sistema de actualización

Fecha: 2026-02-26
Estado: Implementado (MVP endurecido)

## Checklist de seguridad

- [x] **Integridad mínima SHA-256 obligatoria** para artefactos descargados (`require_sha256=true` por defecto).
- [x] **Verificación criptográfica de firma** soportada (`hmac-sha256` / `ed25519`).
- [x] **Modo estricto de firma configurable** (`require_signature=true` fuerza firma válida).
- [x] **Protección anti-downgrade/replay** con `highest_seen_version`.
- [x] **Timeouts y retries** en cliente de GitHub API y descargas de sidecars/artefactos.
- [x] **Validación robusta de respuestas GitHub API**:
  - HTTPS obligatorio.
  - `owner/repo` saneados.
  - límite de tamaño de payload.
  - JSON válido y estructura esperada.
- [x] **Fallback seguro si GitHub no está disponible**:
  - uso de caché local con TTL.
  - rechazo de caché vencida.
- [x] **Bloqueo de ejecución de artefactos no verificados** antes de lanzar instalador.
- [x] **Rollback en fallo de aplicación** (incluye validación posterior opcional y restauración).
- [x] **Telemetría local mínima sin datos sensibles** (`~/.chemuson/update_logs/events.jsonl`).

## Amenazas y mitigaciones

### 1. Manipulación de artefactos en tránsito

Riesgo:
- atacante altera binario descargado.

Mitigaciones:
- comparación SHA-256 obligatoria por defecto.
- verificación de firma cuando hay clave/verificador configurado.
- opción `require_signature` para modo estricto.

### 2. Downgrade / replay attack

Riesgo:
- feed de releases devuelve una versión remota más antigua para forzar downgrade.

Mitigaciones:
- guardado de `highest_seen_version`.
- bloqueo de respuesta cuando versión observada es menor a la máxima histórica.

### 3. API GitHub malformada / indisponible

Riesgo:
- payload inválido, respuestas inesperadas o caída temporal.

Mitigaciones:
- retries con backoff y timeout controlado.
- límites de tamaño de payload.
- validaciones de esquema básico.
- fallback a caché local vigente.
- rechazo explícito de caché vencida.

### 4. Ejecución de instalador no verificado

Riesgo:
- artefacto cambia en disco después de la descarga inicial.

Mitigaciones:
- re-verificación justo antes de `launch_inno_installer`.
- si falla verificación, se bloquea la ejecución.

### 5. Fuga de datos en telemetría

Riesgo:
- logs locales con rutas privadas o URLs sensibles.

Mitigaciones:
- lista blanca de campos permitidos.
- sanitización que elimina rutas/URLs.
- eventos locales únicamente (sin envío remoto).

## Configuración relevante

Campos de seguridad en `UpdateSettings`:

- `require_sha256` (default: `true`)
- `require_signature` (default: `false`)
- `highest_seen_version` (persistente, anti-downgrade)

## Pruebas automáticas de manipulación

Cobertura principal:

- `tests/test_update_core.py`
  - bloqueo sin checksum cuando SHA-256 es requerido.
  - bloqueo sin firma cuando `require_signature=true`.
  - rollback ante fallo de validación post-update.
  - bloqueo de downgrade/replay.

- `tests/test_update_provider.py`
  - fallback a caché cuando GitHub falla.
  - rechazo de caché vencida.
  - rechazo de `api_base` inseguro.
  - descarte de assets no HTTPS.

- `tests/test_update_security.py`
  - rechazo de descargas por HTTP.
  - retries con recuperación ante fallo transitorio.

- `tests/test_update_telemetry.py`
  - sanitización de campos sensibles en telemetría local.

## Riesgos residuales (MVP)

- Si `require_signature=false`, la seguridad depende principalmente de SHA-256 y origen de checksums.
- El esquema de confianza criptográfica final depende de distribuir y rotar correctamente la clave pública/clave HMAC fuera del código duro.
- No hay pinning de certificado TLS (se usa validación TLS del sistema).
