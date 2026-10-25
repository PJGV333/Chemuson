# Validación de preparación beta.3

## Identidad de entrada y preservación de tags

- Rama creada desde `6aecba5a5585d9560b693186578b76bea9e623b8` = HEAD de `origin/release/v0.3.0-beta.2-prep` = commit objetivo del tag anotado `v0.3.0-beta.2`. HEAD local y remoto coincidían al crear la rama; `git status --short` vacío.
- Tag anotado `v0.3.0-beta.1` → `85fa4ed2bd0afe1d31953440228fb73a69b54e7e`; tag anotado `v0.3.0-beta.2` → `6aecba5a…`. `beta.1` es ancestro de `beta.2`. Ninguno fue movido, borrado ni reescrito.
- Evidencia de la causa raíz en GitHub: run `40160549763` (workflow `release`, push de `v0.3.0-beta.2`, `release_gate`): `Release preflight failed closed: Command '['git', 'rev-parse', 'HEAD']' returned non-zero exit status 128.` a los ~7 s del paso `Validate immutable tag and exact source SHA`. No existe GitHub Release para beta.1 ni beta.2.

## Baseline (antes de cambios de producción)

- `git status --short`: vacío. `python -m compileall src tests tools packaging`: exit 0. `pytest --collect-only -q`: 2639 tests collected (incluye las 9 regresiones nuevas; ver nota en `baseline.md`). Ruff global: exit 1 con la única excepción basal F401 preexistente (`tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`), idéntica a beta.2. `pytest -q` completo: no ejecutado por instrucción explícita.

## Corrección y regresiones

- `REPO_ROOT = RELEASE_DIR.parents[1]`; `grep "parents\["` en `packaging`/`tools`/`src`: ningún otro defecto equivalente (los compañeros indexan desde el archivo, correcto).
- `tests/test_release_tag_validation.py`: 23 pruebas pasan, incluyendo constante vs `git rev-parse --show-toplevel`, `verify_tag_checkout` real sobre checkout con tag anotado, `main()` simulado offline con salidas exactas y dos subprocess fail-closed (toque de token exigido tras pasar todos los gates locales; SHA discordante falla antes de la red).
- Comprobación de mutación: una copia con `parents[2]` reprodujo el fallo exacto de beta.2 (`exit status 128`); la versión corregida supera todos los gates locales y sólo se detiene en la consulta HTTP con token inválido (401 fail-closed), comportamiento oficial esperado.
- Dry-run del validador corregido contra el repositorio real con contexto simulado: coincide el SHA del checkout y del tag; el mismatch de SHA sigue fallando cerrado (inmutabilidad intacta).

## Verificación focal local (sin suite masiva)

- `pytest -q` del conjunto focal de release/preview/version/branding/artifacts/contracts: **81 passed**.
- `tests/architecture` (acotado, igual que el gate): **283 passed in 11.31s**.
- `tests/test_iupac_ui.py` (título de ventana beta.3, PyQt6 offscreen): **5 passed**.
- `python packaging/release/validate_baseline_exceptions.py`: PASS (4 registros exactos).
- `python tools/chemname_acceptance.py --quiet`: PASS.
- `python -m compileall -q src tests tools packaging`: PASS. `git diff --check`: limpio.
- Ruff acotado exacto del workflow de release: **All checks passed!**
- Ruff global final: idéntico al baseline (única F401 basal preexistente, sin cambios).
- `appstreamcli validate` del metainfo: validación correcta (nota informativa `redundante` preexistente por licencia por defecto).
- `openspec validate prepare-chemuson-0-3-0-beta-3 --strict`: valid.
- YAML de `build-preview.yml` y `release.yml` parsea correctamente; contracts estáticos del preview exigen paridad byte a byte de los pasos acotados con `release_gate`.

## Protecciones verificadas

- El ensayo del preflight no publica: sin `git push`, sin actions de release, tag temporal solo local con limpieza `if: always()` y doble guarda `ls-remote`; token `github.token` (`contents: read`); sin secrets en el preview.
- `release.yml` sigue siendo exclusivamente `push.tags: v*`; su preflight conserva parseo estricto de tag, triple identidad SHA (evento/tag/checkout), metadata canónica+AppStream y consulta fail-closed de Release existente. La corrección no debilita ninguna de estas validaciones.
- `preview_report` exige `success` en `release_preflight_rehearsal` y `release_gate_bounded_rehearsal` además de los tres builds: sin ensayos verdes, el Build Preview no es apto para publicación.

## Resultado candidato (Build Preview 40162819727, SHA b16ccbcbc31781dc7f90e6d991844cf556023750)

Run oficial `Build Preview` push-disparado por la rama: https://github.com/PJGV333/Chemuson/actions/runs/40162819727 — conclusión `success`; los 7 jobs en verde (`Validate preview source`, ambos ensayos, los 3 builds y `Preview build report`).

### Ensayo real del preflight en el runner (evidencia literal del log del job `Release preflight rehearsal (non-publishing)`)

- `Path resolution check: REPO_ROOT == git toplevel == /home/runner/work/Chemuson/Chemuson` — la constante corregida coincide con la raíz Git del runner; el mismo script que falló en beta.2 (exit 128 fuera del repo) ahora resuelve correctamente.
- El validador oficial ejecutado bajo evento de tag simulado terminó con exit code 0 y printing final: `Validated v0.3.0-beta.3 (beta) at b16ccbcbc31781dc7f90e6d991844cf556023750; no existing release.` (la consulta de inexistencia de Release usó `github.token` de lectura y devolvió 404).
- Salidas exactas del validador (`GITHUB_OUTPUT`): `tag=v0.3.0-beta.3`, `version=0.3.0-beta.3`, `channel=beta`, `sha=b16ccbcbc31781dc7f90e6d991844cf556023750`; el paso `Verify preflight outputs match the candidate identity` las cotejó con `grep -qxF`.
- Ciclo del tag anotado temporal solo en el runner: guardas `git ls-remote` (inexistencia en remoto antes y después), `git tag -a` local sobre el SHA candidato, y limpieza incondicional: `Deleted tag 'v0.3.0-beta.3' (was 5329f86d)` (objeto anotado, distinto del commit). Nunca se ejecutó ningún push; el tag oficial `v0.3.0-beta.3` sigue sin existir.

### Ensayo de los checks acotados del gate (`Release gate bounded checks rehearsal`)

Comandos byte a byte idénticos a `release.yml` (contract verificado): instalación Qt/requirements/-e ., catálogo de excecciones, ChemName acceptance, `compileall`, Ruff acotado → `All checks passed!`, arquitectura acotada → `283 passed in 22.63s`, contracts acotados → `208 passed in 28.91s` (incluye las 9 regresiones nuevas del preflight).

### Paquetes y provenance (artifacts descargados y verificados localmente)

- Grupos: `chemuson-preview-windows-portable`, `chemuson-preview-windows-installer`, `chemuson-preview-linux-appimage`, `chemuson-preview-linux-flatpak` + `chemuson-preview-build-report`.
- `sha256sum -c` de los 4 paquetes y sus provenance: todas las sumas coinciden.
- Los cuatro `preview-provenance.json` reportan `0.3.0-beta.3`, `git_sha=b16ccbcbc31781dc7f90e6d991844cf556023750` (único SHA fuente), `build_status: success`, `publication: false`.
- `preview-build-report.json`: `complete: true`; estados `success` en `source_validation`, `release_preflight_rehearsal`, `release_gate_bounded_rehearsal`, `windows_portable_and_installer`, `linux_appimage_portable`, `linux_flatpak`.

Sin errores de versionado ni rutas en ningún job. Sin suite masiva de pytest. Tags `v0.3.0-beta.1`/`v0.3.0-beta.2` y sus commits permanecen intactos; no se creó `v0.3.0-beta.3` ni GitHub Release. La creación del tag queda pendiente de la confirmación explícita del propietario.

### Actualización final de evidencia (SHA definitivo)

Este documento se actualizó y se empujó como commit final de la rama; el push re-disparó el Build Preview oficial sobre ese SHA definitivo, repitiendo el ensayo completo del preflight (run registrado en `tasks.md`/informed final). El SHA definitivo para el tag es el HEAD verde resultante de esa ejecución.
