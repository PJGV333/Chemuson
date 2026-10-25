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

## Resultado candidato

Pendiente al capturar esta versión del documento: empujar `release/v0.3.0-beta.3-prep`, ejecutar el Build Preview oficial y documentar aquí el run (ensayos verdes con exit code 0 y salidas exactas, tag temporal borrado, cuatro paquetes, provenance/checksums, SHA único). La creación del tag `v0.3.0-beta.3` queda pendiente de la confirmación explícita del propietario. No se ejecutó pytest completo.
