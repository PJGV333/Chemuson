# Validación de preparación beta.2

## Identidad de entrada y preservación de beta.1

- Rama base verificada: `release/v0.3.0-beta.1-prep` / `origin/release/v0.3.0-beta.1-prep` en `aa2e91c98cf1a16e8299e881955b1614edb2e18d`.
- Tag beta.1 anotado: commit objetivo `85fa4ed2bd0afe1d31953440228fb73a69b54e7e`; objetivo local y remoto coincidían al iniciar.
- Rama beta.2 y tag/release `v0.3.0-beta.2`: ausentes al iniciar.
- Los fixes listados en el diseño son ancestros de la rama base.

## Baseline (antes de cambios de producción)

- Python: `python -m compileall src tests tools packaging` — exit 0.
- Colección: `pytest --collect-only -q` — exit 0; 2630 tests collected.
- Suite completa: `pytest -q` — **no ejecutada**, por restricción explícita.
- Ruff: `ruff check src tests tools packaging --select F401,F811,F821,E722,E741` — exit 1 por la única excepción basal `F401` de `math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; no se modifica fuera de alcance.
- La salida exacta de baseline se conserva en `baseline.md` y `baseline-pytest-collect.txt`.

## Verificación focal y protecciones

- OpenSpec: `openspec validate prepare-chemuson-0-3-0-beta-2 --strict` — PASS.
- Primer foco de versión/contratos/UI: **66 passed in 2.88s**.
- Suite de arquitectura acotada: `timeout 300s python -m pytest -q tests/architecture` — **283 passed in 9.93s**.
- Gate focal de release (misma selección de tests del `release.yml`, sin pytest monolítico): **199 passed in 20.05s**.
- `python tools/chemname_acceptance.py --quiet --report /tmp/chemuson-beta2-chemname-source-acceptance.json` y `python packaging/release/validate_baseline_exceptions.py` — PASS; 4 casos ChemName exactos validados.
- `python -m compileall -q src tests tools packaging`, `git diff --check`, `appstreamcli validate packaging/flatpak/io.github.PJGV333.Chemuson.metainfo.xml` y Ruff scoped exacto del workflow de release — PASS. El parser de política acepta `v0.3.0-beta.2` como canal `beta`, prerelease verdadero.
- Ruff global final conserva exactamente el baseline: exit 1 por un único F401 `math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; salida idéntica byte a byte a `baseline-ruff.txt`. No se corrigió fuera de alcance.
- `pytest -q` completo — no ejecutado por instrucción explícita.
- El ruleset `Protect official version tags v*` (ID `25801369`) sigue `active`, aplica a `refs/tags/v*` y restringe creación/actualización/borrado; su único bypass actor es `PJGV333` (User ID `34524163`).
- `beta` y `flatpak-beta`: reviewer requerido `PJGV333` (ID `34524163`), `prevent_self_review=false`, `can_admins_bypass=false`, y deployment policy sólo permite tags `v*-beta.*`. No se inspeccionaron valores de secrets.
- `release.yml` sigue siendo exclusivamente `push.tags: v*`; conserva `contents: read` por defecto, gate de tag/SHA, checkout de SHA validado, builds de los cuatro paquetes, provenance/checksums, `overwrite_files: false`, y escritura sólo en jobs oficiales protegidos. No tiene `workflow_dispatch`. `build-preview.yml` conserva `contents: read` y el filtro `release/**-prep`; no publica.
- La API de GitHub responde 404 para una Release beta.1, coherente con que no hay publicación oficial. Al iniciar, no existían branch/tag/release beta.2. El tag beta.1 local/remoto aún apunta al mismo commit y no fue modificado.

## Resultado candidato

Pendiente al capturar esta versión del documento: empujar la rama de preparación y ejecutar un Build Preview nuevo desde el SHA candidato exacto; comprobar estado `complete=true`, `publication=false`, los cuatro paquetes, checksums y provenance descargados. Aun si todo pasa, esos artifacts son previews temporales, no paquetes oficiales. La creación del tag y la publicación quedan pendientes de la confirmación del propietario. No se ejecutó pytest completo.
