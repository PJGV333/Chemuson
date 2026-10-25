# Diseño: beta.3 con preflight probado en Actions antes del tag

## Base y aislamiento

- Rama base: `release/v0.3.0-beta.3-prep` creada desde `6aecba5a5585d9560b693186578b76bea9e623b8`, que es simultáneamente HEAD de `origin/release/v0.3.0-beta.2-prep`, el commit objetivo del tag anotado `v0.3.0-beta.2` y HEAD local verificado.
- `v0.3.0-beta.1` (anotado → `85fa4ed2bd0afe1d31953440228fb73a69b54e7e`) y `v0.3.0-beta.2` (anotado → `6aecba5a…`) quedan intactos; `beta.1` es ancestro de `beta.2` y todo el trabajo previo permanece en la historia.
- Reglas de publicación inmutables: no se crea el tag beta.3 hasta que el propietario revise el SHA candidato y la evidencia del ensayo.

## Causa raíz y corrección

`RELEASE_DIR = Path(__file__).resolve().parent` es un directorio; `RELEASE_DIR.parents[2]` sube: `packaging` → raíz → **padre de la raíz**. En el runner `/home/runner/work/Chemuson/Chemuson` esto resuelve `/home/runner/work/Chemuson`, sin `.git`, y `_git_output` lanza `CalledProcessError` (128). Corrección mínima: `REPO_ROOT = RELEASE_DIR.parents[1]`. Los demás scripts de `packaging/release` usan `Path(__file__).resolve().parents[2]` desde el archivo, que ya es correcto; la auditoría (`grep 'parents\['`) no encontró otros defectos equivalentes.

## Regresiones

1. `test_module_repo_root_matches_git_toplevel`: la constante del módulo coincide con `git rev-parse --show-toplevel` ejecutado desde `packaging/release`.
2. `test_verify_tag_checkout_runs_real_git_logic_on_rehearsal_checkout`: repositorio Git temporal con tag anotado; ejecuta `verify_tag_checkout` real (positivo y negativos: SHA discordante, `ref_type` de rama, tag inexistente).
3. `test_full_preflight_main_succeeds_on_simulated_tag_run`: `main()` real en el checkout de ensayo con la consulta a GitHub sustituida; valida las salidas `tag/version/channel/sha`.
4. Dos pruebas de subprocess sin red sobre el script completo: todos los gates locales pasan y el fallo queda exactamente en la exigencia del token; y un SHA discordante falla cerrado antes de cualquier red.
5. Comprobación de mutación realizada durante la preparación: reintroducir `parents[2]` en una copia reprodujo el error exacto de beta.2 (`exit status 128`), confirmando que las regresiones detectan el defecto.

## Ensayo del preflight en el Build Preview

Job `preflight_rehearsal` en `build-preview.yml`:

- `ubuntu-latest`, `actions/checkout@v4` con `fetch-depth: 0` sobre `needs.prepare.outputs.sha` (exactamente las condiciones del `release_gate` oficial), Python 3.11.
- Guardas `git ls-remote` que abortan si el tag ya existe en el remoto (el ensayo nunca puede correr contra un tag publicado).
- Tag anotado temporal `v<versión>` creado **solo dentro del runner** apuntando al SHA candidato.
- Simulación completa del contexto de evento de tag mediante `export` en el paso: `GITHUB_REF_NAME`, `GITHUB_REF_TYPE=tag`, `GITHUB_SHA` exacto, `GITHUB_EVENT_PATH` (JSON sin `deleted`), `GITHUB_OUTPUT` volcado a archivo de la `RUNNER_TEMP`.
- Ejecución del código oficial `python packaging/release/validate_release_tag.py` con `RELEASE_PREFLIGHT_TOKEN: github.token` (contents: read), tras comprobar explícitamente `REPO_ROOT == git rev-parse --show-toplevel` del runner.
- Verificación de exit code 0 implícita (`set -euo pipefail`) y de las salidas esperadas con `grep -qxF` (`tag=`, `version=`, `channel=` derivado de la versión, `sha=` exacto).
- Limpieza `if: always()`: `git tag -d` y comprobación de ausencia local y remota. Nunca se ejecuta ningún push.
- El job entra al `needs` y a los estados del `preview_report`: cualquier resultado distinto de `success` marca el preview completo como incompleto (fail-closed).

Job `release_gate_rehearsal`: replica los pasos acotados del `release_gate` posterior al preflight (instalación Qt/requirements/-e ., catálogo de excepciones, ChemName acceptance, compileall, Ruff acotado, arquitectura 5m, contracts acotados 8m) en el mismo runner. El contract estático compara los comandos byte a byte (normalizados) contra `release.yml` para impedir deriva. Su resultado también gatea el `preview_report`.

Los contracts estáticos del preview se ajustan: se permite el uso controlado de `git tag -a`/`git tag -d` y `refs/tags/` (lectura) dentro del ensayo, manteniendo prohibidos `git push`, `gh release`, `softprops/action-gh-release`, `gh-pages`, la invocación de scripts de publicación y los canales públicos; el token sigue sin inyectar secrets y el literal `GITHUB_TOKEN` no aparece.

## Versionado y documentación

`set_version.py` sincroniza `_version.py` y AppStream a `0.3.0-beta.3` (fecha 2026-10-25). Nuevas notas `0.3.0-beta.3-public.md` (contenido equivalente a beta.2 más la nota de fiabilidad de publicación) y documento interno `0.3.0-beta.3.md` con la evidencia del fallo beta.2. `release.yml` sólo cambia su `body_path`. `VERSIONING_POLICY.md` y `PREVIEW_BUILDS.md` documentan la recuperación y los ensayos. Los contratos literales de versión/branding (`test_branding`, título de ventana en `test_iupac_ui`, `body_path` en `test_release_workflow_contract`) se actualizan a beta.3.

## Riesgos de regresión

- El ensayo podría introducir semántica de publicación en el preview: mitigado por contratos byte-a-byte, guardas `ls-remote`, limpieza incondicional y prohibición global de `git push`/actions de release.
- Paridad gate↔ensayo: mitigada por comparación automatizada de comandos en `test_release_gate_rehearsal_mirrors_official_gate_commands_exactly`.
- Ningún cambio en producción, química, GUI, Clean2D, ChemName ni RDKit; la suite masiva no se ejecuta por instrucción explícita.
