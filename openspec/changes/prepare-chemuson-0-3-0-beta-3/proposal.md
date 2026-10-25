# Propuesta: preparar ChemUSON 0.3.0-beta.3

## Motivación

El tag anotado `v0.3.0-beta.2` sí inició el workflow oficial, pero `release_gate` falló cerrado a los 10 segundos: `Release preflight failed closed: Command '['git', 'rev-parse', 'HEAD']' returned non-zero exit status 128.` (run `40160549763`). La causa raíz está confirmada en el código: `packaging/release/validate_release_tag.py` calculaba `REPO_ROOT = RELEASE_DIR.parents[2]` desde el **directorio** `packaging/release`, resolviendo un nivel **fuera** del checkout; en el runner el `git rev-parse` se ejecutaba fuera de cualquier repositorio Git. Los scripts compañeros calculan `parents[2]` desde el **archivo** (correcto), y ninguna prueba auxili tocaba la constante del validador, por lo que el Build Preview y los contracts previos nunca detectaron el defecto antes de gastar el tag.

La recuperación segura mantiene la numeración monotónica: preparar `0.3.0-beta.3` sobre la candidata beta.2 sin tocar los tags beta.1/beta.2, corrigiendo el defecto y — lo esencial — probando el preflight oficial de extremo a extremo **dentro de GitHub Actions** antes de crear el siguiente tag.

## Alcance

- Rama aislada `release/v0.3.0-beta.3-prep` creada desde el SHA exacto de la candidata beta.2 (`6aecba5a5585d9560b693186578b76bea9e623b8`), verificando ancestros y estado Git.
- Corrección localizada de `REPO_ROOT` en `validate_release_tag.py`; auditoría de errores equivalentes en scripts relacionados.
- Regresiones: `REPO_ROOT` contra `git rev-parse --show-toplevel`; ejecución de la lógica real de `verify_tag_checkout` y del script completo en un checkout de ensayo con tag anotado.
- Job de ensayo en `build-preview.yml` que reproduce el preflight oficial en un runner con checkout completo del SHA candidato: tag anotado temporal solo local, evento de tag simulado, token de lectura, exit code 0, salidas esperadas y borrado garantizado del tag temporal. Fail-closed: sin ensayo aprobado, el preview no es apto para publicación.
- Ensayo de las comprobaciones acotadas del job `release_gate` en condiciones equivalentes a su runner, con paridad byte a byte verificada por contracts.
- Actualización únicamente de versión canónica, AppStream, notas, `body_path` del workflow, contratos de pruebas afectados, OpenSpec y documentación de recuperación/política.

## Fuera de alcance

- Crear, borrar, mover o recrear `v0.3.0-beta.1` o `v0.3.0-beta.2`; crear `v0.3.0-beta.3`; crear GitHub Release o publicar el remoto Flatpak.
- Cambiar funcionalidad de ChemUSON: química, geometría, Clean2D, ChemName, RDKit, GUI, persistencia o dependencias.
- Ejecutar la suite masiva de pytest.
- Debilitar la validación SHA o los controles de inmutabilidad del preflight.

## Módulos probables afectados

`packaging/release/validate_release_tag.py`, `.github/workflows/build-preview.yml`, `.github/workflows/release.yml` (sólo `body_path`), `tests/test_release_tag_validation.py`, `tests/test_preview_workflow_contract.py`, `tests/test_release_workflow_contract.py`, `tests/test_branding.py`, `tests/test_iupac_ui.py`, `src/chemuson/_version.py`, metainfo AppStream, documentación `docs/release/*`.
