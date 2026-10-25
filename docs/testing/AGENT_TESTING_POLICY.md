# Política de pruebas y compilación de ChemUSON

Fecha: 2026-10-22; reconciliada con AGENTS.md el 2026-10-25. Esta política separa tres actividades distintas y **no elimina
ninguna prueba existente**: la colección completa de `tests/` se conserva como
herramienta para los agentes y para auditorías extraordinarias.

## A. Desarrollo y diagnóstico local (agentes)

- Las pruebas durante la implementación deben ser **proporcionales al cambio** y
  comprobar su funcionamiento real. Ejemplo: al reparar el worker RDKit
  congelado se probaron el mecanismo de transporte, sus errores y su
  interacción con la GUI (`tests/test_rdkit_live_stream.py`), no la suite
  completa de ~2.600 tests ajenos.
- Los scripts exploratorios temporales se eliminan al terminar o se conservan
  como evidencia **fuera del código de producción** (p. ej. `/tmp/…`).
- Las regresiones realmente útiles y reutilizables se quedan en `tests/`.
- No borrar ni debilitar la colección existente de pruebas.

Directorios de referencia para pruebas focalizadas habituales:

| Cambio | Pruebas sugeridas |
| --- | --- |
| Worker/empaquetado RDKit | `tests/test_rdkit_live_stream.py`, `tests/packaging/test_rdkit.py`, `tests/packaging/test_rdlib.py` |
| Visor 3D / geometría | `tests/viewport3d/`, `tests/geometry3d/` |
| GUI de compchem/apagado | `tests/test_compchem_dock*.py`, suites de apagado asíncrono |
| Arquitectura/importaciones | `tests/architecture/` |
| Empaquetado/workflows | `tests/test_appimage_validation.py`, `tests/test_preview_workflow_contract.py`, `tests/test_release_workflow_contract.py`, `tests/test_linux_distribution_scripts.py` |

## B. Verificación local previa al push

Antes de enviar cambios que afecten a producción, ejecutar localmente según el
alcance:

```bash
chem/bin/python -m pytest tests/<concreto> -q          # pruebas focalizadas
chem/bin/python -m pytest tests/architecture -q        # límites arquitectónicos
python -m compileall src tests tools packaging         # compilación de módulos
ruff check src tests tools packaging --select F401,F811,F821,E722,E741
openspec validate <cambio-activo> --strict             # si hay OpenSpec
git diff --check                                       # higiene de diffs
```

- Si el cambio afecta a la GUI, añadir un smoke funcional (arranque offscreen
  o, para paquetes, el smoke gráfico 2D→3D real descrito en
  `docs/release/PREVIEW_BUILDS.md`).
- **La ejecución masiva de pytest no es obligatoria en cada modificación.** La
  suite completa queda como auditoría extraordinaria cuando el alcance o el
  riesgo lo justifican (ver sección C y ejecución manual abajo).
- Registrar siempre: comandos exactos, resultados, versión del intérprete y de
  dependencias clave, y **SHA exacto del código comprobado**. No atribuir a
  GitHub pruebas realizadas solo en el equipo local.

## C. Compilación de paquetes en GitHub Actions

- `test.yml` ya **no se ejecuta automáticamente en cada push ni PR**: la suite
  completa (plan + 8 shards + summary + smokes Windows/Flatpak) se dispara solo
  manualmente mediante `workflow_dispatch` (pestaña Actions → *test* → *Run
  workflow*). Sigue siendo una auditoría completa, falla-cerrada y sin
  enmascarar: el mismo plan externo verifica que la colección coincide con
  `tests/ci/expected_pytest_nodeids.txt`.
- `Build Preview` (`build-preview.yml`) y `release.yml` construyen paquetes desde
  código previamente validado y **no ejecutan la suite completa**. Incluyen
  únicamente código de producción, dependencias y recursos necesarios; no
  empaquetan `tests/`, fixtures, cachés ni evidencias de OpenSpec.
- Las validaciones de paquete son ligeras, independientes del empaquetado y
  posteriores a la construcción del binario: identidad de código (versión
  canónica y SHA exacto), integridad (SHA-256, manifiestos de procedencia),
  dependencias nativas (iconos, worker RDKit, ChemName) y smoke funcional
  mínimo, incluido el **smoke esencial de generación real 2D→3D**
  (`geometry_flow_smoke` en `packaging/release/validate_appimage.py`), que es la
  comprobación que detectó el fallo del worker congelado.
- No existen dependencias `needs:` entre los jobs de `test.yml` y los de los
  workflows de construcción.
- No se eliminan módulos de producción por confundirlos con pruebas, ni se
  desactivan controles de seguridad/actualización para facilitar una build.

## Auditoría extraordinaria de la suite completa

- **En GitHub:** Actions → workflow *test* → `workflow_dispatch` sobre la
  referencia a auditar. Conserva shards, timeouts acotados y verificación de
  cobertura total de node IDs.
- **Local:** `python tools/ci_pytest_shards.py plan/run/verify` con los mismos
  argumentos que `test.yml` (ver `tests/test_test_workflow_contract.py` para el
  contrato exacto). Requiere el entorno de pruebas del repositorio.

## Directorios de salida de build

- `dist/` y `build/` siguen siendo rutas de construcción locales ignoradas por
  Git.
- `CHEMUSON_OUTPUT_DIR=/ruta/externa` reubica todas las salidas de las
  compilaciones locales (`build_linux_all.sh`, `build_appimage.sh`,
  `build_flatpak.sh` y su PyInstaller): cuelga `dist/`, `dist-appimage/`,
  `dist-flatpak/` y `dist/pyinstaller-work/` de ese directorio. Las variables
  específicas (`DIST_DIR`, `APPIMAGE_OUT_DIR`, `FLATPAK_OUT_DIR`) y los
  argumentos posicionales mantienen prioridad.
- Los paquetes generados nunca se trackean: `dist-appimage-*/`,
  `dist-preview*/`, `*.AppImage`, `*.flatpak`, `*.zsync` e instaladores están
  en `.gitignore`. Los manifiestos de procedencia y checksums ya confirmados
  como evidencia histórica se conservan trackeados deliberadamente
  (p. ej. `dist-appimage-diagnostic/{checksums.sha256,diagnostic-provenance.json}`)
  y no se reescriben.
- Limpiar una salida no afecta al código de producción: basta con borrar el
  directorio de salida correspondiente.

## Reconciliación normativa con `AGENTS.md` (2026-10-25)

La armonización de `AGENTS.md` §1.2 y §1.3 fue **aprobada por el propietario** y queda aplicada en esta revisión. Ambos documentos prescriben ahora el mismo procedimiento:

- Baseline y verificaciones proporcionales al alcance, con `git status --short`, SHA, `compileall`, arquitectura y Ruff focal cuando se modifique código. En cambios exclusivamente documentales, comprobar documentación, referencias y `git diff --check`.
- `pytest --collect-only -q` y el manifiesto de node IDs se verifican cuando cambian la colección, los tests o el plan de shards, no indiscriminadamente en cambios de documentación.
- Tests relevantes y controles acotados durante el desarrollo; la suite completa se reserva a auditoría extraordinaria justificada por el riesgo/alcance y **autorizada expresamente por el propietario**.
- No eliminar, desactivar ni debilitar pruebas. `test.yml` sigue disponible por `workflow_dispatch`. Builds y Release conservan sus gates fail-closed, checksums, proveniencia y smoke de 2D→3D.
- Una prueba omitida por entorno, un warning de baseline o un fallo histórico **no equivalen a PASS**. Registrar comandos, entorno, SHA y resultados reales.

Esta revisión sólo aclara reglas de trabajo; no cambia código, workflows, límites, manifiestos ni resultados históricos de pruebas.
