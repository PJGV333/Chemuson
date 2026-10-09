# Validación — estabilización de pytest y CI

## Estado de aceptación

**NOT READY para aceptación de campaña.** La validación local completa está verde, pero no se pudo publicar `fix/ci-pytest-stabilization` ni observar un run de Actions sobre este cambio. `gh auth status` indica que no hay sesión GitHub. El push normal falló antes de transmitir commits:

```text
fatal: could not read Username for 'https://github.com': terminal prompts disabled
```

`git ls-remote --heads origin fix/ci-pytest-stabilization` no devuelve una rama remota. No hubo merge, PR, tag, release ni publicación de artefactos. Además, el `Inspector x=-3` reportado por CI no se reproduce localmente; su causa remota sigue abierta y no se tocó producción/UI.

## Git y alcance

- Rama: `fix/ci-pytest-stabilization`.
- Base exacta: `fix/chemio-stereo-roundtrip` / `f5f82a63c84c56c9dd3ee07e7c7f56e7eb700bef`.
- Commits de campaña:
  - `2103d11` — OpenSpec independiente y baseline.
  - `7860ccc` — regresiones aisladas, settings NativeFormat, probes y geometría estable.
  - `5ef9e18` — manifiesto, supervisor/verificador de shards y workflow.
- Archivos modificados:
  - `.github/workflows/test.yml`
  - `openspec/changes/stabilize-pytest-ci-reproducibility/{.openspec.yaml,baseline.md,design.md,proposal.md,tasks.md,validation.md}` y `specs/reliable-pytest-ci/spec.md`
  - `tests/ci/expected_pytest_nodeids.txt`
  - `tests/test_clean2d_engine_candidates.py`
  - `tests/test_compchem3d_dock.py`
  - `tests/test_molecular_assistant_ui.py`
  - `tests/test_rdkit_packaged_worker.py`
  - `tests/test_side_panel.py`
  - `tests/test_test_workflow_contract.py`
  - `tools/ci_pytest_shards.py`
- Sin cambios en `src/`, catálogo de arquitectura, lógica química, release workflow ni comportamiento de producción.

## Nodos originales y diagnóstico

Se ejecutaron primero en procesos nuevos con Python 3.11, entorno Qt offscreen y timeout externo de 45 s:

| Nodo | Baseline aislado | Resultado/diagnóstico final |
|---|---|---|
| Clean2D `test_generate_candidates_attempts_rdkit_for_cyclic_graphs` | Fallaba: `rdkit_isolated` no está en la lista devuelta. | PASS. Spy confirma una llamada real al backend, guarda el candidato previo a deduplicación y comprueba que su hash `e9d64e43026270e943aee2b9a5f314c34ee84aa7` coincide con `simple_aromatic_template`; no se exigen dos candidatos para la misma geometría. |
| CompChem `test_compchem_controller_generates_async_with_fake_backend` | PASS aislado. | PASS con observación de `worker.finished → QThread.finished → job_finished`, sin ampliar el timeout de 3 s; también pasó en el shard completo. |
| Molecular Assistant `test_dialog_is_modeless_masks_key_and_rejects_missing_request_fields` | PASS aislado. Un probe con preferencia `resolution_method=ai` precargada antes de crear el diálogo reprodujo el fallo con el fixture antiguo. | PASS. El fixture redirige y limpia `QSettings.NativeFormat/UserScope`, restaura la ruta y verifica `ai_reference`, identidad habilitada y permiso externo deshabilitado. |
| Worker smoke: `test_smoke_allows_qt_runtime_hooks_without_a_gui_application` | PASS aislado; fallaba al colectarlo junto a tests GUI. | PASS en subprocess propio; mantiene intacta la detección. |
| Worker smoke: `test_smoke_detects_active_qt_application_and_chemuson_gui` | PASS aislado; fallaba al colectarlo junto a tests GUI. | PASS en subprocess propio; mantiene intacta la detección. |
| Panel `test_primary_tabs_have_complete_labels_padding_and_separation` | PASS aislado en Python 3.11. | PASS en shard completo; espera muestras geométricas estables y conserva las cotas originales. No se redujo tamaño, padding, gaps ni aserción de clipping. |

El grupo GUI + ambos probes worker reprodujo antes de los cambios **1 passed, 2 failed**: durante colección ya había 110 módulos `chemuson.gui`. Tras aislar los probes, el mismo grupo dio **3 passed**. El helper de geometría inicial, limitado por reloj a 1 s, falló bajo el backlog de eventos de un shard aunque su última medida era `Inspector x=2`; se cambió a un número finito de observaciones Qt estables con diagnóstico de las muestras recientes. La prueba conserva un límite externo del shard de 300 s.

## Colección completa en shards

- Colecta final Python 3.11: **2.130 tests collected**.
- Manifiesto versionado, ordenado y sin duplicados: **2.130 node IDs**; SHA-256 del conjunto `9eb563b221efcc98be35cad5de67d807a0051756e7a9a58d2d447fe83dea48f7`.
- Plan determinista: 8 shards, distribución `[267, 267, 266, 266, 266, 266, 266, 266]`.
- Verificación final: **2.130/2.130 node IDs**, cada uno asignado y reconocido en JUnit exactamente una vez en la matriz de shards; ningún reporte ausente o duplicado.
- Resultado total: **2.111 passed, 19 skipped preexistentes, 0 failed**. No se añadieron skips, exclusiones ni `continue-on-error`.

| Shard | Casos JUnit | Passed | Skipped | Tiempo externo |
|---:|---:|---:|---:|---:|
| 0 | 267 | 260 | 7 | 127.642 s |
| 1 | 267 | 265 | 2 | 56.551 s |
| 2 | 266 | 266 | 0 | 148.016 s |
| 3 | 266 | 261 | 5 | 60.442 s |
| 4 | 266 | 262 | 4 | 180.759 s |
| 5 | 266 | 265 | 1 | 99.384 s |
| 6 | 266 | 266 | 0 | 60.609 s |
| 7 | 266 | 266 | 0 | 29.347 s |

Todos los shards terminaron con código 0, sin timeout ni señal. El máximo fue **180.759 s**, menor que el límite externo de 300 s. Un probe intencional con límite de 1 s sobre un test que espera 1.05 s produjo `timed_out=true`, salida `SIGTERM` y el ID del último test; el proceso fue terminado. Una verificación con artefactos de shard ausentes falló explícitamente con `missing shard result reports`.

El primer pase del verificador detectó además que los IDs parametrizados con IPv6 contienen `::` dentro de corchetes. Se corrigió el parser de IDs de JUnit y se añadió un contrato de regresión; la campaña completa posterior verificó todos los IDs correctamente.

## Resto de verificaciones

- `pytest -q tests/architecture` (externo ≤5 min): **280 passed**.
- `python -m compileall -q src tests tools packaging`: PASS con Python 3.11 y Python del sistema.
- Ruff focal sobre todos los archivos afectados: PASS.
- Ruff global `src tests tools packaging --select F401,F811,F821,E722,E741`: sólo el F401 preexistente `import math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; se deja intacto por estar fuera de alcance.
- `pytest --collect-only -q` externo: **2.130 collected**.
- Contratos estáticos de workflow/manifiesto: PASS dentro de los shards.
- `openspec validate stabilize-pytest-ci-reproducibility --strict`: PASS.
- `git diff --check`: PASS.
- No se ejecutó el pytest monolítico.

## CI remoto pendiente

El run upstream #343 confirmó el SHA base y que falló el job `pytest`; `windows-smoke` y `flatpak-smoke` pasaron entonces. Los logs de Actions accesibles sin autenticación sólo exponen una anotación genérica de exit 1; su descarga respondió 403. Los smoke jobs se mantuvieron en el workflow nuevo, pero no hubo ejecución remota de aquel cambio porque el push no fue autenticado. Por ello los resultados locales no se presentan como evidencia de CI completo ni se declara la campaña aceptada.

## Microcorrección final — pestañas del panel lateral

### Evidencia y causa raíz

El run #344 (`37984526481`, SHA `c01a1f626b51da019d702d87517bc6355ffcd1c5`) confirma que el layout ya estaba estable: a 1440 × 900, viewport 308 px, tira 313 px, rango de scroll 0..5 y offset 5; `Inspector` quedó en x=-3 con ancho 55 px. Sólo falló `test_primary_tabs_have_complete_labels_padding_and_separation`; el plan, los otros siete shards y ambos smoke jobs pasaron. Los logs detallados muestran 264 passed, 1 skipped y 1 failed en shard 5. La suma verificada sobre los ocho artefactos JUnit es **2.109 passed, 20 skipped, 1 failed** (2.130 casos).

La fórmula sumaba 2 px redundantes a cada uno de los cinco botones, además de los 3 px laterales contractuales. Retirarlos reduce la tira exactamente 10 px: con las dimensiones de Actions, 313 → 303 px frente a viewport 308 px; desaparece el rango de desplazamiento que permitía a `ensureWidgetVisible()` empujar el primer botón fuera de vista. No fue necesario cambiar `set_active()` ni `ensureWidgetVisible()`.

### Corrección y prueba

- `side_panel.py`: ancho fijo igual a texto + `2 * sideTabPadX`; se conservan font-size 10 px, padding mínimo 3 px y gap 3 px.
- `test_side_panel.py`: en 1440 × 900 y 980 × 600, espera geometría estable y comprueba las cinco pestañas sin clipping ni rango/offset de scroll, al activar cada una.
- Medición local antes/después en viewport 308 px: tira 305 → 295 px; `scroll=(0,0,0)` en ambas ventanas y en los cinco estados activos. Los anchos locales pasan a 52/55/66/49/57 px, cada uno exactamente 6 px por encima de su etiqueta: 3 px laterales por lado.
- Baseline previo a editar: test afectado **1 passed**, archivo lateral **10 passed**, colecta completa **2.130 tests**, `compileall` y Ruff focal PASS; strict OpenSpec PASS. Ruff global conserva sólo `F401 math` fuera de alcance. La interpretación inicial de Python 3.11 no tenía pytest; el Python 3.14 del sistema sin `chem/lib` no tenía RDKit. Se usó el Python 3.14.7 del sistema con `chem/lib/python3.14/site-packages` del checkout (pytest 9.1.1, PyQt6 6.11.0, RDKit 2026.03.6); colecta completa: 2.130 tests. No se ejecutó la suite monolítica.
- Test afectado: **1 passed** (timeout externo 60 s); `tests/test_side_panel.py`: **10 passed** (120 s); grupo `test_side_panel.py`, `test_ui_theme_foundation.py`, `test_ui_polish.py`: **75 passed** (300 s); arquitectura: **280 passed** (300 s).
- `compileall`: PASS; Ruff focal: PASS; OpenSpec strict: PASS; `git diff --check`: PASS. Ruff global sigue mostrando únicamente el `F401 math` preexistente y fuera de alcance.

**CI posterior a esta corrección: pendiente de push y ejecución. Estado: NOT READY hasta verificar en Actions el plan, ocho shards, resumen, Windows smoke, Flatpak smoke y cobertura exacta de 2.130 IDs.**
