# Validación — pulido final del panel lateral

## Baseline

- Base: `release/v0.3.0-beta.1-prep` / `8f20e675d8be8632c5d750e63699830d6d1223bb`; árbol inicial limpio.
- La suite completa baseline, ejecutada una vez bajo `timeout 600s`, concluyó `1 failed, 2126 passed, 20 skipped in 595.74s`. Falló `tests/test_compchem3d_dock.py::test_compchem_controller_generates_async_with_fake_backend` por el orden de señales `worker.finished`/`thread.finished`/`job_finished`. Es un fallo preexistente, no relacionado y fuera de alcance; no se corrigió ni se repitió la suite monolítica.
- Ruff global baseline conservó únicamente el F401 conocido `math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; fuera de alcance.

## Verificación focalizada

Todos los comandos pytest tuvieron timeout externo real menor a 10 minutos:

| Comando | Resultado | Tiempo |
|---|---:|---:|
| `pytest -q tests/test_side_panel.py tests/test_spectroscopy_dock.py` | 19 passed | 5.93 s (tras el ajuste de anchuras) |
| `pytest -q tests/test_ui_polish.py` | 38 passed | 9.77 s |
| `pytest -q tests/test_ui_tool_rail.py` | 64 passed | 35.81 s |
| `pytest -q tests/test_ui_app_bar_tabs.py` | 24 passed | 4.30 s |
| `pytest -q tests/test_platform_settings.py` | 12 passed | 0.19 s |
| `pytest -q tests/test_ui_theme_foundation.py` | 28 passed | 0.94 s |
| `pytest -q tests/test_main_window_tabs.py` | 16 passed | 6.08 s |
| `pytest -q tests/architecture` | 280 passed | 9.83 s |
| `pytest -q tests/test_test_workflow_contract.py` | 4 passed | 0.19 s |
| `pytest --collect-only -q` | 2152 collected | 0.85 s total |

Se revisó el diff del manifiesto de IDs: cinco casos nuevos de interfaz, una prueba renombrada por su nueva expectativa y ningún ID eliminado fuera de ese rename. `tests/ci/expected_pytest_nodeids.txt` se actualizó con `tools/ci_pytest_shards.py update-expected --confirm`; el contrato CI pasó.

- `python -m compileall -q src tests tools packaging`: PASS (0.04 s).
- Ruff focalizado (`F401,F811,F821,E722,E741`) sobre Python cambiado: PASS (0.01 s). No se tocó el F401 baseline global.
- `openspec validate fix-beta-properties-panel-polish --strict`: PASS.
- `openspec validate prepare-chemuson-0-3-0-beta-1 --strict`: PASS.
- `git diff --check`: PASS.
- Se validaron ambos temas y selección/alineación/legibilidad del texto renderizado; los tests prueban que la copia TSV conserva exactamente los textos fuente sin HTML ni sustituciones Unicode.

## Alcance y aceptación manual

- La tabla cambia sólo el render de fórmula, etiquetas isotópicas NMR y aductos MS; magnitudes y TSV no cambian.
- El divisor QSplitter se arrastró en ambas direcciones en el test, se ejercitaron 980×600 y 1600×900, ancho mínimo, rail intacto, espacio del lienzo al ocultar y restauración de ancho/página.
- El inicio oculto se probó con `visible=true` persistido, manteniendo página activa y sin alterar esa preferencia; las siete acciones Ver se probaron.
- El onboarding del primer arranque se probó con panel oculto y el tercer spotlight en el botón visible del menú AppBar, sin apertura automática.
- Cada pestaña conserva texto, tooltip y nombre accesible completos incluso con viewport angosto.
- **Caveat Qt explícito:** un probe fuente/offscreen con timeout externo de 30 s reajustó el splitter en 340/320/300 px y `window.close()` devolvió `async_shutdown_complete=True`, pero el proceso terminó con SIGSEGV 139 durante el teardown del intérprete Qt. Un segundo intento acotado añadió limpieza explícita del UndoStack y `processEvents()` después del cierre; terminó igual. No quedó core versionado y no se continuó instrumentando. La prueba focalizada pytest de splitter/hide/show y los demás shards Qt indicados arriba pasaron. No se determinó la causa y no se presenta como crash del paquete ni como corregido; coincide con la limitación de teardown que el usuario pidió mantener explícita.
- La aceptación manual de Windows/Linux/AppImage sigue `NOT TESTED`; la automatización local no sustituye la instalación ni la prueba personal solicitada por el propietario. No se creó tag/Release ni se integró a la rama prep.

## CI y entrega remota

- Primer CI tras push, [run 38647067861](https://github.com/PJGV333/Chemuson/actions/runs/38647067861), SHA `4f8c4326de02d5776ed06fd96476c357284f9f04`: falló sólo la nueva aserción de legibilidad en el shard 7; en Ubuntu runner el valor MS ocupaba 176 px y la columna Stretch/Stretch ofrecía 167 px. Se corrigió el reparto de columnas en `ChemicalPropertiesDock` (Propiedad `ResizeToContents`, Valor `Stretch`) y se repitió localmente el shard focalizado: 19 passed.
- CI tras corregir las anchuras, [run 38647689407](https://github.com/PJGV333/Chemuson/actions/runs/38647689407), SHA `6287bf519f4ac0c581c5807045b2913bee7e0419`: **success**; 8/8 shards, 2152/2152 node IDs verificados (269 por shard), pytest-plan, Windows smoke e inspección Flatpak pasaron.
- CI del HEAD final, [run 38648382216](https://github.com/PJGV333/Chemuson/actions/runs/38648382216), SHA `ee75534ba9a02fffeeae001912501db0da4ddfc4`: **success**; 8/8 shards, los mismos 2152 node IDs, pytest-plan, Windows smoke e inspección Flatpak pasaron.
- Rama publicada por push normal: `fix/beta-properties-panel-polish` — https://github.com/PJGV333/Chemuson/tree/fix/beta-properties-panel-polish.
- Commits: `95bec64` código/pruebas; `4f8c432` OpenSpec/notas/matriz; `6287bf5` ajuste de ancho encontrado por CI; `ee75534` registro final de validación. No se abrió PR, no se integró a prep y no se creó/publicó una beta o tag.
