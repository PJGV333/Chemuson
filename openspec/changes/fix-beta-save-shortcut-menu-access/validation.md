# Validación — guardado, cierre y barra de menús

Fecha: 2026-10-15. Rama: `fix/beta-save-shortcut-menu-access`.

## Rama y causa raíz

- Base solicitada: `release/v0.3.0-beta.1-prep` y `origin/release/v0.3.0-beta.1-prep`, ambos en `9bc8192db3835262191ddb5247a8ff87f4d9169f` al crear la rama correctiva.
- El diagnóstico reproducible y outputs baseline están en [`baseline.md`](baseline.md).
- **Causa confirmada:** `action_save` tenía shortcuts Ctrl+S y Ctrl+G y estaba conectada al handler correcto, pero sólo asociada a QMenu/toolbar ocultos. En modo compacto Qt no activaba ese QAction. Asociar el mismo QAction a QMainWindow registra ambos atajos sin duplicarlos y conserva `_on_file_save` → `FileController.save_file`.
- `FileController` ya actualizaba `.cmsn` existente, solicitaba ruta a documentos nuevos y sólo marcaba clean tras escribir con éxito. No se añadió otro flujo ni se modificó su comportamiento; pruebas de interacción ahora lo prueban con bytes/JSON reales.
- El flujo existente de cierre ya condicionaba cerrar tras Save a que el canvas objetivo quedara clean; los tests confirman esa guarda y que cancelar/fallar mantiene la pestaña dirty íntegra.

## Cambios

- `src/chemuson/gui/actions/file_actions.py`: asociación de la acción Guardar original a la ventana para shortcuts.
- `src/chemuson/gui/shell/assembly.py`: aplica la visibilidad configurada al QMenuBar existente.
- `src/chemuson/gui/main_window.py`, `src/chemuson/gui/main_window_ui_builder.py`: opción checkable Ver → Mostrar barra de menús y persistencia/reaplicación.
- `src/chemuson/platform/settings.py`: preferencia `UiPreferences.show_menu_bar` (`ui/show_menu_bar`, default compacto), preservada también al escribir preferencias de tema.
- `tests/test_beta_save_shortcut_menu_access.py`: 15 casos Qt con QTest, QSettings/archivos temporales, guardado efectivo, pestañas múltiples, errores/cancelaciones, modos de menús, temas y tamaños.
- `tests/test_platform_settings.py`: cobertura de normalización/round-trip de la preferencia.
- `tests/ci/expected_pytest_nodeids.txt`: sólo 15 node IDs nuevos, revisados; 2169 IDs actuales, 0 retirados.
- `openspec/changes/fix-beta-save-shortcut-menu-access/`: proposal, design, tareas, spec delta, baseline y este informe.

No se tocaron Clean2D, ChemName, ChemIO, química/estereoquímica ni persistencia interna CMSN. No se modificaron release workflows, beta, `main`, stable o tags.

## Verificación automatizada

- Regresión nueva: **15 passed**; grupo de interfaz/settings/pestañas/hamburguesa: **131 passed**; guardado MOL/document controller: **20 passed**.
- Plan completo actualizado: **2169 node IDs / 8 shards** (`[272, 271, 271, 271, 271, 271, 271, 271]`). Los ocho procesos se limitaron externamente a 300 s cada uno; el más largo duró **185.875 s**. JUnit reconciliado: **2149 passed, 20 skipped, 0 failed, 0 errors**. `tools/ci_pytest_shards.py verify` confirmó cobertura única completa de 2169/2169.
- Arquitectura: **280 passed**.
- `tests/test_test_workflow_contract.py`: **4 passed**; el plan de shards verificó el manifiesto exacto.
- `python -m compileall -q src tests tools packaging`: PASS.
- Ruff focalizado (`F401,F811,F821,E722,E741`) sobre todos los Python modificados: PASS. Ruff global baseline anterior aún reporta el F401 histórico `math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; no se alteró.
- `openspec validate fix-beta-save-shortcut-menu-access --strict`: PASS.
- `openspec validate prepare-chemuson-0-3-0-beta-1 --strict`: PASS.
- `git diff --check`: PASS.
- El baseline monolítico `python -m pytest -q` anterior a cambios alcanzó el límite externo de 480 s; no se declaró PASS ni se repitió. La suite completa posterior pasó en ocho shards acotados.

## Estado de entrega

- Commit de implementación: `24b94f920ed6fa77afc711b9af624a1f89757946` (`fix(gui): restore save shortcuts and optional menus`), descendiente directo de la base indicada arriba.
- Rama `fix/beta-save-shortcut-menu-access` publicada por push normal en `origin`; HEAD remoto verificado en ese SHA. La rama beta sigue en la base original `9bc8192db3835262191ddb5247a8ff87f4d9169f`.
- CI GitHub Actions **success** para el SHA de implementación: [run 38769058477](https://github.com/PJGV333/Chemuson/actions/runs/38769058477). Pasaron pytest-plan, los ocho shards y pytest-summary (reconciliación total), windows-smoke y flatpak-smoke. Hubo únicamente avisos deprecativos de Node.js 20 y migración futura de `ubuntu-latest`; no bloquearon la ejecución.
- No se hizo merge a beta, creado tag/release ni iniciado publicación de artefactos. La aceptación visual interactiva manual sigue pendiente del propietario; las pruebas Qt no la sustituyen.

## Comprobación manual breve

1. Abre la ventana con el modo compacto, añade un átomo y pulsa Ctrl+S: el diálogo debe pedir una ruta `.cmsn`. Modifica de nuevo y pulsa Ctrl+S: debe actualizar el mismo archivo sin diálogo. Repite desde **Archivo → Guardar** y verifica el contenido al reabrir.
2. En **Ver → Mostrar barra de menús**, alterna a clásico: Archivo, Editar, Ver, Estructura, Reacción y Ayuda deben permanecer visibles. Reinicia y confirma que la preferencia se conserva; vuelve a compacto con la misma opción desde Ver. Prueba claro/oscuro y una ventana estrecha/ancha.
3. Con dos documentos, modifica el no activo y pulsa su X. Prueba Guardar, Descartar y Cancelar; Cancelar y un fallo/cancelación del diálogo de Guardar deben dejarlo abierto y dirty. Confirma que Guardar escribe su `.cmsn` antes de cerrarlo y que la pestaña activa no objetivo permanece intacta.