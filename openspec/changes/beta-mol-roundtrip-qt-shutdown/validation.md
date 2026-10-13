# Validación — `beta-mol-roundtrip-qt-shutdown`

Fecha: 2026-10-13. Base verificada: `c75026d5656972e0ee3f002f615d80060113bc68`; rama: `fix/beta-mol-roundtrip-qt-shutdown`.

Commits funcionales: `3471197` (ChemIO/MOL) y `258c67e` (teardown Qt). El informe OpenSpec forma un commit documental independiente.

## Evidencia focalizada

Todos los comandos pytest fueron limitados externamente a ≤300 s; no se ejecutó la suite monolítica.

- `PYTHONPATH="$PWD/chem/lib/python3.14/site-packages:$PWD/src" timeout --signal=TERM --kill-after=5s 300s /usr/bin/pytest -q tests/test_mol_interoperability_p0.py tests/test_molblock_normalization.py tests/test_rdkit_roundtrip.py tests/test_smiles_stereo_import.py tests/test_cml_io.py`: **49 passed, 4 skipped in 18.77 s**. Incluye importación real por `FileController`, rutas de exportación ChemUSON/RDKit, escala, no mutación, S/R, fallback, CML, y casos existentes de carga/isótopo/E-Z y stereo no asignado.
- `PYTHONPATH="$PWD/chem/lib/python3.14/site-packages:$PWD/src" timeout --signal=TERM --kill-after=5s 300s /usr/bin/pytest -q tests/test_p0_qt_shutdown_subprocess.py tests/test_gui_async_worker_shutdown.py tests/test_molecular_assistant_recovery.py tests/test_main_window_tabs.py tests/test_autosave_manager.py tests/test_persistence_contract.py tests/test_ui_app_bar_tabs.py`: **88 passed in 26.40 s**.
- `tests/test_drag_move_undo.py tests/test_molecular_assistant_lifecycle.py`: **12 passed in 3.12 s**; orden inverso: **12 passed in 3.21 s**. Ambos procesos acotados a 150 s. No apareció el SIGSEGV de orden de tests histórico.
- Repetición final tras los últimos ajustes de pruebas: MOL P0 **13 passed in 4.11 s**; shutdown P0 **2 passed in 4.89 s**.

## Verificación de proyecto

- `pytest -q tests/architecture`: **280 passed in 10.58 s**.
- `PYTHONPATH="$PWD/chem/lib/python3.14/site-packages:$PWD/src" pytest --collect-only -q`: **2145 tests collected in 0.94 s**, desde 2130 en el baseline (+15 nuevos).
- `python -m compileall -q src tests tools packaging`: PASS.
- `openspec validate beta-mol-roundtrip-qt-shutdown --strict`: PASS.
- Ruff focal sobre todos los Python cambiados con `--select F401,F811,F821,E722,E741`: PASS.
- Ruff global con los mismos selectores mantiene únicamente el hallazgo baseline no relacionado `F401 math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; no se modificó ese archivo.
- `git diff --check`: PASS.
- Contrato `.cmsn` sin cambios; `tests/test_persistence_contract.py` y autosave/recovery incluidos en el shard Qt anterior.

## Teardown Qt: causa y estado

La reproducción previa al fix de `exported` terminó `SIGABRT` (`-6`) por excepción PyQt, no SIGSEGV. GDB mostró `QUndoStack::~QUndoStack → clear → indexChanged → PyQtSlotProxy → ChemusonWindow._on_active_undo_index_changed`; el slot de la ventana accedía al canvas tras su destrucción. `_begin_async_worker_shutdown` ahora desconecta las señales del canvas activo sólo después de que el cierre es aprobado. La nueva matriz de subprocess (normal, exportado, modificado, canvas cerrado, async cleanup y autosave recuperado) pasa y verifica destrucción/threads.

Esto corrige la causa confirmada de la reproducción source; no prueba que sea idéntica al SIGSEGV informado en el AppImage antiguo. **El SIGSEGV exacto de AppImage sigue sin reproducirse y no se declara resuelto**. El retest manual del propietario sigue pendiente; las C superpuestas no reaparecieron en el caso de recuperación controlada.

## CI y paquetes

El push de la rama fix ejecutará CI por shards; sus resultados se añadirán después de observar el SHA remoto. No se construyó un AppImage local. Build Preview sólo admite `release/**-prep`, por lo que esta rama `fix/**` no es elegible para generar el paquete; no se modificó ni publicó otra rama para sortear esa política. Un AppImage actualizado sólo se generará mediante Build Preview en el flujo de preparación autorizado. DATA-11 **no** se declara completo.
