# Baseline — estabilización de pytest CI

## Git

- Rama de campaña creada: `fix/ci-pytest-stabilization`.
- Base/upstream solicitado: `fix/chemio-stereo-roundtrip` en `f5f82a63c84c56c9dd3ee07e7c7f56e7eb700bef`.
- `HEAD` inicial y `origin/fix/chemio-stereo-roundtrip`: ambos `f5f82a63c84c56c9dd3ee07e7c7f56e7eb700bef`.
- `git status --short`: salida vacía antes de crear este OpenSpec.
- No se cambiaron ramas protegidas ni se transportaron ediciones de ChemName.

## Python / colecta / estáticos

Entorno de compatibilidad con CI instalado en `/tmp/chemuson-ci311`: Python 3.11.16, pytest 9.1.1, PyQt6 6.11.0 / Qt 6.11.2 y RDKit 2026.9.1.

- `python -m compileall -q src tests tools packaging`: **PASS**, tanto Python 3.14 del sistema como Python 3.11.
- `timeout 180s /tmp/chemuson-ci311/bin/python -m pytest --collect-only -q`: **2.128 tests collected in 2.42s**.
- `ruff check src tests tools packaging --select F401,F811,F821,E722,E741`: **1 error baseline**, ajeno al alcance: `F401 math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`.
- Suite monolítica: **NO EJECUTADA** (restricción de campaña y contaminación/abortos Qt históricos).

## GitHub Actions histórico

- Run #343, ID `37968649720`, URL https://github.com/PJGV333/Chemuson/actions/runs/37968649720: rama `fix/chemio-stereo-roundtrip`, SHA exacto `f5f82a63c84c56c9dd3ee07e7c7f56e7eb700bef`, pytest falló; `windows-smoke` y `flatpak-smoke` pasaron.
- El usuario reporta seis fallos: `tests/test_clean2d_engine_candidates.py::test_generate_candidates_attempts_rdkit_for_cyclic_graphs`, `tests/test_compchem3d_dock.py::test_compchem_controller_generates_async_with_fake_backend`, `tests/test_molecular_assistant_ui.py::test_dialog_is_modeless_masks_key_and_rejects_missing_request_fields`, dos pruebas `test_smoke_*` de `tests/test_rdkit_packaged_worker.py`, y `tests/test_side_panel.py::test_primary_tabs_have_complete_labels_padding_and_separation`.
- Run #342, ID `37936302839`, URL https://github.com/PJGV333/Chemuson/actions/runs/37936302839: rama `chemname/iupac-robustness`, SHA `054f6c7fc378975289c15e850a950a257942773b`, terminó con fallo. El usuario confirma que comparte estos seis y añade los tres fallos estéreo anteriores a la corrección ChemIO.
- La API anónima devuelve para #343 sólo la anotación genérica `Process completed with exit code 1`; descargar los logs de Actions respondió 403. No se infieren del API mensajes de test que no están disponibles.

## Nodos focalizados aislados (Python 3.11)

Cada nodo se ejecutó en proceso pytest nuevo con `QT_QPA_PLATFORM=offscreen`, `XDG_CONFIG_HOME` temporal y `timeout` externo de 45 s:

| Nodo | Resultado baseline |
|---|---|
| Clean2D `test_generate_candidates_attempts_rdkit_for_cyclic_graphs` | **FAIL reproducido**: salida `{clean2d_v2,current,rdkit_direct,simple_aromatic_template}`; el test exige `rdkit_isolated` tras deduplicación. |
| CompChem `test_compchem_controller_generates_async_with_fake_backend` | PASS (1 en 0.20 s). |
| Diálogo Molecular `test_dialog_is_modeless_masks_key_and_rejects_missing_request_fields` | PASS (1 en 0.43 s). |
| Worker empaquetado `test_smoke_allows_qt_runtime_hooks_without_a_gui_application` | PASS aislado (1 en 0.04 s). |
| Worker empaquetado `test_smoke_detects_active_qt_application_and_chemuson_gui` | PASS aislado (1 en 0.04 s). |
| Panel lateral `test_primary_tabs_have_complete_labels_padding_and_separation` | PASS aislado (1 en 0.43 s). |

## Reproducciones de orden/estado

- Grupo de colección `test_ai_action_is_discoverable_in_structure_menu_and_command_palette` + ambos nodos `test_smoke_*`: **1 passed, 2 failed**. La colección previa de `test_molecular_assistant_ui.py` importó 110 módulos `chemuson.gui`; los contratos del smoke inspeccionan el `sys.modules` compartido. Los mismos dos nodos pasan solos.
- Probe de QSettings en proceso aislado: tras crear `QApplication` y escribir `ai/molecular_assistant/resolution_method=ai` en NativeFormat, el test del diálogo falla en `identity_verification_enabled is True`. Confirma que cambiar sólo `XDG_CONFIG_HOME` después de crear QApplication no cambia la ruta ya cacheada. El origen externo preciso del valor en Actions permanece sin logs; el fixture sí carece del aislamiento NativeFormat necesario.
- Probe CompChem instrumentado de baseline: orden observado `worker.finished → thread.finished → job_finished`, un resultado y cero jobs activos; la prueba aislada también pasa. No reproduce una carrera local.
- El panel no reproduce el `Inspector x=-3` reportado por el usuario en CI en Python 3.11 local. Una inspección única de geometría estable a 1440x900 dio viewport 308 px, strip 305 px, scroll 0 y `Inspector x=2`; una prueba aislada en 1440x900 y 980x600 pasa. No hay todavía evidencia local para tocar `side_panel.py` ni para debilitar la aserción.
- Una exploración redundante que creaba repetidamente ventanas por familia de fuente agotó su límite externo de 30 s antes de producir salida; fue terminada y no dejó procesos ni directorios temporales. No se repitió el mismo diagnóstico.
