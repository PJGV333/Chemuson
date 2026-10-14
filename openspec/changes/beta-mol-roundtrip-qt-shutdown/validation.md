# Validación — `beta-mol-roundtrip-qt-shutdown`

Fecha: 2026-10-13. Base verificada: `c75026d5656972e0ee3f002f615d80060113bc68`; rama: `fix/beta-mol-roundtrip-qt-shutdown`.

Commits funcionales: `3471197` (ChemIO/MOL) y `258c67e` (teardown Qt). `dc44424` añadió el OpenSpec inicial; `f198095` actualizó el manifiesto CI de los 15 tests nuevos y registró el estado remoto.

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

- CI inicial del push, run [38483483503](https://github.com/PJGV333/Chemuson/actions/runs/38483483503), SHA `dc44424aef6bb96e25e9513a1e30728120eaa99a`: Windows smoke e instalador **PASS**; Flatpak manifest **PASS**; `pytest-plan` falló antes de ejecutar shards porque el manifiesto versionado esperaba 2130 IDs y la colección actualizó a 2145 (+15 tests de este cambio). El summary falló como consecuencia y ningún shard corrió. No fue un fallo de tests de producto.
- Se regeneró explícitamente `tests/ci/expected_pytest_nodeids.txt` mediante `PYTHONPATH="$PWD/chem/lib/python3.14/site-packages:$PWD/src" /usr/bin/python tools/ci_pytest_shards.py update-expected --confirm`. Verificación local: `tools/ci_pytest_shards.py plan --shards 8 ...`: **2145 node IDs**, distribución `[269, 268, 268, 268, 268, 268, 268, 268]`, SHA-256 `23b4746059ae9fee60abc33381e95f2e1d75fe582a12f693a1d2a3a7ae40fc1e`; `tests/test_test_workflow_contract.py`: **4 passed**. La lista sólo incorpora los 15 node IDs nuevos solicitados; no elimina/oculta fallos.
- CI con manifiesto corregido, run [38483720208](https://github.com/PJGV333/Chemuson/actions/runs/38483720208), SHA `f1980958d80b3195da4b039be15364ebc66ddfa6`: **PASS**. Plan + resumen y los ocho shards (0–7) pasaron; se asignaron/ejecutaron **2145/2145** node IDs, **20 skipped**, **0 failures/errors**, sin timeout. Windows smoke/instalador y Flatpak manifest también pasaron.
- No se construyó un AppImage local. Build Preview sólo admite `release/**-prep`, por lo que esta rama `fix/**` no es elegible para generar el paquete; no se modificó ni publicó otra rama para sortear esa política. Un AppImage actualizado sólo se generará mediante Build Preview en el flujo de preparación autorizado. DATA-11 **no** se declara completo.

## Addendum 2026-10-14 — exportación MOL parcial y bloqueo de integración

### Verificación previa

- `origin/fix/beta-mol-roundtrip-qt-shutdown` se verificó en `370f9f968070b96c1a8497155e735eab34b6eb9f`; `origin/release/v0.3.0-beta.1-prep` en `c75026d5656972e0ee3f002f615d80060113bc68`. El árbol fuente estaba limpio; `git rev-list --left-right --count destino...origen` dio `0 5` (fast-forward posible antes de esta corrección).
- La revisión de `FileController.save_file` confirmó un defecto P0: al elegir MOL sólo serializaba `canvas.graph`, después reasignaba la ruta del documento, limpiaba el undo stack y ejecutaba limpieza/reescritura de autosave. Texto y flechas del canvas no están en el CTAB. El test anterior `test_save_as_mol_adds_extension_once_and_updates_document_after_write` exigía explícitamente esa falsa semántica de guardado.
- Antes de la corrección: compileall PASS; colección **2145 tests**; el Ruff global encontró sólo el F401 baseline `math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; no se corrió pytest monolítico. Los tests MOL/persistencia/autosave previos (40 passed) comprobaban la ruta, pero no detectaban el riesgo porque aceptaban que MOL limpiara el estado del documento.

### Corrección acotada

- MOL ahora se expone como `Archivo > Exportar como > Estructura molecular (MOL)...`; Guardar no acepta `.mol` como documento nativo y un documento abierto desde MOL solicita una ruta nativa al guardar cambios.
- El exportador escribe únicamente el grafo molecular, agrega `.mol` una sola vez y deja intactos la ruta `.cmsn`, el estado dirty, el autosave (ruta/limpieza), recientes y los objetos de texto/flecha. El mensaje indica explícitamente que MOL no contiene texto ni flechas y no equivale a guardar el documento completo.
- Se añadió una regresión con `ChemusonCanvas` real, dos átomos/enlace, una anotación y una flecha; verifica el CTAB sólo molecular y que no cambian ruta/dirty/autosave/recientes. La generación del bloque ocurre antes de crear el archivo, para no dejar un MOL vacío si el escritor falla. Se actualizaron el probe de cierre, contrato OpenSpec y matriz manual.
- `DATA-11` queda **PASS manual sólo para `(2S)-2-fluoropropan-1-ol` (`C[C@H](F)CO`)**, basado en confirmación del propietario: exportó MOL y reimportó en ChemUSON con estructura y estereoquímica representada conservadas. No se proporcionó en este registro SHA de paquete ni adjunto del archivo. No se generaliza; aminoácido, enantiómero opuesto, E/Z y otros casos siguen `NOT TESTED`. Se añade `DATA-13` para el documento con anotación/flecha; aceptación manual pendiente, separada del test automatizado.

### Verificación posterior

- MOL/round-trip (`test_mol_interoperability_p0.py`, `test_molblock_normalization.py`, `test_rdkit_roundtrip.py`, `test_smiles_stereo_import.py`): **45 passed, 4 skipped en 23.35 s**.
- Guardado/persistencia/autosave (`test_main_window_tabs.py`, `test_document_controller.py`, `test_persistence_contract.py`, `test_autosave_manager.py`): **49 passed en 12.42 s**.
- Cierre/Qt (`test_p0_qt_shutdown_subprocess.py`, `test_gui_async_worker_shutdown.py`, `test_molecular_assistant_lifecycle.py`): tras actualizar el probe para llamar al nuevo export explícito, **13 passed en 11.24 s**. La primera ejecución del probe agotó su timeout de 25 s porque aún invocaba Guardar con el antiguo filtro MOL y abría el aviso modal; se corrigió el estímulo y la repetición acotada pasó.
- Contrato del test runner/manifiesto: **4 passed**; colección final **2147 node IDs**, 8 shards balanceados `[269,269,269,268,268,268,268,268]`, hash `2563961b73babdce04fc652dca499d0a783fb2bb0ba0600d568c6b5870112274`.
- Arquitectura **280 passed en 16.32 s**; Ruff focal **PASS**; compileall y `git diff --check` **PASS**. Los dos OpenSpec afectados pasan `openspec validate <change> --strict`. `openspec validate --all --strict` conserva **33 fallos preexistentes** de specs no relacionados por placeholders de `Purpose` y 27 cambios/specs validados; no se alteraron contratos ajenos.
- Ruff global mantiene el F401 baseline citado. No se ejecutó pytest monolítico. La corrección está en el commit local `6df2ac01fabb29122e603278b239412208a2e64b`, pero no se ha hecho push/fast-forward; `gh auth status` reporta sesión ausente y `origin/fix/...` continúa en `370f9f968070b96c1a8497155e735eab34b6eb9f`. El propietario debe publicar la rama fuente por GitHub Desktop antes de volver a validar el fast-forward.
- No hay CI ni Build Preview para el código local nuevo. No se han iniciado tags/releases ni modificado `main`, `gh-pages`, canales públicos, Clean2D, ChemName, centros quirales ni persistencia `.cmsn`. El retest del SIGSEGV AppImage y los demás casos manuales permanecen pendientes.

**Dictamen de integración: BLOCKED** hasta publicar la corrección P0 en la rama fuente, verificar CI en su SHA, volver a comprobar el fast-forward y completar el preview de cuatro formatos.
