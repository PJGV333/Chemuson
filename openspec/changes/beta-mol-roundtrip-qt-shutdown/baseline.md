# Baseline — Microcampaña P0 MOL/Qt

Fecha: 2026-10-13. Rama creada: `fix/beta-mol-roundtrip-qt-shutdown`. Antes de crearla, `git status --short --branch` devolvió `## release/v0.3.0-beta.1-prep...origin/release/v0.3.0-beta.1-prep` y el SHA `c75026d5656972e0ee3f002f615d80060113bc68`, igual al SHA solicitado y al HEAD remoto. El árbol de trabajo estaba limpio. Tras `git switch -c fix/beta-mol-roundtrip-qt-shutdown`, el nuevo branch parte del mismo SHA.

## Baseline obligatoria y focal

Comandos ejecutados antes de editar código; cada pytest tuvo timeout externo `timeout --signal=TERM --kill-after=5s 300s`. No se ejecutó `pytest -q` monolítico por instrucción expresa.

- `python -m compileall src tests tools packaging`: exit 0.
- `pytest --collect-only -q`: exit 0, **2130 tests collected in 0.67s**.
- `PYTHONPATH="$PWD/chem/lib/python3.14/site-packages:$PWD/src" timeout ... /usr/bin/pytest -q tests/test_smiles_stereo_import.py tests/test_rdkit_roundtrip.py tests/test_molblock_normalization.py`: **31 passed, 4 skipped in 13.39s**.
- `PYTHONPATH="$PWD/chem/lib/python3.14/site-packages:$PWD/src" timeout ... /usr/bin/pytest -q tests/test_gui_async_worker_shutdown.py tests/test_molecular_assistant_recovery.py tests/test_autosave_manager.py`: **41 passed in 2.30s**.
- `PYTHONPATH="$PWD/chem/lib/python3.14/site-packages:$PWD/src" timeout ... /usr/bin/pytest -q tests/architecture`: **280 passed in 12.65s**.
- `ruff check src tests tools packaging --select F401,F811,F821,E722,E741`: exit 1 únicamente por el F401 baseline `math` no usado en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`. No se modifica ese test.
- `git status --short` al inicio: salida vacía.

## Reproducción inicial de MOL

Entorno source con PyQt6/RDKit instalados; `QT_QPA_PLATFORM=offscreen`. Una ventana real `ChemusonWindow` inició, insertó `C[C@H](F)CO`, llamó al flujo de Guardar mediante filtro `Archivo MOL (*.mol)` y cerró. Salida relevante:

```text
phase=started
phase=structure
phase=export {"exact_path_exists": true, "mol_path_exists": false, "canvas_path": "/tmp/.../prueba_mol_1", "clean": true}
phase=import {"atoms": 5, "bonds": 4}
phase=close True visible= False
phase=exit
```

El archivo sin sufijo se interpreta si se entrega por ruta/contenido, pero no tiene el sufijo que usan el filtro MOL de apertura y algunos lectores externos. La escritura sí limpia el documento porque no ocurrió error; no es evidencia de una falsa limpieza tras fallo.

RDKit baseline del fixture:

| SMILES | Fórmula | Átomos pesados / enlaces | Mediana interna | Mediana CTAB actual | Canonical isomeric SMILES tras import |
|---|---|---:|---:|---:|---|
| `C[C@H](F)CO` (S) | `C3H7FO` | 5 / 4 | 40.0009 | 40.0010 | `C[C@H](F)CO` |
| `C[C@@H](F)CO` (R) | `C3H7FO` | 5 / 4 | 40.0009 | 40.0010 | `C[C@@H](F)CO` |

El grafo fuente usa escala de lienzo (~40); ambos writers exportan tal cual la misma escala, que es la causa geométrica confirmada. Una exploración del fallback con los mismos enlaces cuña/hash mantuvo el SMILES S/R, pero `molfile_to_molgraph` no volvió a poblar `Atom.stereo_cip` tras analizar cuñas; la ausencia de metadata CIP en el MolGraph importado requiere prueba/corrección explícita aunque la estereoquímica derivada del wedge siga legible.

## Cierre Qt y autosave — reproducción baseline

- Una segunda reproducción source, con ventana mostrada, fixture insertado, guardar tras seleccionar MOL y luego `window.close()`, terminó con código 0; no hubo SIGSEGV ni excepción Python. En un probe de destrucción de canvas se usó inicialmente `processEvents()` sin entregar correctamente `DeferredDelete`; por ello ese probe no acredita destrucción C++ y se repetirá con `QEvent.Type.DeferredDelete` y subproceso.
- Snapshot real por `AutosaveManager` y recuperación a través de `ChemusonWindow._open_autosave_document`: **5 átomos, 4 enlaces, undo stack limpio**. Los tres carbonos quedaron como `is_explicit=false` y sus etiquetas `C` estaban invisibles; F/O fueron las únicas etiquetas visibles. No se reprodujeron etiquetas C superpuestas ni duplicación de átomos/enlaces. Es evidencia source, no retest del AppImage.
- El código existente ya hace confirmaciones antes del shutdown, detiene autosave y espera `QThread` descendientes. La corrección histórica de test `indexChanged` no se toma como causa del crash manual. La observación AppImage `SIGSEGV` sigue sin causa raíz y no se declara resuelta por esta única no-reproducción.

## Herramientas/entorno

- Test runner: `/usr/bin/pytest`; RDKit se hace visible mediante `PYTHONPATH="$PWD/chem/lib/python3.14/site-packages:$PWD/src"`.
- Los casos Qt propensos a bloqueo/crash se ejecutarán mediante subprocess con timeout de Python y timeout externo del comando pytest.
- No se ejecutó un AppImage viejo como evidencia de corrección; el nuevo paquete, de solicitarse/autorizarse, sólo se generará con Build Preview sobre el SHA revisado.

## Reproducción adicional y traza Qt (antes de la corrección de cierre)

La primera prueba sólo confirmó `window.close()` y no forzó `DeferredDelete`; por eso había salido 0 sin ejercitar destructor de `QUndoStack`. Se cambió el estímulo de diagnóstico: después de guardar MOL, se cierra la ventana y se fuerza `QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)` sobre el canvas.

- Comando de reproducción: `ulimit -c 0; PYTHONPATH="$PWD/chem/lib/python3.14/site-packages:$PWD/src" QT_QPA_PLATFORM=offscreen XDG_CONFIG_HOME=$(mktemp -d) timeout --signal=TERM --kill-after=5s 45s gdb --batch --quiet -ex 'set pagination off' -ex run -ex 'thread apply all bt 20' --args /usr/bin/python tests/probes/p0_shutdown_workflow.py exported`.
- Resultado previo al fix: Python/PyQt traceback y salida **SIGABRT (`-6`)**, no SIGSEGV. Excepción: `RuntimeError: wrapped C/C++ object of type ChemusonCanvas has been deleted`, desde `ChemusonWindow._on_active_undo_index_changed` al leer `self.canvas.current_opacity_percent()`.
- Backtrace útil: `pthread_kill → raise → abort → QMessageLogger::fatal → pyqt6_err_print → PyQtSlotProxy::unislot → QUndoStack::indexChanged → QUndoStack::clear → QUndoStack::~QUndoStack → QObjectPrivate::deleteChildren → QWidget::~QWidget → sipQGraphicsView::~sipQGraphicsView`. Los demás hilos visibles eran threads pooled Qt inactivos.
- Causa confirmada para ese proceso: señal `indexChanged` del QUndoStack seguía conectada al slot de la ventana mientras se destruía el canvas. No se atribuye automáticamente el SIGSEGV AppImage a la misma señal, pero sí queda sustentado el defecto de ciclo de vida source relacionado al cierre.
- El core dump fue deshabilitado con `ulimit -c 0`; no quedan procesos/cores de la reproducción. La corrección desconecta señales del canvas activo tras aprobar cierre y antes de liberar widgets; las pruebas Qt vuelven a ejecutarse en subprocess con timeout.
