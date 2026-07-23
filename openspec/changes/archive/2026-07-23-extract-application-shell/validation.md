# Validación

Fecha: 2026-07-23

## Resultado

- `pytest -q tests/test_application_startup_non_qt.py tests/architecture`:
  **210 passed**.
- `python -m compileall -q src tests tools packaging`: correcto.
- `openspec validate extract-application-shell --strict`: válido.
- `openspec validate --all --strict`: **20/20 válidos**.
- `git diff --check`: correcto.
- Ruff global con `F401,F811,F821,E722,E741`: conserva únicamente el F401
  histórico de `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py`;
  los archivos de esta fase no introducen diagnósticos.
- Suite completa: la colección se detiene con los mismos **56 errores** por
  ausencia de `libEGL.so.1`; no aparece un error distinto atribuible al cambio.

## Evidencia arquitectónica

- `ChemusonWindow.__init__` ejecuta `QMainWindow.__init__` y delega exactamente
  una vez en `assemble_application_shell(self)`.
- `main_window.py` ya no construye directamente tabs, docks, toolbars ni
  widgets de status.
- `gui/shell/` pertenece exclusivamente a M08 y no importa
  `gui.main_window`, por lo que no introduce ciclo.
- El catálogo conserva las mismas dependencias, cero excepciones temporales y
  cero ciclos.
- `ChemusonWindow` permanece en su ruta pública histórica.

## Alcance revisado

No se modificaron handlers, canvas, items, controllers, química, persistencia,
formatos, renderizado ni estilos. La comprobación visual y funcional completa
debe realizarse en un entorno con Qt/EGL antes de archivar el cambio.

## Validación manual

El 2026-07-23 el usuario probó la rama
`architecture/phase7-extract-application-shell` en su entorno Qt real y
confirmó que Chemuson inicia y funciona correctamente. La prueba cubrió la
composición visible de menús, toolbars, pestaña inicial, docks, canvas, barra
de estado y cierre normal de la ventana. Con esta aprobación explícita se
cumple la tarea 6.1 y el cambio puede archivarse.
