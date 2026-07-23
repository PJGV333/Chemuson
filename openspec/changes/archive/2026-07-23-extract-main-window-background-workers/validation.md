# Validación

Fecha: 2026-07-23

## Resultado

- Caracterización AST ejecutada directamente: **3/3 contratos pasan**.
- `python -m compileall -q src tests tools packaging`: correcto.
- `openspec validate extract-main-window-background-workers --strict`: válido.
- `openspec validate --all --strict`: **21/21 válidos**.
- `git diff --check`: correcto.
- Prueba manual en Qt real: Chemuson inicia y funciona correctamente desde la
  rama publicada; aprobación explícita del usuario recibida el 2026-07-23.
- Pytest y Ruff no están instalados en la imagen actual; la limitación ya está
  registrada en `baseline.md`.

## Evidencia arquitectónica

- `_DescriptorWorker` y `_NameToStructureWorker` tienen un único propietario
  en `gui/background_workers.py`.
- Las señales, constructores, imports diferidos, timeouts y emisiones se
  trasladaron literalmente.
- `background_workers.py` no importa `main_window`.
- `ChemusonWindow` conserva ambos `QThread`, conexiones, registros de jobs,
  cancelación, progreso y callbacks.
- M08 conserva dependencias, API pública, cero excepciones y cero ciclos.

## Alcance revisado

No se modificaron handlers, UI visible, canvas, items, química, persistencia,
formatos, renderizado ni estilos. La prueba manual en Qt real fue satisfactoria
y autoriza el archivado.
