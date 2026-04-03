# Refactor fase 2

## Objetivo
Reducir acoplamiento de `gui/main_window.py` delegando flujos complejos a controladores.

## Controladores incorporados
- `ExportController`: export PNG/SVG/PDF.
- `Clean2DController`: pipeline clean-2D con fallback acíclico y degradación sin RDKit.
- `TextFormatController`: formato de texto externo temporal y sincronización toolbar.
- `RecoveryController`: metadata/listado/archivo/apertura de autosaves y UI de recovery center.

## Delegaciones aplicadas en main_window
- `_on_export` delega a `ExportController.export`.
- `_run_clean_2d` delega a `Clean2DController.run_clean_2d`.
- Formato de texto externo delega a `TextFormatController`.
- `_read_autosave_metadata`, `_list_autosave_entries`, `_archive_autosave`, `_open_autosave_document`, `_show_recovery_center`, `check_autosaves` delegan a `RecoveryController`.

## Cambios estructurales
- `docs/refactor_arquitectura_2026_04.md` movido a `docs/archive/`.
- `tests/data/orbitals/fit_report/` movido a `tests/archive/orbitals_fit_report/`.
