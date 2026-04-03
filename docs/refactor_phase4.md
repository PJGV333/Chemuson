# Refactor fase 4

## Objetivo
Reducir el rol hub de `main_window.py` mediante la extracción de operaciones de documento.

## Cambios
- Se agrega `DocumentController` con operaciones de archivos recientes, ruta de documento y confirmación de descarte.
- `main_window.py` delega:
  - `_load_recent_files`
  - `_save_recent_files`
  - `_set_canvas_file_path`
  - `_update_tab_title`
  - `_add_recent_file`
  - `_confirm_discard_changes`
- Se divide `project_actions.py` en módulos por dominio:
  - `view_actions.py`
  - `numbering_actions.py`
  - `structure_actions.py`
- `project_actions.py` queda como compositor liviano.
