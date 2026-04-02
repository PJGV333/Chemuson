# Refactor arquitectónico incremental (abril 2026)

## 1) Mapa real de imports internos (actual)

- Total de aristas internas detectadas: **130**
- Top 10 hubs por dependencias salientes:
  - chemuson.gui.main_window: 22
  - chemuson.gui.canvas: 18
  - chemuson.update: 10
  - chemuson.gui.toolbar: 6
  - chemuson.update.core: 6
  - chemuson.gui.items: 5
  - chemuson.chemname.engine: 4
  - chemuson.gui.commands: 4
  - chemuson.gui.composite_diagram_item: 4
  - chemuson.gui.actions: 4

## 2) main_window como hub (antes/después)

### Antes
```mermaid
graph LR
  MW[main_window] --> N1[chemuson.chemio.persistence]
  MW[main_window] --> N2[chemuson.chemio.rdkit_io]
  MW[main_window] --> N3[chemuson.core.model]
  MW[main_window] --> N4[chemuson.gui.canvas]
  MW[main_window] --> N5[chemuson.gui.commands]
  MW[main_window] --> N6[chemuson.gui.composite_diagram_item]
  MW[main_window] --> N7[chemuson.gui.dialogs]
  MW[main_window] --> N8[chemuson.gui.docks]
  MW[main_window] --> N9[chemuson.gui.energy_diagrams]
  MW[main_window] --> N10[chemuson.gui.icons]
  MW[main_window] --> N11[chemuson.gui.orbitals]
  MW[main_window] --> N12[chemuson.gui.periodic_table]
  MW[main_window] --> N13[chemuson.gui.styles]
  MW[main_window] --> N14[chemuson.gui.template_library]
  MW[main_window] --> N15[chemuson.gui.templates]
  MW[main_window] --> N16[chemuson.gui.text_toolbar]
  MW[main_window] --> N17[chemuson.gui.toolbar]
  MW[main_window] --> N18[chemuson.update]
  MW[main_window] --> N19[chemuson.utils]
  MW[main_window] --> N20[chemuson.utils.autosave]
  MW[main_window] --> N21[chemuson.version]
```

### Después
```mermaid
graph LR
  MW[main_window] --> N1[chemuson.chemio.persistence]
  MW[main_window] --> N2[chemuson.chemio.rdkit_io]
  MW[main_window] --> N3[chemuson.core.model]
  MW[main_window] --> N4[chemuson.gui.actions]
  MW[main_window] --> N5[chemuson.gui.canvas]
  MW[main_window] --> N6[chemuson.gui.commands]
  MW[main_window] --> N7[chemuson.gui.composite_diagram_item]
  MW[main_window] --> N8[chemuson.gui.dialogs]
  MW[main_window] --> N9[chemuson.gui.docks]
  MW[main_window] --> N10[chemuson.gui.energy_diagrams]
  MW[main_window] --> N11[chemuson.gui.icons]
  MW[main_window] --> N12[chemuson.gui.orbitals]
  MW[main_window] --> N13[chemuson.gui.periodic_table]
  MW[main_window] --> N14[chemuson.gui.styles]
  MW[main_window] --> N15[chemuson.gui.template_library]
  MW[main_window] --> N16[chemuson.gui.templates]
  MW[main_window] --> N17[chemuson.gui.text_toolbar]
  MW[main_window] --> N18[chemuson.gui.toolbar]
  MW[main_window] --> N19[chemuson.update]
  MW[main_window] --> N20[chemuson.utils]
  MW[main_window] --> N21[chemuson.utils.autosave]
  MW[main_window] --> N22[chemuson.version]
```

- Imports internos directos `main_window.py`: **21 → 22**.
- Cambio principal: externalización de construcción de acciones por dominio (`gui/actions/*`).

## 3) Archivos creados/movidos/eliminados

### Creados
- `src/chemuson/gui/actions/__init__.py`
- `src/chemuson/gui/actions/file_actions.py`
- `src/chemuson/gui/actions/edit_actions.py`
- `src/chemuson/gui/actions/project_actions.py`
- `src/chemuson/gui/actions/update_actions.py`

### Movidos
- Ninguno (refactor incremental sin renombrados masivos).

### Eliminados
- Ninguno en esta iteración.

## 4) Justificación breve por cambio

- Se separa la creación de acciones de GUI en módulos por dominio (`file`, `edit`, `project`, `update`) para reducir responsabilidades de `main_window.py`.
- `main_window.py` queda más orientado a coordinación y wiring general.
- Se mantiene compatibilidad funcional porque los `QAction` conservan nombres, atajos y señales.

## 5) Tests obsoletos eliminados y tests conservados

- Tests eliminados: ninguno en esta iteración (no se identificó evidencia suficiente de obsolescencia segura).
- Tests conservados y ejecutados para validar regresión de capas core/update:
  - `tests/test_core_model.py`
  - `tests/test_update_core.py`
