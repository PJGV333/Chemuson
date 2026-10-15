# Validación

## Baseline (antes de la corrección)

- Rama de trabajo: `fix/beta-save-shortcut-menu-access`; HEAD inicial `2c705ad`
  (árbol limpio, sincronizado con `origin`).
- `python -m compileall -q src tests tools packaging`: PASS.
- `pytest --collect-only -q`: 2174 IDs; el manifiesto revisado pasó a 2200 IDs.
- Focalizado de UI/preferencias (`test_beta_save_shortcut_menu_access.py`,
  `test_ui_app_bar_tabs.py`, `test_ui_tool_rail.py`, `test_platform_settings.py`):
  **120 passed** en 143 s.
- Ruff `--select F401,F811,F821,E722,E741` sobre todo el árbol: único hallazgo
  preexistente `F401 math` en
  `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; no se tocó.

## Diagnóstico reproducido en Qt (offscreen, antes de editar)

Medición directa sobre `ChemusonWindow` con INI temporal:

| Superficie | Baseline | Esperado |
| --- | --- | --- |
| `QMenuBar` (oculta por defecto) | alto 43 px, ítems 50 px | visible y compacta |
| `main_toolbar` | alto 65 px, 9 botones + separadores | 3 botones, fila compacta |
| `app_bar` | 54 px | 54 px (sin cambio) |

Pestaña con un solo documento (baseline):

- limpio: `tabRect=(0,0,148,37)`, bloque lateral `28×24`, X en `x=2`, tamaño real
  de la X `26×26` (QSS `min-width/min-height: 24px` + borde 1 px ganaba a
  `setFixedSize(24, 24)`), punto dirty oculto.
- modificado: `tabRect` idéntico pero la X salta a `x=14` con `26×26` dentro de un
  bloque de 28 px → desborda el borde derecho de la pestaña (14+26 > 28+133).
  Causa raíz: el layout del bloque lateral exigía 40 px
  (2+7+5+24+2) y Qt repartía el hueco disponible, así que el punto dirty robaba
  12 px al botón de cierre.

## Verificación de esta corrección

- **A3 (geometría de pestaña):** tras el bloque lateral de tamaño fijo (40×24) con
  el hueco del punto dirty reservado, la medición Qt da `tabRect=(0,0,149,37)`,
  bloque `40×24`, X `24×24` en `x=14` y holgura constante de 3 px al borde
  derecho, **idénticos en limpio → modificado → limpio**. Píxeles pintados con el
  token `accent` dentro del bloque: 0 en limpio, >0 en modificado, 0 al volver a
  limpio; los píxeles del tinte `icon` (la X) son idénticos en ambos estados.
- **A1 (barra de menús):** `QMenuBar` visible por defecto, altura **30 px**
  (baseline 43 px), seis ítems con geometría válida, sin solapamiento con la app
  bar de 54 px. Los seis `QMenu` del agregador hamburguesa son los mismos objetos
  (`_menu_aggregator` → `menuBar()`), sin QActions duplicadas.
- **A2 (barra de herramientas):** `main_toolbar` con
  `[Nuevo, Abrir, Guardar]` como únicas acciones no separador, altura **37 px**
  (baseline 65 px), tres botones de 31×31 con icono de 18 px. Undo/Redo siguen
  sólo en la app bar; rotaciones, volteos, Clean2D y SMILES siguen en Editar /
  Estructura, sus atajos, el rail y la paleta de comandos.
- **Preferencias:** `UiPreferences.show_menu_bar` pasa a `True`; migración única
  por `ui/shell_defaults_version` (`SHELL_DEFAULTS_VERSION = 1`). Un INI con
  `ui/show_menu_bar=false` y `ui/show_main_toolbar=false` previos abre con ambas
  barras visibles; tras la primera elección explícita, la nueva ventana respeta lo
  persistido.

## Suite completa por shards (timeout externo 300 s por shard)

`tools/ci_pytest_shards.py plan --shards 8` recolectó **2200 node IDs**
(distribución 275×8, sha256 `74f3d24259fc67ed058b5549c238650e18358c557ff9fed28ab4a112d7896802`).
Ocho shards ejecutados y verificados con `verify`:
**2200 unique node IDs, 8 complete successful shards** — 2180 passed, 20 skipped,
sin timeouts ni fallos.

| Shard | Resultado | Tiempo |
| --- | --- | --- |
| 0 | 275 passed | 63.9 s |
| 1 | 273 passed, 2 skipped | 87.0 s |
| 2 | 273 passed, 2 skipped | 89.7 s |
| 3 | 270 passed, 5 skipped | 43.3 s |
| 4 | 271 passed, 4 skipped | 80.5 s |
| 5 | 272 passed, 3 skipped | 94.1 s |
| 6 | 271 passed, 4 skipped | 60.0 s |
| 7 | 275 passed | 63.4 s |

- El fallo intermitente CompChem
  (`tests/test_compchem3d_dock.py::test_compchem_controller_generates_async_with_fake_backend`)
  **no reprodujo** en las dos pasadas completas de ocho shards de esta ronda
  (16 ejecuciones de shard, exit 0). Se deja registrado como flaky preexistente
  sin cambiar el test ni sus contratos, conforme a la instrucción de no relajar
  CompChem.
- Arquitectura: `tests/architecture` → **280 passed**.
- Focalizado de shell/preferencias/pestañas:
  `tests/test_beta_compact_shell_and_tabs.py` (25 pruebas nuevas),
  `tests/test_beta_save_shortcut_menu_access.py`, `tests/test_ui_app_bar_tabs.py`,
  `tests/test_ui_tool_rail.py`, `tests/test_platform_settings.py`,
  `tests/test_ui_polish.py`, `tests/test_ui_theme_foundation.py`,
  `tests/test_command_palette.py`, `tests/test_side_panel.py`,
  `tests/test_branch_rotation_shortcuts.py`: PASS.
- `python -m compileall -q src tests tools packaging`: PASS.
- Ruff focalizado (`F401,F811,F821,E722,E741`) sobre los módulos y tests
  cambiados: PASS. En todo el árbol queda únicamente el `F401 math` baseline.
  El hallazgo `I001` de `tests/conftest.py` es preexistente y no se tocó.
- `git diff --check`: PASS.
- `openspec validate fix-beta-tab-close-and-toolbar-visibility --strict` y
  `openspec validate fix-beta-save-shortcut-menu-access --strict`: PASS.

## Decisiones registradas

- La opción heredada **Mostrar copiar/pegar/zoom en barra superior** se conserva
  como opt-in: por defecto la barra tiene exactamente tres botones, y sólo esa
  opción explícita puede añadir acciones auxiliares. Retirarla habría eliminado
  una funcionalidad existente fuera del alcance pedido; queda marcada para
  decisión del propietario.
- Se añadió el fixture `isolated_shell_settings` en `tests/conftest.py` porque
  varios tests construían `ChemusonWindow` con el INI real del equipo, lo que
  hacía no deterministas las pruebas de defaults.
- El marcador `ui/shell_defaults_version` es la única forma de aplicar los
  defaults nuevos de la beta a configuraciones previas sin borrar después las
  elecciones del usuario.

## Aceptación y publicación

- La aceptación visual interactiva en las plataformas del propietario queda
  pendiente: los renders y la automatización no la sustituyen.
- Commit de implementación: `8a61225`
  (`fix(gui): compact shell bars and stable tab geometry`), publicado por push
  normal únicamente a `fix/beta-save-shortcut-menu-access`.
- No se integra a `release/v0.3.0-beta.1-prep` ni a `main`, y no se crean tags ni
  releases.

## Pasos de aceptación manual sugeridos

1. Abrir ChemUSON con la configuración actual: debe verse la barra de menús
   (`Archivo | Editar | Ver | Estructura | Reacción | Ayuda`) en una fila fina,
   debajo la barra con tres iconos (nuevo/abrir/guardar) y debajo la app bar.
   Si en una ejecución anterior se ocultó la barra de menús, activar
   `Ver → Mostrar barra de menús` una vez.
2. Alternar ambas opciones de Ver, cambiar de tema y reiniciar: cada barra debe
   reaparecer como se eligió, sin afectar al panel lateral ni al rail.
3. Escribir en un documento y observar la pestaña: la X no debe moverse al
   aparecer el punto azul; probar pestañas largas, varias, reordenar y cerrar con
   Guardar/Descartar/Cancelar.
4. Comprobar Undo/Redo desde la app bar y sus atajos, y Clean2D/rotaciones/SMILES
   desde Estructura/Editar y la paleta de comandos (Ctrl+P).
