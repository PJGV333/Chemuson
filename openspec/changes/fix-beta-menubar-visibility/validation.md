# Validación — visibilidad de la fila de menús

## Causa raíz y wiring

En `Ver` coexistían dos acciones checkable: `action_show_menu_bar` («Mostrar barra de menús») que controlaba `QMenuBar` y guardaba `ui/show_menu_bar`, y el `toggleViewAction()` de `main_toolbar` («Mostrar barra de herramientas») que mostraba Nuevo/Abrir/Guardar y guardaba `ui/show_main_toolbar`. El label que el propietario quería usar estaba conectado a la superficie equivocada. La hamburguesa ya abría `_menu_aggregator` con los mismos seis QMenu originales.

La corrección reaprovecha `action_show_menu_bar` y su persistencia existente, cambia su label al solicitado, y retira del menú Ver el toggle de la toolbar. La QToolBar legada se fuerza oculta y se deshabilita su toggle nativo, incluso si el INI mantiene `ui/show_main_toolbar=true`. No se cambia la preferencia de menú, los defaults (`show_menu_bar=True`) ni su versión de migración.

## Baseline

- Rama/HEAD: `fix/beta-final-ui-clean2d` / `94ffdd80080c9904a4a2d472b06a5174f65d8630`, igual a `origin`; árbol limpio antes del cambio.
- `python -m compileall src tests tools packaging` pasó.
- `pytest --collect-only -q`: `2236 tests collected in 0.62s` antes del cambio.
- Baseline `pytest -q`, con timeout externo 600 s: exit 124 al alcanzar 38% de progreso; no se mostró fallo de test. No se repitió el monolítico y no quedaron procesos pytest/RDKit activos.
- Baseline focal UI/settings: 127 passed in 165.55s (`test_beta_compact_shell_and_tabs`, `test_ui_app_bar_tabs`, `test_ui_tool_rail`, `test_platform_settings`). Arquitectura: 280 passed in 9.59s.
- Ruff global acotado: único hallazgo baseline conocido `F401 math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`.
- Mapa previo: `action_show_menu_bar` → `QMenuBar`/`ui/show_menu_bar`; `main_toolbar.toggleViewAction()` → `QToolBar`/`ui/show_main_toolbar`; AppBar hamburguesa → agregador de los seis menús originales.

## Verificación local tras el cambio

- Primera pasada de regresión encontró cuatro aserciones antiguas que esperaban que «Mostrar barra de herramientas» mostrara la QToolBar. Se actualizaron para el nuevo contrato; no fue un fallo de runtime ni de acciones del documento.
- `pytest -q tests/test_beta_compact_shell_and_tabs.py`: 25 passed in 18.82s.
- Suite focalizada final, timeout 600 s: `pytest -q tests/test_beta_compact_shell_and_tabs.py tests/test_ui_app_bar_tabs.py tests/test_ui_tool_rail.py tests/test_platform_settings.py tests/test_beta_save_shortcut_menu_access.py tests/test_side_panel.py`: **161 passed in 309.35s**.
- `pytest -q tests/architecture`: 280 passed in 9.29s. `pytest -q tests/test_test_workflow_contract.py`: 4 passed in 0.06s.
- Ruff sobre todos los archivos GUI/tests modificados: pasó. Ruff global mantiene sólo el F401 baseline indicado arriba.
- `python -m compileall src tests tools packaging`: pasó. `openspec validate fix-beta-menubar-visibility --strict --no-interactive`: válido. `git diff --check`: pasó.
- Colecta final/plan: **2235 IDs**, 8 shards con distribución `[280,280,280,279,279,279,279,279]`, SHA-256 `7356386ef65f4686c066facb08da25d8a5aee2cb118c8e50fecc2e66bc9b6aad`. Se retiraron/reemplazaron node IDs de pruebas que afirmaban la visibilidad de la toolbar; el manifiesto fue actualizado con el updater tras revisar la diferencia.

## CI/aceptación

- GitHub Actions run [`38836149649`](https://github.com/PJGV333/Chemuson/actions/runs/38836149649), HEAD de implementación `d4d2f67bd1b4101b9ebd734fb34d096939da46dd`: **success**. Plan, los ocho shards, verificación completa de IDs/JUnit, `windows-smoke` y `flatpak-smoke` pasaron. Shard local más largo en Actions: 4m30s, dentro del límite.
- Push normal únicamente a `fix/beta-final-ui-clean2d`; no se integró a `release/v0.3.0-beta.1-prep`, no se modificó `main` ni se publicaron tags/releases/paquetes.
- La comprobación manual del propietario sigue pendiente.

## Seguimiento de revisión — popup de Ver accesible

### Cambios

La segunda revisión confirmó que corregir sólo el texto/target no bastaba: el popup mantenía 22 opciones directas y casi ocupaba toda la pantalla disponible. `Ver` ahora presenta siete filas, con «Mostrar barra de herramientas» primera y marcada; Visualización, Zoom, Numeración, Tamaño de lienzo, Paletas y Paneles reúnen las QActions existentes. Se movieron las páginas del panel lateral dentro de Paneles conservando sus mismas acciones y conexiones.

Se eliminó la acción obsoleta «Mostrar copiar/zoom en barra de herramientas» de la creación de acciones, el registro de comandos y el wiring de la QToolBar; Copiar permanece en Editar con `QKeySequence.StandardKey.Copy`, y Zoom en su submenú con los atajos estándar. La QToolBar sigue oculta y su toggle nativo deshabilitado. No se modificaron Clean2D, química, ChemIO, ChemName, canvas, Ctrl+S, undo/redo ni pestañas.

### Verificación local

- Baseline antes de esta iteración: `QScreen.availableGeometry=(0,0,800,800)` en Qt offscreen; `Ver` medía `778×786` px, con 22 opciones directas en columnas.
- `python -m pytest -q tests/test_beta_compact_shell_and_tabs.py tests/test_beta_save_shortcut_menu_access.py`: **48 passed in 69.41s**. Incluye mouse QTest desde la QMenuBar, primer/último renglón dentro del QScreen, cada submenú, el matrix de tema claro/oscuro y ventanas 1600×900, 1280×720 y 1024×768, restauración por hamburguesa y persistencia en dos reinicios.
- `python -m pytest -q tests/test_side_panel.py`: **15 passed in 5.03s**. La primera pasada encontró tres pruebas antiguas que aún buscaban las páginas directamente en `Ver`; se actualizaron para navegar por el submenú `Paneles` y la repetición completa pasó.
- Una primera variante del test geométrico basaba los submenús en hover y ventanas separadas; Qt no abrió `Paneles` de forma estable al repetir la matriz. Se cambió a clic QTest y una misma ventana/QMenu a través de los seis casos; la regresión final pasó.
- `python -m pytest -q tests/test_ui_app_bar_tabs.py tests/test_ui_tool_rail.py tests/test_platform_settings.py tests/test_command_palette.py`: **128 passed in 145.43s**.
- `python -m pytest -q tests/architecture tests/test_test_workflow_contract.py`: **284 passed in 10.70s**.
- `python -m compileall src tests tools packaging`: pasó. `git diff --check`: pasó. OpenSpec estricto: válido.
- Colecta: **2238 IDs**, exactamente tres IDs nuevos de las pruebas de Ver y ninguno retirado. El manifiesto se actualizó con el updater tras revisar esos IDs. Plan de 8 shards: `[280,280,280,280,280,280,279,279]`, SHA-256 `d013b567a97a9f162c44598b632b417f066a004737d7d0161d8f32974bc281c4`.
- Ruff prescrito (`F401,F811,F821,E722,E741`) sobre los archivos modificados pasó. El chequeo global prescrito conserva sólo el baseline conocido `F401 math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`. Un chequeo Ruff adicional sin selector salió con 38 infracciones preexistentes de estilo/importación; la comparación de resultados contra HEAD confirmó cero hallazgos nuevos en los archivos editados, por lo que no se ejecutó autofix ni refactor oportunista.
- Evidencia: [`evidence/view-menu-screen-bounds-qt-offscreen-800x800.png`](evidence/view-menu-screen-bounds-qt-offscreen-800x800.png) muestra una captura compuesta del `QMenu.grab()` real dentro del rectángulo disponible de 800×800; la marca y la primera y última fila son visibles.

### Límite de verificación en pantalla nativa

El backend offscreen no permite configurar `QScreen` a los tres modos solicitados; esos valores se probaron como tamaños de ventana, comparando siempre el popup real contra `availableGeometry`. En una inspección exploratoria separada con XCB, Qt reportó `availableGeometry=(0,0,1429,740)`, DPR `1.34375`, y el menú raíz midió `(347,65,271,271)`; la primera fila `(356,78,253,35)` y el menú completo cabían. No obstante, esa automatización adicional con QTest no abrió de forma fiable todos los submenús en XCB y el proceso de captura terminó con SIGSEGV durante el teardown de Qt. No se trata como prueba aprobada ni se alteró código de lifecycle para ocultarlo; la suite oficial usa offscreen y la validación manual del propietario sigue siendo obligatoria.

### Publicación y CI por shards

- Publicado por push normal únicamente a `fix/beta-final-ui-clean2d`, commit de implementación `9b7bd010643de119e27f0824d600db79da8ee788`; `origin/fix/beta-final-ui-clean2d` coincidió con el HEAD local.
- GitHub Actions [run 38842965169](https://github.com/PJGV333/Chemuson/actions/runs/38842965169), para ese SHA: **success**. El plan recogió los 2238 IDs y distribuyó los ocho shards `[280,280,280,280,280,280,279,279]`; todos los shards y la verificación de IDs/JUnit pasaron. El shard más largo tardó 3m21s. También pasaron `windows-smoke` y `flatpak-smoke`.
- No se integró en `release/v0.3.0-beta.1-prep` ni en `main`; no hubo tag, Release ni publicación de paquetes. La revisión manual del propietario sigue pendiente; no declarar beta lista hasta recibirla.
