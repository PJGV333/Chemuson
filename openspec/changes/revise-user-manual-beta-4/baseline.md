# Baseline y registro de discrepancias — manual de usuario 0.3.0-beta.4

## Estado antes del cambio

- **Rama de trabajo:** `docs/manual-usuario-beta4-campana`, creada desde `main`.
- **HEAD de `main`:** `bc39ef35c5b7baac4bbbc4a7e9a395528f747f99` («docs(manual): redirigir guia historica hacia manual actual beta.4»).
- **Tag de referencia:** `v0.3.0-beta.4`.
- **Equivalencia de código:** `git diff --stat v0.3.0-beta.4..main -- src` no devuelve diferencias. El código de `main` es el de `0.3.0-beta.4`; las diferencias están solo en `README.md`, `README.en.md`, `docs/MANUAL_USUARIO.md` y `docs/manual_usuario.md`.
- **`git status --short` antes del cambio:** árbol limpio.
- **Manual afectado:** `docs/MANUAL_USUARIO.md`, 1201 líneas; su propia introducción declaraba que los capítulos detallados procedían de una auditoría del commit `140a080c336515650bbeea0a4e6ead67a9999b23` y que «no han sido revalidados íntegramente» para `0.3.0-beta.4`.
- **Escala del desfase:** `git log --oneline 140a080c..v0.3.0-beta.4 -- src/chemuson/gui` devuelve 66 commits de interfaz entre la base del manual y la versión publicada.

## Método de verificación

1. Lectura del código de la interfaz en el tag: `main_window_ui_builder.py`, `actions/*.py`, `toolbar.py`, `tool_rail.py`, `canvas/canvas_context_menu.py`, `canvas/*`, `side_panel.py`, `docks.py`, `energy_diagrams.py`, `orbitals.py`, `diagram_presets.py`, `clean2d/`, `chemname/`.
2. Ejecución real de la aplicación con `QT_QPA_PLATFORM=offscreen` sobre ese código, volcando menús, acciones con atajo, paletas del *rail* y diálogos, y capturando los widgets.
3. Los dos procedimientos son reproducibles: el generador de capturas queda versionado en `docs/manual-figures/capturar_figuras.py`.

## Discrepancias encontradas (antes de corregir)

### A. Afirmaciones falsas o insostenibles en el manual

| # | Afirmación del manual | Realidad verificada en 0.3.0-beta.4 |
|---|---|---|
| A1 | §11.5 «Accesos directos: C, N, O, S, P, Cl, Br, I, H» para elegir elementos | Esas letras son **celdas de la paleta Átomos**. Como teclas, `C` activa la herramienta Átomos y `N` la de Flechas. No existe selección de elemento por tecla. |
| A2 | §11.7 y §14.5 «Rotación 3D precisa: solicita parámetros exactos» presentada como herramienta suelta con `Alt+Arrastre` y `Alt+0` | La función existe, pero vive en la **paleta de la herramienta Seleccionar** (`tool_rotate_3d_precise`) y abre el diálogo `TrackballRotationDialog` (Ángulo X/Y, paso, Restablecer). `Alt+Arrastre` y `Alt+0` son gestos de rotación interactiva, no la entrada del diálogo. |
| A3 | §12.7 «varios orbitales f» e híbridos `sp`/`sp²`/`sp³` en la herramienta Orbitales | La paleta de **Orbitales** tiene 23 glifos y **ninguno `f`**; el híbrido disponible es **Lóbulo sp** (un lóbulo), `sp3` y `sp2`. Los subniveles `f` y el híbrido `sp` sí existen en la paleta de **Diagramas de energía** (8 presets). |
| A4 | §12.6 constructores de diagramas de energía listado con ocho entradas | **Electronic Diagrams** ofrece tres constructores: `Atomic diagram...`, `Diatomic MO diagram...`, `Ligand field diagram...`; **Electronic Diagram Presets** despliega **Atomic Presets**, **Diatomic MO Presets** y **Ligand Field Presets**. |
| A5 | §26 tabla de atajos con `Ctrl+Q` (Salir) | `action_quit` usa `QKeySequence.StandardKey.Quit`, que **resuelve a cadena vacía** en esta plataforma: la entrada no muestra atajo. Tampoco existe `Ctrl+W` para cerrar documento. |
| A6 | §26 `Ctrl++` como superíndice alternativo | `Ctrl++` está ocupado por **Zoom +**; el superíndice es `Ctrl+Shift+=`. |
| A7 | §9 «Ayuda → Manual de usuario» | El menú **Ayuda** solo tiene **Guía rápida...**, **Buscar actualizaciones...** y **Acerca de ChemUSON...**. |
| A8 | §16/§17/§18 describiendo flujos propios para paneles, placas y plantillas | Contenido duplicado de §7.4, §12, §13 y §22; las capturas enlazadas eran de `docs/ui-modernization/`, previas a `0.3.0-beta.4`. |
| A9 | §15.4 «matrices comparativas de representación geométrica» | Eran galerías de pruebas internas de desarrollo (`tools/visual/render_*_gallery.py`), no material didáctico. |

### B. Duplicaciones internas

| # | Duplicación | Decisión |
|---|---|---|
| B1 | Procedimiento de plantillas en §7.4, §18.1 y §22 | Un único capítulo (13) con enlaces internos. |
| B2 | Barras flotantes y simbología en §12 y §16 | §16 absorbido en el capítulo de anotaciones (7). |
| B3 | Placas TLC/gel en §12.5 y §17 | Un único capítulo (9). |
| B4 | Gestión de paneles en §12, §18 y §22 | Capítulo 12 limitado a organización, activación y navegación; el contenido de cada vista se explica en el capítulo de su función. |
| B5 | Limpieza 2D en §7.1, §12 y §20 | Único desarrollo en §4.4.1. |
| B6 | Menú contextual repetido en §12, §13 y §19 | Único desarrollo en §11, con submenús por objeto. |

### C. Hallazgos de comportamiento (se documentan, no se modifican)

| # | Hallazgo | Evidencia |
|---|---|---|
| C1 | **`Estructura → Limpiar 2D`** (sin atajo) es conductualmente idéntico a **`Limpiar 2D (1 paso)`** (`Ctrl+K`): `Clean2DController.run_clean_2d` ignora `step_ratio`, `fallback_iterations` y `status_suffix`, y ambas rutas usan `mode="quick"`. | `src/chemuson/clean2d/controller.py` |
| C2 | Resultados de análisis se **insertan como texto en el lienzo** (deshacible) con etiquetas en inglés: `Chemical Formula:`, `Exact Mass:`, `Molecular Weight:`, `Elemental Analysis:`. | `canvas_structure.py::_insert_analysis_text` |
| C3 | **Análisis elemental** tabula ~17 elementos; el resto devuelve `N/D`. Desde el menú contextual inserta texto; el diálogo interactivo solo está en el menú Estructura y usa la fórmula de **todo** el documento. | `canvas_context_menu.py`, `analysis/*` |
| C4 | Menú contextual ofrece **Nombre (IUPAC)**; el submenú **Estructura → Análisis** no. | volcado de ambos menús |
| C5 | **`Nombre a estructura...` consulta PubChem con `allow_network=True` por defecto.** | `gui/background_workers.py:45` |
| C6 | **Ver → Mostrar barra de herramientas** en realidad muestra/oculta la **barra de menús**. | `main_window_ui_builder.py:451` |
| C7 | `.sdf` rechazado al abrir; no hay arrastrar y soltar; **Exportar SMILES...** copia al portapapeles y no escribe archivo. | diálogo de apertura, `file_actions.py` |
| C8 | **CompChem** tiene seis botones (**Generar conformero 3D**, **Optimizar**, **Reset/regenerar**, **Proyectar a 2D**, **Exportar XYZ**, **Exportar input ▾** con ORCA/Gaussian/NWChem), no tres. | `gui/docks.py` (CompChemDock) |
| C9 | **Espectroscopía** no tiene botón «Analizar»: la predicción se recalcula sola y ofrece **Copiar tabla** y **Exportar tabla...**. Cabecera informa `heuristic-v1` y confianza global. | `gui/docks.py::SpectroscopyDock`, `main_window.py::_refresh_spectroscopy_dock` |
| C10 | El **Asistente molecular** está expuesto en el menú Estructura sin ninguna marca de experimental en la interfaz. | volcado del menú Estructura |

### D. Omisiones del manual anterior

Sin documentación de: submenú **Copiar como**, **Redimensionar selección...**, Grosor de trazo, **Ver → Visualización** completo, **Numeración**, **Tamaño de lienzo**, **Paletas**, gestos de creación de enlaces/cadenas/anillos, atajos numéricos, resultados del análisis como texto en el lienzo, política de guardado automático (120 s + 8 s, 5 copias), centros de recuperación con tres pestañas, y la distinción entre pegado de SMILES, Molfile e imagen.

## Decisiones de documentación

1. `Ctrl+K` pasa a presentarse como **método recomendado**; `Ctrl+Shift+K` y `Ctrl+Alt+K` como alternativas con propósito distinto; `Limpiar 2D` se documenta como entrada heredada equivalente (hallazgo C1).
2. No se añade ninguna entrada de menú falsa: se enlaza el manual en la sección Ayuda y se registra la entrada «Ayuda → Manual de usuario» como **propuesta independiente**, sin implementarla.
3. Las figuras se regeneran contra 0.3.0-beta.4. Las galerías de `docs/visual/` dejan de enlazarse desde el manual pero se conservan en el repositorio porque varios cambios de `openspec/changes/` las citan como evidencia.
4. El visor 3D (OpenGL) no se pudo capturar de forma fiable sin pantalla: se documenta sin figura y se declara el hueco.

## Riesgos de regresión

Nulos en ejecución: el cambio no toca `src/`, `tests/`, `packaging/` ni workflows. El único riesgo es documental (enlaces rotos, cifras de capítulo desactualizadas), cubierto con el validador de enlaces e imágenes y `git diff --check`.
