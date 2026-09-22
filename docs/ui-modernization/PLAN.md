# Plan de Modernización de la UI de Chemuson

**Fecha:** 2026-09-21
**Estado:** Propuesta (pendiente de aprobación y de OpenSpec por fase)
**Entregable complementario:** [`mockup-ui.html`](./mockup-ui.html) (maqueta interactiva de la propuesta)

---

## 0. Resumen ejecutivo

- **Framework actual: PyQt6** (no PySide). La dependencia está en `requirements.txt` y `pyproject.toml`, y todos los imports de la GUI son `from PyQt6...`.
- La sensación de "viejo" **no viene de PyQt6**: Qt6 es moderno y permite UIs de primera calidad. El problema está en (a) iconos dibujados a mano con `QPainter` (`gui/icons.py`, 1320 líneas) que se ven inconsistentes, (b) una distribución densa de funciones entre menús, 2 barras laterales de iconos con submenús, una barra superior y 7 docks ocultos por defecto, y (c) un sistema de temas QSS que existe pero no está gobernado por tokens de diseño.
- **Decisión propuesta: NO cambiar de framework.** Se moderniza **dentro de PyQt6 (QtWidgets)** con: un sistema de design tokens, iconografía SVG curada, una barra de aplicación unificada con pestañas de documento, un panel de herramientas unificado (rail + flyouts), paneles laterales organizados en tabs y una paleta de comandos (Ctrl+K).
- El trabajo se ejecuta en **8 fases**, cada una como un cambio OpenSpec independiente con baseline/verificación, de forma que la app sigue funcional y con la suite verde después de cada fase.

---

## 1. Diagnóstico del estado actual

### 1.1 Mapa de la UI actual

| Región | Widget | Contenido |
|---|---|---|
| Barra de menús | `QMenuBar` (`main_window_ui_builder.py`) | Archivo, Editar, Ver, Estructura, Reacción (placeholder), Ayuda |
| Toolbar superior | `QToolBar "Principal"` (`build_main_toolbar`) | Nuevo/Abrir/Guardar · Deshacer/Rehacer · Rotar/voltear · Limpiar 2D · SMILES. Sección auxiliar (copiar/pegar/zoom) oculta por defecto |
| Toolbar izquierda | `ChemusonToolbar` (`toolbar.py`, vertical) | Seleccion▾ · Enlace▾ · Cadena · Anillo▾ · Átomo▾ · Centro de coordinación · Rotación 3D |
| Toolbar derecha | `SymbolPaletteToolbar` (`toolbar.py`, vertical) | Texto▾ · Corchetes▾ · Flechas▾ · Placas▾ · Símbolos▾ · Diagramas de energía▾ · Orbitales▾ |
| Toolbar de texto | `TextFormatToolbar` (superior, contextual) | Fuente, tamaño, negrita, sub/superíndice, alineación, color, opacidad |
| Centro | `QTabWidget` document-mode (`shell/assembly.py`) | Pestañas de documento (cierre, reorden) con `ChemusonCanvas` (QGraphicsView) |
| Docks derecha | 7 × `QDockWidget` (`docks.py`) | Plantillas, Inspector, Validación, Propiedades químicas, Espectroscopía, CompChem 3D, Apariencia — **todos ocultos por defecto** (solo accesibles desde menú Ver) |
| Barra de estado | `QStatusBar` | Nombre IUPAC + carga total |

**Hechos clave sobre el look actual:**

1. **Iconos:** 100 % generados por `QPainter` en `gui/icons.py` (círculos con letras, glifos a mano, flechas, etc.). No hay ningún asset SVG ni asset de iconos para la UI en el repo (los únicos PNG existentes — `assets/baseline/orbitals_palette.png`, `src/repro_v2.png`, `tests/archive/…` — son baselines de regresión, no interfaz). Se regeneran al cambiar de tema (coste de CPU innecesario y look frágil).
2. **Temas:** ya existe `gui/styles.py` con paletas claro/oscuro (slate + acento cyan `#0891B2`/`#22D3EE`) y QSS generado por f-string. Es una base buena pero: los colores están dispersos (también en `icons.py`), no hay tokens de espaciado/tipografía, y varias piezas (docks, tab widget) heredan estilos del sistema.
3. **Descubrimiento:** la mitad de las herramientas (espectroscopía, compchem, validación, propiedades…) vive oculta en el menú *Ver*; las paletas de enlace/anillo/átomo/flechas solo se conocen abriendo cada submenú.
4. **Pestañas:** `QTabWidget` estándar por debajo de la toolbar; sin integración visual con la barra de la app.

### 1.2 Archivos involucrados (estado actual)

```
src/chemuson/gui/
├── main_window.py                  (2164 l)  ventana principal + handlers
├── main_window_ui_builder.py       (536 l)   menús + toolbar principal
├── toolbar.py                      (1517 l)  ChemusonToolbar + SymbolPaletteToolbar
├── text_toolbar.py                 toolbar contextual de texto
├── docks.py                        (940 l)   7 docks
├── styles.py                       (773 l)   paletas + QSS
├── icons.py                        (1320 l)  iconos QPainter
├── shell/assembly.py               (271 l)   composición del shell (región central, docks, toolbars, status)
├── tab_manager.py                  (265 l)   pestañas + autosave
└── canvas/…                        M09: editor (NO se toca la lógica)
```

### 1.3 ¿Por qué no cambiar de framework?

| Opción | Veredicto | Razón |
|---|---|---|
| Migrar a PySide6 | ❌ | Mismo Qt6; solo rompería todos los imports y tests sin ganar nada visual. |
| QtQuick/QML para todo | ❌ | El corazón de la app es un `QGraphicsView` con input quirúrgico (bonds, selección, drag). QML no aporta para ese componente y mixear dos paradigmas duplicaría estado. |
| Frontend web (QtWebEngine) | ❌ (por ahora) | Riesgo alto sobre el subsistema más protegido (AGENTS.md §4), pérdida de precisión de input y de empaquetado actual. Se puede reevaluar **solo para diálogos** en una fase posterior. |
| **Modernizar dentro de PyQt6 (QtWidgets)** | ✅ | QSS + widgets personalizados + iconos SVG cubren el 100 % del objetivo visual con riesgo controlado. |

---

## 2. Diseño objetivo

### 2.1 Nueva estructura de la ventana (racional y ordenada)

```
┌────────────────────────────────────────────────────────────────────────────┐
│  [◆ Chemuson]  [ Sin título * | cafeina.cmsn | + ]   [Buscar o ejecutar ⌘K] [↶ ↷ | 🌓 | ⚙]   ← Barra de aplicación (1 sola)
├──┬──────────────────────────────────────────────────────┬──────────────────┤
│ ▸│                                                      │  Inspector      │
│ S │   Lienzo (QGraphicsView) sobre fondo de trabajo,    │  Validación     │
│ e │   hoja blanca con sombra, rejilla opcional,         │  Propiedades    │
│ l │   pill de zoom abajo a la derecha.                  │  Plantillas     │
│ e │                                                      │  Apariencia / … │
│ c │                                                      │  (docks en tabs)│
│ c │                                                      │                  │
│ i │                                                      │                  │
│ ó │                                                      │                  │
│ n│                                                      │                  │
│ ─ │  Panel de herramientas (rail + flyouts):            │                  │
│ D │  Seleccion · Enlace · Anillo · Cadena · Átomo ·     │                  │
│ i │  Rotación 3D · Texto · Flechas · Corchetes ·        │                  │
│ b │  Símbolos · Diagramas · Orbitales · Placas ·        │                  │
│ u │  Limpiar 2D · Validar                                │                  │
│ j │                                                      │                  │
├──┴──────────────────────────────────────────────────────┴──────────────────┤
│  ✓ Seleccionar (V)                                        C₇H₇NO₂ · IUPAC · Carga 0 · ✓ autosave 14:32 │
└────────────────────────────────────────────────────────────────────────────┘
```

Principios:

1. **Una sola barra de aplicación** marca + pestañas de documento + búsqueda/comandos + acciones globales (deshacer/rehacer, tema, ajustes). El menú `QMenuBar` se mantiene (accesible por teclado y por compatibilidad) pero queda como referencia; lo frecuente vive arriba.
2. **Pestañas de documento integradas en la barra** (no una franja aparte): título + punto de suciedad + cerrar; tab `+` para nuevo.
3. **Un solo panel de herramientas** a la izquierda, agrupado por tarea: *Seleccionar · Dibujar · Anotar · Diagramas · Placas · Acciones*. Cada grupo expone sus opciones como **flyout en cuadrícula** (no como submenús nativos): visible, navegable por teclado y con atajos (`V` seleccionar, `B` enlace, `R` anillo, `C` átomo, `T` texto…).
4. **Paneles de la derecha en tabs** (docks reorganizables): Inspector, Validación, Propiedades, Plantillas, Apariencia como tabs visibles; Espectroscopía, CompChem y el resto accesibles desde un tab `…`. El Inspector se muestra por defecto (hoy todo está oculto).
5. **Barra de estado moderna**: izquierda = herramienta activa + atajo + posición del cursor; derecha = fórmula + IUPAC + carga + indicador de autosave.
6. **Paleta de comandos (Ctrl+K)**: búsqueda sobre todas las acciones (archivos, vista, estructura, análisis, plantillas, docks). Es la mayor ganancia de "práctico" por esfuerzo.
7. **Estado vacío del lienzo**: mensaje de onboarding ("Dibuja con la herramienta Enlace o importa un SMILES") con acciones clicables; desaparece al primer trazo.

### 2.2 Sistema de design tokens

Nuevo subpaquete `src/chemuson/gui/theme/` (se registrará en `architecture/modules.yml` como parte del OpenSpec; dentro del alcance del módulo M08 de GUI):

| Token | Claro | Oscuro | Uso |
|---|---|---|---|
| `bg-app` | `#F1F5F9` | `#0B1120` | fondo de ventana |
| `bg-surface` | `#FFFFFF` | `#0F172A` | paneles, docks, menú |
| `bg-surface-2` | `#F8FAFC` | `#1E293B` | filas alternas, hover suave |
| `border` / `border-strong` | `#E2E8F0` / `#CBD5E1` | `#2B3A55` / `#334155` | separadores |
| `text-1/2/3` | `#0F172A`/`#475569`/`#94A3B8` | `#F1F5F9`/`#CBD5E1`/`#64748B` | jerarquía tipográfica |
| `accent` | `#0E7490` | `#22D3EE` | acción principal, herramienta activa |
| `accent-soft` | `#ECFEFF` | `rgba(34,211,238,.12)` | relleno de selección activa |
| `danger/warn/ok` | `#DC2626`/`#D97706`/`#059669` | variantes claras | severidades de validación |

- **Tipografía:** familia del sistema (Inter en Linux/macOS, Segoe UI en Windows), 13 px base, 12 px en tablas/badges, `tabular-nums` en la barra de estado.
- **Espaciado:** grilla de 8 px (márgenes 8/12/16), radio 8–10 px en superficies, 6 px en botones.
- **Elevación:** máx. 2 niveles (hoja del lienzo y flyouts/paleta de comandos) con sombras suaves; el resto se separa con bordes.
- `styles.py` se reescribe para **consumir los tokens** (función única `get_main_stylesheet(theme)`); `icons.py` deja de contener colores propios.
- Modo **seguir sistema**: Qt6 ya expone el modo del sistema; se agrega la opción a Preferencias (hoy hay claro/oscuro manual).

### 2.3 Iconografía

- **Nuevo set SVG** (licencia libre: Tabler/Lucide/Phosphor o equivalentes CC0, o set propio) embebido como recursos: `src/chemuson/gui/theme/icons/*.svg`, rejilla 24 px, trazo 1.75 px, esquinas redondeadas, `currentColor` donde sea posible.
- `IconProvider` (`theme/icon_provider.py`): carga SVG → `QIcon` con caché por (nombre, tamaño, tema); tamaños 16/20/24/28 px; soporte HiDPI.
- **`icons.py` se mantiene como fachada compatible** (`draw_generic_icon`, `draw_atom_icon`, …) delegando en el provider; así ningún caller de M08/M09 se rompe y la migración es incremental.
- **Fichas de elemento químico** (C, N, O, …) y anillos: se generan como SVG dinámico (plantilla parametrizada: símbolo + color CPK actualizado), con el mismo lenguaje visual que el set (trazo fino, sin círculo relleno saturado).
- Se elimina el redibujado al cambiar de tema (los SVG se teñen una vez por tema y se cachean).

### 2.4 Paleta de comandos

`CommandPalette` (QWidget modal, 560 px): input con filtro + lista de acciones con icono, título, atajo y sección. Fuentes de acciones:

- Todas las `QAction` del menubar (registradas con `data-section`).
- Plantillas del `TemplateBrowserService`.
- Docks (mostrar/ocultar).
- Vistas/tema (claro/oscuro/seguir sistema), zoom, limpiar 2D, validación, análisis.

---

## 3. Plan de ejecución paso a paso

> Cada fase = **1 cambio OpenSpec** (`openspec/changes/2026-MM-DD-<slug>`) con `proposal.md`, `design.md`, `tasks.md` y `specs/`. Antes de cada fase: baseline (`git status`, `compileall`, `pytest --collect-only`, `pytest -q`, ruff F401/F811/F821/E722/E741) registrado en `baseline.md`. Al final de cada fase: suite completa + smoke Qt offscreen.

### Fase 0 — Preparación y decisiones de diseño (esfuerzo: S)

1. Crear el OpenSpec `moderna-ui-foundation` (o uno por fase) y aprobar con el equipo.
2. Capturar baseline completa y **screenshots de referencia de la UI actual** (claro/oscuro) en `docs/ui-modernization/baseline/` para comparar después.
3. Decidir licencia/origen del set de iconos y meter los SVG en el repo (asset, sin dependencias nuevas).
4. Decidir nombres públicos: mantener las señales `tool_changed(str)` y los `tool_id` existentes (contrato con M09) — **no renombrar nada**.
5. Definir el registro en `architecture/modules.yml`: subpaquete `gui/theme` dentro de M08 (o nuevo id si el comité prefiere), sin dependencias circulares nuevas.

**Criterio de aceptación:** OpenSpec validado (`openspec validate --strict`), baseline guardada, set de iconos en el repo, sin cambios de código aún.

### Fase 1 — Design tokens y QSS (esfuerzo: M)

1. `theme/tokens.py`: dict de tokens claro/oscuro (la tabla §2.2).
2. `theme/qss.py`: generadores de hoja de estilo a partir de tokens (`get_main_stylesheet`, `get_tool_palette_stylesheet`, `get_dialog_stylesheet`); migrar las reglas actuales de `styles.py` (que pasa a ser fachada de compatibilidad) y añadir: pestañas, dock titles, tooltips, estados vacíos.
3. `theme/__init__.py`: `apply_theme(app_or_window, theme_name)` + `set_theme_from_system()`.
4. Preferencias: opción "Seguir sistema" (persistida vía `platform.settings`, M21 — sin cambiar su API).
5. Tests: fixture de tokens (claro/oscuro se resuelven a colores válidos), smoke Qt offscreen con ambos temas.

**Criterio de aceptación:** la app luce los dos temas sin colores hardcoded fuera de tokens; suite verde; `styles.py` solo re-exporta.

### Fase 2 — Sistema de iconos SVG (esfuerzo: M)

1. Inventario de iconos actuales (funciones de `icons.py` × call sites) → mapa 1:1 a nombres SVG.
2. `theme/icon_provider.py` + carpeta `theme/icons/*.svg` (~60–80 iconos).
3. Reescribir `icons.py` como fachada: mismas firmas, delega en el provider; `draw_atom_icon`/`draw_ring_icon`/`draw_charge_icon` generan SVG dinámico con las mismas paletas CPK (mejoradas) y mismo resultado conceptual.
4. Reemplazo de `QIcon.fromTheme` + pixel-loop de tint (lento) por tint de SVG (recoloreo por path, sin iterar píxeles).
5. Tests: provider cachea (misma instancia por clave), tamaños 16/20/24/28 no nulos, fachada devuelve `QIcon` no nulos para todos los nombres del inventario.

**Criterio de aceptación:** ningún icono `QPainter` manual queda fuera de la fachada; cambio de tema no redibuja (mismo `QIcon` cacheado); suite verde.

### Fase 3 — Barra de aplicación y pestañas de documento (esfuerzo: L)

1. `gui/document_tabs.py`: `DocumentTabBar(QWidget)` — tabs con icono de documento, título truncado, punto de suciedad (sincronizado con `CanvasTabManager.update_tab_title`), botón cerrar, tab `+` "Nuevo"; drag-reorder preservado.
2. `gui/app_bar.py`: barra superior compuesta (logo, `DocumentTabBar`, pill de búsqueda, undo/redo, tema, ajustes); reemplaza visualmente la toolbar "Principal" manteniendo las mismas `QAction` (mismos atajos de teclado).
3. `shell/assembly.py`: composición nueva (app bar en lugar de menubar+toolbar superior); el `QMenuBar` se conserva pero puede ocultarse opcionalmente (preference) — la accesibilidad por teclado (`Alt`) sigue funcionando.
4. `main_window_ui_builder.py`: migrar acciones de "toolbar principal" a la app bar; sección auxiliar (copiar/pegar/zoom) se integra siempre visible.
5. Tests: smoke offscreen — la ventana expone app bar, tabs, acciones; abrir/cerrar/reordenar pestañas; suciedad refleja en el título.

**Criterio de aceptación:** toda acción previa sigue alcanzable (menú o app bar o Ctrl+K), atajos intactos, suite verde + smoke.

### Fase 4 — Panel de herramientas unificado (esfuerzo: L — la fase más visible)

1. `gui/tool_panel.py`: rail vertical (56–64 px) con grupos (Seleccionar, Dibujar, Anotar, Diagramas, Placas, Acciones) y botones con icono SVG; herramienta activa resaltada con `accent-soft` + borde `accent`.
2. `gui/flyout.py`: flyout reutilizable (cuadrícula 3–4 columnas, ítem = icono + etiqueta, estado activo, búsqueda opcional, atajos) que reemplaza los `QMenu` de paleta de `ChemusonToolbar`/`SymbolPaletteToolbar`.
3. Migrar paletas existentes **conservando `tool_id` y señales**: enlace (11 estilos), anillo (benceno + 3–12 + personalizado), átomo (10 + tabla periódica), flechas (15), corchetes (10), símbolos (cargas/radicales/pares de electrones), diagramas de energía, orbitales, placas (TLC/gel).
4. Atajos de herramienta nuevos y documentados: `V` seleccionar, `B` enlace, `R` anillo, `C` átomo, `T` texto, `N` flecha, `G` corchetes, `E` diagramas de energía, `O` orbitales.
5. `toolbar.py` queda como fachada (o se elimina si todos los call sites migran, registrado en el OpenSpec) — M09 y canvas **no cambian** (solo reciben `tool_id`).
6. Tests: cada `tool_id` existente tiene un camino de activación; flyout muestra todas las entradas de cada paleta; selección exclusiva (una herramienta a la vez); smoke.

**Criterio de aceptación:** paridad 1:1 de herramientas con la UI actual (checklist del inventario §Fase 2), sin `tool_id` huérfano, suite verde.

### Fase 5 — Paneles laterales (tabs) y barra de estado (esfuerzo: M)

1. Contenedor de paneles: `gui/side_panel.py` — `QTabWidget` dockable que aloja los docks existentes **sin reescribir su contenido** (los widgets de `docks.py` se reutilizan tal cual; solo cambia el contenedor y el título/ícono del tab).
2. Tabs por defecto: Inspector (visible), Validación, Propiedades, Plantillas, Apariencia + `…` (Espectroscopía, CompChem). Persistencia de tab activo y orden en `platform.settings`.
3. `gui/status_bar.py`: barra de estado con widget permanente izquierdo (herramienta + atajo + cursor) y derecho (fórmula, IUPAC, carga, autosave) — mismo contenido de hoy, mejor jerarquía visual.
4. Estado vacío del lienzo: overlay decorativo en escena vacía con 2 acciones (Dibujar enlace / Importar SMILES) que no interfieren con hit-testing (se agrega vía escena pública existente o item nuevo catalogado; M09 no modifica su lógica).
5. Tests: docks abren/cierran desde tabs y desde menú Ver (back-compat); contenido de Inspector/Validación idéntico al de hoy (tests existentes de docks pasan sin cambios); smoke.

**Criterio de aceptación:** los 7 docks siguen funcionales; los tests de docks actuales pasan sin modificación; suite verde.

### Fase 6 — Paleta de comandos Ctrl+K (esfuerzo: M)

1. `gui/command_palette.py`: registro de acciones (`register(action, section, keywords)`) alimentado desde `main_window_ui_builder` + `TemplateBrowserService` + docks.
2. Filtro: substring por defecto + ranking simple por prefijo; navegación ↑↓, Enter ejecuta, Esc cierra; último resultado usado como atajo rápido.
3. Integración: pill de búsqueda en la app bar y atajo global `Ctrl+K`.
4. Tests: registro/filtro/ejecución de acciones; suite verde.

**Criterio de aceptación:** ≥ 60 acciones buscables; todos los docks y las exportaciones alcanzables desde la paleta; smoke.

### Fase 7 — Pulido (esfuerzo: M)

1. Tooltips uniformes (nombre + atajo) en rail, app bar y flyouts; estado "deshacer" deshabilitado visible.
2. HiDPI: verificación de iconos y QSS a 125/150/200 % (tests offscreen con `QT_SCALE_FACTOR=2`).
3. Oscuro pulido: contraste AA en texto de secundarios; canvas sigue siendo hoja blanca (consistente con export PNG).
4. Onboarding: overlay de bienvenida 1ª vez (3 puntos: rail, canvas, panel derecho) con "No mostrar de nuevo".
5. Alineación fina del mockup vs implementación (comparar con `docs/ui-modernization/mockup-ui.html`).
6. Actualizar `manual_usuario.md` con capturas nuevas (claro y oscuro).

**Criterio de aceptación:** checklist visual firmado por el equipo; manual actualizado; suite verde.

### Fase 8 — QA final y release (esfuerzo: S)

1. Regresión completa: `pytest -q` (1492 passed/55 skipped como referencia), `ruff` scoped, smoke Qt offscreen, AppImage local manual.
2. Comparativa de screenshots antes/después (contra la baseline de Fase 0) en `docs/ui-modernization/after/`.
3. Notas de release 0.4.0-dev: "Nueva interfaz" + guía de atajos.
4. Archivar los OpenSpecs de las fases 1–7.

**Criterio de aceptación:** release candidate con la nueva UI, sin regresiones de comportamiento (canvas, Clean2D, nomenclatura, persistencia intactos).

### Dependencias y orden

```
F0 (preparación) → F1 (tokens/QSS) → F2 (iconos) ─┐
                       │                          ├→ F3 (app bar/tabs)
                       └──────────────────────────┤      │
                                                  └→ F4 (tool panel) ─ F5 (side panel/status)
                                                                            │
                                                                    F6 (Ctrl+K) → F7 (pulido) → F8 (QA)
```

F2 y F3 son independientes entre sí (pueden paralelizarse con dos desarrolladores). F4 depende de F2 (iconos). F5 depende de F1. F6 depende de F3.

### Estimación global

| Fase | Esfuerzo | Riesgo |
|---|---|---|
| 0 | S | — |
| 1 | M | bajo |
| 2 | M | bajo |
| 3 | L | medio (layout central) |
| 4 | L | medio (contrato `tool_id` con M09) |
| 5 | M | bajo (reutilización de docks) |
| 6 | M | bajo |
| 7 | M | bajo |
| 8 | S | — |

---

## 4. Riesgos y mitigaciones

| Riesgo | Mitigación |
|---|---|
| Romper el contrato `tool_id`/señales con M09 (canvas) | Fase 4 conserva los ids y señales; checklist de paridad 1:1; tests de activación por id. |
| QSS no soporta todo el CSS (transiciones, flex) | El diseño objetivo evita dependencias de CSS avanzado; flyouts y tabs son widgets propios pequeños. |
| Licencia del set de iconos | Solo sets MIT/CC0/ISC (Tabler=MIT, Lucide=ISC, Phosphor=MIT) o set propio; README de assets con atribución. |
| Docks existentes pierden comportamientos al moverse a tabs | Reutilizar los widgets de `docks.py` sin reescribirlos; tests de docks actuales como puerta. |
| Cambio de tema lento (redibujado) | Caché de `QIcon` por (nombre, tamaño, tema); test que verifica instancia cacheada. |
| Alcance scope-creep (mejoras de canvas) | Regla dura: ninguna fase toca `canvas/` (M09), `clean2d/`, `chemname/`, `chemio/persistence.py` (AGENTS.md §2.1 y §4). |

## 5. Guardarraíles (AGENTS.md) — qué NO cambia

- **Lógica del canvas/editor** (`src/chemuson/gui/canvas/`, `editor2d/`, M09): solo se le *habla* con los mismos `tool_id`.
- **Jerarquía de mixins y orden de despacho de eventos de Qt** en `main_window.py`.
- **Clean2D, ChemName, persistencia `.cmsn`**: intactos.
- **Nuevas dependencias externas:** cero (SVGs y QSS son assets propios).
- **`architecture/modules.yml`**: se actualiza (registro de `gui/theme`) dentro del OpenSpec de cada fase; excepciones `temporary_exceptions` solo si fueran estrictamente necesarias, documentadas.
- **Prohibido:** refactor oportunista fuera del alcance de cada fase; ocultar fallos de tests.

## 6. Criterio de éxito final

1. La app se ve **moderna y coherente** en claro/oscuro (tokens + iconos SVG uniformes), comparable al mockup.
2. Toda función de la UI actual sigue disponible **con el mismo o menor número de clics** (checklist de paridad).
3. Las acciones frecuentes están a ≤ 2 interacciones: herramienta en el rail, paneles en tabs, resto en Ctrl+K.
4. Suite completa verde + smoke Qt; **cero regresiones** en canvas, Clean2D, nomenclatura y persistencia.
5. Documentación (manual + atajos) actualizada y screenshots antes/después en `docs/ui-modernization/`.
